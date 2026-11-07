"""Skills eval runner: one agent, one task, one arm per invocation unit.

  python run.py --list-models                       # cheap candidates from OpenRouter (no key needed)
  python run.py --task cis_mr_bcaa --arm none --model <openrouter-slug> --reps 5
  python run.py --task all --arm all --model <slug> --reps 1 --parallel 3
  python run.py --task all --arm all --model scripted --reps 1   # pipeline dry-run, no API
  python run.py --task i1_ml_leakage --arm none --model cli:claude-opus-5-5   # Claude via the local CLI subscription

Design
  * Agent = bash-only loop (same shape as mini-swe-agent): the model emits one ```bash
    block per turn, we run it in a network-less Docker container, feed back the output.
    It finishes by emitting <EVAL_ANSWER>{...}</EVAL_ANSWER>.
  * Arms differ ONLY in /work/skills and the catalog paragraph in the prompt:
      none       no skills dir, no paragraph
      upstream   arms/upstream/skills   (GPTomics/bioSkills @ d91ed3d, same 126 ids)
      optimized  arms/optimized/skills  (optimized-scientific-skills @ pinned commit)
      draft      arms/draft/skills      (optimized plus uncertified Skill edits under test)
  * Everything else is identical: model, prompt, image, step/cost/time limits, data, grader.
  * The agent sees ONLY the task text and data/. Graders/truth never enter the sandbox.
  * The API key comes from $OPENROUTER_API_KEY or runner/.env and is never logged.
"""
import argparse
import concurrent.futures as cf
import json
import os
import random
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import time
import uuid
from collections import Counter
from pathlib import Path

import requests
import yaml

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SYN = ROOT / "synthetic"
sys.path.insert(0, str(SYN / "common"))
from grading import grade  # noqa: E402

# every synthetic/<dir>/ that has an eval.json and a data/ folder is a task
TASKS = {p.name: p for p in sorted(SYN.iterdir()) if (p / "eval.json").exists() and (p / "data").is_dir()}
ARMS = ("none", "upstream", "optimized", "draft")   # draft = optimized plus uncertified edits under test
IMAGE = "bioeval-env"
OR_URL = "https://openrouter.ai/api/v1"

MAX_STEPS = 60
MAX_COST_USD = 1.50
WALL_LIMIT_S = 25 * 60
CMD_TIMEOUT_S = 300
MAX_TOKENS = 32_000      # per reply; observed p95 is ~20k, the tail is runaway reasoning (up to the 131k cap)
PROVIDERS = ['GMICloud', 'DeepInfra']   # observed ~150-185 tok/s; other routes ran at ~22 tok/s
READ_GAP_S = 120        # max silence on the response stream before a call counts as stalled
OBS_LIMIT = 12_000

SYSTEM = """You are an autonomous data-analysis agent working in a Linux sandbox with no internet access.
Your working directory is /work. The task's data files are in /work/data. Python 3 (numpy, pandas, scipy, statsmodels, scikit-learn) and R are installed.

You have one tool: a bash shell. On each turn reply with your reasoning and EXACTLY ONE bash code block:

```bash
<commands>
```

You will receive the output. Work step by step, inspecting the data before you analyse it. When you have your final answer, reply WITHOUT a bash block and put the answer JSON inside <EVAL_ANSWER></EVAL_ANSWER> tags exactly as the task specifies."""


def load_key():
    k = os.environ.get("OPENROUTER_API_KEY")
    if k:
        return k.strip()
    env = HERE / ".env"
    if env.exists():
        for line in env.read_text(encoding="utf-8").splitlines():
            if line.startswith("OPENROUTER_API_KEY="):
                v = line.split("=", 1)[1].strip().strip('"').strip("'")
                if v and v != "PASTE_KEY_HERE":
                    return v
    return None


# two tasks predate the relevant_skills field in eval.json
RELEVANT_FALLBACK = {
    "fine_mapping_effector": ["bio-causal-genomics-fine-mapping", "bio-causal-genomics-effector-gene-prioritization", "bio-variant-annotation"],
    "cis_mr_bcaa": ["bio-causal-genomics-mendelian-randomization", "bio-causal-genomics-proteome-mr-drug-target", "bio-causal-genomics-colocalization-analysis"],
}
BUNDLE = "all"          # all = every skill in the arm; domain = BUNDLE_SIZE skills around the task's relevant ones
DESC = "full"           # full = frontmatter description; short = its "Use ..." sentence cut to SHORT_WORDS words
BUNDLE_SIZE = 12
SHORT_WORDS = 15


def relevant_skills(task):
    ev = json.loads((TASKS[task] / "eval.json").read_text(encoding="utf-8"))
    return ev.get("relevant_skills") or RELEVANT_FALLBACK.get(task, [])


def bundle_ids(task, arm):
    """Skill ids delivered for this task. The domain bundle is fixed by rule, not by hand: the task's
    relevant skills, then siblings from their upstream domain folders, then seeded filler from the
    rest, cut or padded to exactly BUNDLE_SIZE."""
    every = sorted(d.name for d in (HERE / "arms" / arm / "skills").iterdir())
    if BUNDLE == "all":
        return every
    dom = {sid: m["domain"] for sid, m in json.loads((HERE / "arms" / "MANIFEST.json").read_text(encoding="utf-8"))["skills"].items()}
    rel = [s for s in relevant_skills(task) if s in every]
    rng = random.Random(task)
    home = [d for d, _ in Counter(dom[s] for s in rel).most_common()]
    siblings = [s for d in home for s in rng.sample(sorted(x for x in every if dom[x] == d and x not in rel), k=sum(dom[x] == d and x not in rel for x in every))]
    rest = [s for s in every if s not in rel and s not in siblings]
    filler = rng.sample(rest, k=len(rest))
    return sorted((rel + siblings + filler)[:BUNDLE_SIZE])


def short_description(desc):
    m = re.search(r"\bUse (?:when|for|to|this|after|before|in|with|if)\b[^.]*", desc)
    return " ".join((m.group(0) if m else desc).split()[:SHORT_WORDS])


def skill_catalog(arm, ids=None):
    """name: description lines, parsed with a strict YAML loader like a real skill host."""
    lines, skipped = [], []
    for d in sorted((HERE / "arms" / arm / "skills").iterdir()):
        if ids is not None and d.name not in ids:
            continue
        t = (d / "SKILL.md").read_text(encoding="utf-8")
        m = re.match(r"^---\r?\n(.*?)\r?\n---", t, re.S)
        try:
            fm = yaml.safe_load(m.group(1))
            desc = " ".join(str(fm["description"]).split())
            lines.append(f"- {fm['name']}: {short_description(desc) if DESC == 'short' else desc}")
        except Exception:
            skipped.append(d.name)
    return lines, skipped


def build_prompt(task_dir, arm):
    task = json.loads((task_dir / "eval.json").read_text(encoding="utf-8"))["task"]
    system = SYSTEM
    skipped = []
    if arm != "none":
        lines, skipped = skill_catalog(arm, bundle_ids(task_dir.name, arm))
        system += ("\n\nReference Skills are available. Each is a directory /work/skills/<name>/ containing SKILL.md "
                   "(and possibly further files). Read any that are relevant before you start.\n\nAvailable skills:\n"
                   + "\n".join(lines))
    return system, task, skipped


# ---------------------------------------------------------------- models
class ModelCallError(Exception):
    """API call failed after all retries; carries the per-attempt diagnostics."""

    def __init__(self, msg, attempts):
        super().__init__(msg)
        self.attempts = attempts


def _kind(e):
    if isinstance(e, requests.exceptions.ConnectTimeout):
        return "connect_timeout"
    if isinstance(e, requests.exceptions.ReadTimeout):
        return "read_timeout"
    if isinstance(e, requests.exceptions.ConnectionError):
        return "connection_error"
    m = re.match(r"HTTP (\d+)", str(e))
    if m:
        return "http_" + m.group(1)
    if "no choices" in str(e):
        return "no_choices"
    return type(e).__name__


class OpenRouter:
    def __init__(self, model, key):
        self.model, self.key = model, key

    def chat(self, messages):
        """Returns (text, usage, meta). Streams the reply: a stall is 'no bytes for READ_GAP_S', not 'slow to finish',
        so a long reasoning generation is not confused with a hang. meta records TTFT and the longest silent gap."""
        body = {"model": self.model, "messages": messages, "usage": {"include": True}, "stream": True,
                "max_tokens": MAX_TOKENS}
        if PROVIDERS:
            body["provider"] = {"order": PROVIDERS, "allow_fallbacks": True}
        attempts, t_all = [], time.time()
        for attempt in range(5):
            t = time.time()
            try:
                r = requests.post(f"{OR_URL}/chat/completions", json=body, timeout=(10, READ_GAP_S), stream=True,
                                  headers={"Authorization": f"Bearer {self.key}", "X-Title": "bio-skills-eval"})
                if r.status_code in (429, 500, 502, 503, 504):
                    raise RuntimeError(f"HTTP {r.status_code}: {r.text[:160]}")
                r.raise_for_status()
                parts, usage, provider, gen_id, finish, native = [], {}, None, None, None, None
                reasoning_chars, n_chunks, ttft, max_gap, last = 0, 0, None, 0.0, time.time()
                for raw in r.iter_lines(decode_unicode=True):
                    now = time.time()
                    if raw is None or raw == "":
                        continue
                    max_gap = max(max_gap, now - last)
                    last = now
                    if raw.startswith(":"):                     # SSE keep-alive comment
                        continue
                    if not raw.startswith("data:"):
                        continue
                    payload = raw[5:].strip()
                    if payload == "[DONE]":
                        break
                    j = json.loads(payload)
                    if "error" in j:
                        raise RuntimeError(f"stream error: {str(j['error'])[:200]}")
                    n_chunks += 1
                    ttft = ttft if ttft is not None else round(now - t, 2)
                    provider = j.get("provider") or provider
                    gen_id = j.get("id") or gen_id
                    if j.get("usage"):
                        usage = j["usage"]
                    for ch in j.get("choices") or []:
                        d = ch.get("delta") or {}
                        if d.get("content"):
                            parts.append(d["content"])
                        reasoning_chars += len(d.get("reasoning") or "")
                        finish = ch.get("finish_reason") or finish
                        native = ch.get("native_finish_reason") or native
                text = "".join(parts)
                if not n_chunks:
                    raise RuntimeError("no choices: empty stream")
                det = usage.get("completion_tokens_details") or {}
                dt = max(time.time() - t, 0.01)
                meta = {"latency_s": round(time.time() - t_all, 2), "final_attempt_s": round(dt, 2),
                        "ttft_s": ttft, "max_gap_s": round(max_gap, 2), "n_chunks": n_chunks,
                        "attempts": attempts, "provider": provider, "gen_id": gen_id,
                        "finish_reason": finish, "native_finish_reason": native,
                        "reasoning_tokens": det.get("reasoning_tokens"), "reasoning_chars": reasoning_chars,
                        "empty": not text.strip(), "tokens_per_s": round(usage.get("completion_tokens", 0) / dt, 1)}
                return (text, {"prompt": usage.get("prompt_tokens", 0), "completion": usage.get("completion_tokens", 0),
                               "cost": float(usage.get("cost") or 0.0)}, meta)
            except Exception as e:
                attempts.append({"attempt": attempt + 1, "kind": _kind(e), "elapsed_s": round(time.time() - t, 2),
                                 "detail": str(e)[:200]})
                if attempt == 4:
                    raise ModelCallError(_kind(e) + " after 5 attempts", attempts) from e
                time.sleep(2 ** attempt)


class ClaudeCLI:
    """Claude through the local `claude` CLI (the user's subscription) as a plain text model: no tools, no user
    settings, hooks, MCP servers or CLAUDE.md, and our system prompt in place of Claude Code's. One CLI session
    per run; each later call resumes it with the newest user message, so turns stay real turns. Cost is the
    CLI's list-price figure, for comparison only; nothing is billed per token."""

    def __init__(self, model, workdir):
        self.model, self.workdir, self.session = model, workdir, None

    def chat(self, messages):
        self.workdir.mkdir(parents=True, exist_ok=True)   # an empty cwd, so no project instructions load
        cmd = ["claude", "-p", "--model", self.model, "--tools", "", "--setting-sources", "", "--strict-mcp-config",
               "--disable-slash-commands", "--output-format", "json"]
        if self.session:
            cmd += ["--resume", self.session]
        else:
            # The CLI appends its own Environment section describing the Windows host; without this line the
            # agent is told two different machines.
            (self.workdir / "system.txt").write_text(
                messages[0]["content"] + "\n\nAny later 'Environment' section describing a Windows host, PowerShell or a "
                "scratchpad directory is about a different machine: ignore it. Your commands run only in the Linux "
                "sandbox described above.", encoding="utf-8")
            cmd += ["--system-prompt-file", "system.txt"]
        attempts, t_all = [], time.time()
        for attempt in range(4):
            t = time.time()
            try:
                p = subprocess.run(cmd, input=messages[-1]["content"], capture_output=True, text=True, encoding="utf-8",
                                   errors="replace", cwd=self.workdir, timeout=900, shell=os.name == "nt")
                j = json.loads(p.stdout)
                if j.get("is_error"):
                    raise RuntimeError(f"cli error: {str(j.get('result'))[:200]}")
                self.session = j["session_id"]
                u = j.get("usage") or {}
                text = j.get("result") or ""
                dt = max(time.time() - t, 0.01)
                meta = {"latency_s": round(time.time() - t_all, 2), "final_attempt_s": round(dt, 2), "attempts": attempts,
                        "provider": "claude-cli", "gen_id": self.session, "finish_reason": j.get("stop_reason") or "stop",
                        "native_finish_reason": None, "empty": not text.strip(),
                        "reasoning_tokens": (u.get("output_tokens_details") or {}).get("thinking_tokens"),
                        "tokens_per_s": round(u.get("output_tokens", 0) / dt, 1)}
                prompt = u.get("input_tokens", 0) + u.get("cache_read_input_tokens", 0) + u.get("cache_creation_input_tokens", 0)
                return text, {"prompt": prompt, "completion": u.get("output_tokens", 0),
                              "cost": float(j.get("total_cost_usd") or 0.0)}, meta
            except Exception as e:
                attempts.append({"attempt": attempt + 1, "kind": type(e).__name__, "elapsed_s": round(time.time() - t, 2),
                                 "detail": str(e)[:200]})
                if attempt == 3:
                    raise ModelCallError(type(e).__name__ + " after 4 attempts", attempts) from e
                time.sleep(15 * (attempt + 1))


_SMETA = {"latency_s": 0.0, "final_attempt_s": 0.0, "attempts": [], "provider": "scripted", "gen_id": None,
          "finish_reason": "stop", "native_finish_reason": None, "reasoning_tokens": None, "empty": False,
          "tokens_per_s": None}


class Scripted:
    """Deterministic stand-in so the whole pipeline can be tested with no API."""

    def __init__(self, task_name, has_skills):
        self.task, self.has_skills, self.n = task_name, has_skills, 0

    def chat(self, messages):
        self.n += 1
        if self.n == 1:
            return "Look around.\n```bash\nls -R /work | head -40\n```", {"prompt": 0, "completion": 0, "cost": 0.0}, _SMETA
        if self.n == 2 and self.has_skills:
            rel = json.loads((TASKS[self.task] / "eval.json").read_text(encoding="utf-8")).get("relevant_skills") or []
            sk = rel[0] if rel else "bio-causal-genomics-fine-mapping"
            return f"Read the relevant skill.\n```bash\nhead -20 /work/skills/{sk}/SKILL.md\n```", {"prompt": 0, "completion": 0, "cost": 0.0}, _SMETA
        sys.path.insert(0, str(TASKS[self.task]))
        import importlib
        mod = importlib.import_module("solvers")
        fn = getattr(mod, "integrated", None) or getattr(mod, "correct", None) or next(
            getattr(mod, n) for n in sorted(dir(mod)) if n.startswith("correct"))
        ans = fn()
        ans = ans[0] if isinstance(ans, tuple) else ans
        if not isinstance(ans, dict) and hasattr(mod, "ans"):
            ans = mod.ans(ans)
        sys.path.pop(0)
        sys.modules.pop("solvers", None)
        return f"Done.\n<EVAL_ANSWER>\n{json.dumps(ans)}\n</EVAL_ANSWER>", {"prompt": 0, "completion": 0, "cost": 0.0}, _SMETA


# ---------------------------------------------------------------- sandbox
def dk(*args, timeout=None, input=None):
    return subprocess.run(["docker", *args], capture_output=True, text=True, timeout=timeout,
                          encoding="utf-8", errors="replace", input=input)


def clip(s):
    if len(s) <= OBS_LIMIT:
        return s
    h = OBS_LIMIT // 2
    return s[:h] + f"\n...[{len(s) - OBS_LIMIT} chars elided]...\n" + s[-h:]


BASH_RE = re.compile(r"```(?:bash|sh)\s*\n(.*?)```", re.S)
# Some models answer in their native tool-call markup instead of a fenced block; take its command the same way.
INVOKE_RE = re.compile(r'invoke name="bash">.*?parameter name="command"[^>]*>(.*?)</[^>]*parameter>', re.S)
ANS_RE = re.compile(r"<EVAL_ANSWER>\s*(.*?)\s*</EVAL_ANSWER>", re.S)
SKILL_REF = re.compile(r"skills/([A-Za-z0-9_.-]+)")


POSTMORTEM = ("The task is over; no further commands will run. In at most 5 short lines, document anything that "
              "went wrong or that you were unsure about: errors or tool problems you hit, ambiguities in the task or "
              "the data, methods you considered but did not use and why, and how confident you are in your answer. "
              "Plain text only.")


def diagnose(events, container_start_s, total_s, err_attempts):
    """Summarise per-step events into the numbers that explain hangs and flakiness."""
    calls = [e for e in events if "model_s" in e]
    lat = sorted(e["model_s"] for e in calls)

    def pct(p):
        return lat[min(len(lat) - 1, int(p * len(lat)))] if lat else None
    retry_kinds = {}
    for e in events:
        for a in e.get("attempts", []):
            retry_kinds[a["kind"]] = retry_kinds.get(a["kind"], 0) + 1
    model_total = sum(e["model_s"] for e in calls)
    exec_total = sum(e.get("exec_s", 0) for e in events)
    return {
        "container_start_s": container_start_s, "n_calls": len(calls),
        "model_s_total": round(model_total, 1), "exec_s_total": round(exec_total, 1),
        "other_s": round(max(total_s - model_total - exec_total - container_start_s, 0), 1),
        "latency_p50": pct(0.5), "latency_p95": pct(0.95), "latency_max": lat[-1] if lat else None,
        "slow_calls_over_60s": [e["step"] for e in calls if e["model_s"] > 60],
        "stalled_calls_gap_over_60s": [e["step"] for e in calls if (e.get("max_gap_s") or 0) > 60],
        "retry_kinds": retry_kinds, "n_retries": sum(retry_kinds.values()),
        "failed_call": err_attempts, "cmd_timeouts": [e["step"] for e in events if e.get("cmd_timeout")],
        "slowest_exec_s": max([e.get("exec_s", 0) for e in events] or [0]),
        "empty_replies": [e["step"] for e in calls if e.get("empty")],
        "truncated_by_length": [e["step"] for e in calls if e.get("finish_reason") == "length"],
        "providers": sorted({e["provider"] for e in calls if e.get("provider")}),
        "reasoning_tokens": sum(e.get("reasoning_tokens") or 0 for e in calls),
        "ttft_max_s": max([e.get("ttft_s") or 0 for e in calls] or [0]),
        "max_silent_gap_s": max([e.get("max_gap_s") or 0 for e in calls] or [0]),
        "completion_tokens_max_step": max([e.get("completion_tokens", 0) for e in calls] or [0]),
    }


def failure_class(passed, reason, graders):
    """Mechanical first-pass label; the trap analysis and the agent's own note refine it."""
    if passed:
        return None
    if reason != "answered":
        return reason.split(":")[0]
    failed = [g["type"] for g in graders if not g["passed"]]
    return "wrong_answer:" + "+".join(failed) if failed else "wrong_answer"


def run_one(task, arm, model_name, rep, out_root, key):
    task_dir = TASKS[task]
    run_id = f"{task}__{arm}__{model_name.replace('/', '_').replace(':', '-')}__r{rep}__{uuid.uuid4().hex[:6]}"
    run_dir = out_root / run_id
    ws = run_dir / "ws"
    ws.mkdir(parents=True)
    shutil.copytree(task_dir / "data", ws / "data")
    if arm != "none":
        for sid in bundle_ids(task, arm):
            shutil.copytree(HERE / "arms" / arm / "skills" / sid, ws / "skills" / sid)
    system, task_text, skipped = build_prompt(task_dir, arm)
    if model_name == "scripted":
        model = Scripted(task, arm != "none")
    elif model_name.startswith("cli:"):
        model = ClaudeCLI(model_name[4:], Path(tempfile.gettempdir()) / "bioeval-cli" / run_id)
    else:
        model = OpenRouter(model_name, key)

    def status(phase, step=0):
        (run_dir / "status.json").write_text(json.dumps({"phase": phase, "step": step, "t": time.time()}), encoding="utf-8")

    def event(ev):
        with open(run_dir / "events.jsonl", "a", encoding="utf-8") as f:
            f.write(json.dumps(ev) + "\n")

    events, err_attempts, self_report = [], None, None
    status("container_start")
    t_c = time.time()
    started = dk("run", "-d", "--rm", "--network", "none", "--cpus", "2", "--memory", "8g",
                 "-w", "/work", IMAGE, "sleep", "86400",
                 timeout=120)   # a wedged Docker leaves the container in Created forever; fail loudly instead
    cid = started.stdout.strip()
    if not cid:
        raise RuntimeError(f"could not start container: {started.stderr.strip()[:300]}")
    # Copy in rather than bind-mount: Docker Desktop asks the user to approve every new host directory.
    copied = dk("cp", f"{ws}{os.sep}.", f"{cid}:/work", timeout=300)
    if copied.returncode != 0:
        dk("rm", "-f", cid)
        raise RuntimeError(f"could not copy the workspace into the container: {copied.stderr.strip()[:300]}")
    container_start_s = round(time.time() - t_c, 2)
    msgs = [{"role": "system", "content": system}, {"role": "user", "content": task_text}]
    traj = [{"role": "system", "content": system}, {"role": "user", "content": task_text}]
    touched, first_touch, cost, ptok, ctok = {}, None, 0.0, 0, 0
    answer, reason, bad_fmt, t0, step = None, "step_limit", 0, time.time(), 0
    try:
        for step in range(1, MAX_STEPS + 1):
            if time.time() - t0 > WALL_LIMIT_S:
                reason = "wall_limit"; break
            if cost > MAX_COST_USD:
                reason = "cost_limit"; break
            status("model_call", step)
            t_call = time.time()
            ev = {"step": step, "t": round(t_call - t0, 1)}
            try:
                text, u, meta = model.chat(msgs)
            except ModelCallError as e:
                ev.update({"model_call_failed": True, "attempts": e.attempts, "latency_s": round(time.time() - t_call, 1)})
                events.append(ev)
                event(ev)
                raise
            ev.update({"model_s": meta["latency_s"], "attempts": meta["attempts"], "provider": meta["provider"],
                       "finish_reason": meta["finish_reason"], "reasoning_tokens": meta["reasoning_tokens"],
                       "completion_tokens": u["completion"], "prompt_tokens": u["prompt"], "empty": meta["empty"],
                       "tokens_per_s": meta["tokens_per_s"], "gen_id": meta["gen_id"],
                       "ttft_s": meta.get("ttft_s"), "max_gap_s": meta.get("max_gap_s")})
            cost += u["cost"]; ptok += u["prompt"]; ctok += u["completion"]
            msgs.append({"role": "assistant", "content": text})
            traj.append({"role": "assistant", "content": text, "step": step, "usage": u})
            # One action per turn, as in mini-swe-agent: if the reply contains a bash block, run the
            # first one and ignore everything after it (models sometimes hallucinate whole
            # transcripts, tool output and a made-up answer, in a single reply). Only a reply with
            # NO bash block can finish, and then the LAST <EVAL_ANSWER> block counts.
            b = BASH_RE.search(text) or INVOKE_RE.search(text)
            if not b:
                answers = ANS_RE.findall(text)
                if answers:
                    try:
                        answer = json.loads(answers[-1])
                        reason = "answered"
                    except Exception:
                        answer, reason = None, "answer_not_json"
                    events.append(ev)
                    event(ev)
                    break
            if not b:
                bad_fmt += 1
                if ev.get("finish_reason") == "length" or ev.get("empty"):
                    obs = ("Your last reply was cut off or empty because it was too long. Think much less: reply with a short "
                           "note and exactly one ```bash block, or give the final <EVAL_ANSWER>.")
                else:
                    obs = "Format error: reply with exactly one ```bash block, or give the final <EVAL_ANSWER>."
                if bad_fmt >= 4:
                    reason = "format_errors"
                    events.append(ev)
                    event(ev)
                    break
            else:
                cmd = b.group(1)
                status("exec", step)
                t_ex = time.time()
                for sid in SKILL_REF.findall(cmd):
                    touched[sid] = touched.get(sid, 0) + 1
                    first_touch = first_touch or step
                try:
                    p = dk("exec", cid, "timeout", str(CMD_TIMEOUT_S), "bash", "-c", cmd, timeout=CMD_TIMEOUT_S + 30)
                    obs = clip((p.stdout or "") + (p.stderr or "")) or "(no output)"
                    if p.returncode == 124:
                        obs += f"\n[command timed out after {CMD_TIMEOUT_S}s]"
                except subprocess.TimeoutExpired:
                    obs = "[command timed out]"
                ev["exec_s"] = round(time.time() - t_ex, 1)
                ev["cmd_timeout"] = "timed out" in obs
                ev["obs_chars"] = len(obs)
            events.append(ev)
            event(ev)
            msgs.append({"role": "user", "content": obs})
            traj.append({"role": "user", "content": obs, "step": step})
            (run_dir / "trajectory.json").write_text(json.dumps(traj, indent=1), encoding="utf-8")
    except Exception as e:  # API / infra failure: record, do not score as a model failure silently
        reason = f"error: {type(e).__name__}: {str(e)[:200]}"
        err_attempts = getattr(e, "attempts", None)
    finally:
        status("finished", step)
        dk("rm", "-f", cid)

    if not reason.startswith("error") and model_name != "scripted":
        # Ask for a brief note on failure modes WITHOUT saying whether the answer was right, so nothing leaks.
        try:
            status("self_report", step)
            note, u2, _ = model.chat(msgs + [{"role": "user", "content": POSTMORTEM}])
            cost += u2["cost"]
            self_report = note.strip()[:1500]
        except Exception as e:
            self_report = f"(self-report failed: {type(e).__name__})"

    graded = grade(task_dir / "eval.json", answer) if isinstance(answer, dict) else {"passed": False, "graders": []}
    result = {
        "run_id": run_id, "task": task, "arm": arm, "model": model_name, "rep": rep,
        "passed": graded["passed"], "graders": [{"type": g["type"], "passed": g["passed"]} for g in graded["graders"]],
        "answer": answer, "finish": reason, "steps": step, "seconds": round(time.time() - t0, 1),
        "cost_usd": round(cost, 5), "prompt_tokens": ptok, "completion_tokens": ctok,
        "skills_touched": touched, "first_skill_step": first_touch,
        "catalog_skipped": skipped, "infra_error": reason.startswith("error"),
        "bundle": BUNDLE, "desc": DESC, "skills_delivered": [] if arm == "none" else bundle_ids(task, arm),
        "config": {"max_tokens": MAX_TOKENS, "providers": PROVIDERS, "max_steps": MAX_STEPS, "wall_limit_s": WALL_LIMIT_S},
        "failure_class": failure_class(graded["passed"], reason, graded["graders"]),
        "self_report": self_report,
        "diagnostics": diagnose(events, container_start_s, time.time() - t0, err_attempts),
    }
    (run_dir / "result.json").write_text(json.dumps(result, indent=1), encoding="utf-8")
    (run_dir / "trajectory.json").write_text(json.dumps(traj, indent=1), encoding="utf-8")
    shutil.rmtree(ws, ignore_errors=True)        # keep logs, drop the working copy of data/skills
    return result


def list_models():
    j = requests.get(f"{OR_URL}/models", timeout=60).json()["data"]
    rows = []
    for m in j:
        p = m.get("pricing") or {}
        try:
            pin, pout = float(p.get("prompt", 0)) * 1e6, float(p.get("completion", 0)) * 1e6
        except Exception:
            continue
        rows.append((pin + pout, m["id"], m.get("context_length"), pin, pout,
                     "tools" if "tools" in (m.get("supported_parameters") or []) else ""))
    for _, i, ctx, pin, pout, tl in sorted(rows):
        print(f"{i:60s} ctx={ctx!s:>8s}  in=${pin:6.2f}/M out=${pout:6.2f}/M {tl}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--list-models", action="store_true")
    ap.add_argument("--task", default="all")
    ap.add_argument("--arm", default="all")
    ap.add_argument("--model", default="scripted")
    ap.add_argument("--reps", type=int, default=1)
    ap.add_argument("--rep-start", type=int, default=1, help="first repeat number (to fill gaps)")
    ap.add_argument("--parallel", type=int, default=2)
    ap.add_argument("--resume", action="store_true", help="skip task/arm/rep cells already completed in --out")
    ap.add_argument("--max-tokens", type=int, default=MAX_TOKENS)
    ap.add_argument("--providers", default=",".join(PROVIDERS), help="comma list in preference order; empty = let OpenRouter route")
    ap.add_argument("--bundle", choices=("all", "domain"), default=BUNDLE)
    ap.add_argument("--desc", choices=("full", "short"), default=DESC)
    ap.add_argument("--out", default=str(HERE / "runs"))
    a = ap.parse_args()
    if a.list_models:
        return list_models()
    globals()["MAX_TOKENS"] = a.max_tokens
    globals()["BUNDLE"], globals()["DESC"] = a.bundle, a.desc
    globals()["PROVIDERS"] = [p for p in a.providers.split(",") if p]
    tasks = list(TASKS) if a.task == "all" else a.task.split(",")
    for t in tasks:
        if t not in TASKS:
            sys.exit(f"unknown task {t!r}; known: {sorted(TASKS)}")
    arms = list(ARMS) if a.arm == "all" else [a.arm]
    key = None
    if a.model != "scripted" and not a.model.startswith("cli:"):
        key = load_key()
        if not key:
            sys.exit("No OpenRouter key: set OPENROUTER_API_KEY or put it in runner/.env")
    out = Path(a.out).resolve()
    out.mkdir(exist_ok=True)
    jobs = [(t, arm, r) for r in range(a.rep_start, a.rep_start + a.reps) for t in tasks for arm in arms]
    if a.resume:                      # skip cells that already have a clean result; park half-finished run dirs
        done = set()
        for f in out.glob("*/result.json"):
            x = json.loads(f.read_text(encoding="utf-8"))
            if x["model"] == a.model and not x.get("infra_error"):
                done.add((x["task"], x["arm"], x["rep"]))
        park = out / "_incomplete"
        for d in out.iterdir():
            if d.is_dir() and d.name != "_incomplete" and not (d / "result.json").exists():
                park.mkdir(exist_ok=True)
                shutil.move(str(d), str(park / d.name))
        jobs = [j for j in jobs if j not in done]
        print(f"resume: {len(done)} cells already done, {len(jobs)} to run", flush=True)
    lock = threading.Lock()
    with cf.ThreadPoolExecutor(a.parallel) as ex:
        futs = {ex.submit(run_one, t, arm, a.model, r, out, key): (t, arm, r) for t, arm, r in jobs}
        for f in cf.as_completed(futs):
            t, arm, r = futs[f]
            try:
                res = f.result()
                with lock:
                    print(f"{t:22s} {arm:9s} r{r} pass={res['passed']!s:5s} finish={res['finish']:14s} steps={res['steps']:2d} "
                          f"${res['cost_usd']:.3f} skills={list(res['skills_touched'])[:3]}", flush=True)
            except Exception as e:
                print(f"{t} {arm} r{r} CRASH {type(e).__name__}: {e}", flush=True)


if __name__ == "__main__":
    main()
