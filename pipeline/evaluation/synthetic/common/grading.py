"""Grade an answer against an eval.json using LatchBio's own graders.

The `latch_eval_tools` package __init__ pulls in the Latch SDK, which we do not
have and cannot log into. The graders themselves are pure Python, so we register
a bare parent package and import only `latch_eval_tools.graders`.

Requires: ../../latch-eval-tools (cloned from github.com/latchbio/latch-eval-tools)
"""
import json
import re
import sys
import types
from pathlib import Path

_SRC = Path(__file__).resolve().parents[2] / "latch-eval-tools" / "src" / "latch_eval_tools"


def _load_graders():
    if "latch_eval_tools" not in sys.modules:
        pkg = types.ModuleType("latch_eval_tools")
        pkg.__path__ = [str(_SRC)]
        sys.modules["latch_eval_tools"] = pkg
    from latch_eval_tools.graders import get_grader
    return get_grader


def extract_answer(text: str) -> dict:
    """Pull the JSON between <EVAL_ANSWER> tags, as the benchmark prompt demands."""
    m = re.search(r"<EVAL_ANSWER>\s*(.*?)\s*</EVAL_ANSWER>", text, re.S)
    if not m:
        raise ValueError("no <EVAL_ANSWER> block")
    return json.loads(m.group(1))


def grade(eval_path, answer: dict) -> dict:
    """Return {'passed': bool, 'graders': [{'type', 'passed', 'reasoning'}]}.

    Passed means every grader passed. Per-grader results are kept so partial
    credit can be reported separately.
    """
    get_grader = _load_graders()
    spec = json.loads(Path(eval_path).read_text(encoding="utf-8"))
    results = []
    for g in spec["graders"]:
        r = get_grader(g["type"]).evaluate_answer(answer, g["config"])
        results.append({"type": g["type"], "passed": bool(r.passed), "reasoning": r.reasoning})
    return {"passed": all(x["passed"] for x in results), "graders": results}
