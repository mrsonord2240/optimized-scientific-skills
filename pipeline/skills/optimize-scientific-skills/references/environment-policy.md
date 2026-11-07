# Environment policy

Use the `science` WSL distro first. Preserve its boundary:

- `/mnt/openscience` is the only Windows filesystem mount;
- Windows interop remains disabled;
- use the existing `science` user and micromamba-based environments;
- never widen visibility or enable interop for convenience.

Stage inputs and outputs under `F:\OpenScience`. Record the distro,
environment, relevant versions, launch command, mounts, and checked smoke
result. Prefer an isolated environment when one Skill would downgrade or
materially alter a shared environment.

Use Docker when it is the cleaner supported route; record image version or
digest, mounts, command, resource assumptions, and output checks. Use native
Windows only when the actual tool or integration requires it, and record why
WSL and Docker are unsuitable and how to repeat the native run.

From a Windows agent shell, save the commands as a script under
`F:\OpenScience` with the Write tool and launch it from Git Bash with output
to a log file; a command passed inline or left attached to the terminal hangs:

```bash
MSYS_NO_PATHCONV=1 timeout 250 wsl -d science -- bash -l /mnt/openscience/<path>/run.sh > run.log 2>&1 < /dev/null; cat run.log
```

`wsl: Failed to start the systemd user session for 'sci'` is a harmless
warning; the command still runs. An ecosystem whose environment is a native
Windows one (the staging `TOOLS.md` says so) runs there.

Inside the script call an environment's interpreter by its full path
(`<env>/bin/python`, `<env>/bin/Rscript`); piping `micromamba run` scrambles
the output order.

Install or update dependencies only in a tooling phase. Coordinate shared
mutations with a lock, snapshot relevant versions, run from saved scripts when
practical, kill only run-owned PIDs, and keep heavy artifacts outside the Skill
product and records repository.

Private, authenticated, paid, license-bound, registration-gated, unavailable,
or resource-infeasible requirements remain explicit blockers. Record the exact
error, affected workflow, evidence obtained, user action, and rerun steps.
Never silently replace the advertised workflow with a materially different
public alternative. A heavy optional surface left untooled by policy (defined
in the tooling worker Skill) is not a blocker: the Skill labels it as not
executed and no user action is requested.
