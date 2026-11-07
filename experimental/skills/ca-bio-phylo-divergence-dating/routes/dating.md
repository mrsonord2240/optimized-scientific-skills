# Run configured prior and posterior analyses
~~~bash
mcmctree mcmctree_prior.ctl
~~~
Supply a reviewed MCMCTree control file with usedata=0 in its own run directory. The source pipeline also uses usedata=3 to generate out.BV, then usedata=2 in.BV for the posterior in a separate directory. scripts/mcmctree_setup.py accepts an alignment and calibrated tree and launches that entire pipeline; it creates a temporary run directory and uses fixed default priors, so inspect and adapt settings before use. A single input is not a setup-only CLI. For BEAST use reviewed prior-only and data XMLs in separate output locations. See [calibrations and workflows](../references/methods.md).
Done when effective prior, posterior, calibration assumptions, units and independent convergence diagnostics are reviewed.
