"""Score one physical capture with Scanpy Scrublet without dropping cells."""
import argparse
import scanpy as sc

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("input")
p.add_argument("output")
p.add_argument("--expected-rate", type=float, required=True)
p.add_argument("--seed", type=int, default=0)
args = p.parse_args()
if not 0 < args.expected_rate < 1:
    p.error("expected rate must be between 0 and 1")
adata = sc.read_h5ad(args.input)
sc.pp.scrublet(adata, expected_doublet_rate=args.expected_rate, random_state=args.seed)
adata.write_h5ad(args.output)
