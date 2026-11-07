## Minimal SVDQuartets and BPP Commands

Neither PAUP* (SVDQuartets) nor BPP is installed here -- PAUP* is commercial/registration-gated freeware, not conda/pip-installable, so this is not a version gap to fix but a licensing one. Both blocks below are checked against the PAUP*/SVDQuartets command reference (phylosolutions.com SVDQuartets tutorial) and BPP's own shipped example control files (`bpp/bpp` GitHub, `examples/frogs/A00.bpp.ctl` and `A10.bpp.ctl`, BPP 4.x) -- not executed on this machine.

**SVDQuartets:** map tips to species with a `taxpartition`, then run `svdq` against it.

```
#NEXUS
begin sets;
  taxpartition species = sp1: ind1a ind1b, sp2: ind2a ind2b, sp3: ind3a ind3b, sp4: ind4a ind4b;
end;

begin paup;
  log file=svdq.log;
  svdq taxpartition=species evalq=all bootstrap=standard nreps=100 treeFile=svdq.tre;
end;
```
Run with `paup svdquartets.nex -L svdq.log` (batch) or `exe svdquartets.nex` at the interactive `paup>` prompt. `taxpartition` groups multiple individuals per species tip (RADseq/UCE/SNP matrices); a single-sequence-per-species matrix still works with singleton partitions (`sp1: ind1a`).

**BPP A00** (fixed species tree, parameter estimation only -- `speciesdelimitation = 0`):
```
seed = -1
seqfile = loci.txt
Imapfile = species.imap
jobname = out

speciesdelimitation = 0
speciestree = 0
species&tree = 4  A  B  C  D
                  2  2  2  2
                 ((A, B), (C, D));

phase = 0 0 0 0
usedata = 1
nloci = 20
cleandata = 0

thetaprior = gamma 2 2000
tauprior = gamma 2 1000

finetune = 1
print = 1 0 0 0
burnin = 8000
sampfreq = 2
nsample = 100000
```
`species&tree` lines are: species count and names, max sequences per species at any locus, then the fixed guide tree. `Imapfile` is two columns, one row per sequence label (`<sequence-name> <species-name>`). **BPP A10** (species delimitation on this same guide tree) changes only `speciesdelimitation = 1 0 2` (rjMCMC algorithm 0, finetune 2) with `speciestree = 0` kept -- everything else above is unchanged (checked against `A10.bpp.ctl`).
