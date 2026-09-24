# Similarity Searching - Usage Guide

## Overview
Find structurally similar molecules using fingerprint-based Tanimoto similarity. Includes library searching, clustering, and maximum common substructure analysis.

## Prerequisites
```bash
pip install rdkit
pip install numpy pandas
```

## Quick Start
Tell your AI agent what you want to do:
- "Find compounds similar to my lead molecule"
- "Search my library for molecules with Tanimoto > 0.7 to this query"
- "Cluster my compound library by structural similarity"
- "Find the maximum common substructure among these molecules"

## Example Prompts

### Similarity Search
> "Find all compounds in my library with Tanimoto similarity > 0.7 to aspirin."

> "Search for analogs of this SMILES using ECFP4 fingerprints."

### Clustering
> "Cluster my compound library at 70% similarity using Butina clustering."

> "Group similar molecules together with a 0.3 distance cutoff."

### Scaffold Analysis
> "Find the maximum common substructure for my hit series."

> "What structural features are shared by these active compounds?"

## What the Agent Will Do
1. Generate fingerprints for query and library
2. Calculate pairwise Tanimoto similarities
3. Filter/rank by similarity threshold
4. Perform clustering if requested
5. Find MCS for shared scaffold analysis

## Tips
- For ECFP4, >0.85 and >0.7 are repository starting heuristics for close and related analogs; calibrate them on the project's known series and fingerprint settings
- Butina distance cutoff = 1 - centroid-neighbor similarity threshold (cutoff 0.3 assigns neighbors at >=0.7 similarity to a selected centroid; it does not guarantee >=0.7 for every member pair)
- BulkTanimotoSimilarity is faster for large libraries
- ECFP4 (radius=2) is the most common fingerprint for similarity
- Tversky is a feature-overlap ranker, not proof of substructure containment; confirm retained candidates with a SMARTS match
- Set a finite MCS timeout and inspect `result.canceled`: a tiny MCS with `canceled=False` needs pre-clustering, not a longer timeout
- Calibrate each fingerprint's cutoff to a retained-pair percentile on the actual library rather than reusing a raw score
- Strip salts before comparison and use chirality-aware Morgan fingerprints when stereoisomers matter; use standardized canonical SMILES or InChIKey for identity

## Related Skills
- molecular-descriptors - Generate fingerprints for similarity
- substructure-search - Pattern-based searching
- molecular-io - Load molecules for searching
