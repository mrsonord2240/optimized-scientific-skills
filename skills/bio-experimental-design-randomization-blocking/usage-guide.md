# Randomization and Blocking Usage Guide

## Overview

This guide covers the structural design of biological experiments: deciding the experimental unit, randomizing treatments, replicating the right thing, and removing known nuisance variation by blocking. The central goal is that the analysis model mirrors how the experiment was actually run, so that the resulting inference is valid by construction. The single most common failure it prevents is pseudoreplication: counting cells, wells, or technical aliquots as independent replicates when the true sample size is the number of animals, donors, or cages.

## Prerequisites

```r
# R packages
install.packages(c('designit', 'lme4', 'lmerTest', 'pwr'))
```

Familiarity with the difference between a fixed effect (a question) and a random effect (a level of the randomization/sampling structure) is helpful but not required; the agent will explain the mapping.

## Quick Start

Tell your AI agent what you want to do:
- "What is the experimental unit in my study, and what is my real n?"
- "I have 10,000 cells from 3 mice per group; can I test on cells?"
- "Help me randomize 24 samples to 3 processing days without confounding"
- "Should this be a factorial, split-plot, or nested design?"
- "Write the mixed-model random-effects structure that matches my design"

## Example Prompts

### Experimental Unit and Pseudoreplication

> "I measured a marker in 200 cells from each of 4 control and 4 treated animals. A reviewer says my n is 4, not 1600. How should I analyze this?"

> "Treatment is delivered through the drinking water of co-housed cages. What is my experimental unit?"

### Randomization and Blocking

> "I can only process 8 of my 24 samples per day over 3 days. How do I assign and randomize them so that processing day is not confounded with condition?"

> "Set up a randomized complete block design for 2 conditions across 3 litters and give me the model formula to analyze it."

### Design Structure

> "My incubator can only hold one temperature at a time, but I want to test temperature and genotype. Is this a split-plot, and how does that change the analysis?"

> "I want to test genotype and drug together; design a factorial and tell me how to read the interaction."

## Related Skills

See SKILL.md's Related Skills section.
