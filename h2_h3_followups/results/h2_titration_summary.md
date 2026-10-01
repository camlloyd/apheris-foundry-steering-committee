# H2 titration — result against the pre-registered criterion

Criterion: at some depth, >= 2/5 seeds give mean pLDDT >= 70 AND H12-out with margin >= 1 A.

| depth | structures | folded (pLDDT>=70) | mean pLDDT | H12-out structures | seeds with >=1 H12-out | criterion |
|---|---|---|---|---|---|---|
| 8 | 50 | 0 | 46.5 | 0 | 0 ([]) | **not met** |
| 32 | 50 | 0 | 44.5 | 0 | 0 ([]) | **not met** |
| 128 | 50 | 50 | 89.6 | 0 | 0 ([]) | **not met** |
| 512 | 50 | 50 | 89.3 | 0 | 0 ([]) | **not met** |

## Outcome

**H2 not supported — declared negative applies.** On RORγ LBD, fold stability and H12-state accessibility could not be separated by MSA depth: the state was inaccessible at every depth that preserved the fold. This closes the mechanism the literature suggests should work, and is the stronger of the two negatives.

100 structure(s) failed the pLDDT floor and are recorded as `unscorable` in titration_results.csv — listed, not dropped.
