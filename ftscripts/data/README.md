## Silva prevalence: human-gut v2.0.2

`silva_prevalence_human_gut_v2.0.2.tsv` has one row for each of the
4,744 MGnify species representatives. Of these, 3,027 map to Table S5 of
[Silva et al.](https://doi.org/10.1016/j.chom.2026.01.013).
The mapping uses Table S5 `Genome` -> MGnify `Genome_accession` ->
`Species_rep`. All 3,027 mapped entries resolve to different representatives.

Prevalence is the percentage of healthy samples in which the species was
detected in each continent. The global value is the arithmetic mean of the
six continent percentages, as reported by Silva; it is not abundance and is
not the fraction of all pooled samples worldwide.

The `global_mean_ge_10_percent` through `global_mean_ge_100_percent`
columns flag nested subsets with global mean **>=** the cutoff.
Blank values indicate absent data in Table S5, not zero prevalence or a
failed cutoff. `prevalence_status` identifies the 1,717 such representatives.
Per-continent prevalence columns preserve the values reported in the paper.
The companion JSON records source fingerprints and threshold denominators.

## Silva core: human-gut v2.0.2

`silva_core_human_gut_v2.0.2.tsv` contains the 31 Table S5 entries explicitly
marked core, mapped to the current representative. Names from the paper
are retained separately from MGnify taxonomy. In this catalogue, this panel
matches the global mean >=70% subset of the prevalence table.
