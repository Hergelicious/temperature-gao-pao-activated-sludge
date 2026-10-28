# Temperature and the GAO:PAO balance in full-scale activated sludge

Code and data to reproduce all statistical results in the manuscript
"Temperature shifts the balance from polyphosphate- to glycogen-accumulating organisms in full-scale activated sludge worldwide".

## Contents

```
analysis.py          discovery cohort: all results from the MiDAS global genome catalog (run first)
validation_midas4.py validation cohort: all results from the MiDAS 4 amplicon survey (run second)
requirements.txt     Python dependencies
results/             output tables written by analysis.py
```

## Input data (`data/raw/`)

The input files are not redistributed here. Download them from the original publications and place them in `data/raw/` with their original file names (details in `data/raw/README.md`).

| File | Source | Used for |
| --- | --- | --- |
| `media-2.xlsx` | Supplementary Tables of Liu et al. (2026), The MiDAS global genome catalog, bioRxiv, https://doi.org/10.64898/2026.07.10.737647 | Plant metadata (S1), sequencing yield (S2), MAG taxonomy and abundance (S4, S5), trait conservation (S13) |
| `41467_2021_22203_MOESM5_ESM.xlsx` | Supplementary Data 3 of Singleton et al. (2021), Nat. Commun. 12, 2009, https://doi.org/10.1038/s41467-021-22203-2 | Taxonomy of Danish HQ MAGs |
| `41467_2021_22203_MOESM15_ESM.xlsx` | Supplementary Data 13 of Singleton et al. (2021) | MAG abundances in 23 Danish plants, 2016–2018 |
| `midas4/` (six files) | `Dueholm2021a_data.zip` from Dueholm et al. (2022), MiDAS 4, Nat. Commun. 13, 1908; figshare https://doi.org/10.6084/m9.figshare.16566408 | ASV tables (V1–V3 and V4), MiDAS 4.8 classifications, sequencing and plant metadata for the validation cohort |

## Running

```
pip install -r requirements.txt
# download the input files into data/raw/ first
python analysis.py
python validation_midas4.py
```

Runtime is a few minutes for each script (Python 3.10 or later; tested with Python 3.12, pandas 3.0, statsmodels 0.15).

## Outputs (`results/`)

| File | Content |
| --- | --- |
| `primary_results.csv` | Confirmatory models H1a, H1b, H2, H3: temperature effect per 10 °C with 95% and 90% CI, p, equivalence margin Δ and TOST p (H1), Holm-adjusted p |
| `h4_shift_vs_yearly_variation.csv` | H4: predicted shift across the interquartile temperature range and its ratio to Danish year-to-year variation |
| `h5_ebpr_subset.csv` | H5: EBPR and non-EBPR subsets, temperature × EBPR interaction, post hoc power and minimum detectable effect |
| `sensitivity_analyses.csv` | Sensitivity analyses S1–S9 (leave-one-out rows show the weakest estimate and the range) |
| `sensitivity_leave_one_out_detail.csv` | Every leave-one-continent-out (S5) and leave-one-genus-out (S6) fit |
| `trait_scan.csv` | Exploratory scan of 23 community-weighted traits with Benjamini–Hochberg q-values |
| `genus_slopes.csv` | Genus-level temperature effects (CLR) for PAO, GAO, storage-niche and nitrite-oxidising genera |
| `danish_baseline_samples.csv` | Outcome values for each Danish sample (23 plants, Q3 2016–2018) |
| `danish_baseline_delta.csv` | Δ, the median within-plant standard deviation across years, per outcome |
| `plant_dataset.csv` | Plant-level analysis dataset: metadata, covariates, outcomes and community-weighted traits |
| `validation_results.csv` | Validation analyses V1–V5 in the MiDAS 4 survey: primary replication, EBPR subsets and interaction, sensitivity analyses, Nitrotoga:Nitrospira, catalog plants re-measured by amplicon, power |
| `validation_genus_slopes.csv` | Genus-level temperature effects (CLR) for the guild genera in the independent validation plants |
| `validation_method_concordance.csv` | Correlation between amplicon- and genome-based GAO:PAO balances in the plants present in both datasets |
| `validation_plant_dataset.csv` | Validation plant-level dataset: metadata, covariates, balances, genus CLR values and cohort (catalog or independent) |

## Variable definitions

- `T10`: reported plant temperature / 10 (°C). Effects are per 10 °C.
- `EBPR`, `Nrem`, `AS`: indicators for phosphorus removal, nitrogen removal and activated sludge plant type; `lind`: log(1 + industrial share of influent, %).
- `H1a`: log-odds of the storage-niche share of the trait-annotated core pool.
- `H1b`: log total nitrite-oxidiser relative abundance.
- `H2`: GAO:PAO balance, mean CLR of GAO genera minus mean CLR of PAO genera.
- `H3`: log ratio of *Nitrotoga* to *Nitrospira_A* + *Nitrospira_D* (validation: *Nitrotoga* to *Nitrospira*).
- Validation guilds (MiDAS 4.8 names): PAO = Ca_Accumulibacter, Tetrasphaera, Dechloromonas; GAO = Ca_Competibacter, Ca_Contendobacter, Defluviicoccus, Propionivibrio. Libraries with fewer than 10,000 reads are discarded and replicate libraries of a plant are averaged.

All models are ordinary least squares with HC3 standard errors and the covariates EBPR, nitrogen removal, plant type and industrial load, unless stated otherwise.