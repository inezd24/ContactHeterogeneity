# GraphComparison

Code accompanying the manuscript '*Contact heterogeneity in network reconstruction: evaluating synthetic graph generators for rabies transmission modeling*', which compares five synthetic network generators (SBM, DCSBM, ERM, NCRG, SENCA from Laager et al., 2018) against empirical dog-contact networks from four locations, and evaluating how well each generator reproduces rabies (SEIR) and generic (SIS) outbreak dynamics on those networks. This manuscript has been submitted to *PLOS Computational Biology*. DOI:

**Code author:** Inez Derkx (inez.derkx@swisstph.ch)

**Date:** September 2026

## Status of this code

This repository contains the analysis scripts exactly as they were run on my HPC cluster to produce the results reported in the manuscript. The R code is presented here as it was run, only file locations have been changed. In comparison to my local versions, the following comments should be noted:

1. All analyses for the original manuscript were run on University of Basel's sciCORE HPC using SLURM. Both the .R and .sh files are included in this repository.
2. The scripts contain **hardcoded, machine-specific paths** to my HPC account (see [Adapting paths to your own system](#adapting-paths-to-your-own-system)). These must be edited by anyone who wants to actually *execute* the pipeline elsewhere. Editing these paths does not change any computation; it only tells R and SLURM where to read/write files. I intentionally left this for the user to adapt to keep the code as identical to its original state.
3. All .R scripts contain there own notes in the header. Please make sure to read these prior to execution. 

## Repository structure

```
GraphComparison/
├── README.md
├── LICENSE
├── .gitignore
├── renv.lock                                 # Information about the full project environment
├── Derkx2026_publication.Rproj               # This project holds all scripts and the environment 
├── data/
│   └── README.md                             # What input data the pipeline expects (not included; see Data availability)
└── src/
    ├── run_pipeline.sh                       # This runs all scripts (01-06) in a single pipeline
    ├── 00_constructNetworksFunction.R   
    ├── 01_NetworkGeneratorComparison.R/.sh
    ├── 01a_Run_NCRG.R/.sh
    ├── 02_collectGraphs.R
    ├── 03a_SEIR_Randomized.R/.sh
    ├── 03b_SEIR_Factorial.R/.sh
    ├── 04a_SIS_Randomized.R/.sh
    ├── 04b_SIS_Factorial.R/.sh
    ├── 05_Graphs_Other_Locations.R/.sh
    ├── 06_collectGraphsOtherLocations.R
    ├── manuscript_plots.R
    └── manuscript_SI_plots.R
```

Each `.R`/`.sh` pair is kept **in the same folder** on purpose: the SLURM scripts call `Rscript <script_name>.R` using a bare filename, which only resolves correctly if the shell script is submitted (`sbatch`) from inside the same directory the `.R` file lives in. This is why they are not separated into an 'src/' and 'slurm/' folder. 

## Pipeline / execution order

Step numbers match the file name prefixes as provided. Scripts without a `.sh` file were run as single (non-array) jobs, usually in a medium-sized session on a single node with 16 cores. You are welcome to adapt memory use, which is especially important when loading the larger data files for plotting. 

| Step | Script(s) | What it does | SLURM array |
|---|---|---|---|
| — | `00_constructNetworksFunction.R` (+ `.RData`) | Reference/master copy of the shared helper functions (`empirical_net`, `construct_network`, `generate_sbm_graph`, `generate_dcsbm_graph`, `grid_optimization`, `construct_five_networks`, `compare_ks_distances`, etc.). Not run directly — these functions are copy-pasted into the scripts below so each SLURM array task is self-contained without a `source()` call. You are welcome to adapt this, hence they are included as such. | — |
| 0 | `run_pipeline.sh` | This .sh script runs the whole pipeline (01-06, but not 00 as this is sourced separately in the R scripts) in a single go | — |
| 1 | `01_NetworkGeneratorComparison.R/.sh` | For the primary ("Chad") location: generates and compares the five synthetic network types against the empirical network via grid-search parameter optimization and K-S distance. | 1–5000 |
| 1a | `01a_Run_NCRG.R/.sh` | Reruns/patches the NCRG (Newman) generator for each seed produced in Step 1. | 1–5000 |
| 2 | `02_collectGraphs.R` | Aggregates the per-seed outputs of Steps 1 and 1a into the master ensemble of vetted graphs and summary tables used downstream. | — (single run) |
| 3a | `03a_SEIR_Randomized.R/.sh` | Rabies SEIR outbreak simulations under randomly sampled parameters, submitted separately per location (Romana/Sabaneta/Habi/Hepang) via a loop of `sbatch` calls inside the `.sh` file. | 0–9999 per location |
| 3b | `03b_SEIR_Factorial.R/.sh` | SEIR simulations under a full factorial parameter design. | 0–9999 |
| 4a | `04a_SIS_Randomized.R/.sh` | SIS outbreak simulations, randomly sampled parameters. | 0–999 |
| 4b | `04b_SIS_Factorial.R/.sh` | SIS simulations, factorial design. | 0–9999 |
| 6 | `05_Graphs_Other_Locations.R/.sh` | Repeats the Step 1 network-generation comparison for the three additional empirical locations. | 1–5000 |
| 7 | `06_collectGraphsOtherLocations.R` | Aggregates Step 5's outputs across the additional locations (companion to Step 2). | — (single run) |
| — | `manuscript_plots.R` | Generates the main-text manuscript figures from the pipeline outputs above. | — |
| — | `manuscript_SI_plots.R` | Generates the Supplementary Information figures. | — |

## Software environment

All scripts were run under:

- **R 4.2.1** (`R/4.2.1-foss-2022a` environment module)
- Cluster: scicore (sciCORE HPC, University of Basel), SLURM workload manager
- You can find the corresponding renv.lock file in this directory with information about all packages used in the .Rproject.
- If you clone this repository, you can simply run `renv::restore()` or `renv::restore(lockfile = "path/to/renv.lock")` using the `renv` R package to recreate the environment.

## Adapting paths to your own system

Every script sets an absolute working directory before doing anything else, and several also hardcode an HPC library path or input-data path. To run these scripts on a different machine, update the following (line numbers refer to the files as provided):

| File | Line(s) | What's hardcoded |
|---|---|---|
| `01_NetworkGeneratorComparison.R` | 20, 33 | `.libPaths(...)`; `LOCAL_ROOT_DIR` |
| `01a_Run_NCRG.R` | 20, 33 | `.libPaths(...)`; `LOCAL_ROOT_DIR` |
| `02_collectGraphs.R` | 30 | `LOCAL_ROOT_DIR` |
| `03a_SEIR_Randomized.R` | 26, 268–272 | `LOCAL_ROOT_DIR`; absolute `GRAPH_FILE`/`EMP_FILE` paths |
| `03b_SEIR_Factorial.R` | 25 | `LOCAL_ROOT_DIR` (`GRAPH_FILE` is passed in from the `.sh` script instead) |
| `04a_SIS_Randomized.R` | 26, 257 | `LOCAL_ROOT_DIR`; a hardcoded `vetted_graphs_...rds` fallback path |
| `04b_SIS_Factorial.R` | 25 | `LOCAL_ROOT_DIR` (`GRAPH_FILE` passed in from the `.sh` script) |
| `05_Graphs_Other_Locations.R` | 20, 35 | `.libPaths(...)`; `LOCAL_ROOT_DIR` |
| `06_collectGraphsOtherLocations.R` | 3 | `file_dir` |
| `manuscript_plots.R` | 37, 46 | `LOCAL_ROOT_DIR`; `files_directory` |
| `manuscript_SI_plots.R` | 38 | `LOCAL_ROOT_DIR` |

In addition, several scripts (`01_NetworkGeneratorComparison.R`, `01a_Run_NCRG.R`, `manuscript_plots.R`, `manuscript_SI_plots.R`) read a base empirical dataset from `~/BaseData/Chad/ServerData.rds` or the equivalent for the other three locations — see [`data/README.md`](data/README.md). The `.sh` files also hardcode `#SBATCH --output`/`--error` log paths under `/scicore/.../GraphComparison/SLURM/`, which you'll want to point at your own cluster account.

## Data availability

The raw/input data (dog-contact survey and GPS-derived edgelists underlying `ServerData.rds` for each of the four locations) are **not included in this repository**. Add a pointer here to wherever that data is deposited (e.g. Dryad/Zenodo) once available, to match the manuscript's Data Availability Statement. See [`data/README.md`](data/README.md) for the specific files the pipeline expects.

## Reproducibility notes

- Random seeds are handled explicitly throughout (`MASTER_SEED`, per-task `base_seed <- 1000 + task_id * 10000`, etc.), so individual array tasks are independently reproducible given the same seed and R/package versions.

## License

Released under the MIT License — see [`LICENSE`](LICENSE).

## Citation

If you use this code, please cite: [to be adapted once versionized]

> Derkx I. et al. (in review). *Contact heterogeneity in network reconstruction: evaluating synthetic graph generators for rabies transmission modeling*. PLOS Computational Biology*. [DOI once available]

## Contact

Inez Derkx — inez.derkx@swisstph.ch
If you encounter any issues with the code, please feel free to write me an email or open an issue on this repository. 
If you're interested in this work and would like to collaborate on a similar project, feel free to reach out. 
