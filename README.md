# From Atoms to Neuronal Spikes: A Multiscale Study of AMPA Receptors

[![DOI](https://img.shields.io/badge/DOI-10.1021%2Facs.jctc.5c01793-blue.svg)](https://pubs.acs.org/doi/10.1021/acs.jctc.5c01793)

**Authors:** Ana Damjanovic, Vincenzo Carnevale, Thorsten Hater, Nauman Sultan, Giulia Rossetti, Sandra Diaz-Pier, Paolo Carloni  
**Affiliations:** Johns Hopkins University | National Institutes of Health | Temple University | Forschungszentrum Jülich GmbH | University Hospital Aachen

## Overview
This repository contains the system preparation files, simulation inputs, analysis scripts, and plotting code (for MD simulations) used in the publication: **"From Atoms to Neuronal Spikes: A Multiscale Simulation Framework"** (*Journal of Chemical Theory and Computation*). 

The computational pipeline outlines the workflow to simulate the wild-type AMPA receptor (AMPAR) and its Q586 mutants (Q586E, Q586G, Q586R). It bridges atomistic properties to macroscopic observables by simulating systems under an applied electric field to measure ion conductance, pore geometry, and electrostatic profiles.

## Repository Structure
```text
AMPAR-Multiscale-MD/
├── 1_system_setup/               # PDB and PSF topology files for WT and mutants
│   ├── Q586E_0
│   ├── Q586E_1
│   ├── Q586E_2
│   ├── Q586G
│   ├── Q586R
│   └── wild_type
├── 2_simulation_input/           # AMBER input files for NPT and NVT production
│   ├── amber_equilibration_NPT
│   └── amber_production
├── 3_analysis_scripts/           # Trajectory analysis and conductance calculations
└── 4_plot_scripts/               # Python routines for manuscript figures
```

## System Setup (`1_system_setup`)
Contains the starting configurations and CHARMM force field topologies required to build the simulation systems. All structures are pre-aligned to the Z-axis to facilitate downstream flux and pore analyses.

* **`wild_type/`**: Contains the base structural files (`.pdb`, `.psf`) for the native GLN 586 receptor.
* **`Q586E_*/`, `Q586G/`, `Q586R/`**: Contains the structural files for the respective selectivity filter mutants, including independent replicates (e.g., `_0`, `_1`, `_2`).

## Simulation Inputs (`2_simulation_input`)
Contains the configuration files to execute the Molecular Dynamics pipeline using the AMBER engine.

* **`amber_equilibration_NPT/`**: Includes `preproduction_NPT.in`, which equilibrates the membrane system using semi-isotropic pressure scaling to stabilize the area-per-lipid.
* **`amber_production/`**: Includes `production_ef.in`, the production run utilizing the NVT ensemble (to avoid barostat artifacts) while applying an external electric field (`efz = -0.11`) to drive ion permeation.

## Analysis Scripts (`3_analysis_scripts`)
The parsing utilities and analytical scripts responsible for extracting biophysical properties from the raw MD trajectories. 

* **`calc_ampar_1d_profiles.sh`**: A shell script utilizing VMD's PMEpot to compute the 3D electrostatic potential and map it to a 1D profile alongside ion density.
* **`calc_ampar_conductance.py`**: Computes macroscopic ion channel current (pA) and conductance (pS) from tracked permeation events and applied voltage.
* **`calc_ampar_ion_flux.py`**: Quantifies the rate, direction, and specific pathway of ion permeation across the selectivity filter and lower gate.
* **`calc_ampar_pore_radius.py`**: A Python wrapper for HOLE2 that computes the physical dimensions and bottleneck radius of the conductive pore over time.
* **`calc_ampar_rmsd_box.py`**: Tracks the Root Mean Square Deviation (RMSD) of the protein backbone and simulation box dimensions to verify system stability.

## Figure Generation (`4_plot_scripts`)
Python plotting routines using `matplotlib` and `seaborn` to generate the final visualizations for the manuscript.

* **`plot_ampar_1d_profiles.py`**: Generates a two-panel publication plot comparing the 1D electrostatic potential and corresponding K+/Cl- ion densities.
* **`plot_ampar_rmsd.py`**: Plots the combined backbone RMSD trajectories across multiple simulated systems to compare structural stability over time.

## Usage
1. Clone the repository to your local machine or cluster.
2. Ensure you have the required dependencies installed (AMBER, VMD, HOLE2, Python 3 with `numpy`, `pandas`, `mdtraj`, `MDAnalysis`, `matplotlib`, `seaborn`).
3. Execute the workflow sequentially from directories 1 through 4. 

*(Note: Raw MD trajectory `.nc` files are not hosted directly on this repository due to size limits. Execution of the scripts in step 3 requires the base trajectory data.)*

*(Note: The AMBER `.in` files in step 2 are configured for 1 ns blocks. For adequate sampling of permeation events, submit these as continuous loops on your HPC cluster).*

## Citation
If you use these scripts, parameters, or methods in your research, please cite our paper:
```bibtex
@article{Damjanovic_2026,
  author = {Damjanovic, Ana and Carnevale, Vincenzo and Hater, Thorsten and Sultan, Nauman and Rossetti, Giulia and Diaz-Pier, Sandra and Carloni, Paolo},
  title = {From Atoms to Neuronal Spikes: A Multiscale Simulation Framework},
  journal = {Journal of Chemical Theory and Computation},
  year = {2026},
  volume = {22},
  number = {2},
  pages = {783-793},
  doi = {10.1021/acs.jctc.5c01793}
}
```

## License
This project is licensed under the MIT License - see the `LICENSE` file for details.