# Lipid Tracking in Membrane Simulations

A Python project developed for a computational biology examination
to analyse lipid motion in a membrane molecular-dynamics trajectory
provided as a multi-model PDB file.

## Overview

The script uses MDAnalysis to select membrane components by residue
name and track their centres of mass over time.

Implemented routines include:

- Tracking of CHL, DAPC, and DSPC residues.
- Interactive XY trajectory animations for selected lipid types.
- Selection of upper, lower, or both membrane leaflets.
- An XZ trajectory animation for CHL.
- Calculation and CSV export of lateral squared displacements.
- Plotting and linear fitting of displacement data.

## Requirements

- Python
- MDAnalysis
- NumPy
- Matplotlib

```bash
pip install MDAnalysis numpy matplotlib
```

## Usage

```bash
python LipidTracking.py data
python LipidTracking.py plots
python LipidTracking.py xy
python LipidTracking.py xz
```

The `data`, `xy`, and `xz` routines prompt for the trajectory path.
The XY routine additionally requests the residue name and leaflet.

The input must provide multiple trajectory frames, the residue names
`CHL`, `DAPC`, and `DSPC`, and suitable atomic masses for
centre-of-mass calculations.

## Output

- `CHL.csv`, `DAPC.csv`, `DSPC.csv` — displacement-analysis arrays.
- `MSD.png` — displacement plots and linear fits.
- Interactive Matplotlib windows displaying lipid trajectories.

## Project status

This repository documents an examination project, not a validated
general-purpose analysis package.

The original code assumes 201 trajectory frames and a 100 ps frame
spacing in its plotting routine. The displacement calculation currently
evaluates only a single-frame lag; a complete lag-dependent MSD analysis
requires extending this loop.

Before quantitative use, the plotting code, CSV-to-residue mapping,
diffusion-coefficient calculation, and periodic-boundary handling
require correction or verification.

Leaflet membership is assigned using the mean Z coordinate of the
selected lipid type at the initially loaded frame and is kept fixed
during animation.
