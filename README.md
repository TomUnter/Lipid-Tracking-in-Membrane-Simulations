# Lipid Tracking and Membrane Diffusion Analysis

A Python project developed for a computational biology examination
to analyse lipid motion in a membrane molecular-dynamics trajectory.

The project combines molecular tracking, lateral mean-squared
displacement (MSD) calculations, visualisation, and comparison
with GROMACS-generated MSD data.

## Features

- Selection of CHL, DAPC, and DSPC residues using MDAnalysis.
- Tracking of residue centres of mass.
- Lateral MSD calculation in the membrane XY plane.
- Averaging over multiple time origins for each lag time.
- Subtraction of the all-atom centre of mass to account for drift.
- Per-molecule CSV output and lipid-type averages.
- Linear fitting and comparison with GROMACS MSD output.
- Interactive trajectory animations.

## Files

| File | Purpose |
|------|---------|
| `MDLipid.py` | Full lag-time MSD calculation and CSV export. |
| `Auswertung.py` | CHL MSD plotting, linear fitting, and comparison with GROMACS output. |
| `LipidTracking.py` | Command-line version with analysis and animation routines; its current MSD loop evaluates only one lag. |
| `Complete.py` | Combined development version of the analysis and animation routines. |
| `Animation.py` | Early standalone animation draft. |

## Requirements

Python packages:

```bash
pip install numpy matplotlib MDAnalysis
```

The early `Animation.py` draft also imports pandas.

## Main analysis workflow

Place the multi-model PDB trajectory in the working directory as
`mem.pdb`, then run:

```bash
python MDLipid.py
```

This generates:

- `CHL.csv`
- `DAPC.csv`
- `DSPC.csv`

Each CSV stores lag times along rows and individual molecules
along columns. Coordinates are converted from angstroms to
nanometres before calculating squared displacements.

To compare the CHL results with GROMACS output, also provide
`msd_CHOL.xvg` and run:

```bash
python Auswertung.py
```

Trajectory and reference-output files must be supplied separately.

## Assumptions and project status

This repository contains the original examination scripts and
development versions, not a validated general-purpose package.

The MSD calculation assumes 201 trajectory frames. The plotting
scripts assume a frame spacing of 100 ps.

Before quantitative interpretation, the diffusion-coefficient
conversion, fitting interval, periodic-boundary treatment, and
consistency with the GROMACS reference analysis require verification.

Some development scripts contain known syntax or implementation
issues and need correction before execution.
