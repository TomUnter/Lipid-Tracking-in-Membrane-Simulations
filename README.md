# lipid-msd

Lateral diffusion of lipids in a membrane simulation: per-lipid mean squared displacement (MSD),
diffusion coefficients from a linear fit, and animated centre-of-mass tracks.
Developed from my master's thesis project (Python, MDAnalysis, NumPy).

## What it does

1. Reads a multi-frame membrane trajectory (e.g. the PDB written by GROMACS).
2. Computes the centre of mass (COM) of every lipid residue (mass-weighted), taking periodic
   boundaries into account.
3. Removes the drift of the whole membrane.
4. Computes the time-averaged lateral (xy) MSD for every lipid and lag time.
5. Fits MSD = 4·D·t and reports D in 10⁻⁵ cm²/s. Optionally compares to `gmx msd`.
6. Makes an animation of the lipid tracks (xy or xz).

## Installation

```
pip install -r requirements.txt
```

## Input requirements

- Multi-frame file with **box information in every frame** (`CRYST1` in each `MODEL` for PDB files).
- **Orthorhombic box** (right angles).
- Residue names as in the topology (default: `CHL`, `DAPC`, `DSPC`; set with `--lipids`).
- The **time between frames in ps** (`--dt`, required). PDB files do not store it.

If the file also contains water or ions, restrict the membrane selection with
`--selection "resname CHL or resname DAPC or resname DSPC"`. The drift is calculated from this selection.

## Usage

```
# 1. MSD per lipid (reads the trajectory once, can take a while for large files)
python lipid_msd.py msd mem.pdb --dt 100 --lipids CHL DAPC DSPC --out results

# quick test on the first 30 frames
python lipid_msd.py msd mem.pdb --dt 100 --max-frames 30 --out test_run

# 2. plots + diffusion coefficients
python lipid_msd.py plot results --fit 0.1 0.5

# compare with GROMACS (gmx msd ... -lateral z)
python lipid_msd.py plot results --gromacs CHL=msd_CHOL.xvg --gromacs-dim 2

# 3. animation (GIF) of the tracks
python lipid_msd.py animate results --lipid CHL --plane xy
```

Only the `msd` step needs MDAnalysis and the trajectory. `plot` and `animate` work on the files in
`results/`.

## Output (`results/`)

| File | Content |
|---|---|
| `tracks.npz` | unwrapped, drift-corrected COM tracks (nm) of all residues |
| `msd_<lipid>.npz` / `.csv` | lags (ps) and MSD per lipid (npz), mean ± SEM over lipids (csv) |
| `msd.png`, `diffusion.csv` | plots and fitted D per lipid type |
| `<lipid>_<plane>.gif` | animation |

## Method

- **COM of a lipid:** `r_COM = Σ mᵢ rᵢ / Σ mᵢ`. Before that, all atoms of a residue are shifted to the
  periodic image closest to its first atom, so lipids that cross the box boundary are not torn apart.
- **Unwrapping in time:** frame-to-frame displacements are folded back into [−L/2, L/2] and summed.
  This assumes that no lipid moves more than half a box between two frames.
- **Drift:** the mass-weighted COM of all residues is subtracted in every frame.
- **MSD:** `MSD(τ) = ⟨ |r(t+τ) − r(t)|² ⟩`, averaged over all start times t (sliding window, only x and y).
  The largest lag is half the trajectory by default, so that every lag is averaged over many windows.
- **Diffusion coefficient:** in two dimensions `MSD = 4·D·τ`, so `D = slope / 4`.
  Unit conversion: 1 nm²/ps = 10⁻² cm²/s = 10³ · 10⁻⁵ cm²/s.
- **Fit range:** 10–50 % of the largest lag by default (`--fit`). Very short lags show sub-diffusive
  behaviour, and the longest lags are poorly sampled.

## Validation

`python test_lipid_msd.py` runs tests on synthetic random walks with a known D (no trajectory needed):

- free diffusion: recovered D within 5 %
- drift: a common drift inflates D, `remove_drift` corrects it
- periodic boundaries: two-atom lipids split by the box, D recovered after COM and unwrapping step,
  while the naive wrapped coordinates give a wrong D (> 20 % off)
- unit conversion and leaflet selection

I have **not** tested the trajectory reader (`load_tracks`) against a real file in this form. Before using the
numbers, run it on a few frames and compare D to `gmx msd`. Statistical uncertainty is considerable with few
lipids or short runs (about 10 % with 60 lipids and 400 frames in my synthetic check); the SEM band in the plot
shows the spread between lipids, but windows overlap, so it underestimates the true error.

## Limitations

- Orthorhombic boxes only; leaflets are assigned by the sign of z in the first frame, which assumes that the
  bilayer is not split by the z boundary (check that the printed upper/lower counts are about equal).
- Only lateral (xy) diffusion; no finite-size correction (Saffman–Delbrück / Yeh–Hummer).
- Time-averaged MSD of single lipids; no anomalous diffusion exponent.

## Changes compared to my first scripts

- D is now `slope / 4` (2D); before, the raw slope was used.
- Files for CHL, DAPC and DSPC were mixed up when plotting; fixed.
- Periodic boundaries are handled explicitly, and drift is removed on unwrapped tracks.
- Number of frames and `dt` are no longer hard-coded; one animation function instead of three copies.
- Computation reads the trajectory once and uses array slicing instead of re-reading frames for every lag.
