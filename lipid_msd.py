#!/usr/bin/env python3
"""
lipid_msd.py - lateral mean squared displacement (MSD) and diffusion coefficients
of lipids in a membrane simulation.

Subcommands
    msd      read a trajectory, compute per-lipid MSD, save results
    plot     plot MSD curves, fit the diffusion coefficient D
    animate  animated trajectories of lipid centres of mass (GIF)

Example
    python lipid_msd.py msd mem.pdb --dt 100 --lipids CHL DAPC DSPC --out results
    python lipid_msd.py plot results --fit 0.1 0.5
    python lipid_msd.py animate results --lipid CHL --plane xy

Units: lengths in nm, times in ps, D reported in 1e-5 cm^2/s.
"""

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np

# 1 nm^2/ps = 1e-14 cm^2 / 1e-12 s = 1e-2 cm^2/s = 1e3 * (1e-5 cm^2/s)
NM2_PS_TO_1E5_CM2_S = 1.0e3


# ---------------------------------------------------------------------------
# Core numerics (pure NumPy, no MDAnalysis needed; covered by test_lipid_msd.py)
# ---------------------------------------------------------------------------
def residue_coms(positions, masses, res_idx, n_res, box):
    """Mass-weighted centre of mass of every residue.

    Residues split by the periodic boundary are made whole first (every atom is
    shifted to the periodic image closest to the first atom of its residue).

    positions : (n_atoms, 3)   res_idx : (n_atoms,) integers 0..n_res-1
    masses    : (n_atoms,)     box     : (3,) orthorhombic box lengths
    returns   : (n_res, 3) in the units of `positions`
    """
    box = np.asarray(box, dtype=float)
    _, first = np.unique(res_idx, return_index=True)
    ref = positions[first][res_idx]
    d = positions - ref
    d -= box * np.round(d / box)
    whole = ref + d
    m = np.bincount(res_idx, weights=masses, minlength=n_res)
    com = np.stack(
        [np.bincount(res_idx, weights=masses * whole[:, k], minlength=n_res) for k in range(3)],
        axis=1,
    )
    return com / m[:, None]


def unwrap_time(com, boxes):
    """Remove periodic jumps from centre-of-mass tracks.

    com   : (frames, residues, 3), boxes : (frames, 3)
    Assumes that no lipid moves more than half a box length between two frames.
    """
    d = np.diff(com, axis=0)
    b = np.asarray(boxes)[1:, None, :]
    d -= b * np.round(d / b)
    out = np.empty_like(com)
    out[0] = com[0]
    out[1:] = com[0] + np.cumsum(d, axis=0)
    return out


def remove_drift(com, weights):
    """Subtract the (weighted) centre of mass of all residues in every frame."""
    w = np.asarray(weights, dtype=float)
    w = w / w.sum()
    return com - (com * w[None, :, None]).sum(axis=1, keepdims=True)


def compute_msd(xy, max_lag):
    """Time-averaged MSD for lags 0..max_lag (in frames).

    xy : (frames, residues, 2)
    returns (max_lag + 1, residues); row 0 is zero.
    """
    n_frames = xy.shape[0]
    if max_lag >= n_frames:
        raise ValueError(f"max_lag ({max_lag}) must be smaller than the number of frames ({n_frames})")
    msd = np.zeros((max_lag + 1, xy.shape[1]))
    for tau in range(1, max_lag + 1):
        disp = xy[tau:] - xy[:-tau]
        msd[tau] = (disp ** 2).sum(axis=-1).mean(axis=0)
    return msd


def fit_diffusion(lags_ps, msd_mean, frac=(0.1, 0.5), dim=2):
    """Linear fit MSD = 2*dim*D*t + c over a fraction of the lag range.

    Returns (D in 1e-5 cm^2/s, slope in nm^2/ps, intercept in nm^2).
    dim = 2 for lateral (xy) diffusion: D = slope / 4.
    """
    lags_ps = np.asarray(lags_ps, dtype=float)
    t_max = lags_ps[-1]
    sel = (lags_ps >= frac[0] * t_max) & (lags_ps <= frac[1] * t_max)
    if sel.sum() < 3:
        raise ValueError("fit range contains fewer than 3 points")
    slope, intercept = np.polyfit(lags_ps[sel], np.asarray(msd_mean)[sel], 1)
    return slope / (2 * dim) * NM2_PS_TO_1E5_CM2_S, slope, intercept


def leaflet_mask(z0, leaflet):
    """Select the upper/lower leaflet from centred z positions of the first frame."""
    z0 = np.asarray(z0)
    if leaflet == "upper":
        return z0 > 0
    if leaflet == "lower":
        return z0 <= 0
    return np.ones(z0.shape, dtype=bool)


# ---------------------------------------------------------------------------
# Trajectory I/O (needs MDAnalysis)
# ---------------------------------------------------------------------------
def load_tracks(path, selection="all", max_frames=None):
    """Read a trajectory and return unwrapped, drift-corrected residue COM tracks in nm."""
    try:
        import MDAnalysis as mda
    except ImportError:
        sys.exit("MDAnalysis is required: pip install MDAnalysis")

    u = mda.Universe(str(path))
    ag = u.select_atoms(selection)
    if len(ag) == 0:
        sys.exit(f"Selection '{selection}' matched no atoms.")

    uniq, res_idx = np.unique(ag.resindices, return_inverse=True)
    n_res = len(uniq)
    resnames = np.array(u.residues[uniq].resnames, dtype=str)

    masses = ag.masses.astype(float)
    if np.any(masses <= 0):
        warnings.warn("Some atoms have no mass (guessing failed); using equal weights.")
        masses = np.ones_like(masses)

    coms, boxes = [], []
    for ts in u.trajectory:
        if max_frames is not None and ts.frame >= max_frames:
            break
        dims = ts.dimensions
        if dims is None or np.any(np.asarray(dims[:3]) <= 0):
            sys.exit("No box information (CRYST1) found; it is needed for periodic boundaries.")
        if not np.allclose(dims[3:], 90.0, atol=0.5):
            sys.exit("Only orthorhombic boxes are supported.")
        coms.append(residue_coms(ag.positions, masses, res_idx, n_res, dims[:3]))
        boxes.append(np.array(dims[:3], dtype=float))
        if (ts.frame + 1) % 20 == 0:
            print(f"read frame {ts.frame + 1}", end="\r", flush=True)
    print()

    com = np.array(coms) / 10.0      # Angstrom -> nm
    boxes = np.array(boxes) / 10.0
    com = unwrap_time(com, boxes)
    res_mass = np.bincount(res_idx, weights=masses, minlength=n_res)
    com = remove_drift(com, res_mass)
    upper = com[0, :, 2] > 0
    return com, resnames, upper


# ---------------------------------------------------------------------------
# Subcommands
# ---------------------------------------------------------------------------
def cmd_msd(args):
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    com, resnames, upper = load_tracks(args.trajectory, args.selection, args.max_frames)
    n_frames = com.shape[0]
    print(f"{n_frames} frames, {com.shape[1]} residues; upper leaflet: {upper.sum()}, lower: {(~upper).sum()}")
    np.savez(out / "tracks.npz", com=com, resnames=resnames, upper=upper, dt=args.dt)

    max_lag = args.max_lag if args.max_lag else n_frames // 2
    lags = np.arange(max_lag + 1) * args.dt
    leaf = leaflet_mask(com[0, :, 2], args.leaflet)

    for lipid in args.lipids:
        mask = (resnames == lipid) & leaf
        if not mask.any():
            print(f"[warning] no residues named {lipid} found - skipped")
            continue
        msd = compute_msd(com[:, mask, :2], max_lag)
        np.savez(out / f"msd_{lipid}.npz", lags_ps=lags, msd=msd)
        mean = msd.mean(axis=1)
        sem = msd.std(axis=1, ddof=1) / np.sqrt(msd.shape[1]) if msd.shape[1] > 1 else np.zeros_like(mean)
        np.savetxt(
            out / f"msd_{lipid}.csv",
            np.column_stack([lags, mean, sem]),
            delimiter=",",
            header="lag_ps,msd_mean_nm2,msd_sem_nm2",
            comments="",
        )
        D, _, _ = fit_diffusion(lags, mean)
        print(f"{lipid}: {mask.sum()} lipids, D = {D:.2f} x 1e-5 cm^2/s (fit 10-50 % of the lag range)")
    print(f"Results written to {out.resolve()}")


def load_xvg(path):
    data = np.loadtxt(path, comments=("#", "@"))
    return data[:, 0], data[:, 1]


def cmd_plot(args):
    import matplotlib.pyplot as plt

    out = Path(args.results)
    lipids = args.lipids or [f.stem[4:] for f in sorted(out.glob("msd_*.npz"))]
    if not lipids:
        sys.exit("No msd_*.npz files found. Run the 'msd' subcommand first.")
    gromacs = dict(item.split("=", 1) for item in (args.gromacs or []))

    fig, axes = plt.subplots(len(lipids), 1, figsize=(6, 3.3 * len(lipids)), squeeze=False)
    rows = []
    for ax, lipid in zip(axes[:, 0], lipids):
        data = np.load(out / f"msd_{lipid}.npz")
        lags, msd = data["lags_ps"], data["msd"]
        n = msd.shape[1]
        mean = msd.mean(axis=1)
        sem = msd.std(axis=1, ddof=1) / np.sqrt(n) if n > 1 else np.zeros_like(mean)
        D, slope, icpt = fit_diffusion(lags, mean, tuple(args.fit))

        ax.plot(lags, mean, lw=2, label=f"{lipid} MSD (n = {n})")
        ax.fill_between(lags, mean - sem, mean + sem, alpha=0.25)
        ax.plot(lags, slope * lags + icpt, "--", label=f"fit: D = {D:.2f}e-5 cm$^2$/s")

        if lipid in gromacs:
            t, m = load_xvg(gromacs[lipid])
            Dg, _, _ = fit_diffusion(t, m, tuple(args.fit), dim=args.gromacs_dim)
            ax.plot(t, m, lw=1, label=f"GROMACS: D = {Dg:.2f}e-5 cm$^2$/s")
        ax.set_xlabel("lag time $\\tau$ (ps)")
        ax.set_ylabel("MSD (nm$^2$)")
        ax.legend()
        rows.append((lipid, n, D, args.fit[0] * lags[-1], args.fit[1] * lags[-1]))

    fig.tight_layout()
    fig.savefig(out / "msd.png", dpi=200)
    np.savetxt(
        out / "diffusion.csv",
        np.array([(r[1], r[2], r[3], r[4]) for r in rows]),
        delimiter=",",
        header="n_lipids,D_1e-5_cm2_per_s,fit_from_ps,fit_to_ps",
        comments="",
    )
    for lipid, n, D, lo, hi in rows:
        print(f"{lipid:6s} n = {n:4d}   D = {D:6.2f} x 1e-5 cm^2/s   (fit {lo:.0f}-{hi:.0f} ps)")
    print(f"Figure saved to {(out / 'msd.png').resolve()}")
    if args.show:
        plt.show()


def cmd_animate(args):
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation, PillowWriter

    out = Path(args.results)
    data = np.load(out / "tracks.npz")
    com, resnames, upper = data["com"], data["resnames"], data["upper"]
    mask = (resnames == args.lipid) & leaflet_mask(com[0, :, 2], args.leaflet)
    if not mask.any():
        sys.exit(f"No residues named {args.lipid} found in the selected leaflet.")

    axes_idx = {"xy": (0, 1), "xz": (0, 2)}[args.plane]
    tr = com[:, mask][:, :, list(axes_idx)]
    frames = list(range(0, tr.shape[0], args.step))

    fig, ax = plt.subplots(figsize=(5, 5))
    pad = 0.05 * (tr.max() - tr.min())
    ax.set_xlim(tr[..., 0].min() - pad, tr[..., 0].max() + pad)
    ax.set_ylim(tr[..., 1].min() - pad, tr[..., 1].max() + pad)
    ax.set_xlabel(f"{args.plane[0]} (nm)")
    ax.set_ylabel(f"{args.plane[1]} (nm)")
    lines = [ax.plot([], [], lw=0.6)[0] for _ in range(tr.shape[1])]
    dots = ax.scatter([], [], s=8, c="k")

    def update(i):
        for k, line in enumerate(lines):
            line.set_data(tr[: i + 1, k, 0], tr[: i + 1, k, 1])
        dots.set_offsets(tr[i])
        ax.set_title(f"{args.lipid}, {args.plane} movement (frame {i})")
        return lines + [dots]

    anim = FuncAnimation(fig, update, frames=frames, interval=100)
    target = out / f"{args.lipid}_{args.plane}.gif"
    anim.save(target, writer=PillowWriter(fps=args.fps))
    print(f"Animation saved to {target.resolve()}")
    if args.show:
        plt.show()


# ---------------------------------------------------------------------------
def build_parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="command", required=True)

    m = sub.add_parser("msd", help="compute MSD per lipid from a trajectory")
    m.add_argument("trajectory", help="multi-frame PDB/XTC/... (needs box information)")
    m.add_argument("--dt", type=float, required=True, help="time between frames in ps (required)")
    m.add_argument("--lipids", nargs="+", default=["CHL", "DAPC", "DSPC"], help="residue names")
    m.add_argument("--selection", default="all", help="MDAnalysis selection for the membrane (default: all)")
    m.add_argument("--leaflet", choices=["both", "upper", "lower"], default="both")
    m.add_argument("--max-lag", type=int, default=None, help="largest lag in frames (default: half the frames)")
    m.add_argument("--max-frames", type=int, default=None, help="only use the first N frames (for quick tests)")
    m.add_argument("--out", default="results", help="output folder")
    m.set_defaults(func=cmd_msd)

    q = sub.add_parser("plot", help="plot MSD curves and fit D")
    q.add_argument("results", help="folder written by the msd subcommand")
    q.add_argument("--lipids", nargs="+", default=None)
    q.add_argument("--fit", nargs=2, type=float, default=[0.1, 0.5], metavar=("LO", "HI"),
                   help="fit range as fractions of the maximum lag (default 0.1 0.5)")
    q.add_argument("--gromacs", nargs="+", default=None, metavar="LIPID=FILE.xvg",
                   help="overlay gmx msd output, e.g. CHL=msd_CHOL.xvg")
    q.add_argument("--gromacs-dim", type=int, choices=[2, 3], default=2,
                   help="2 if gmx msd was run with -lateral z, 3 for the default (default 2)")
    q.add_argument("--show", action="store_true", help="open the figure window")
    q.set_defaults(func=cmd_plot)

    a = sub.add_parser("animate", help="animated centre-of-mass trajectories (GIF)")
    a.add_argument("results", help="folder written by the msd subcommand")
    a.add_argument("--lipid", required=True)
    a.add_argument("--plane", choices=["xy", "xz"], default="xy")
    a.add_argument("--leaflet", choices=["both", "upper", "lower"], default="both")
    a.add_argument("--step", type=int, default=2, help="use every N-th frame")
    a.add_argument("--fps", type=int, default=10)
    a.add_argument("--show", action="store_true")
    a.set_defaults(func=cmd_animate)
    return p


def main():
    args = build_parser().parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
