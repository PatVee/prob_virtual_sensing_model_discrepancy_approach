"""Rebuild the KW51 probabilistic-ensemble manuscript figures in Python.

Python counterpart of ``plot_kw51_manuscript_figures.m``. It operates on the
frozen MATLAB run products and does not rerun the FE inference.

Figures written
---------------
kw_fe_mesh.pdf/png
kw_inputs.pdf/png
kw_posterior_load_benchmark.pdf/png
kw_posterior_load_benchmark_17passage_summary.pdf/png
kw_ensemble_validation_summary.pdf/png
kw_prediction_histories_trackB.pdf/png
kw_prediction_histories_trackA.pdf/png

The FE mesh panel is a Python structural schematic generated from the released
KW51 coordinate system and the seven sensor stations. All run-dependent plots
are reconstructed directly from the saved inference/validation products.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional, Sequence, Tuple, Union
import csv
import re

import numpy as np
from scipy.io import loadmat
from scipy import sparse
import matplotlib.pyplot as plt

from plot_kw51_selected_event import plot_kw51_selected_event, resolve_run_folder

PathLike = Union[str, Path]


def _loadmat(path: PathLike) -> Dict[str, Any]:
    return loadmat(str(path), simplify_cells=True)


def _set_manuscript_font() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "font.size": 9.5,
            "axes.titlesize": 10,
            "axes.labelsize": 10,
            "legend.fontsize": 8,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def _save(fig: plt.Figure, folder: Path, stem: str, visible: bool) -> Tuple[Path, Path]:
    pdf = folder / f"{stem}.pdf"
    png = folder / f"{stem}.png"
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(png, dpi=220, bbox_inches="tight")
    if not visible:
        plt.close(fig)
    return pdf, png


def _format_axis(ax: plt.Axes) -> None:
    ax.grid(True, alpha=0.25)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(direction="out")


def _as_inputs(value: Any) -> list[dict]:
    if isinstance(value, dict):
        return [value]
    if isinstance(value, list):
        return value
    if isinstance(value, np.ndarray):
        return list(value.reshape(-1))
    return []




def _resolve_input_file(run_folder: Path, item: dict) -> Path:
    """Resolve a catalog input path, repairing stale absolute paths if needed."""
    raw = Path(str(item.get("input_file", "")))
    if raw.is_file():
        return raw
    event = str(item["event_id"])
    local = run_folder / "inputs" / f"{event}_input.mat"
    if local.is_file():
        return local
    raise FileNotFoundError(f"Could not resolve saved input file for {event}")

def _resolve_package_root(run_folder: Path, package_root: Optional[PathLike]) -> Path:
    if package_root is not None:
        p = Path(package_root).expanduser().resolve()
        if not (p / "fe_models" / "model_15.mat").is_file():
            raise FileNotFoundError(f"model_15.mat was not found below package root {p}")
        return p

    search = [run_folder, *run_folder.parents]
    for p in search:
        if (p / "fe_models" / "model_15.mat").is_file():
            return p
        # Normal layout: <root>/KW51ProbabilisticResults/run_...
        if (p.parent / "fe_models" / "model_15.mat").is_file():
            return p.parent
    raise FileNotFoundError(
        "Could not locate fe_models/model_15.mat automatically. Pass package_root=... "
        "to plot_kw51_manuscript_figures()."
    )


def _sensor_locations() -> list[dict]:
    return [
        {"key": "R1", "channel": "sgBD1011A", "role": "reference", "x": 32.5915, "y": 0.0, "z": 0.2719},
        {"key": "R2", "channel": "sgBD1415A", "role": "reference", "x": 48.3245, "y": 0.0, "z": 0.2719},
        {"key": "R3", "channel": "sgBD2728A", "role": "reference", "x": 99.3080, "y": 0.0, "z": 0.2719},
        {"key": "T1", "channel": "sgBD1718C", "role": "target_strain", "x": 60.4245, "y": 11.75, "z": 0.2719},
        {"key": "T2", "channel": "sgBD2324C", "role": "target_strain", "x": 83.7245, "y": 11.75, "z": 0.2719},
        {"key": "D1", "channel": "dABLA", "role": "target_disp", "x": 115.2570, "y": 0.0, "z": -0.15},
        {"key": "D2", "channel": "dABLC", "role": "target_disp", "x": 115.2570, "y": 11.75, "z": -0.15},
    ]


def _make_fe_mesh(folder: Path, visible: bool) -> None:
    """Create the Python structural schematic used as the FE-location figure.

    The full 6079-node/4857-element geometry is generated natively by the
    MATLAB Stabil model. Python does not execute Stabil here; this panel uses
    the released bridge coordinate system, bowstring arch equation, transverse
    deck grid, representative diagonals, and exact sensor operator stations.
    """
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

    sensors = _sensor_locations()
    x = np.linspace(0.0, 116.316, 120)
    z_arch = -0.0056679878 * x**2 + 0.659277672 * x + 0.341695
    z_arch = np.maximum(z_arch, 0.2719)

    fig = plt.figure(figsize=(11.7, 6.1), constrained_layout=True)
    gs = fig.add_gridspec(2, 1, height_ratios=[2.1, 1.0])
    ax = fig.add_subplot(gs[0], projection="3d")

    # Main girders and arches.
    for y in (0.0, 11.75):
        ax.plot(x, np.full_like(x, y), np.full_like(x, 0.2719), linewidth=1.0)
        ax.plot(x, np.full_like(x, y), z_arch, linewidth=1.0)

    # Representative transverse deck mesh.
    for xx in np.linspace(0, 116.316, 48):
        ax.plot([xx, xx], [0, 11.75], [0.2719, 0.2719], linewidth=0.35, alpha=0.55)
    for yy in np.linspace(0, 11.75, 9):
        ax.plot(x, np.full_like(x, yy), np.full_like(x, 0.2719), linewidth=0.25, alpha=0.45)

    # Representative bowstring diagonals.
    lower = np.array([11.258, 22.758, 34.558, 46.358, 58.158, 69.958, 81.758, 93.558, 105.058])
    upper = np.array([17.158, 28.658, 40.458, 52.258, 64.058, 75.858, 87.658, 99.458])
    for y in (0.0, 11.75):
        for a, b in zip(lower[:-1], upper):
            za = 0.2719
            zb = max(-0.0056679878 * b**2 + 0.659277672 * b + 0.341695, 0.2719)
            ax.plot([a, b], [y, y], [za, zb], linewidth=0.8)
        for a, b in zip(upper, lower[1:]):
            za = max(-0.0056679878 * a**2 + 0.659277672 * a + 0.341695, 0.2719)
            ax.plot([a, b], [y, y], [za, 0.2719], linewidth=0.8)

    for s in sensors:
        marker = "s" if s["role"] == "reference" else ("o" if s["role"] == "target_strain" else "D")
        ax.scatter(s["x"], s["y"], s["z"], marker=marker, s=34, edgecolors="black", linewidths=0.8)

    ax.set_xlabel("Longitudinal x (m)")
    ax.set_ylabel("y (m)")
    ax.set_zlabel("z (m)")
    ax.set_title("(a) FE geometry with measurement locations", fontweight="normal")
    ax.view_init(elev=19, azim=-57)
    ax.grid(True, alpha=0.22)

    ax2 = fig.add_subplot(gs[1])
    ax2.plot([0, 116.316], [1, 1], "k-", linewidth=1.0)
    ax2.plot([0, 116.316], [0, 0], "k-", linewidth=1.0)
    ax2.text(-3, 1, "A", ha="right", va="center")
    ax2.text(-3, 0, "C", ha="right", va="center")
    for s in sensors:
        yy = 1.0 - s["y"] / 11.75
        marker = "s" if s["role"] == "reference" else ("o" if s["role"] == "target_strain" else "D")
        ax2.scatter(s["x"], yy, marker=marker, s=42, edgecolors="black", linewidths=0.8)
        dy = 0.18 if yy >= 0.5 else -0.18
        ax2.text(s["x"], yy + dy, f"{s['key']}: {s['channel']}", ha="center", va="center", fontsize=8.2)
    ax2.set_xlim(-4, 120)
    ax2.set_ylim(-0.45, 1.45)
    ax2.set_yticks([])
    ax2.set_xlabel("Longitudinal FE coordinate x (m)")
    ax2.set_title("(b) Girder-axis locations and sensor identities", fontweight="normal")
    ax2.spines[["top", "right", "left"]].set_visible(False)

    # Dummy legend handles.
    ax2.scatter([], [], marker="s", edgecolors="black", label="Reference strain")
    ax2.scatter([], [], marker="o", facecolors="none", edgecolors="black", label="Withheld strain")
    ax2.scatter([], [], marker="D", facecolors="none", edgecolors="black", label="Withheld displacement")
    ax2.legend(loc="lower center", bbox_to_anchor=(0.5, -0.48), ncol=3, frameon=False)

    _save(fig, folder, "kw_fe_mesh", visible)


def _column_quantiles(X: np.ndarray, probs: Sequence[float]) -> np.ndarray:
    return np.quantile(np.asarray(X, dtype=float), probs, axis=0)


def _short_event_label(event: str) -> str:
    m = re.match(r"traindata_(\d{8})_(\d{6})$", event)
    if not m:
        return event
    d, t = m.groups()
    return f"{d[6:8]}/{d[4:6]} {t[0:2]}:{t[2:4]}"


def _make_inputs_figure(run_folder: Path, inputs: list[dict], folder: Path, visible: bool) -> None:
    event_id = "traindata_20191217_170709"
    match = next((x for x in inputs if str(x["event_id"]) == event_id), None)
    if match is None:
        raise RuntimeError(f"Manuscript input event {event_id} is absent from this run.")
    E = _loadmat(_resolve_input_file(run_folder, match))
    D = _loadmat(run_folder / "draws" / f"{event_id}_ensemble_draws.mat")

    fig, axes = plt.subplots(2, 2, figsize=(11.7, 7.1), constrained_layout=True)
    ax = axes[0, 0]
    y_ref = np.asarray(E["y_ref"], dtype=float)
    names = [str(x) for x in np.asarray(E["reference_names"]).reshape(-1)]
    for j in range(y_ref.shape[1]):
        ax.plot(np.asarray(E["time_s"]).reshape(-1), y_ref[:, j], linewidth=1.0, label=names[j])
    ax.set_xlabel("Record time (s)")
    ax.set_ylabel("Strain (µε)")
    ax.set_title("(a) Reference strain histories", fontweight="normal")
    ax.legend(loc="best")
    _format_axis(ax)

    ax = axes[0, 1]
    tau = np.asarray(E["rail_offsets_s"], dtype=float).reshape(-1)
    for value in tau:
        ax.plot([value, value], [0, 1], linewidth=0.85)
    ax.set_xlim(min(0, tau.min() - 0.2), tau.max() + 0.2)
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("Time from first rail peak (s)")
    ax.set_ylabel("Peak indicator")
    ax.set_title(f"(b) Rail-derived candidate timings ({len(tau)} peaks)", fontweight="normal")
    _format_axis(ax)

    ax = axes[1, 0]
    loads = np.asarray(D["loads_100kN"], dtype=float) * 100.0
    med, lo, hi = _column_quantiles(loads, [0.50, 0.05, 0.95])
    xx = np.arange(1, loads.shape[1] + 1)
    ax.fill_between(xx, lo, hi, alpha=0.25, label="5th–95th percentiles")
    ax.plot(xx, med, "-o", markersize=3.0, linewidth=1.0, label="Median")
    ax.set_xlabel("Candidate axle number")
    ax.set_ylabel("Load (kN)")
    ax.set_title(f"(c) Prior load ensemble: {loads.shape[0]} candidates", fontweight="normal")
    ax.legend(loc="best")
    _format_axis(ax)

    ax = axes[1, 1]
    track_colors = {"A": "tab:blue", "B": "tab:orange"}
    for item in inputs:
        Ein = _loadmat(_resolve_input_file(run_folder, item))
        event = str(Ein["event_id"])
        vpath = run_folder / "validation" / f"{event}_validation.mat"
        if not vpath.is_file():
            continue
        V = _loadmat(vpath)
        speed = float(np.asarray(V["kinematics"]).reshape(-1)[0]) * 3.6
        peaks = len(np.asarray(Ein["rail_offsets_s"]).reshape(-1))
        track = str(Ein["track"])
        ax.plot(speed, peaks, "o", markersize=4, color=track_colors.get(track, "0.4"), label=None)
    ax.plot([], [], "o", color=track_colors["A"], label="Track A")
    ax.plot([], [], "o", color=track_colors["B"], label="Track B")
    ax.set_xlabel("Reference-inferred entry speed (km/h)")
    ax.set_ylabel("Detected rail peaks")
    ax.set_title(f"(d) All {len(inputs)} accepted passages", fontweight="normal")
    ax.legend(loc="best")
    _format_axis(ax)

    _save(fig, folder, "kw_inputs", visible)


def _smooth_basis(t: np.ndarray, max_frequency: float) -> np.ndarray:
    t = np.asarray(t, dtype=float).reshape(-1)
    duration = t[-1] - t[0]
    if duration <= 0:
        raise ValueError("Time vector must have positive duration.")
    nfreq = int(np.floor(duration * max_frequency))
    if nfreq == 0:
        return np.ones((len(t), 1))
    freq = np.arange(1, nfreq + 1, dtype=float) / duration
    return np.column_stack(
        [
            np.ones_like(t),
            np.cos(2 * np.pi * np.outer(t, freq)),
            np.sin(2 * np.pi * np.outer(t, freq)),
        ]
    ) / np.sqrt(1.0 + nfreq)


def _kw51_basis_refs(B: Dict[str, Any], E: Dict[str, Any], t: np.ndarray, kin: np.ndarray) -> np.ndarray:
    """Exact Python translation of kw51_basis for reference outputs 1:3."""
    t = np.asarray(t, dtype=float).reshape(-1)
    tau = np.asarray(E["rail_offsets_s"], dtype=float).reshape(-1)
    v, entry, direction, acc = np.asarray(kin, dtype=float).reshape(-1)
    dt = t - entry
    vr = v + acc * (float(E["rail_first_time_s"]) - entry)
    x = v * dt[:, None] + 0.5 * acc * dt[:, None] ** 2 - (vr * tau[None, :] + 0.5 * acc * tau[None, :] ** 2)
    if direction < 0:
        x = 116.316 - x
    speed2 = (v + acc * dt) ** 2
    track_index = 0 if str(E["track"]) == "A" else 1
    grid = np.asarray(B["xgrid"], dtype=float).reshape(-1)
    dx = grid[1] - grid[0]
    valid = (x >= grid[0]) & (x <= grid[-1])
    q = (x - grid[0]) / dx
    idx = np.clip(np.floor(q).astype(int), 0, len(grid) - 2)
    a = q - idx

    H = np.zeros((len(t), 3, len(tau)), dtype=float)
    H0 = np.asarray(B["H0"], dtype=float)
    H2 = np.asarray(B["H2"], dtype=float)
    for j in range(3):
        h0 = H0[:, j, track_index]
        h2 = H2[:, j, track_index]
        u0 = h0[idx]
        d0 = h0[idx + 1] - h0[idx]
        u2 = h2[idx]
        d2 = h2[idx + 1] - h2[idx]
        z = ((u0 + a * d0) + speed2[:, None] * (u2 + a * d2)) * valid
        H[:, j, :] = z * 1e5
    return H


def _posterior_weights(
    B: Dict[str, Any],
    E: Dict[str, Any],
    V: Dict[str, Any],
    c: Dict[str, Any],
    D: Dict[str, Any],
) -> Tuple[np.ndarray, float]:
    """Replay the saved ensemble conditioning exactly for the three references."""
    t = np.asarray(E["time_s"], dtype=float).reshape(-1)
    kin = np.asarray(V["kinematics"], dtype=float).reshape(-1)
    H = _kw51_basis_refs(B, E, t, kin)
    T, R, n_axles = H.shape

    loads = np.asarray(D["loads_100kN"], dtype=float)
    kappa = np.asarray(D["stiffness_states"], dtype=float).reshape(-1)
    N = len(kappa)
    nominal = np.zeros((N, T, R), dtype=float)
    for r in range(R):
        nominal[:, :, r] = (loads @ H[:, r, :].T) / kappa[:, None]

    pairloads = np.asarray(D["paired_loads_100kN"], dtype=float)
    M = pairloads.shape[0]
    nominal_pair = np.zeros((M, T, R), dtype=float)
    pair_kappa = float(c["pair_nominal_stiffness"])
    for r in range(R):
        nominal_pair[:, :, r] = (pairloads @ H[:, r, :].T) / pair_kappa

    psi = _smooth_basis(t, float(c["support_frequency_limit_Hz"]))
    factors = np.asarray(D["paired_stiffness_factors"], dtype=float).reshape(-1)
    recovery = np.asarray(D["paired_recovery"], dtype=float)[:, :3]
    kins = np.asarray(D["paired_kinematics"], dtype=float)
    forcecoef = np.asarray(D["support_force_coefficients"], dtype=float)
    Graw = B["G"]
    Garr = Graw.toarray() if sparse.issparse(Graw) else np.asarray(Graw, dtype=float)
    if Garr.shape[1] != 2:
        raise RuntimeError("Run FE bank is not the two-coordinate blueprint bank.")
    Gref = Garr[:3, :]

    delta = np.zeros((M, T, R), dtype=float)
    for m in range(M):
        forces = psi @ forcecoef[m] * float(c["support_force_sd_N"])
        # Blueprint BC: Brussels anchors longitudinal rigid motion; do not
        # remove the mean of the two Leuven support-force histories.
        Hp = _kw51_basis_refs(B, E, t, kins[m])
        yp = np.zeros((T, R), dtype=float)
        for r in range(R):
            yp[:, r] = Hp[:, r, :] @ pairloads[m]
        scale = pair_kappa * factors[m]
        yp = yp / scale + (forces @ Gref.T) / scale
        yp *= recovery[m][None, :]
        delta[m] = yp - nominal_pair[m]

    mu = delta.mean(axis=0)
    dc = delta - mu[None, :, :]
    y = np.asarray(E["y_ref"], dtype=float)
    ll = np.zeros(N, dtype=float)
    stride = int(c["likelihood_time_stride"])
    likelihood_power = float(c["likelihood_power"])
    ref_sigma2 = float(c["reference_noise_sigma_microstrain"]) ** 2
    cov_ridge = float(c["cov_ridge"])

    for it in range(0, T, stride):
        Dref = dc[:, it, :]
        C = (Dref.T @ Dref) / (M - 1)
        C = 0.5 * (C + C.T)
        ridge = cov_ridge * max(np.trace(C) / R, 1e-12) + 1e-12
        S = C + np.eye(R) * (ref_sigma2 + ridge)
        Sinv = np.linalg.solve(S, np.eye(R))
        innovation = y[it, :] - nominal[:, it, :] - mu[it, :]
        ll -= 0.5 * likelihood_power * np.sum((innovation @ Sinv) * innovation, axis=1)

    weights = np.exp(ll - ll.max())
    weights /= weights.sum()
    ess = 1.0 / np.sum(weights**2)
    return weights, ess


def _weighted_summary(X: np.ndarray, w: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    X = np.asarray(X, dtype=float)
    w = np.asarray(w, dtype=float).reshape(-1)
    w /= w.sum()
    mean = w @ X
    lo = np.empty(X.shape[1])
    hi = np.empty(X.shape[1])
    for j in range(X.shape[1]):
        order = np.argsort(X[:, j])
        v = X[order, j]
        cw = np.cumsum(w[order])
        lo[j] = np.interp(0.025, np.r_[0.0, cw], np.r_[v[0], v])
        hi[j] = np.interp(0.975, np.r_[0.0, cw], np.r_[v[0], v])
    return mean, lo, hi


def _make_load_benchmark_figures(
    run_folder: Path,
    bank_file: Path,
    inputs: list[dict],
    folder: Path,
    visible: bool,
) -> Tuple[list[str], list[list[Any]]]:
    c = _loadmat(run_folder / "settings.mat")
    B = _loadmat(bank_file)

    benchmark_axle = np.array([132.0, 131.4, 116.8, 116.4, 94.7, 94.6, 94.6, 94.7, 116.4, 116.8, 131.4, 132.0])
    benchmark_bogie = benchmark_axle.reshape(6, 2).sum(axis=1)
    benchmark_total = benchmark_axle.sum()

    rows: list[list[Any]] = []
    example = None
    for item in sorted(inputs, key=lambda x: str(x["event_id"])):
        E = _loadmat(_resolve_input_file(run_folder, item))
        offsets = np.asarray(E["rail_offsets_s"]).reshape(-1)
        if len(offsets) != 12:
            continue
        event = str(E["event_id"])
        D = _loadmat(run_folder / "draws" / f"{event}_ensemble_draws.mat")
        V = _loadmat(run_folder / "validation" / f"{event}_validation.mat")
        weights, ess = _posterior_weights(B, E, V, c, D)
        saved_ess = float(V["posterior_ess"])
        if abs(ess - saved_ess) >= 1e-7:
            raise RuntimeError(f"Posterior replay mismatch for {event}: {ess} versus {saved_ess}")

        axle = np.asarray(D["loads_100kN"], dtype=float) * 100.0
        am, al, ah = _weighted_summary(axle, weights)
        bogie = axle.reshape(axle.shape[0], 6, 2).sum(axis=2)
        bm, bl, bh = _weighted_summary(bogie, weights)
        total = axle.sum(axis=1, keepdims=True)
        tm, tl, th = _weighted_summary(total, weights)
        tm, tl, th = float(tm[0]), float(tl[0]), float(th[0])

        rows.append(
            [
                event,
                str(V["track"]),
                ess,
                tm,
                tl,
                th,
                benchmark_total,
                float(np.mean(np.abs(am - benchmark_axle))),
                float(np.mean(np.abs(bm - benchmark_bogie))),
            ]
        )
        if example is None:
            example = {
                "event": event,
                "track": str(V["track"]),
                "am": am,
                "al": al,
                "ah": ah,
                "bm": bm,
                "bl": bl,
                "bh": bh,
            }

    if example is None:
        raise RuntimeError("No accepted 12-peak passages found.")

    fig, axes = plt.subplots(1, 2, figsize=(11.7, 5.0), constrained_layout=True)
    ax = axes[0]
    x = np.arange(1, 13)
    ax.errorbar(x, example["am"], yerr=np.vstack([example["am"] - example["al"], example["ah"] - example["am"]]), fmt="o", capsize=4, linewidth=1.0, label="Posterior mean and 95% interval")
    ax.plot(x, benchmark_axle, "s-", linewidth=1.0, label="Measured empty-train benchmark")
    ax.set_xticks(x)
    ax.set_xlabel("Axle number")
    ax.set_ylabel("Load (kN)")
    ax.set_title("(a) Axle-load comparison", fontweight="normal")
    _format_axis(ax)

    ax = axes[1]
    x = np.arange(1, 7)
    ax.errorbar(x, example["bm"], yerr=np.vstack([example["bm"] - example["bl"], example["bh"] - example["bm"]]), fmt="o", capsize=4, linewidth=1.0, label="Posterior mean and 95% interval")
    ax.plot(x, benchmark_bogie, "s-", linewidth=1.0, label="Measured empty-train benchmark")
    ax.set_xticks(x)
    ax.set_xlabel("Bogie number")
    ax.set_ylabel("Load (kN)")
    ax.set_title("(b) Bogie-load comparison", fontweight="normal")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.18), frameon=False)
    _format_axis(ax)
    fig.suptitle(f"{example['event']} | Track {example['track']} | 12 detected rail peaks", fontsize=10.5, fontweight="normal")
    _save(fig, folder, "kw_posterior_load_benchmark", visible)

    event_names = [r[0] for r in rows]
    total_mean = np.array([r[3] for r in rows], dtype=float)
    total_lo = np.array([r[4] for r in rows], dtype=float)
    total_hi = np.array([r[5] for r in rows], dtype=float)
    bogie_mae = np.array([r[8] for r in rows], dtype=float)

    fig, axes = plt.subplots(1, 2, figsize=(11.7, 5.2), constrained_layout=True)
    x = np.arange(1, len(rows) + 1)
    labels = [_short_event_label(x) for x in event_names]
    ax = axes[0]
    ax.errorbar(x, total_mean, yerr=np.vstack([total_mean - total_lo, total_hi - total_mean]), fmt="o", capsize=3, linewidth=0.9)
    ax.axhline(benchmark_total, linestyle="--", linewidth=1.2)
    ax.set_xticks(x, labels, rotation=55, ha="right")
    ax.set_ylabel("Total train load (kN)")
    ax.set_xlabel("Passage")
    ax.set_title("(a) Posterior total load", fontweight="normal")
    _format_axis(ax)

    ax = axes[1]
    ax.bar(x, bogie_mae)
    ax.set_xticks(x, labels, rotation=55, ha="right")
    ax.set_ylabel("Mean absolute bogie-load difference (kN)")
    ax.set_xlabel("Passage")
    ax.set_title("(b) Bogie-load difference from benchmark", fontweight="normal")
    _format_axis(ax)
    _save(fig, folder, "kw_posterior_load_benchmark_17passage_summary", visible)

    headers = [
        "event_id",
        "track",
        "posterior_ESS",
        "posterior_total_mean_kN",
        "posterior_total_q025_kN",
        "posterior_total_q975_kN",
        "benchmark_total_kN",
        "axle_mean_MAE_kN",
        "bogie_mean_MAE_kN",
    ]
    return headers, rows


def _reporting_split(summary: np.ndarray) -> str:
    values = {str(x) for x in summary[:, 0]}
    if "evaluation" in values:
        return "evaluation"
    if "development" in values:
        return "development"
    raise RuntimeError("No recognised saved reporting split found.")


def _summary_row(summary: np.ndarray, split: str, target: str, sigma: float) -> Optional[np.ndarray]:
    for row in summary:
        if str(row[0]) == split and str(row[1]) == target and abs(float(row[2]) - sigma) < 1e-12:
            return row
    return None


def _make_validation_summary(summary: np.ndarray, c: Dict[str, Any], folder: Path, visible: bool) -> None:
    targets = [str(x) for x in np.asarray(c["protected_names"]).reshape(-1)]
    split = _reporting_split(summary)
    sigma = float(c["translation_sd_mm"])
    values = np.full((4, 3), np.nan)
    for j, target in enumerate(targets):
        row = _summary_row(summary, split, target, sigma)
        if row is not None:
            values[j] = [float(row[4]), 100.0 * float(row[7]), float(row[8])]

    fig, axes = plt.subplots(1, 3, figsize=(11.7, 4.5), constrained_layout=True)
    titles = ["(a) Median prediction error", "(b) Nominal interval coverage", "(c) Interval width"]
    ylabels = ["Range-normalised RMSE", "Mean coverage (%)", "Median width / measured range"]
    x = np.arange(1, 5)
    for p, ax in enumerate(axes):
        ax.bar(x, values[:, p])
        if p == 0:
            ax.axhline(0.15, color="black", linestyle=":", linewidth=1.1)
        elif p == 1:
            ax.axhline(95, color="black", linestyle="--", linewidth=1.0)
            ax.set_ylim(0, 110)
        ax.set_xticks(x, targets, rotation=25, ha="right")
        ax.set_ylabel(ylabels[p])
        ax.set_title(titles[p], fontweight="normal")
        _format_axis(ax)
    _save(fig, folder, "kw_ensemble_validation_summary", visible)


def plot_kw51_manuscript_figures(
    run: PathLike,
    *,
    visible: bool = False,
    output_folder: Optional[PathLike] = None,
    package_root: Optional[PathLike] = None,
) -> Path:
    """Rebuild all manuscript figures from a completed MATLAB run."""
    _set_manuscript_font()
    run_folder = resolve_run_folder(run)
    # package_root is retained for call compatibility; all run-dependent
    # posterior replay uses the exact FE bank frozen inside the run folder.
    bank_file = run_folder / "fe_bank" / "model_15.mat"
    if not bank_file.is_file():
        raise FileNotFoundError(f"Missing run-local blueprint FE bank: {bank_file}")
    folder = Path(output_folder).expanduser().resolve() if output_folder else run_folder / "manuscript_figures_python"
    folder.mkdir(parents=True, exist_ok=True)

    required = [run_folder / "validation_results.mat", run_folder / "settings.mat", run_folder / "input_catalog.mat"]
    for path in required:
        if not path.is_file():
            raise FileNotFoundError(f"Required run product is missing: {path}")

    c = _loadmat(run_folder / "settings.mat")
    R = _loadmat(run_folder / "validation_results.mat")
    I = _loadmat(run_folder / "input_catalog.mat")
    inputs = _as_inputs(I["inputs"])

    _make_fe_mesh(folder, visible)
    _make_inputs_figure(run_folder, inputs, folder, visible)
    headers, rows = _make_load_benchmark_figures(run_folder, bank_file, inputs, folder, visible)
    _make_validation_summary(np.asarray(R["summary"], dtype=object), c, folder, visible)

    # Exact two passage histories used in the manuscript.
    plot_kw51_selected_event(
        run_folder,
        "traindata_20191201_080952",
        visible=visible,
        output_folder=folder,
        output_name="kw_prediction_histories_trackB",
    )
    plot_kw51_selected_event(
        run_folder,
        "traindata_20191201_171118",
        visible=visible,
        output_folder=folder,
        output_name="kw_prediction_histories_trackA",
    )

    with (folder / "kw_load_benchmark_summary.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        writer.writerows(rows)

    print(f"\nAll Python manuscript figures written to:\n{folder}")
    return folder


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", help="run_* folder or KW51ProbabilisticResults folder")
    parser.add_argument("--package-root", default=None, help="KW51 probabilistic package root if auto-detection fails")
    parser.add_argument("--output-folder", default=None)
    parser.add_argument("--show", action="store_true", help="leave figures open after saving")
    args = parser.parse_args()

    plot_kw51_manuscript_figures(
        args.run,
        visible=args.show,
        output_folder=args.output_folder,
        package_root=args.package_root,
    )
