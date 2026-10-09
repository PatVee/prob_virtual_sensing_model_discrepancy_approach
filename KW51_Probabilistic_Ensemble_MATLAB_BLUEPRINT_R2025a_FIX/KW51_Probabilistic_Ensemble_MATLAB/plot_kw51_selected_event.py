"""Plot the four withheld KW51 responses for one accepted passage.

Python counterpart of ``plot_kw51_selected_event.m``.

The preferred selector is an integer 1, 2, ..., N. Accepted event IDs are
sorted chronologically, so 1 is the earliest accepted passage in the run.
A full event ID such as ``traindata_20191201_171118`` is also accepted.

The function reads saved validation products only. It does not refit the FE
model, loads, discrepancy model, or posterior weights.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional, Tuple, Union

import numpy as np
from scipy.io import loadmat
import matplotlib.pyplot as plt

PathLike = Union[str, Path]


def _loadmat(path: PathLike) -> Dict[str, Any]:
    return loadmat(str(path), simplify_cells=True)


def resolve_run_folder(run: PathLike) -> Path:
    """Resolve either a run folder or KW51ProbabilisticResults folder.

    If an output root is supplied, the lexically latest ``run_*`` folder is
    used. This mirrors the usual workflow where each inference creates a new
    timestamped run folder.
    """
    p = Path(run).expanduser().resolve()
    if not p.exists():
        raise FileNotFoundError(f"Path does not exist: {p}")

    if (p / "validation_results.mat").is_file() and (p / "input_catalog.mat").is_file():
        return p

    candidates = sorted(q for q in p.glob("run_*") if q.is_dir())
    if not candidates:
        raise FileNotFoundError(
            f"Could not find a KW51 run folder in {p}. Expected validation_results.mat "
            "or one or more run_* directories."
        )
    return candidates[-1]


def _accepted_event_ids(run_folder: Path) -> list[str]:
    catalog = _loadmat(run_folder / "input_catalog.mat")
    inputs = catalog.get("inputs", [])
    if isinstance(inputs, dict):
        inputs = [inputs]
    names = sorted(str(item["event_id"]) for item in inputs)
    if not names:
        raise RuntimeError("No accepted events found in input_catalog.mat")
    return names


def resolve_event_selector(run_folder: Path, event_selector: Union[int, str]) -> Tuple[str, int, int]:
    names = _accepted_event_ids(run_folder)
    n = len(names)

    if isinstance(event_selector, (int, np.integer)):
        idx = int(event_selector)
        if not 1 <= idx <= n:
            raise ValueError(f"Event number must be between 1 and {n}.")
        event_id = names[idx - 1]
        event_number = idx
    elif isinstance(event_selector, str):
        if event_selector not in names:
            raise ValueError(f"Event ID {event_selector!r} is not an accepted event in this run.")
        event_id = event_selector
        event_number = names.index(event_id) + 1
    else:
        raise TypeError("event_selector must be an integer event number or a full event ID string.")

    print(f"Selected event {event_number} of {n}: {event_id}")
    return event_id, event_number, n


def _object_vector_item(value: Any, index: int) -> np.ndarray:
    """Return one MATLAB cell entry as a 1-D float array."""
    if isinstance(value, np.ndarray) and value.dtype == object:
        return np.asarray(value.flat[index], dtype=float).reshape(-1)
    if isinstance(value, (list, tuple)):
        return np.asarray(value[index], dtype=float).reshape(-1)
    arr = np.asarray(value)
    if arr.ndim >= 2 and arr.shape[0] > index:
        return np.asarray(arr[index], dtype=float).reshape(-1)
    return np.asarray([], dtype=float)


def _score_lookup(
    validation_results: Dict[str, Any],
    event_id: str,
    target: str,
    sigma: float,
) -> Tuple[float, float]:
    scores = validation_results.get("scores")
    headers = validation_results.get("score_headers")
    if scores is None or headers is None:
        return np.nan, np.nan

    headers = [str(x) for x in np.asarray(headers).reshape(-1)]
    col = {name: i for i, name in enumerate(headers)}
    required = {"event", "target", "translation_sd_mm", "NRMSE", "coverage95"}
    if not required.issubset(col):
        return np.nan, np.nan

    scores = np.asarray(scores, dtype=object)
    for row in scores:
        if (
            str(row[col["event"]]) == event_id
            and str(row[col["target"]]) == target
            and abs(float(row[col["translation_sd_mm"]]) - sigma) < 1e-12
        ):
            return float(row[col["NRMSE"]]), float(row[col["coverage95"]])
    return np.nan, np.nan


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


def plot_kw51_selected_event(
    run: PathLike,
    event_selector: Union[int, str],
    *,
    visible: bool = True,
    output_folder: Optional[PathLike] = None,
    output_name: Optional[str] = None,
) -> Dict[str, Any]:
    """Plot withheld strains/displacements for one accepted passage.

    Parameters
    ----------
    run
        Either the timestamped ``run_*`` folder or its parent
        ``KW51ProbabilisticResults`` folder.
    event_selector
        Integer 1..N in chronological accepted-event order, or a full event ID.
    visible
        Leave the matplotlib window open after saving when True.
    output_folder
        Defaults to ``<run>/selected_event_figures``.
    output_name
        Optional filename stem, without extension.
    """
    _set_manuscript_font()
    run_folder = resolve_run_folder(run)
    event_id, event_number, event_count = resolve_event_selector(run_folder, event_selector)

    output_dir = Path(output_folder).expanduser().resolve() if output_folder else run_folder / "selected_event_figures"
    output_dir.mkdir(parents=True, exist_ok=True)

    settings_path = run_folder / "settings.mat"
    validation_path = run_folder / "validation" / f"{event_id}_validation.mat"
    results_path = run_folder / "validation_results.mat"
    if not settings_path.is_file():
        raise FileNotFoundError(f"Missing settings.mat in {run_folder}")
    if not validation_path.is_file():
        raise FileNotFoundError(f"No saved validation file for event {event_id}")

    c = _loadmat(settings_path)
    V = _loadmat(validation_path)
    R = _loadmat(results_path) if results_path.is_file() else {}

    targets = [str(x) for x in np.asarray(c["protected_names"]).reshape(-1)]
    t = np.asarray(V["time_s"], dtype=float).reshape(-1)
    sigmas = np.asarray(V["translation_sigma_values"], dtype=float).reshape(-1)
    default_sigma = float(c["translation_sd_mm"])
    idx_default = int(np.argmin(np.abs(sigmas - default_sigma)))
    if abs(sigmas[idx_default] - default_sigma) > 1e-12:
        raise RuntimeError("Default translation scale is absent from validation file.")

    filtered_mean = np.asarray(V["filtered_mean"], dtype=float)
    lower95 = np.asarray(V["lower95"], dtype=float)
    upper95 = np.asarray(V["upper95"], dtype=float)
    if lower95.ndim == 2:
        lower95 = lower95[:, :, None]
        upper95 = upper95[:, :, None]

    fig, axes = plt.subplots(2, 2, figsize=(11.7, 7.15), constrained_layout=True)
    axes = axes.ravel()

    for j, ax in enumerate(axes):
        ax.fill_between(
            t,
            lower95[:, j, idx_default],
            upper95[:, j, idx_default],
            alpha=0.28,
            label="Nominal 95% interval",
        )

        ax.plot(t, filtered_mean[:, j], linewidth=1.45, label="Posterior mean")

        tx = _object_vector_item(V.get("native_time", []), j)
        yy = _object_vector_item(V.get("native_truth", []), j)
        if tx.size and yy.size:
            if j <= 1:
                ax.plot(tx, yy, "k--", linewidth=1.0, label="Withheld measurement")
            else:
                ax.plot(tx, yy, "k.--", linewidth=0.9, markersize=4, label="Withheld measurement")

        kin = np.asarray(V["kinematics"], dtype=float).reshape(-1)
        rail_offsets = np.asarray(V["rail_offsets_s"], dtype=float).reshape(-1)
        start_time = kin[1] - 1.0
        finish_time = kin[1] + 116.316 / kin[0] + np.max(rail_offsets) + 2.0
        ax.set_xlim(start_time, finish_time)
        ax.grid(True, alpha=0.25)
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_title(f"({chr(ord('a') + j)}) {targets[j]}", fontweight="normal")
        ax.set_xlabel("Time from record start (s)")
        ax.set_ylabel("Strain (µε)" if j <= 1 else "Displacement (mm)")

        nr, coverage = _score_lookup(R, event_id, targets[j], default_sigma)
        if np.isfinite(nr):
            ax.text(
                0.02,
                0.97,
                f"NRMSE = {nr:.3f}; coverage = {100.0 * coverage:.1f}%",
                transform=ax.transAxes,
                va="top",
                fontsize=8.5,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.8, "pad": 2},
            )
        if j == 2:
            ax.legend(loc="best")

    fig.suptitle(
        f"Event {event_number}/{event_count}: {event_id} | Track {V['track']} | KW51 blueprint BC",
        fontsize=11,
        fontweight="normal",
    )

    stem = output_name or f"kw_selected_event_{event_number:02d}_{event_id}"
    pdf_path = output_dir / f"{stem}.pdf"
    png_path = output_dir / f"{stem}.png"
    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(png_path, dpi=220, bbox_inches="tight")

    if visible:
        plt.show(block=False)
    else:
        plt.close(fig)

    print(f"Selected-event figure written to:\n  {pdf_path}\n  {png_path}")
    return {"pdf": pdf_path, "png": png_path, "figure": fig, "event_id": event_id}


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", help="run_* folder or KW51ProbabilisticResults folder")
    parser.add_argument("event", help="1..N or full event ID")
    parser.add_argument("--hide", action="store_true", help="close figure after saving")
    parser.add_argument("--output-folder", default=None)
    parser.add_argument("--output-name", default=None)
    args = parser.parse_args()

    try:
        selector: Union[int, str] = int(args.event)
    except ValueError:
        selector = args.event

    plot_kw51_selected_event(
        args.run,
        selector,
        visible=not args.hide,
        output_folder=args.output_folder,
        output_name=args.output_name,
    )
