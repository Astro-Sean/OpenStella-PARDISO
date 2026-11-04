#!/usr/bin/env python3
"""OpenStella-PARDISO validation plots.

Compares the PARDISO parallel build against the original Y12M serial
build on an identical STELLA model.  Reads the ``.swd`` structure
output (which also carries per-zone composition) and the ``.tt``
photometry from each run directory and produces:

* ``validation_profiles.png`` -- density, velocity and gas temperature
  radial profiles at three common epochs, all four runs overlaid.
* ``validation_lightcurve.png`` -- the complete 15-band light curve
  for the 1-CPU PARDISO run.
* ``validation_results.txt`` -- wall times, exit status and the
  maximum profile differences between each PARDISO run and Y12M.
* ``validation_runtime_vs_model_time.png`` -- wall-clock time since
  the STELLA solver start versus model simulation time for each run.

The .swd column layout (strad/lbalsw.trf, FORMAT 181/281):

    t[d] zone log10(Mext/Msun) log10(R/cm) v(1e8 cm/s) log10(Tgas)
    log10(Trad) log10(rho_int) log10(P) log10(Qvis) log10(E)
    L/1e40  kappa  n_bar(cm-3)  n_e(cm-3)  X(H..Ni)  (15 fractions)

Only the zone-1 and zone-Nzon rows of each epoch block carry the
timestamp; interior rows write 0 in the time column.
"""

import argparse
import glob
import os
import re
import subprocess
import sys
from datetime import datetime

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import AutoMinorLocator

AVOGADRO = 6.02214076e23

# Elements follow DataChames in strad/tt4strad.f.
ELEMENTS = [
    "H", "He", "C", "N", "O", "Ne", "Na", "Mg",
    "Al", "Si", "S", "Ar", "Ca", "Fe", "Ni",
]

RUN_STYLE = {
    "Y12M": {"color": "k", "lw": 2.0, "ls": "-"},
    "P1": {"color": "tab:blue", "lw": 1.2, "ls": "--"},
    "P2": {"color": "tab:orange", "lw": 1.2, "ls": "-."},
    "P3": {"color": "tab:red", "lw": 1.2, "ls": ":"},
}

RUN_LABEL = {
    "Y12M": "Y12M serial",
    "P1": "PARDISO 1 CPU",
    "P2": "PARDISO 2 CPU",
    "P3": "PARDISO 3 CPU",
}


def style_axes(ax, minor_x=True, minor_y=True, grid_axis="both"):
    """Apply the shared validation tick and grid style."""
    if minor_x and ax.get_xscale() == "linear":
        ax.xaxis.set_minor_locator(AutoMinorLocator())
    if minor_y and ax.get_yscale() == "linear":
        ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.grid(axis=grid_axis, alpha=0.3)
    ax.grid(axis=grid_axis, alpha=0.12, which="minor")


def parse_swd(path):
    """Return {epoch_time: (radius_cm, v_kms, rho, tgas, trad)} dict."""
    blocks = {}
    cur_t = None
    cur_rows = []
    for line in open(path):
        f = line.split()
        if len(f) != 30:
            continue
        zone = int(f[1])
        t = float(f[0])
        if zone == 1 and cur_rows:
            # Flush the previous block before starting a new one.
            blocks[cur_t] = np.array(cur_rows)
            cur_rows = []
        if t != 0.0:
            cur_t = t
        # cols: 3 logR, 4 v(1e8), 5 logTg, 6 logTrad, 13 n_bar
        cur_rows.append(
            (10.0 ** float(f[3]), float(f[4]) * 1e3,
             float(f[13]) / AVOGADRO, 10.0 ** float(f[5]),
             10.0 ** float(f[6]))
        )
    if cur_rows:
        blocks[cur_t] = np.array(cur_rows)
    return blocks


def parse_tt(path):
    """Return {band: (time, mag)} from a .tt lightcurve file.

    The file carries ~90 lines of run preamble (frequency grids,
    parameter echoes) before the photometry block; the real header is
    the last line starting with "time", followed by a units row.
    """
    with open(path) as fh:
        lines = [l for l in fh if l.strip()]
    header = None
    for i, l in enumerate(lines):
        if l.split()[0] == "time":
            header = i  # keep the last match; preamble may echo "time"
    if header is None:
        raise ValueError(f"no header row in {path}")
    cols = lines[header].split()
    rows = []
    for l in lines[header + 1:]:
        f = l.split()
        try:
            rows.append([float(x) for x in f])
        except ValueError:
            continue  # units row and any trailing text
    data = np.array([r for r in rows if len(r) == len(cols)])
    t = data[:, 0]
    out = {}
    for j, name in enumerate(cols[1:], start=1):
        out[name] = (t, data[:, j])
    return out


def common_epochs(blocks_by_run, n=3):
    """Pick n epochs present in every run (nearest-time match)."""
    ref = sorted(blocks_by_run["Y12M"])
    usable = []
    for t in ref:
        if all(
            min(abs(np.array(list(b)) - t)) < 0.05
            for b in blocks_by_run.values()
        ):
            usable.append(t)
    if len(usable) <= n:
        return usable
    idx = np.linspace(0, len(usable) - 1, n).astype(int)
    return [usable[i] for i in idx]


def nearest_block(blocks, t):
    return blocks[min(blocks, key=lambda k: abs(k - t))]


def plot_profiles(blocks_by_run, epochs, outpath):
    fig, axes = plt.subplots(
        3, len(epochs), figsize=(4.2 * len(epochs), 10),
        sharex=True, squeeze=False,
    )
    quantities = [
        ("density", 2, r"$\rho$ (g cm$^{-3}$)", "log"),
        ("velocity", 1, "v (km s$^{-1}$)", "log"),
        ("temperature", 3, r"$T_{\rm gas}$ (K)", "log"),
    ]
    # One x range and one y range per metric across all epochs, so the
    # panels show the structure evolving on fixed axes.
    xlim = [np.inf, -np.inf]
    ylims = [[np.inf, -np.inf] for _ in quantities]
    for t in epochs:
        for blocks in blocks_by_run.values():
            b = nearest_block(blocks, t)
            xlim[0] = min(xlim[0], b[:, 0].min())
            xlim[1] = max(xlim[1], b[:, 0].max())
            for row, (_, idx, _, _) in enumerate(quantities):
                v = b[:, idx]
                v = v[np.isfinite(v) & (v > 0.0)]
                if v.size:
                    ylims[row][0] = min(ylims[row][0], v.min())
                    ylims[row][1] = max(ylims[row][1], v.max())

    def pad_log(lim):
        lo, hi = np.log10(lim[0]), np.log10(lim[1])
        pad = 0.06 * (hi - lo)
        return [10.0 ** (lo - pad), 10.0 ** (hi + pad)]

    xlim = pad_log(xlim)
    ylims = [pad_log(lim) for lim in ylims]
    for col, t in enumerate(epochs):
        for row, (_, idx, label, scale) in enumerate(quantities):
            ax = axes[row][col]
            for run, blocks in blocks_by_run.items():
                b = nearest_block(blocks, t)
                style = RUN_STYLE[run]
                ax.plot(
                    b[:, 0], b[:, idx], label=RUN_LABEL[run], **style
                )
            ax.set_xscale("log")
            ax.set_yscale(scale)
            ax.set_xlim(xlim)
            ax.set_ylim(ylims[row])
            ax.set_xlabel("radius (cm)")
            if col == 0:
                ax.set_ylabel(label)
            if row == 0:
                ax.set_title(f"t = {t:.1f} d")
            if row == 0 and col == 0:
                ax.legend(fontsize=8, loc="lower left")
            style_axes(ax)
    fig.suptitle(
        "Type IIn validation: Y12M vs PARDISO structure agreement",
        y=0.995,
    )
    fig.tight_layout()
    fig.savefig(outpath, dpi=160)
    plt.close(fig)


def plot_lightcurve(tt, outpath):
    bands = [
        ("MU", "Johnson U"), ("MB", "Johnson B"), ("MV", "Johnson V"),
        ("MR", "Cousins R"), ("MI", "Cousins I"),
        ("Mu", "SDSS u"), ("Mg", "SDSS g"), ("Mr", "SDSS r"),
        ("Mi", "SDSS i"), ("Mz", "SDSS z"),
        ("MJ", "2MASS J"), ("MH", "2MASS H"), ("MKs", "2MASS Ks"),
        ("Mc", "ATLAS c"), ("Mo", "ATLAS o"),
    ]
    cmap = plt.cm.turbo(np.linspace(0.05, 0.95, len(bands)))
    fig, ax = plt.subplots(figsize=(9, 6))
    for (key, label), c in zip(bands, cmap):
        if key not in tt:
            continue
        t, m = tt[key]
        ax.plot(t, m, label=label, color=c, lw=1.3)
    ax.invert_yaxis()
    ax.set_xlabel("time since shock breakout (days)")
    ax.set_ylabel("absolute magnitude")
    ax.set_title(
        "Type IIn -- PARDISO 1 CPU: full 15-band STELLA light curve"
    )
    ax.legend(ncol=3, fontsize=8, loc="upper center",
              bbox_to_anchor=(0.5, -0.12))
    style_axes(ax)
    fig.tight_layout()
    fig.savefig(outpath, dpi=160, bbox_inches="tight")
    plt.close(fig)


def profile_metrics(blocks_by_run, epochs):
    """Max |dlog10| between each PARDISO run and Y12M, on a common
    radius grid, for rho, v and T at each epoch."""
    lines = []
    for t in epochs:
        ref = nearest_block(blocks_by_run["Y12M"], t)
        r_ref = ref[:, 0]
        lines.append(f"epoch t = {t:.2f} d")
        for run in ("P1", "P2", "P3"):
            if run not in blocks_by_run:
                continue
            b = nearest_block(blocks_by_run[run], t)
            lo = max(r_ref.min(), b[:, 0].min())
            hi = min(r_ref.max(), b[:, 0].max())
            grid = np.logspace(np.log10(lo), np.log10(hi), 800)
            res = []
            for idx, name in ((2, "rho"), (1, "v"), (3, "Tgas")):
                a = np.interp(
                    np.log10(grid), np.log10(r_ref), np.log10(ref[:, idx])
                )
                c = np.interp(
                    np.log10(grid), np.log10(b[:, 0]), np.log10(b[:, idx])
                )
                res.append(f"{name}: {np.max(np.abs(a - c)):.4f} dex")
            lines.append(f"  {run} vs Y12M  max|dlog10|  " + "  ".join(res))
    return lines


def job_elapsed(jobid):
    try:
        out = subprocess.run(
            [
                "sacct", "-j", str(jobid), "-n", "-X",
                "--format=JobID,JobName,Elapsed,ExitCode,State",
            ],
            capture_output=True, text=True, timeout=30,
        ).stdout.strip()
        return out or "  (sacct unavailable)"
    except Exception:
        return "  (sacct unavailable)"


def elapsed_minutes(jobid):
    """Wall time in minutes from sacct, or None if unavailable."""
    line = job_elapsed(jobid).splitlines()
    if not line:
        return None
    fields = line[0].split()
    if len(fields) < 3:
        return None
    try:
        days, sep, hms = fields[2].partition("-")
        if not sep:
            days, hms = "0", days
        h, m, s = (int(x) for x in hms.split(":"))
        return int(days) * 1440 + h * 60 + m + s / 60.0
    except (ValueError, IndexError):
        return None


def swd_epoch_times(path):
    """Return the zone-1 timestamp from every .swd block in order."""
    times = []
    with open(path, errors="ignore") as fh:
        for line in fh:
            f = line.split()
            if len(f) == 30 and int(f[1]) == 1:
                t = float(f[0])
                if t != 0.0:
                    times.append(t)
    return times


def clock_seconds(day, clock):
    """Convert REAL_CLOCK date and HHMMSS.sss fields to seconds."""
    value = float(clock)
    hour = int(value // 10000)
    minute = int((value - hour * 10000) // 100)
    second = value - hour * 10000 - minute * 100
    day_number = datetime.strptime(day, "%Y%m%d").toordinal()
    return day_number * 86400 + hour * 3600 + minute * 60 + second


def find_job_log(model_dir, name, jobid):
    """Locate the SLURM .out file for one validation run."""
    patterns = [
        os.path.join(model_dir, "stella", f"*{jobid}*.out"),
        os.path.join(model_dir, "stella", f"stella_{name}_*.out"),
        os.path.join(model_dir, "stella", "*.out"),
    ]
    for pattern in patterns:
        matches = sorted(glob.glob(pattern))
        if matches:
            return matches[0]
    return None


def parse_runtime_progress(log_path, swd_path):
    """Match wall-clock stamps to .swd epochs via the next ipout line.

    Each DATE_TIME/REAL_CLOCK stamp immediately precedes an ``ipout`` call
    when it marks model progress.  The ``ipout`` calls occur in the same
    order as the zone-1 rows in the .swd file, so their positions provide
    the model time.  Stamps after the last ``ipout`` are finalization work
    and use the last model epoch.
    """
    sim_times = swd_epoch_times(swd_path)
    clock_events = []
    ipout_lines = []
    day = None
    clock_re = re.compile(r"REAL_CLOCK\((\d)\)=\s*([0-9.+Ee-]+)")

    with open(log_path, errors="ignore") as fh:
        for lineno, line in enumerate(fh, start=1):
            match = clock_re.search(line)
            if match:
                if match.group(1) == "1":
                    day = match.group(2)
                elif match.group(1) == "2" and day:
                    clock_events.append((lineno, clock_seconds(day, match.group(2))))
                    day = None
            if "ipout fullgp" in line:
                ipout_lines.append(lineno)

    if len(sim_times) < 2 or len(clock_events) < 2 or not ipout_lines:
        return None

    pairs = []
    next_output = 0
    for lineno, wall_time in clock_events:
        while next_output < len(ipout_lines) and ipout_lines[next_output] < lineno:
            next_output += 1
        epoch_index = min(next_output, len(sim_times) - 1)
        pairs.append((sim_times[epoch_index], wall_time))

    start = pairs[0][1]
    sim_days = np.array([p[0] for p in pairs])
    elapsed = np.array([(p[1] - start) / 60.0 for p in pairs])
    return sim_days, elapsed


def plot_runtime_progress(progress_by_run, outpath):
    """Plot wall-clock runtime as a function of model simulation time."""
    order = [k for k in ("Y12M", "P1", "P2", "P3") if k in progress_by_run]
    fig, ax = plt.subplots(figsize=(8, 5))
    end_offsets = {
        "Y12M": (5, 0), "P1": (5, 8), "P2": (5, -8), "P3": (5, 0),
    }
    for run in order:
        sim_days, elapsed = progress_by_run[run]
        ax.plot(sim_days, elapsed, label=RUN_LABEL[run], **RUN_STYLE[run])
        ax.annotate(
            f"{elapsed[-1]:.0f} min",
            (sim_days[-1], elapsed[-1]),
            xytext=end_offsets[run], textcoords="offset points",
            color=RUN_STYLE[run]["color"], fontsize=8, va="center",
        )
    all_days = np.concatenate([progress_by_run[k][0] for k in order])
    all_elapsed = np.concatenate([progress_by_run[k][1] for k in order])
    day_span = all_days.max() - all_days.min()
    ax.set_xlim(
        all_days.min() - 0.02 * day_span,
        all_days.max() + 0.14 * day_span,
    )
    ax.set_ylim(0, all_elapsed.max() * 1.10)
    ax.set_xlabel("Simulation time until shock breakout (days)")
    ax.set_ylabel("wall time since STELLA start (min)")
    ax.set_title("Type IIn: model progress versus execution time")
    ax.legend(fontsize=9, loc="upper left")
    style_axes(ax)
    fig.tight_layout()
    fig.savefig(outpath, dpi=160)
    plt.close(fig)


def plot_runtime(minutes_by_run, outpath):
    """Wall-time bar chart with per-run speedup over the Y12M baseline."""
    order = [k for k in ("Y12M", "P1", "P2", "P3") if k in minutes_by_run]
    labels = [RUN_LABEL[k] for k in order]
    vals = [minutes_by_run[k] for k in order]
    base = minutes_by_run["Y12M"]
    colors = ["0.4"] + ["tab:blue", "tab:orange", "tab:red"][: len(order) - 1]
    fig, ax = plt.subplots(figsize=(7, 4.2))
    bars = ax.bar(labels, vals, color=colors)
    for bar, k, v in zip(bars, order, vals):
        tag = f"{v:.0f} min"
        if k != "Y12M":
            tag += f"\n({base / v:.2f}x)"
        ax.text(
            bar.get_x() + bar.get_width() / 2, v + base * 0.01, tag,
            ha="center", va="bottom", fontsize=9,
        )
    ax.set_ylabel("wall time (min)")
    ax.set_ylim(0, base * 1.18)
    ax.set_title("Type IIn (500 zones, 60 d): Y12M vs PARDISO wall time")
    style_axes(ax, minor_x=False, grid_axis="y")
    fig.tight_layout()
    fig.savefig(outpath, dpi=160)
    plt.close(fig)


def find_outputs(model_dir, name):
    """Locate the .swd and .tt a run produced."""
    swd = glob.glob(os.path.join(model_dir, "stella", "**", f"{name}.swd"),
                    recursive=True)
    swd += glob.glob(os.path.join(model_dir, "**", f"{name}.swd"),
                     recursive=True)
    tt = glob.glob(os.path.join(model_dir, "**", f"{name}.tt"),
                   recursive=True)
    swd = sorted(set(swd), key=len)
    tt = sorted(set(tt), key=len)
    return (swd[0] if swd else None, tt[0] if tt else None)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--runs", nargs=4, metavar=("Y12M_DIR", "P1_DIR", "P2_DIR", "P3_DIR"),
        default=[
            "TypeIIn_Y12M", "TypeIIn_P1", "TypeIIn_P2", "TypeIIn_P3",
        ],
        help="model directories, in Y12M/P1/P2/P3 order",
    )
    ap.add_argument("--run-root", default=".")
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--jobids", nargs=4, default=None,
                    help="SLURM job ids for the timing table")
    ap.add_argument("--profile-runs", default=None,
                    help="comma list of run keys for the profile "
                         "figure and metrics (e.g. Y12M,P1,P3). "
                         "Default: all runs with a .swd")
    args = ap.parse_args()

    keys = ["Y12M", "P1", "P2", "P3"]
    blocks_by_run = {}
    tt_by_run = {}
    swd_paths = {}
    model_dirs = {}
    for key, run in zip(keys, args.runs):
        mdir = run if os.path.isabs(run) else os.path.join(
            args.run_root, run
        )
        model_dirs[key] = mdir
        name = os.path.basename(mdir.rstrip("/"))
        swd, tt = find_outputs(mdir, name)
        swd_paths[key] = swd
        if swd:
            blocks_by_run[key] = parse_swd(swd)
            print(f"{key}: {len(blocks_by_run[key])} epochs from {swd}")
        else:
            print(f"{key}: no .swd found under {mdir}")
        if tt:
            tt_by_run[key] = parse_tt(tt)

    os.makedirs(args.outdir, exist_ok=True)

    if args.profile_runs:
        wanted = [k.strip() for k in args.profile_runs.split(",")]
        prof_runs = {k: b for k, b in blocks_by_run.items()
                     if k in wanted}
    else:
        prof_runs = blocks_by_run
    if "Y12M" in prof_runs and len(prof_runs) >= 2:
        epochs = common_epochs(prof_runs)
        print("epochs:", [f"{t:.2f}" for t in epochs])
        plot_profiles(
            prof_runs, epochs,
            os.path.join(args.outdir, "validation_profiles.png"),
        )
        metrics = profile_metrics(prof_runs, epochs)
    else:
        epochs = []
        metrics = ["  insufficient runs for comparison"]

    if "P1" in tt_by_run:
        plot_lightcurve(
            tt_by_run["P1"],
            os.path.join(args.outdir, "validation_lightcurve.png"),
        )

    report = ["Type IIn validation -- OpenStella Y12M vs PARDISO", ""]
    report.append("inputs (.swd):")
    for k in keys:
        report.append(f"  {k}: {swd_paths.get(k)}")
    if args.jobids:
        report.append("")
        report.append("SLURM jobs:")
        minutes = {}
        for k, j in zip(keys, args.jobids):
            report.append(f"  {k} (job {j}):")
            report.append("    " + job_elapsed(j).strip())
            m = elapsed_minutes(j)
            if m is not None:
                minutes[k] = m
        if "Y12M" in minutes and len(minutes) >= 2:
            plot_runtime(
                minutes,
                os.path.join(args.outdir, "validation_runtime.png"),
            )
        progress = {}
        for k, jobid in zip(keys, args.jobids):
            mdir = model_dirs[k]
            name = os.path.basename(mdir.rstrip("/"))
            log_path = find_job_log(mdir, name, jobid)
            swd_path = swd_paths.get(k)
            if not log_path or not swd_path:
                continue
            series = parse_runtime_progress(log_path, swd_path)
            if series is not None:
                progress[k] = series
        if len(progress) >= 2:
            plot_runtime_progress(
                progress,
                os.path.join(args.outdir, "validation_runtime_vs_model_time.png"),
            )
    report.append("")
    report.append("profile agreement (max abs difference in dex):")
    report += metrics
    text = "\n".join(report) + "\n"
    with open(os.path.join(args.outdir, "validation_results.txt"), "w") as fh:
        fh.write(text)
    print(text)


if __name__ == "__main__":
    sys.exit(main())
