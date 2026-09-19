#!/usr/bin/env python3
"""Standalone shell+wind (Type IIn) model generator for OpenStella.

Builds the initial conditions for a supernova shell colliding with a
steady 1/r^2 circumstellar wind - no inner ejecta, no radioactive
nickel.  Writes every input file the run chain needs:

    <outdir>/modmake/<model>.hyd   structure (ASCII)
    <outdir>/modmake/<model>.abn   composition (ASCII)
    <outdir>/modmake/<model>.xni   Ni56 (Fortran record, zeroed)
    <outdir>/modmake/<model>.dat   STELLA run parameters
    <outdir>/eve/<model>.eve       eve2 parameter file
    <outdir>/eve/eve.1             eve2 index
    <outdir>/vladsf/ronfict.1      xronfict input
    <outdir>/strad/strad.1         STELLA index

Only numpy is required.  Physics matches the SN2024pba pipeline
(prepare_SN_B_model.py --no-ejecta) but has no pipeline dependencies.
"""

import argparse
import os
import struct

import numpy as np

MSUN = 1.98847e33
YR_TO_S = 3.15576e7

# Asplund+2009 solar, restricted to the six elements STELLA tracks
# nontrivially in .abn.  Sum with H/He gives Z = 0.0134.
SOLAR = {
    "fH": 0.7381,
    "fHe": 0.2485,
    "fO": 0.00573,
    "fSi": 0.00067,
    "fCa": 0.00006,
    "fFe": 0.00129,
}


def create_shell(mass_msun, v_kms, nzon_shell, t_start_days,
                 n_rho=8.0, v0_kms=3000.0, n_T=0.4,
                 T_interface_K=2000.0, T_max_K=65000.0):
    """Dessart 2010jl shell: flat core plus a steep power-law envelope.

    rho = rho0 for v < v0, rho = rho0*(v0/v)^n_rho beyond; rho0 is set
    by mass conservation.  r_outer = v_shell * t_start because the shell
    is already expanding homologously when STELLA takes over - starting
    at t=0 would give diffusion times the code cannot evolve.
    """
    v_shell = v_kms * 1e5
    v0 = v0_kms * 1e5
    r_outer = v_shell * t_start_days * 86400.0
    frac_v0 = v0_kms / v_kms
    r0 = frac_v0 * r_outer
    r_inner = max(1e8, r0 * 0.01)

    if n_rho > 3:
        factor = 1 + (3.0 / (n_rho - 3)) * (1 - frac_v0 ** (n_rho - 3))
    else:
        factor = 1 + (3.0 / (3 - n_rho)) * (frac_v0 ** (n_rho - 3) - 1)
    rho_0 = mass_msun * MSUN / ((4 / 3) * np.pi * r0 ** 3 * factor)
    rho_edge = rho_0 * frac_v0 ** n_rho

    r = np.logspace(np.log10(r_inner), np.log10(r_outer), nzon_shell)
    rho = np.where(r <= r0, rho_0, rho_0 * (r0 / r) ** n_rho)
    v = v_shell * r / r_outer
    T = np.clip(T_interface_K * (rho / rho_edge) ** n_T,
                T_interface_K, T_max_K)

    print(f"  shell: M={mass_msun} Msun v_outer={v_kms} km/s "
          f"r=[{r_inner:.3e},{r_outer:.3e}] cm")
    print(f"    rho0={rho_0:.3e} rho_edge={rho_edge:.3e} "
          f"contrast={rho_0 / rho_edge:.0f}")
    t_diff = 0.34 * rho_0 * r_outer ** 2 / 3e10
    print(f"    diffusion time ~{t_diff / 86400:.0f} d "
          f"(must be < ~10x target)")

    out = {"r": r, "rho": rho, "v": v, "T": T,
           "r_inner": r_inner, "r_outer": r_outer,
           "r0": r0, "rho_0": rho_0}
    for k, val in SOLAR.items():
        out[k] = np.full(nzon_shell, val)
    return out


def create_wind(mdot_msun_yr, v_wind_kms, T_wind_K,
                r_inner_cm, r_outer_cm, nzon_wind):
    """Steady 1/r^2 wind, constant velocity, isothermal."""
    mdot = mdot_msun_yr * MSUN / YR_TO_S
    v_wind = v_wind_kms * 1e5
    r = np.logspace(np.log10(r_inner_cm), np.log10(r_outer_cm), nzon_wind)
    rho = mdot / (4 * np.pi * r ** 2 * v_wind)
    m_wind = mdot * (r_outer_cm - r_inner_cm) / v_wind
    tau = 0.34 * mdot / (4 * np.pi * v_wind * r_inner_cm)
    print(f"  wind: mdot={mdot_msun_yr} Msun/yr v={v_wind_kms} km/s "
          f"M={m_wind / MSUN:.4f} Msun tau~{tau:.1f}")

    out = {"r": r, "rho": rho, "v": np.full(nzon_wind, v_wind),
           "T": np.full(nzon_wind, T_wind_K),
           "r_inner": r_inner_cm, "r_outer": r_outer_cm}
    for k, val in SOLAR.items():
        out[k] = np.full(nzon_wind, val)
    return out


def cap_zones(nzon, span_dex, min_dlogr, max_dlogr, label):
    """Keep dlog10(r) inside [min_dlogr, max_dlogr] for a segment.

    Too fine a spacing stalls the radiation diffusion step; too coarse
    a jump at a segment boundary makes the stiff Jacobian singular.
    """
    if span_dex <= 0 or nzon <= 0:
        return nzon
    dlogr = span_dex / nzon
    if dlogr < min_dlogr:
        nzon = max(int(span_dex / min_dlogr), 2)
        print(f"  {label}: capped to {nzon} zones (min dlogr)")
        dlogr = span_dex / nzon
    if max_dlogr is not None and dlogr > max_dlogr:
        nzon = max(int(np.ceil(span_dex / max_dlogr)), 2)
        print(f"  {label}: grown to {nzon} zones (max dlogr)")
    return nzon


def shell_grid(r_inner, r0, r_outer, nzon_shell, core_fraction=0.25,
               min_dlogr=0.005):
    """Split shell zones between the flat core and the steep envelope.

    A uniform log grid wastes zones in the core where nothing happens;
    the envelope carries the shock physics and needs the resolution.
    """
    if r0 <= r_inner or r0 >= r_outer:
        return np.logspace(np.log10(r_inner), np.log10(r_outer),
                           nzon_shell)
    nzon_core = max(10, int(nzon_shell * core_fraction))
    nzon_env = nzon_shell - nzon_core
    env_span = np.log10(r_outer) - np.log10(r0)
    max_env = max(int(env_span / min_dlogr), 10)
    if nzon_env > max_env:
        nzon_core += nzon_env - max_env
        nzon_env = max_env
    r_core = np.logspace(np.log10(r_inner), np.log10(r0),
                         nzon_core + 1)[:nzon_core]
    r_env = np.logspace(np.log10(r0), np.log10(r_outer),
                        nzon_env + 1)[1:]
    return np.concatenate([r_core, r_env])


def combine(shell, wind, nzon_total, sw_transition_dex=0.5,
            min_dlogr=0.005):
    """Join shell + buffer + wind into one grid of exactly nzon_total.

    The buffer is a log-linear density ramp; without it the 4-5 dex
    jump at the contact discontinuity breaks the solver.  Zones the
    wind cannot hold under min_dlogr are redistributed to the shell.
    """
    max_dlogr = 2 * min_dlogr
    log_r_min = np.log10(shell["r_inner"])
    log_r_shell = np.log10(shell["r_outer"])
    log_r_max = np.log10(wind["r_outer"])

    mdot = wind["rho"][0] * 4 * np.pi * wind["r_inner"] ** 2 * wind["v"][0]
    v_wind, T_wind = wind["v"][0], wind["T"][0]

    sw_transition_dex = min(sw_transition_dex,
                            (log_r_max - log_r_shell) / 3.0)
    r_buf_end = shell["r_outer"] * 10 ** sw_transition_dex
    rho_buf_start = shell["rho"][-1]
    rho_buf_end = mdot / (4 * np.pi * r_buf_end ** 2 * v_wind)
    contrast = np.log10(rho_buf_start / max(rho_buf_end, 1e-300))

    nzon_buf = max(20, int(contrast / 16 * sw_transition_dex / 0.02))
    nzon_buf = min(nzon_buf, 60)
    nzon_buf = cap_zones(nzon_buf, sw_transition_dex,
                         min_dlogr, max_dlogr, "SW buffer")

    nzon_shell = max(100, int(nzon_total * 0.40))
    nzon_wind = max(60, nzon_total - nzon_shell - nzon_buf)
    wind_span = log_r_max - np.log10(r_buf_end)
    nzon_wind = cap_zones(nzon_wind, wind_span,
                          min_dlogr, max_dlogr, "wind")

    freed = nzon_total - nzon_shell - nzon_buf - nzon_wind
    if freed < 0:
        excess = -freed
        take = min(excess, nzon_shell - 100)
        nzon_shell -= take
        excess -= take
        take = min(excess, nzon_wind - 60)
        nzon_wind -= take
        excess -= take
        freed = 0
        if excess > 0:
            print(f"  WARNING: {excess} zones unplaceable")
    if freed > 0:
        shell_span = log_r_shell - log_r_min
        room = max(int(shell_span / min_dlogr), nzon_shell) - nzon_shell
        give = min(freed, room)
        nzon_shell += give
        freed -= give
        print(f"  redistributed {give} freed zones to shell")
        if freed > 0:
            print(f"  WARNING: {freed} zones unplaceable")

    n = nzon_shell + nzon_buf + nzon_wind
    print(f"  zones: {nzon_shell} shell + {nzon_buf} buffer + "
          f"{nzon_wind} wind = {n}")

    # Shell values interpolated onto the split core/envelope grid.
    r_sh = shell_grid(shell["r_inner"], shell["r0"], shell["r_outer"],
                      nzon_shell, min_dlogr=min_dlogr)
    ls = np.log10(shell["r"])
    rho_sh = 10 ** np.interp(np.log10(r_sh), ls,
                             np.log10(np.maximum(shell["rho"], 1e-300)))
    v_sh = np.interp(np.log10(r_sh), ls, shell["v"])
    T_sh = np.interp(np.log10(r_sh), ls, shell["T"])

    # Renormalize the interpolated shell to conserve its mass exactly.
    n_rho = 8.0
    frac_v0 = shell["r0"] / shell["r_outer"]
    factor = 1 + (3 / (n_rho - 3)) * (1 - frac_v0 ** (n_rho - 3))
    m_target = (4 / 3) * np.pi * shell["r0"] ** 3 * shell["rho_0"] * factor
    r_lo = np.concatenate([[r_sh[0] * 0.5], r_sh[:-1]])
    m_interp = ((4 / 3) * np.pi * (r_sh ** 3 - r_lo ** 3) * rho_sh).sum()
    rho_sh *= m_target / m_interp
    print(f"  shell mass: target={m_target / MSUN:.5f} "
          f"renorm={m_target / m_interp:.4f}")

    r_buf = np.logspace(log_r_shell, np.log10(r_buf_end),
                        nzon_buf + 1)[1:]
    rho_buf = 10 ** np.linspace(np.log10(rho_buf_start),
                                np.log10(max(rho_buf_end, 1e-300)),
                                nzon_buf)
    v_buf = np.linspace(shell["v"][-1], v_wind, nzon_buf)
    T_buf = np.linspace(shell["T"][-1], T_wind, nzon_buf)

    r_w = np.logspace(np.log10(r_buf_end), log_r_max, nzon_wind + 1)[1:]
    rho_w = mdot / (4 * np.pi * r_w ** 2 * v_wind)

    out = {
        "r": np.concatenate([r_sh, r_buf, r_w]),
        "rho": np.concatenate([rho_sh, rho_buf, rho_w]),
        "v": np.concatenate([v_sh, v_buf,
                             np.full(nzon_wind, v_wind)]),
        "T": np.concatenate([T_sh, T_buf,
                             np.full(nzon_wind, T_wind)]),
    }

    # Smoothstep blend across the buffer: a hard metallicity step makes
    # a Rosseland opacity jump that can stall the diffusion timestep.
    u = np.linspace(0.0, 1.0, nzon_buf)
    w = u * u * (3.0 - 2.0 * u)
    for key in SOLAR:
        vs, vw = float(shell[key][0]), float(wind[key][0])
        out[key] = np.concatenate([
            np.full(nzon_shell, vs),
            vs * (1 - w) + vw * w,
            np.full(nzon_wind, vw),
        ])
    s = np.maximum(sum(out[k] for k in SOLAR), 1e-10)
    for key in SOLAR:
        out[key] /= s

    grad = np.abs(np.diff(np.log10(np.maximum(out["rho"], 1e-300)))
                  / np.diff(np.log10(out["r"])))
    print(f"  max |dlog rho/dlog r| = {grad.max():.1f} "
          f"(interior {grad[:-5].max():.1f})")
    if len(out["r"]) != nzon_total:
        print(f"  WARNING: {len(out['r'])} zones != nzon={nzon_total}")
    return out


def write_hyd(path, g, t_start_days, mass_cut_msun=0.0):
    """header: t_start nzon mass_cut rcen rho_cen;
    rows: km dMr[Msun] r[cm] rho T u Mr[Msun] 0"""
    r, rho = g["r"], g["rho"]
    nzon = len(r)
    r_lo = np.concatenate([[r[0] * 0.5], r[:-1]])
    dm = np.maximum((4 / 3) * np.pi * (r ** 3 - r_lo ** 3) * rho, 0.0)
    mr = np.cumsum(dm)
    with open(path, "w") as f:
        f.write(f"{t_start_days:12.3E}{nzon:6d}{mass_cut_msun:13.5E}"
                f"{r[0] * 0.5:13.5E}{rho[0]:13.5E}\n")
        for k in range(nzon):
            f.write(f" {k + 1:4d}{dm[k] / MSUN:12.4E}{r[k]:15.7E}"
                    f"{rho[k]:12.4E}{g['T'][k]:12.4E}{g['v'][k]:12.4E}"
                    f"{mr[k] / MSUN + mass_cut_msun:12.4E}{0.0:12.4E}\n")
    print(f"  wrote {path} ({nzon} zones, "
          f"M={mr[-1] / MSUN:.4f} Msun)")
    return mr[-1] / MSUN


def write_abn(path, g):
    """km 0 0 0 H He C N O Ne Na Mg Al Si S Ar Ca FePeak Ni58 Ni56"""
    with open(path, "w") as f:
        for k in range(len(g["r"])):
            f.write(f"{k + 1:4d}{0.0:10.3E}{0.0:10.3E}{0.0:10.3E}"
                    f"{g['fH'][k]:10.3E}{g['fHe'][k]:10.3E}"
                    f"{0.0:10.3E}{0.0:10.3E}{g['fO'][k]:10.3E}"
                    f"{0.0:10.3E}{0.0:10.3E}{0.0:10.3E}{0.0:10.3E}"
                    f"{g['fSi'][k]:10.3E}{0.0:10.3E}{0.0:10.3E}"
                    f"{g['fCa'][k]:10.3E}{g['fFe'][k]:10.3E}"
                    f"{0.0:10.3E}{0.0:10.3E}\n")
    print(f"  wrote {path}")


def write_xni(path, nzon):
    """Zeroed Ni56 record: 2 + 2*nzon REAL*4 in a Fortran record."""
    n = 2 + 2 * nzon
    floats = [0.0] * n
    floats[0] = 5.605194e-42
    floats[1] = -1.417409e-34
    data = struct.pack("f" * n, *floats)
    with open(path, "wb") as f:
        f.write(struct.pack("i", len(data)))
        f.write(data)
        f.write(struct.pack("i", len(data)))
    print(f"  wrote {path} (zero Ni56)")


def write_dat(path, target_days, r_max_cm, v_max_cm_s,
              epoch_interval_days=1.0):
    """STELLA run parameters; EKO=0 keeps the .hyd velocity field.

    TcurB is the observer-frame target plus a light-travel-time offset
    across the grid - STELLA stops on its internal clock.
    """
    ltt_days = r_max_cm / 3e10 / 86400.0
    tcurb = target_days + ltt_days
    epochs = [round(ltt_days + i * epoch_interval_days, 4)
              for i in range(1, int(target_days / epoch_interval_days) + 1)]
    with open(path, "w") as f:
        f.write(
            "   EPS   (STIFF ACCURACY)   HMIN (MIN.STEP(SEC)) "
            "HMAX (MAX.STEP(SEC))\n"
            "   0.003                    1e-14               3600\n"
            "   METH  ADAMS(1)-GEAR(2)-GEAR+BGH(3)   JSTART  "
            "(0-START_ELSE_ORDER) MAXORD KNadap\n"
            "    3                        0                    "
            "        4        -4\n"
            "   abs(NSTA)>abs(NSTB)=>reads last   TcurA     TcurB\n"
            f"   -9000     -1                      0.0000   {tcurb:.4f}\n"
            "   NSTMAX  NDebug  NOUT    IOUT (0 -full <0 print zones)"
            "   MBATCH\n"
            "   2000000   100000  100      10                   "
            "          900\n"
            "   AMNI (CONTAM. MASS)     XMNI (56NI MASS)\n"
            "   0D0                  0D0\n"
            "   AMHT(Solar)  EBurst(1e50) tBurst(s)  tstart -- "
            "Heated Core, Energy & time\n"
            "   0d0      0d0      1d0      0d0\n"
            "   EKO (kin.energy,1e50) mass fract.tri. u profile  "
            "us (out +1)\n"
            "   0d0                      0.d0  9.d-1  1.d0          +1.\n"
            "   THRMAT   CRAP     CONV     EDTM     CHNCND  Givdtl\n"
            "   1.D-30   1.D0     F        T        T       T\n"
            "   FLOOR(R) FLOOR(V) FLOOR(T) FLOOR(RADIAT)    "
            "FLOOR(UCONV)\n"
            "   1.D-12   1.D-03   0.03D+00   0.01D+00           1.D-04\n"
            "   Wacc(R)  Wacc(V)  Wacc(T)  Wacc(RADIAT)\n"
            "   1d-0    1d-0    0.5D-00   0.5D-00\n"
            "   FitTau   TAUTOL   Rvis   AQ    BQ       DRT     NRT"
            "     SCAT\n"
            "   5.00D0   1.30D0     1.D0   10.0D0  1.0D+00   1.    "
            "  1       T\n"
            "    NTO then    DATA TO\n"
            f"    {len(epochs)}\n"
        )
        for ep in epochs:
            f.write(f"    {ep:.4f}\n")
    print(f"  wrote {path} (TcurB={tcurb:.2f} d, "
          f"{len(epochs)} epochs)")


def write_eve(path, model, total_mass_msun):
    """eve2 parameter file; Kadapt=5 reads the structure from .hyd."""
    with open(path, "w") as f:
        f.write(
            " ANTROP             DENSCE          DENSCE-old\n"
            "  1.00               1.d-17              1.\n"
            " Rbeg    CRHO    RWIND   Tpwind   pwind\n"
            "1.d0     1.d14   -3.D3    2.5D3   0.d0\n"
            " Rce  -- Core Radius (Solar units)\n"
            "0.d0\n"
            " TPROC -- CENTRAL TEMPERATURE\n"
            "1.d-3\n"
            " VSTRE -- FITTING POINT\n"
            "0.55\n"
            " ALPHA\n"
            "0.3333333333333\n"
            " BMK       BM1   AMeveNi\n"
            f"{total_mass_msun:.4f}d0  0.d0   0.\n"
            " STEP\n"
            "-0.98\n"
            " expfac   velofac   rhofac   xnifac\n"
            "1.d0     1.0d0     1.d0     1.d0\n"
            " Kadapt   1-FPRho   2-Polytropes  3-Isothermal  "
            "4-foreign  5-homo  6-woo87a\n"
            "5\n"
            " US\n"
            "+1.\n"
            " TC  ISOTP. CENTRE   TS   SURFACE TEMPERATURE\n"
            "0.3                 0.3\n"
            " AMISO  GAMT\n"
            "1.D0   1.6666666667\n"
            " AS       ZS\n"
            "12.      6.\n"
            " XH       XHE        Depl\n"
            "0.7      0.296      0.0d0\n"
            " AVERG    IWOO\n"
            "28.       2\n"
            " AMZinn   AMZout   AMHein   AMHeout   BXY\n"
            "0.0      0.0      0.0      0.0       0.0001\n"
            " IPRINT\n"
            "2\n"
            " ITER1   ITER2\n"
            "100     100\n"
            " EPS1    EPSA    EPSD\n"
            "-14.0    -14.0   -10.0\n"
            " ULGCAP  ULGEPS\n"
            "0.0     12.0\n"
            " H\n"
            "0.00001\n"
            " DM1         DM2         PLATO        CENTR   DMID    "
            "MIDDLE\n"
            "1.D-02      1.D-02      0.4D0        0.4D0   1.D-2    "
            "4.550\n"
        )
    print(f"  wrote {path} (BMK={total_mass_msun:.4f} Msun)")


def main():
    ap = argparse.ArgumentParser(
        description="Shell+wind Type IIn initial conditions")
    ap.add_argument("--model-name", default="TypeIIn")
    ap.add_argument("--output-dir", default="TypeIIn_run")
    ap.add_argument("--shell-mass", type=float, default=7.0,
                    help="shell mass, Msun")
    ap.add_argument("--shell-velocity", type=float, default=10000.0,
                    help="outer shell velocity, km/s")
    ap.add_argument("--wind-mdot", type=float, default=0.01,
                    help="wind mass-loss rate, Msun/yr")
    ap.add_argument("--wind-velocity", type=float, default=100.0,
                    help="wind velocity, km/s")
    ap.add_argument("--wind-temperature", type=float, default=2000.0,
                    help="wind temperature, K")
    ap.add_argument("--wind-r-outer", type=float, default=1e16,
                    help="wind outer radius, cm")
    ap.add_argument("--nzon", type=int, default=500,
                    help="total zones (must match the build NZON)")
    ap.add_argument("--t-start", type=float, default=5.0,
                    help="start time since explosion, days")
    ap.add_argument("--target-days", type=float, default=60.0,
                    help="evolution target, days")
    args = ap.parse_args()

    modmake = os.path.join(args.output_dir, "modmake")
    eve = os.path.join(args.output_dir, "eve")
    vladsf = os.path.join(args.output_dir, "vladsf")
    strad = os.path.join(args.output_dir, "strad")
    for d in (modmake, eve, vladsf, strad):
        os.makedirs(d, exist_ok=True)

    m = args.model_name
    print(f"Building {m}: {args.shell_mass} Msun shell at "
          f"{args.shell_velocity} km/s into "
          f"{args.wind_mdot} Msun/yr wind")
    shell = create_shell(args.shell_mass, args.shell_velocity,
                         int(args.nzon * 0.6), args.t_start,
                         T_interface_K=args.wind_temperature)
    wind = create_wind(args.wind_mdot, args.wind_velocity,
                       args.wind_temperature,
                       shell["r_outer"] * 1.01, args.wind_r_outer,
                       int(args.nzon * 0.4))
    g = combine(shell, wind, args.nzon)

    total_mass = write_hyd(os.path.join(modmake, f"{m}.hyd"),
                           g, args.t_start)
    write_abn(os.path.join(modmake, f"{m}.abn"), g)
    write_xni(os.path.join(modmake, f"{m}.xni"), args.nzon)
    write_dat(os.path.join(modmake, f"{m}.dat"),
              args.target_days, g["r"].max(), g["v"].max())
    write_eve(os.path.join(eve, f"{m}.eve"), m, total_mass)

    with open(os.path.join(eve, "eve.1"), "w") as f:
        f.write("Data             Result          Model                "
                "Opacity            Rho               ExtMod           "
                "          CompMod                    Ni56 file\n"
                f"'{m}.eve'  'res.{m}'  '{m}.mod'  'DummyOpafile'  "
                f"'{m}.rho'  '../modmake/{m}.hyd'  "
                f"'../modmake/{m}.abn'  '{m}.xni'\n")

    with open(os.path.join(vladsf, "ronfict.1"), "w") as f:
        f.write("  Model                          NiDist               "
                "          results          opafile\n"
                f"'../eve/{m}.mod'  '../eve/{m}.xni'  'ronfict.res'  "
                f"'{m}.0'\n")

    with open(os.path.join(strad, "strad.1"), "w") as f:
        f.write("Run             Results             Model            "
                "Nickel             Opacity\n"
                f"{m}  {m}.res  {m}  {m}  {m}.1\n")

    print(f"\nDone.  Run with: bash run_model.sh {args.output_dir} "
          f"[threads]")


if __name__ == "__main__":
    main()
