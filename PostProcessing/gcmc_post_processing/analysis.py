# =============================================================================
# Reyhaneh — analysis.py
# Focus: ORIENTATIONAL ORDER (2D). Translational parts are disabled for now.
# Dependencies: numpy, freud (tested with freud >= 2.13).
# =============================================================================

from __future__ import annotations
import numpy as np
from dataclasses import dataclass

try:
    import freud
except Exception as e:
    raise ImportError("freud is required for this analysis module.") from e


# =============================================================================
# Small utilities
# =============================================================================

def _ensure_box(box_length):
    """Return (Lx, Ly) and a freud.box.Box."""
    if isinstance(box_length, (float, int)):
        Lx, Ly = float(box_length), float(box_length)
    else:
        Lx, Ly = float(box_length[0]), float(box_length[1])
    return (Lx, Ly), freud.box.Box(Lx=Lx, Ly=Ly)


def _to_3d(positions: np.ndarray) -> np.ndarray:
    """Ensure (N,3) array for freud."""
    if positions.ndim != 2:
        raise ValueError("positions must be (N,2) or (N,3).")
    if positions.shape[1] == 3:
        return positions.astype(np.float64, copy=False)
    if positions.shape[1] == 2:
        out = np.zeros((positions.shape[0], 3), dtype=np.float64)
        out[:, :2] = positions
        return out
    raise ValueError("positions must be (N,2) or (N,3).")


# =============================================================================
# RDF (kept minimal, used occasionally for diagnostics)
# =============================================================================

def compute_rdf(all_positions, box_length, r_max=None, dr=0.02):
    """
    Radial distribution function g(r) averaged over frames.

    Parameters
    ----------
    all_positions : list[np.ndarray]  # each (N,2 or 3)
    box_length : float or (2,)
    r_max : float or None
    dr : float

    Returns
    -------
    r_centers : (nbins,)
    g_avg     : (nbins,)
    """
    (Lx, Ly), box = _ensure_box(box_length)
    if r_max is None:
        r_max = 0.49 * min(Lx, Ly)
    nbins = int(max(10, np.floor(r_max / max(dr, 1e-6))))
    rdf = freud.density.RDF(bins=nbins, r_max=float(r_max))
    accum = []
    for pos in all_positions:
        pts = _to_3d(pos)
        rdf.compute((box, pts))
        accum.append(np.array(rdf.rdf, copy=True))
    g_avg = np.mean(accum, axis=0) if accum else np.zeros(nbins)
    edges = np.linspace(0.0, r_max, nbins + 1)
    r_centers = 0.5 * (edges[:-1] + edges[1:])
    return r_centers, g_avg


def average_rdf_over_trajectory(all_positions, box_length, dr=0.02, rcutoff=None):
    """
    Alias for compute_rdf with the same output signature used elsewhere.
    """
    r_max = rcutoff if rcutoff is not None else None
    return compute_rdf(all_positions, box_length, r_max=r_max, dr=dr)


# =============================================================================
# Voronoi neighbors
# =============================================================================

def voronoi_neighborhood_list(positions: np.ndarray, box_length):
    """
    Return freud Voronoi neighbor list for a single frame.
    """
    (Lx, Ly), box = _ensure_box(box_length)
    pts = _to_3d(positions)
    voro = freud.locality.Voronoi()
    voro.compute((box, pts))
    return voro.nlist


# =============================================================================
# Local density (simple 2D histogram; used for quick maps)
# =============================================================================

def compute_local_density(all_positions, box_length, nbins=(64, 64)):
    """
    2D number-density map averaged over frames using histogram2d.

    Returns
    -------
    x_edges, y_edges, rho_map  # rho_map shape = (ny, nx)
    """
    (Lx, Ly), _ = _ensure_box(box_length)
    nx, ny = int(nbins[0]), int(nbins[1])
    acc = np.zeros((ny, nx), dtype=np.float64)
    T = len(all_positions)
    for pos in all_positions:
        p = np.asarray(pos)
        x = np.mod(p[:, 0], Lx)
        y = np.mod(p[:, 1], Ly)
        H, xedges, yedges = np.histogram2d(y, x, bins=[ny, nx], range=[[0, Ly], [0, Lx]])
        acc += H
    rho = acc / max(T, 1) / ((Lx / nx) * (Ly / ny))
    return xedges, yedges, rho


# =============================================================================
# Hexatic order: Ψ6 series, Binder, χ6
# =============================================================================

def compute_psi6_series(all_positions, box_length, order_k=6):
    """
    Per-frame global hexatic order parameter Ψ_k (complex).

    Returns
    -------
    psi6_series     : (T,) complex
    psi6_abs_mean   : float
    psi6_abs_stderr : float
    """
    (Lx, Ly), box = _ensure_box(box_length)
    psi_list = []
    for pos in all_positions:
        pts = _to_3d(pos)
        voro = freud.locality.Voronoi()
        voro.compute((box, pts))
        op = freud.order.Hexatic(k=order_k)
        op.compute((box, pts), neighbors=voro.nlist)
        psi_list.append(np.mean(op.particle_order))
    psi6_series = np.asarray(psi_list, dtype=np.complex128)
    abs_vals = np.abs(psi6_series)
    n = len(abs_vals)
    mean = float(np.mean(abs_vals)) if n else 0.0
    stderr = float(np.std(abs_vals, ddof=1) / np.sqrt(n)) if n > 1 else 0.0
    return psi6_series, mean, stderr


def binder_and_susceptibility(psi6_series, N_particles):
    """
    Binder U4^(6) and susceptibility χ6 from Ψ6 series.
    """
    psi = np.asarray(psi6_series, dtype=np.complex128)
    m = np.mean(psi)
    m2 = np.mean(np.abs(psi) ** 2)
    m4 = np.mean(np.abs(psi) ** 4)
    U4 = float(1.0 - m4 / max(3.0 * m2 * m2, 1e-16))
    chi6 = float(N_particles * (m2 - np.abs(m) ** 2))
    return U4, chi6


# =============================================================================
# Bond-orientational correlation g6(r)
# =============================================================================

def compute_g6_avg(all_positions, box_length, r_max, nbins=120, return_per_frame=False):
    """
    Average g6(r) over frames using freud.density.CorrelationFunction on ψ6.

    Returns
    -------
    if return_per_frame=False:
        r, g6_avg
    else:
        r, g6_avg, g6_per_frame  # (T, nbins)
    """
    (Lx, Ly), box = _ensure_box(box_length)
    r_max = float(min(r_max, 0.49 * min(Lx, Ly)))
    nbins = int(nbins)
    cf = freud.density.CorrelationFunction(bins=nbins, r_max=r_max)
    per_frame = []
    for pos in all_positions:
        pts = _to_3d(pos)
        voro = freud.locality.Voronoi()
        voro.compute((box, pts))
        op = freud.order.Hexatic(k=6)
        op.compute((box, pts), neighbors=voro.nlist)
        vals = op.particle_order
        cf.compute(system=(box, pts), values=vals, query_points=pts, query_values=vals)
        per_frame.append(np.maximum(np.real(cf.correlation), 1e-16))
    per_frame = np.asarray(per_frame, dtype=np.float64) if len(per_frame) else np.zeros((0, nbins))
    g6_avg = np.mean(per_frame, axis=0) if per_frame.size else np.zeros(nbins)
    edges = np.linspace(0.0, r_max, nbins + 1)
    r = 0.5 * (edges[:-1] + edges[1:])
    if return_per_frame:
        return r, g6_avg, per_frame
    return r, g6_avg


# =============================================================================
# Simple eta6 fit (log–log)
# =============================================================================

def fit_eta6_from_g6(r, g6, rmin=None, rmax=None):
    """
    Fit g6(r) ~ r^{-eta6} on a specified window (log–log).
    """
    r = np.asarray(r, float)
    y = np.asarray(g6, float)
    m = (r > 0) & np.isfinite(y) & (y > 0)
    if rmin is not None:
        m &= (r >= float(rmin))
    if rmax is not None:
        m &= (r <= float(rmax))
    r_fit = r[m]; y_fit = y[m]
    if len(r_fit) < 6:
        return np.nan
    A = np.vstack([np.log(r_fit), np.ones_like(r_fit)]).T
    slope, _ = np.linalg.lstsq(A, np.log(y_fit), rcond=None)[0]
    return -float(slope)


# =============================================================================
# Automatic model selection (power vs exp) — simple
# =============================================================================

@dataclass
class G6Fit:
    kind: str       # 'power' | 'exp' | 'insufficient'
    eta6: float     # valid if 'power'
    r2_power: float
    r2_exp: float
    rmin: float
    rmax: float

def fit_g6_auto(r, g6, rmin=None, rmax=None, start_bin=10, r2_margin=0.02):
    """
    Simple selector between power law and exponential on a fixed window.
    """
    r = np.asarray(r, float)
    y = np.maximum(np.asarray(g6, float), 1e-16)
    if rmin is None or rmax is None:
        # default: drop first 'start_bin' points; use tail
        idx0 = min(start_bin, len(r))
        rr, yy = r[idx0:], y[idx0:]
    else:
        mask = (r >= rmin) & (r <= rmax)
        rr, yy = r[mask], y[mask]
    if len(rr) < 10:
        return G6Fit("insufficient", np.nan, 0.0, 0.0, np.nan, np.nan)

    # power: log y vs log r
    A = np.vstack([np.log(rr), np.ones_like(rr)]).T
    sp, ip = np.linalg.lstsq(A, np.log(yy), rcond=None)[0]
    yhat_p = np.exp(ip) * rr ** sp
    r2p = 1.0 - np.sum((yy - yhat_p) ** 2) / max(np.sum((yy - yy.mean()) ** 2), 1e-16)

    # exp: log y vs r
    B = np.vstack([rr, np.ones_like(rr)]).T
    se, ie = np.linalg.lstsq(B, np.log(yy), rcond=None)[0]
    yhat_e = np.exp(ie + se * rr)
    r2e = 1.0 - np.sum((np.log(yy) - np.log(yhat_e)) ** 2) / max(np.sum((np.log(yy) - np.log(yy).mean()) ** 2), 1e-16)

    if r2p >= r2e + r2_margin:
        return G6Fit("power", -float(sp), float(r2p), float(r2e), float(rr.min()), float(rr.max()))
    if r2e >= r2p + r2_margin:
        return G6Fit("exp", np.nan, float(r2p), float(r2e), float(rr.min()), float(rr.max()))
    return G6Fit("insufficient", np.nan, float(r2p), float(r2e), float(rr.min()), float(rr.max()))


# =============================================================================
# Robust g6 fitting: CONST vs POWER vs EXP with AIC + bootstrap
# =============================================================================

@dataclass
class G6ModelFit:
    model: str             # 'const' | 'power' | 'exp' | 'undecided'
    c_mean: float
    eta6_mean: float
    eta6_std: float
    aic_const: float
    aic_power: float
    aic_exp: float
    best_aic: float
    delta_to_next: float
    rmin: float
    rmax: float
    bins_used: int
    rel_var: float
    tail_mean: float
    notes: str

def _weighted_linfit(x, y, w):
    X = np.vstack([x, np.ones_like(x)]).T
    W = np.diag(w)
    XtW = X.T @ W
    beta = np.linalg.lstsq(XtW @ X, XtW @ y, rcond=None)[0]
    yhat = X @ beta
    rss = float(np.sum(w * (y - yhat) ** 2))
    return float(beta[0]), float(beta[1]), rss

def _weighted_const(y, w):
    c = float(np.sum(w * y) / np.sum(w))
    rss = float(np.sum(w * (y - c) ** 2))
    return c, rss

def _aic_from_rss(rss, n, k):
    rss = max(rss, 1e-30)
    n = max(n, 1)
    return float(n * np.log(rss / n) + 2 * k)

def fit_g6_models(r, g6_avg, g6_per_frame, start_bin=10, n_boot=200, rng_seed=0):
    """
    Fit CONST vs POWER vs EXP on a tail window and compare by AIC.
    - Weights from per-frame SE
    - Bootstrap η6 by resampling frames
    """
    r = np.asarray(r, float)
    y = np.maximum(np.asarray(g6_avg, float), 1e-16)
    idx0 = min(start_bin, len(r))
    rw, yw = r[idx0:], y[idx0:]
    n = len(rw)
    if n < 12:
        return G6ModelFit('undecided', np.nan, np.nan, np.nan, np.inf, np.inf, np.inf,
                          np.inf, 0.0, np.nan, np.nan, n, np.nan,
                          float(np.mean(y[-max(3, len(y)//5):])) if len(y) else np.nan,
                          "insufficient bins")

    # weights
    if g6_per_frame is not None and g6_per_frame.shape[0] >= 2:
        pf_tail = np.maximum(g6_per_frame[:, idx0:], 1e-16)  # (T, n)
        se = np.std(pf_tail, axis=0, ddof=1) / np.sqrt(g6_per_frame.shape[0])
        w = 1.0 / np.maximum(se, 1e-10) ** 2
    else:
        w = np.ones(n, float)

    rel_var = float((np.max(yw) - np.min(yw)) / max(np.mean(yw), 1e-16))
    tail_mean = float(np.mean(y[-max(3, len(y)//5):])) if len(y) else np.nan

    # const
    c_mean, rss_c = _weighted_const(yw, w)
    aic_c = _aic_from_rss(rss_c, n, k=1)

    # power
    sp, ip, rss_p = _weighted_linfit(np.log(rw), np.log(yw), w)
    eta_hat = -sp
    aic_p = _aic_from_rss(rss_p, n, k=2)

    # exp
    se, ie, rss_e = _weighted_linfit(rw, np.log(yw), w)
    aic_e = _aic_from_rss(rss_e, n, k=2)

    aics = np.array([aic_c, aic_p, aic_e])
    labels = np.array(["const", "power", "exp"])
    order = np.argsort(aics)
    best = labels[order[0]]
    best_aic = float(aics[order[0]])
    delta_to_next = float(aics[order[1]] - aics[order[0]])

    # bootstrap η6
    rng = np.random.default_rng(rng_seed)
    etas = []
    if g6_per_frame is not None and g6_per_frame.shape[0] >= 5:
        T = g6_per_frame.shape[0]
        for _ in range(int(n_boot)):
            pick = rng.integers(0, T, size=T)
            g_bs = np.maximum(g6_per_frame[pick][:, idx0:].mean(axis=0), 1e-16)
            sb, ib, _ = _weighted_linfit(np.log(rw), np.log(g_bs), np.ones_like(g_bs))
            etas.append(-sb)
    eta_mean = float(np.mean(etas)) if len(etas) else float(eta_hat)
    eta_std = float(np.std(etas, ddof=1)) if len(etas) > 1 else np.nan

    return G6ModelFit(best, float(c_mean), eta_mean, eta_std,
                      float(aic_c), float(aic_p), float(aic_e),
                      best_aic, delta_to_next,
                      float(rw.min()), float(rw.max()), int(n),
                      rel_var, tail_mean,
                      "const vs power vs exp; weighted; bootstrap eta6")


# =============================================================================
# Phase classifiers (orientational)
# =============================================================================

@dataclass
class PhaseResult:
    phase: str     # 'fluid' | 'hexatic' | 'solid' | 'transition'
    reason: str
    eta6: float
    psi6_abs_mean: float
    U4: float
    chi6: float
    extras: dict

def classify_orientational_phase(frames, box_length, n_particles,
                                 rmax_frac=0.45, nbins=120,
                                 eta_hex_upper=0.25, eta_solid_tol=0.05,
                                 tail_min=0.20, psi6_solid_min=0.50,
                                 r2_margin=0.02):
    """
    Simple classifier using fit_g6_auto + tail level + <|Ψ6|>.
    """
    (Lx, Ly), _ = _ensure_box(box_length)
    r, g6 = compute_g6_avg(frames, box_length, r_max=rmax_frac * min(Lx, Ly), nbins=nbins)
    psi6_series, psi6_abs_mean, _ = compute_psi6_series(frames, box_length, order_k=6)
    U4, chi6 = binder_and_susceptibility(psi6_series, n_particles)

    fit = fit_g6_auto(r, g6, start_bin=10, r2_margin=r2_margin)
    tail_mean = float(np.mean(g6[-max(3, len(g6)//5):])) if len(g6) else 0.0

    if fit.kind == "exp":
        return PhaseResult("fluid", "exp decay fits g6(r)", np.nan, psi6_abs_mean, U4, chi6,
                           {"tail_mean": tail_mean, "fit": fit})

    if fit.kind == "power":
        if (fit.eta6 <= eta_solid_tol) and (tail_mean >= tail_min) and (psi6_abs_mean >= psi6_solid_min):
            return PhaseResult("solid", "g6 ~ const; high tail; high |Ψ6|", fit.eta6, psi6_abs_mean, U4, chi6,
                               {"tail_mean": tail_mean, "fit": fit})
        if (fit.eta6 > eta_solid_tol) and (fit.eta6 < eta_hex_upper):
            return PhaseResult("hexatic", "algebraic with 0<η6<1/4", fit.eta6, psi6_abs_mean, U4, chi6,
                               {"tail_mean": tail_mean, "fit": fit})
        return PhaseResult("transition", "algebraic but near thresholds", fit.eta6, psi6_abs_mean, U4, chi6,
                           {"tail_mean": tail_mean, "fit": fit})

    return PhaseResult("transition", "insufficient or ambiguous", np.nan, psi6_abs_mean, U4, chi6,
                       {"tail_mean": tail_mean, "fit": fit})


@dataclass
class PhaseRules:
    phase: str
    reason: str
    eta6_mean: float
    eta6_std: float
    psi6_abs_mean: float
    U4: float
    chi6: float
    fit: G6ModelFit

def classify_orientational_phase_rules(frames, box_length, n_particles,
                                       rmax_frac=0.45, nbins=120,
                                       start_bin=25,
                                       n_boot=200, rng_seed=0,
                                       tail_liquid_cut=2e-1,
                                       const_rel_tol=0.01,
                                       aic_margin=4.0,
                                       eta_hex=0.25,
                                       eta_margin=0.02,
                                       psi6_solid_min=0.50):
    """
    Robust rule-based classifier (const/power/exp + thresholds).
    """
    (Lx, Ly), _ = _ensure_box(box_length)
    r, g6_avg, g6_pf = compute_g6_avg(frames, box_length, r_max=rmax_frac * min(Lx, Ly),
                                      nbins=nbins, return_per_frame=True)
    psi6_series, psi6_abs_mean, _ = compute_psi6_series(frames, box_length, order_k=6)
    U4, chi6 = binder_and_susceptibility(psi6_series, n_particles)

    fit = fit_g6_models(r, g6_avg, g6_pf, start_bin=start_bin, n_boot=n_boot, rng_seed=rng_seed)

    if fit.tail_mean < float(tail_liquid_cut):
        return PhaseRules("fluid", f"g6 tail < {tail_liquid_cut}", fit.eta6_mean, fit.eta6_std,
                          psi6_abs_mean, U4, chi6, fit)

    if (fit.rel_var <= float(const_rel_tol)) or (fit.model == "const" and fit.delta_to_next >= float(aic_margin)):
        if psi6_abs_mean >= float(psi6_solid_min):
            return PhaseRules("solid", "g6 nearly constant and |Psi6| high", fit.eta6_mean, fit.eta6_std,
                              psi6_abs_mean, U4, chi6, fit)
        return PhaseRules("transition", "g6 nearly constant but |Psi6| not high", fit.eta6_mean, fit.eta6_std,
                          psi6_abs_mean, U4, chi6, fit)

    if (fit.model == "exp") and (fit.delta_to_next >= float(aic_margin)):
        return PhaseRules("fluid", "exp better than power/const by AIC", fit.eta6_mean, fit.eta6_std,
                          psi6_abs_mean, U4, chi6, fit)

    if (fit.model == "power") and (fit.delta_to_next >= float(aic_margin)):
        low = fit.eta6_mean - (fit.eta6_std if np.isfinite(fit.eta6_std) else 0.0)
        high = fit.eta6_mean + (fit.eta6_std if np.isfinite(fit.eta6_std) else 0.0)
        if (abs(fit.eta6_mean - eta_hex) <= float(eta_margin)) or (low < eta_hex < high):
            return PhaseRules("transition", "η6 near 1/4 (hex–liquid boundary)", fit.eta6_mean, fit.eta6_std,
                              psi6_abs_mean, U4, chi6, fit)
        if fit.eta6_mean < eta_hex - float(eta_margin):
            return PhaseRules("hexatic", "algebraic with η6 < 1/4", fit.eta6_mean, fit.eta6_std,
                              psi6_abs_mean, U4, chi6, fit)
        if fit.eta6_mean > eta_hex + float(eta_margin):
            return PhaseRules("transition", "algebraic with η6 > 1/4 (liquid-side)", fit.eta6_mean, fit.eta6_std,
                              psi6_abs_mean, U4, chi6, fit)

    return PhaseRules("transition", "models too close or low confidence", fit.eta6_mean, fit.eta6_std,
                      psi6_abs_mean, U4, chi6, fit)


# =============================================================================
# Scans across state points (density or mu)
# =============================================================================

def scan_by_density_or_mu(dataset: dict, box_length, n_particles,
                          key="density", rmax_frac=0.45, nbins=120):
    """
    Compute <|Ψ6|>, U4, χ6 and g6(r) for each state point.

    Returns
    -------
    table : np.ndarray (structured with columns: key, psi6_mean_abs, psi6_stderr_abs, U4, chi6)
    g6_curves : dict[key -> (r, g6)]
    """
    (Lx, Ly), _ = _ensure_box(box_length)
    rows = []
    curves = {}
    for x in sorted(dataset.keys()):
        frames = dataset[x]
        psi6_series, psi_mean, psi_err = compute_psi6_series(frames, box_length, order_k=6)
        U4, chi6 = binder_and_susceptibility(psi6_series, n_particles)
        r, g6 = compute_g6_avg(frames, box_length, r_max=rmax_frac * min(Lx, Ly), nbins=nbins)
        rows.append((float(x), psi_mean, psi_err, float(U4), float(chi6)))
        curves[float(x)] = (r, g6)
    table = np.array(rows, dtype=[(key, 'f8'),
                                  ('psi6_mean_abs', 'f8'),
                                  ('psi6_stderr_abs', 'f8'),
                                  ('U4', 'f8'),
                                  ('chi6', 'f8')])
    return table, curves


def scan_and_classify_orientational(dataset, box_length, n_particles, key="density",
                                    rmax_frac=0.45, nbins=120, **kwargs):
    rows = []
    for x in sorted(dataset.keys()):
        pr = classify_orientational_phase(dataset[x], box_length, n_particles,
                                          rmax_frac=rmax_frac, nbins=nbins, **kwargs)
        rows.append({
            key: float(x),
            "phase": pr.phase,
            "reason": pr.reason,
            "eta6": pr.eta6,
            "psi6_mean_abs": pr.psi6_abs_mean,
            "U4": pr.U4,
            "chi6": pr.chi6,
        })
    # boundaries (adjacent changes)
    bounds = []
    rs = sorted(rows, key=lambda d: d[key])
    for i in range(len(rs) - 1):
        if rs[i]["phase"] != rs[i + 1]["phase"]:
            bounds.append({"between": (rs[i][key], rs[i + 1][key]),
                           "at": 0.5 * (rs[i][key] + rs[i + 1][key]),
                           "from": rs[i]["phase"], "to": rs[i + 1]["phase"]})
    return rs, bounds


def scan_and_classify_orientational_robust(dataset, box_length, n_particles, key="density",
                                           rmax_frac=0.45, nbins=120, **kwargs):
    rows = []
    for x in sorted(dataset.keys()):
        pr = classify_orientational_phase_rules(dataset[x], box_length, n_particles,
                                                rmax_frac=rmax_frac, nbins=nbins, **kwargs)
        rows.append({
            key: float(x),
            "phase": pr.phase,
            "reason": pr.reason,
            "eta6_mean": pr.eta6_mean,
            "eta6_std": pr.eta6_std,
            "psi6_mean_abs": pr.psi6_abs_mean,
            "U4": pr.U4,
            "chi6": pr.chi6,
            "fit_model": pr.fit.model,
            "delta_to_next": pr.fit.delta_to_next,
            "fit_rmin": pr.fit.rmin,
            "fit_rmax": pr.fit.rmax,
            "bins_used": pr.fit.bins_used,
            "rel_var": pr.fit.rel_var,
            "tail_mean_g6": pr.fit.tail_mean,
        })
    bounds = []
    rs = sorted(rows, key=lambda d: d[key])
    for i in range(len(rs) - 1):
        if rs[i]["phase"] != rs[i + 1]["phase"]:
            bounds.append({"between": (rs[i][key], rs[i + 1][key]),
                           "at": 0.5 * (rs[i][key] + rs[i + 1][key]),
                           "from": rs[i]["phase"], "to": rs[i + 1]["phase"]})
    return rs, bounds


# =============================================================================
# Legacy/compat stubs (translational order — disabled now)
# =============================================================================

def compute_gG(*args, **kwargs):
    """
    Positional correlation at reciprocal vector G.
    Disabled for now. We will reintroduce after orientational pipeline is stable.
    """
    raise NotImplementedError("Translational order (compute_gG) is disabled for now.")


def sub_system_translational(*args, **kwargs):
    """
    Subsystem translational analysis.
    Disabled for now. We will reintroduce after orientational pipeline is stable.
    """
    raise NotImplementedError("Translational analysis is disabled for now.")


# =============================================================================
# Compatibility shim: compute_psi and compute_psi_density (old names)
# =============================================================================

def compute_psi(all_positions, box_length, order_number=6):
    """
    Back-compat wrapper: returns <|ψ_k|> over frames and stderr.
    (ψ_k here refers to local particle-wise order; we report frame-averaged |ψ_k|.)
    """
    (Lx, Ly), box = _ensure_box(box_length)
    per_frame_abs_mean = []
    for pos in all_positions:
        pts = _to_3d(pos)
        voro = freud.locality.Voronoi()
        voro.compute((box, pts))
        op = freud.order.Hexatic(k=order_number)
        op.compute((box, pts), neighbors=voro.nlist)
        per_frame_abs_mean.append(np.mean(np.abs(op.particle_order)))
    arr = np.asarray(per_frame_abs_mean)
    n = len(arr)
    return float(np.mean(arr)), float(np.std(arr, ddof=1) / np.sqrt(n)) if n > 1 else 0.0


def compute_psi_density(all_positions, box_length, nbins=(64, 64)):
    """
    Produce 2D maps of:
      - number density ρ(x,y)
      - mean |ψ6|(x,y) (binned per pixel)
    Returns
    -------
    x_edges, y_edges, rho_map, psi_abs_map
    """
    (Lx, Ly), _ = _ensure_box(box_length)
    nx, ny = int(nbins[0]), int(nbins[1])
    rho_acc = np.zeros((ny, nx), dtype=np.float64)
    psi_acc = np.zeros((ny, nx), dtype=np.float64)
    cnt_acc = np.zeros((ny, nx), dtype=np.float64)

    for pos in all_positions:
        pts2 = np.asarray(pos)
        pts = _to_3d(pos)
        # ψ6 per particle
        _, box = _ensure_box(box_length)
        voro = freud.locality.Voronoi()
        voro.compute((box, pts))
        op = freud.order.Hexatic(k=6)
        op.compute((box, pts), neighbors=voro.nlist)
        psi_abs = np.abs(op.particle_order)

        # bin
        x = np.mod(pts2[:, 0], Lx); y = np.mod(pts2[:, 1], Ly)
        H, yedges, xedges = np.histogram2d(y, x, bins=[ny, nx], range=[[0, Ly], [0, Lx]])
        rho_acc += H

        # weighted sum of |ψ6|
        for i in range(pts2.shape[0]):
            xi = int(min(nx - 1, np.floor(x[i] / (Lx / nx))))
            yi = int(min(ny - 1, np.floor(y[i] / (Ly / ny))))
            psi_acc[yi, xi] += psi_abs[i]
            cnt_acc[yi, xi] += 1.0

    # number density per area
    rho_map = rho_acc / ((Lx / nx) * (Ly / ny)) / max(len(all_positions), 1)
    psi_map = np.divide(psi_acc, np.maximum(cnt_acc, 1.0), out=np.zeros_like(psi_acc), where=cnt_acc > 0)
    return xedges, yedges, rho_map, psi_map

# -----------------------------------------------------------------------------
# Back-compat: legacy compute_g6 API (for __init__ imports and old scripts)
# -----------------------------------------------------------------------------
def compute_g6(all_positions, box_length, r_max=None, nbins=120, return_per_frame=False):
    """
    Legacy wrapper around compute_g6_avg.

    Parameters
    ----------
    all_positions : list[np.ndarray]
    box_length : float or (2,)
    r_max : float or None
        If None, uses 0.45 * min(Lx, Ly).
    nbins : int
    return_per_frame : bool
        If True, also returns per-frame g6 curves (T x nbins).

    Returns
    -------
    r, g6_avg              (if return_per_frame=False)
    r, g6_avg, g6_per_frame (if return_per_frame=True)
    """
    # normalize box + default r_max
    if isinstance(box_length, (float, int)):
        Lx, Ly = float(box_length), float(box_length)
    else:
        Lx, Ly = float(box_length[0]), float(box_length[1])
    if r_max is None:
        r_max = 0.45 * min(Lx, Ly)

    return compute_g6_avg(
        all_positions,
        (Lx, Ly),
        r_max=float(r_max),
        nbins=int(nbins),
        return_per_frame=bool(return_per_frame),
    )
