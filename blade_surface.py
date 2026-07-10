#!/usr/bin/env python3
# -----------------------------------------------------------------------------
#  Author: JOO Seunghyun <chlrh45351@gmail.com>  (2026-07)
#  blade_surface.py — structured multiblock PLOT3D surface of a rotor/wing
#  blade skin (closed and watertight) for pyHyp hyperbolic extrusion.
#
#  Solver-agnostic: outputs plain formatted PLOT3D; the marched volume can be
#  used with any structured/unstructured CFD solver.
#
#  AXIS CONVENTION (rotor frame)
#      x : rotor axis (rotor wake direction +x); section thickness axis,
#          suction side toward -x
#      y : spanwise, root -> TIP (+y)
#      z : chordwise, LE -> TE (the leading edge FACES -z)
#      twist : applied about the local QUARTER CHORD (0.25 c),
#              positive = nose-up (LE rotates toward -x)
#
#  INPUT: a single plain-text .dat file — keyword lines, then a SECTIONS
#  table (see examples/*.dat):
#
#      R          1.143          # tip radius [m]
#      nChord     200            # chordwise points per side
#      ...
#      SECTIONS
#      # r/R   chord[m]  twist[deg]  LE_z[m]  airfoil
#      0.19    0.1905    8.0         0.0      naca0012
#      1.00    0.1905    8.0         0.0      naca0012
#
#  Per-section planform: chord, twist and the fore/aft LE position LE_z
#  (in METERS, along +z) vary linearly between stations; airfoil
#  ordinates ("nacaXXXX" or a Selig .dat path, relative to the input file)
#  are blended linearly.
#
#  TOPOLOGY (pyHyp-compatible watertight skin, no degenerate edges):
#      block 1 : main O-grid   (i = airfoil perimeter, j = span)
#      block 2 : tip  cap — Coons-TFI H-grid whose 4 edges are exact subsets
#                of the main end perimeter (LE nose arc / lower / TE seal /
#                upper); nTE must be ODD (nose-arc symmetry)
#      block 3 : root cap
#  The outward orientation of the whole quilt is checked numerically and all
#  blocks are flipped together if needed, so the march always goes outward.
#
#  Usage:  blade_surface.py <input.dat> <out.fmt>
# -----------------------------------------------------------------------------
import sys, os, math
import numpy as np


# ---------------------------------------------------------------------------
# input file
# ---------------------------------------------------------------------------
SURF_KEYS = {"nChord", "nTE", "teCut", "dLE_c", "dTE_c",
             "nSpan", "dRootFrac", "dTipFrac", "closedTips", "closedSock",
             "rootCut",
             "datSmooth", "capDome"}
MARCH_KEYS = {"firstLayer", "nLayers", "marchDist", "splay", "volSmoothIter",
              "volBlend", "volCoef", "cMax", "epsE", "epsI", "theta",
              "nConstantStart", "matchFactor", "autoMatch", "maxRatio"}
BG_KEYS = {"bgSpacing", "bgGrowth", "bgQuarter", "bgCyl", "bgXmin", "bgXmax", "bgYmin", "bgYmax",
           "bgZmin", "bgZmax", "refXmin", "refXmax", "refYmin", "refYmax",
           "refZmin", "refZmax"}


def _num(tok):
    try:
        return int(tok)
    except ValueError:
        try:
            return float(tok)
        except ValueError:
            return tok


def read_input(path):
    """Parse the .dat input -> config dict (planform/surface/march)."""
    cfg = {"planform": {"sections": []}, "surface": {}, "march": {},
           "background": {}}
    in_sections = False
    for raw in open(path):
        line = raw.split("#")[0].strip()
        if not line:
            continue
        if line.upper() == "SECTIONS":
            in_sections = True
            continue
        p = line.split()
        if not in_sections:
            if len(p) < 2:
                continue
            k, v = p[0], _num(p[1])
            if k in SURF_KEYS:
                cfg["surface"][k] = v
            elif k in MARCH_KEYS:
                cfg["march"][k] = v
            elif k in BG_KEYS:
                cfg["background"][k] = v
            elif k == "R":
                cfg["planform"]["R"] = float(v)
            elif k == "nBlades":
                cfg["planform"]["nBlades"] = int(v)
            else:
                sys.stderr.write("[blade_surface] WARNING: unknown key '%s'\n" % k)
        else:
            if len(p) < 5:
                raise SystemExit("SECTIONS row needs: r/R chord twist LE_z airfoil")
            cfg["planform"]["sections"].append({
                "rR": float(p[0]), "chord": float(p[1]),
                "twistDeg": float(p[2]), "LE_z": float(p[3]),
                "airfoil": p[4]})
    if "R" not in cfg["planform"]:
        raise SystemExit("input must define R")
    if len(cfg["planform"]["sections"]) < 1:
        raise SystemExit("input must contain a SECTIONS table")
    return cfg


# ---------------------------------------------------------------------------
# airfoil ordinates
# ---------------------------------------------------------------------------
def airfoil_ul(name, x, base=None, smooth=5):
    """(yu, yl) ordinates at chord fractions x.  For tabulated (.dat) data,
    'smooth' endpoint-preserving Laplacian passes are applied to the raw
    ordinates first: digitised coordinates carry point-to-point noise that a
    cubic spline reproduces faithfully, and the hyperbolic march inverts in
    the resulting micro-concavities (seen with UIUC sc1095.dat)."""
    nm = str(name).strip().lower()
    if nm.startswith("naca") and len(nm) == 8 and nm[4:].isdigit():
        code = nm[4:]
        m = int(code[0])/100.0; p = int(code[1])/10.0; t = int(code[2:])/100.0
        yt = 5*t*(0.2969*np.sqrt(x) - 0.1260*x - 0.3516*x*x
                  + 0.2843*x**3 - 0.1036*x**4)
        if p == 0:
            yc = np.zeros_like(x)
        else:
            yc = np.where(x < p, m/p**2*(2*p*x - x*x),
                          m/(1-p)**2*((1-2*p) + 2*p*x - x*x))
        return yc + yt, yc - yt
    path = name if os.path.isabs(name) else os.path.join(base or ".", name)
    pts = []
    for line in open(path):
        p2 = line.split()
        if len(p2) == 2:
            try:
                pts.append((float(p2[0]), float(p2[1])))
            except ValueError:
                continue
    P = np.array(pts)
    # Format auto-detection (UIUC database mixes the two):
    #   Lednicer: first numeric line is the POINT COUNTS (values > 1.1),
    #             then upper LE->TE, then lower LE->TE
    #   Selig   : one continuous loop TE -> LE -> TE
    if P[0, 0] > 1.1 and P[0, 1] > 1.1:
        nup = int(round(P[0, 0])); nlo = int(round(P[0, 1]))
        up = P[1:1+nup]
        lo = P[1+nup:1+nup+nlo]
    else:
        xmin = P[:, 0].min()
        ile = int(np.argmin(P[:, 0]))
        up = P[:ile+1][::-1]
        lo = P[ile:]
    for S in (up, lo):
        S[:, 0] -= S[:, 0].min()
    scale = max(up[:, 0].max(), lo[:, 0].max())
    up = up/scale; lo = lo/scale
    for S in (up, lo):
        for _ in range(int(smooth)):
            S[1:-1, 1] += 0.25*(S[:-2, 1] + S[2:, 1] - 2*S[1:-1, 1])
    # cubic-spline resampling in u = sqrt(x): tabulated data is far coarser
    # than the LE point spacing used here, and LINEAR interpolation leaves a
    # faceted nose that collapses the hyperbolic march (y ~ sqrt(x) at the
    # nose becomes linear in u, so the spline stays smooth and bounded)
    return (_cubic(np.sqrt(np.abs(up[:, 0])), up[:, 1], np.sqrt(x)),
            _cubic(np.sqrt(np.abs(lo[:, 0])), lo[:, 1], np.sqrt(x)))


def _cubic(xd, yd, xq):
    """Natural cubic spline through (xd, yd), evaluated at xq (numpy only).
    xd must be strictly increasing; duplicates are collapsed."""
    xd = np.asarray(xd, float); yd = np.asarray(yd, float)
    keep = np.concatenate([[True], np.diff(xd) > 1e-12])
    xd = xd[keep]; yd = yd[keep]
    n = len(xd)
    if n < 3:
        return np.interp(xq, xd, yd)
    h = np.diff(xd)
    # solve for second derivatives M (natural BCs) via Thomas algorithm
    a = h[:-1].copy(); b = 2.0*(h[:-1] + h[1:]); c = h[1:].copy()
    d = 6.0*((yd[2:] - yd[1:-1])/h[1:] - (yd[1:-1] - yd[:-2])/h[:-1])
    for i in range(1, n-2):
        w = a[i]/b[i-1]
        b[i] -= w*c[i-1]
        d[i] -= w*d[i-1]
    M = np.zeros(n)
    M[n-2] = d[-1]/b[-1]
    for i in range(n-3, 0, -1):
        M[i] = (d[i-1] - c[i-1]*M[i+1])/b[i-1]
    j = np.clip(np.searchsorted(xd, xq) - 1, 0, n-2)
    t = xq - xd[j]
    return (yd[j] + t*((yd[j+1] - yd[j])/h[j] - h[j]*(2*M[j] + M[j+1])/6.0)
            + t*t*M[j]/2.0 + t*t*t*(M[j+1] - M[j])/(6.0*h[j]))


# ---------------------------------------------------------------------------
# point distributions
# ---------------------------------------------------------------------------
def tanh_dist(n, s0, s1):
    """Vinokur two-sided tanh distribution on [0,1]; s0/s1 = end spacings
    (normalised); None => free end / uniform."""
    if s0 is None and s1 is None:
        return np.linspace(0.0, 1.0, n)
    if s0 is None:
        return 1.0 - tanh_dist(n, s1, None)[::-1]
    if s1 is None:
        s1 = 2.0/(n - 1)
    A = math.sqrt(s1/s0)
    B = 1.0/((n - 1)*math.sqrt(s0*s1))
    if B > 1.0:
        de = 1.0
        for _ in range(200):
            f = math.sinh(de)/de - B
            fp = (math.cosh(de)*de - math.sinh(de))/de**2
            step = f/fp
            de -= step
            if abs(step) < 1e-12:
                break
        t = np.linspace(0.0, 1.0, n)
        u = 0.5*(1.0 + np.tanh(de*(t - 0.5))/math.tanh(de/2))
    else:
        u = np.linspace(0.0, 1.0, n)
    x = u/(A + (1 - A)*u)
    x = np.clip(x, 0.0, 1.0); x[0] = 0.0; x[-1] = 1.0
    return x


# ---------------------------------------------------------------------------
# planform interpolation
# ---------------------------------------------------------------------------
class Planform:
    def __init__(self, cfg, base):
        pf = cfg["planform"]
        self.R = float(pf["R"])
        secs = sorted(pf["sections"], key=lambda s: float(s["rR"]))
        if len(secs) == 1:
            secs = [dict(secs[0]), dict(secs[0])]
            secs[1]["rR"] = 1.0
        self.rR = np.array([float(s["rR"]) for s in secs])
        self.chord = np.array([float(s["chord"]) for s in secs])
        self.twist = np.array([float(s["twistDeg"]) for s in secs])
        self.lez = np.array([float(s["LE_z"]) for s in secs])
        self.airfoils = [s["airfoil"] for s in secs]
        self.base = base
        self.smooth = int(cfg.get("surface", {}).get("datSmooth", 5))

    def at(self, rR, xc):
        rR = float(np.clip(rR, self.rR[0], self.rR[-1]))
        i = int(np.searchsorted(self.rR, rR, side="right") - 1)
        i = min(max(i, 0), len(self.rR) - 2)
        w = ((rR - self.rR[i])/(self.rR[i+1] - self.rR[i])
             if self.rR[i+1] > self.rR[i] else 0.0)
        chord = (1-w)*self.chord[i] + w*self.chord[i+1]
        twist = (1-w)*self.twist[i] + w*self.twist[i+1]
        lez = (1-w)*self.lez[i] + w*self.lez[i+1]   # [m]
        yu0, yl0 = airfoil_ul(self.airfoils[i], xc, self.base, self.smooth)
        if self.airfoils[i+1] == self.airfoils[i] or w == 0.0:
            yu, yl = yu0, yl0
        else:
            yu1, yl1 = airfoil_ul(self.airfoils[i+1], xc, self.base,
                                  self.smooth)
            yu = (1-w)*yu0 + w*yu1
            yl = (1-w)*yl0 + w*yl1
        return chord, twist, lez, yu, yl


# ---------------------------------------------------------------------------
# perimeter (closed O-loop) in (chord-fraction, thickness) space
# ---------------------------------------------------------------------------
def chord_stations(surf, rle=0.0159):
    nchord = int(surf.get("nChord", 160))
    nte = int(surf.get("nTE", 7))
    if (nte + 2) % 2 == 0:
        raise SystemExit("nTE must be ODD (cap nose-arc symmetry)")
    tecut = float(surf.get("teCut", 0.97))
    k = (nte + 1)//2
    dle = surf.get("dLE_c", None)
    dte = surf.get("dTE_c", 0.003)
    if dle is None:
        # The cap nose arc spans k points each side of the LE; the march folds
        # there unless the arc stays well inside the nose radius.  Scale the
        # LE spacing with the (estimated) nose radius so arc/r_LE ~ 0.16, the
        # ratio validated on NACA0012 (r_LE = 0.0159 c -> dLE ~ 6.2e-4 c).
        dle = min(0.16*rle, 2.5e-3)/max(1, k)
        sys.stderr.write("[blade_surface] r_LE ~ %.4g c -> dLE_c = %.3g\n"
                         % (rle, dle))
    xc = tecut*tanh_dist(nchord, float(dle)/tecut, float(dte)/tecut)
    return xc, nchord, nte, tecut


def perimeter_2d(xc, yu, yl, nte):
    """Closed O-loop: TE_lower -> LE -> TE_upper -> blunt seal -> close.
    Returns (chord-fraction, thickness-fraction) along the loop."""
    C = np.concatenate([xc[::-1], xc[1:]])
    T = np.concatenate([yl[::-1], yu[1:]])
    ca = np.linspace(C[-1], C[0], nte + 2)[1:-1]
    ta = np.linspace(T[-1], T[0], nte + 2)[1:-1]
    C = np.concatenate([C, ca, C[:1]])
    T = np.concatenate([T, ta, T[:1]])
    return C, T


def place(C, T, chord, twistDeg, zle):
    """(chord-frac, thickness) -> (x, z) in the rotor frame:
    z = chordwise (LE at -z side, LE_z offset), x = -thickness (suction -x),
    twist about the local quarter chord, nose-up positive."""
    z0 = zle + C*chord
    x0 = -T*chord
    zp = zle + 0.25*chord
    th = math.radians(twistDeg)
    ct, st = math.cos(th), math.sin(th)
    dx = x0; dz = z0 - zp
    X = dx*ct + dz*st
    Z = -dx*st + dz*ct + zp
    return X, Z


def nose_radius(pf):
    """Smallest leading-edge radius (in chords) over the section airfoils,
    estimated from the thickness at x0:  t ~ 2 sqrt(2 r x)  ->  r = t^2/(8 x)."""
    x0 = 0.002
    x = np.array([x0])
    r = []
    for a in set(pf.airfoils):
        yu, yl = airfoil_ul(a, x, pf.base, pf.smooth)
        r.append(float(yu[0] - yl[0])**2/(8.0*x0))
    return max(1e-4, min(r))


def build(cfg, out, base="."):
    pf = Planform(cfg, base)
    s = cfg.get("surface", {})
    xc, nchord, nte, tecut = chord_stations(s, nose_radius(pf))
    R = pf.R
    rR0 = float(s.get("rootCut", pf.rR[0]))
    nspan = int(s.get("nSpan", 60))
    span = tanh_dist(nspan, s.get("dRootFrac", None), s.get("dTipFrac", None))
    stations = (rR0 + (pf.rR[-1] - rR0)*span)*R

    nper = 2*nchord - 1 + nte + 1
    Xg = np.zeros((nper, nspan)); Yg = np.zeros((nper, nspan)); Zg = np.zeros((nper, nspan))
    sec2d = []
    for j, yst in enumerate(stations):
        chord, twist, zle, yu, yl = pf.at(yst/R, xc)
        C, T = perimeter_2d(xc, yu, yl, nte)
        Xs, Zs = place(C, T, chord, twist, zle)
        Xg[:, j] = Xs; Yg[:, j] = yst; Zg[:, j] = Zs
        sec2d.append((C.copy(), T.copy(), chord, twist, zle))

    # datum: put the ROOT quarter-chord at z = 0 (the twist pivot already
    # sits at x = 0), so the blade pitch/quarter-chord axis passes through
    # the rotor axis; LE_z entries remain RELATIVE planform shifts.
    chord0, _, zle0, _, _ = pf.at(stations[0]/R, xc)
    zshift = zle0 + 0.25*chord0
    Zg -= zshift
    for k in range(len(sec2d)):
        C, T, chord, twist, zle = sec2d[k]
        sec2d[k] = (C, T, chord, twist, zle - zshift)
    if abs(zshift) > 0:
        sys.stderr.write("[blade_surface] z datum: root 0.25c -> z=0 "
                         "(shift %.6g m)\n" % (-zshift))

    blocks = [(Xg, Yg, Zg)]
    if int(s.get("closedTips", s.get("closedSock", 1))):
        dome = float(s.get("capDome", 0.0))
        blocks.append(cap_tfi(sec2d[-1], nchord, nte, stations[-1], True, dome))
        blocks.append(cap_tfi(sec2d[0], nchord, nte, stations[0], False, dome))

    # ---- outward-orientation self check (pyHyp marches along t_i x t_j) ----
    P = np.stack([Xg, Yg, Zg], -1)
    ti = P[1:, :-1] - P[:-1, :-1]
    tj = P[:-1, 1:] - P[:-1, :-1]
    nrm = np.cross(ti, tj)
    ctr = P[:-1, :-1].mean(axis=0, keepdims=True)
    outw = np.einsum('ijk,ijk->ij', nrm, P[:-1, :-1] - ctr)
    if np.median(outw) < 0:
        blocks = [(Bx[::-1].copy(), By[::-1].copy(), Bz[::-1].copy())
                  for (Bx, By, Bz) in blocks]
        sys.stderr.write("[blade_surface] flipped all blocks -> outward march\n")

    with open(out, "w") as f:
        f.write("%d\n" % len(blocks))
        for (A, _, _) in blocks:
            f.write("%d %d %d\n" % (A.shape[0], A.shape[1], 1))
        for (Bx, By, Bz) in blocks:
            ni, nj = Bx.shape
            for A in (Bx, By, Bz):
                for jj in range(nj):
                    for i in range(ni):
                        f.write("%.10g\n" % A[i, jj])
    sys.stderr.write("[blade_surface] wrote %s : main %dx%d%s\n"
                     % (out, nper, nspan,
                        " + 2 TFI caps (watertight)" if len(blocks) > 1 else
                        " (open ends)"))
    return len(blocks)


def cap_tfi(sec, nchord, nte, ystation, isTip, capDome=0.0):
    """Non-degenerate Coons-TFI cap; its boundary points are exact copies of
    the main end perimeter (autoConnect/pointReduce stitchable).
    capDome > 0 bulges the cap INTERIOR outboard into a slightly rounded tip
    (local half-thickness dome; 1.0 ~ a semicircular cross-section).  The
    boundary points are untouched, so the block stitching and the
    directed-edge normal check are unaffected."""
    C, T, chord, twist, zle = sec
    nth = nte + 2
    k = (nth - 1)//2
    LE = nchord - 1
    upTE = 2*nchord - 2
    arc = list(range(LE - k, LE + k + 1))                     # lower -> upper
    lower = list(range(LE - k, -1, -1))                       # arc end -> TE_lo
    upper = list(range(LE + k, upTE + 1))                     # arc end -> TE_up
    seal = [0] + list(range(upTE + nte, upTE, -1)) + [upTE]   # lo -> up
    il = len(lower)
    c2 = np.zeros((il, nth)); t2 = np.zeros((il, nth))
    for i in range(il):
        u = i/(il - 1)
        for jt in range(nth):
            v = jt/(nth - 1)
            c2[i, jt] = ((1-v)*C[lower[i]] + v*C[upper[i]]
                         + (1-u)*C[arc[jt]] + u*C[seal[jt]]
                         - ((1-u)*(1-v)*C[arc[0]] + (1-u)*v*C[arc[-1]]
                            + u*(1-v)*C[seal[0]] + u*v*C[seal[-1]]))
            t2[i, jt] = ((1-v)*T[lower[i]] + v*T[upper[i]]
                         + (1-u)*T[arc[jt]] + u*T[seal[jt]]
                         - ((1-u)*(1-v)*T[arc[0]] + (1-u)*v*T[arc[-1]]
                            + u*(1-v)*T[seal[0]] + u*v*T[seal[-1]]))
    c2[:, 0] = C[lower]; t2[:, 0] = T[lower]
    c2[:, -1] = C[upper]; t2[:, -1] = T[upper]
    c2[0, :] = C[arc];   t2[0, :] = T[arc]
    c2[-1, :] = C[seal]; t2[-1, :] = T[seal]
    Xh, Zh = place(c2, t2, chord, twist, zle)
    Yh = np.full_like(Xh, ystation)
    if capDome > 0.0:
        u = (np.arange(il)/(il - 1.0))[:, None]          # 0 = LE arc, 1 = seal
        w = (2.0*np.arange(nth)/(nth - 1.0) - 1.0)[None, :]   # -1 lower, +1 up
        halft = 0.5*(T[upper] - T[lower])[:, None]*chord      # local half-thick
        dome = capDome*halft*np.sin(np.pi*u)*np.cos(0.5*np.pi*w)
        Yh = Yh + (dome if isTip else -dome)
    if isTip:
        Xh = Xh[::-1].copy(); Yh = Yh[::-1].copy(); Zh = Zh[::-1].copy()
    return Xh, Yh, Zh


if __name__ == "__main__":
    inp = sys.argv[1]
    cfg = read_input(inp)
    out = sys.argv[2] if len(sys.argv) > 2 else "bladeSurf.fmt"
    build(cfg, out, base=os.path.dirname(os.path.abspath(inp)))
