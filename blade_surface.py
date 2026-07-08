#!/usr/bin/env python3
# -----------------------------------------------------------------------------
#  blade_surface.py — structured multiblock PLOT3D surface of a rotor/wing
#  blade skin (closed watertight "sock") for pyHyp hyperbolic extrusion.
#
#  Solver-agnostic: outputs plain formatted PLOT3D; the marched volume can be
#  used with any structured/unstructured CFD solver.
#
#  AXIS CONVENTION (blade frame)
#      x : chordwise, LE -> TE  (section wake direction, +x)
#      y : spanwise,  root -> TIP (+y)
#      z : thickness, suction side up (+z)
#      positive twist : leading edge UP (+z), i.e. nose-up incidence
#      twist is applied about the section LEADING EDGE.
#
#  PLANFORM: arbitrary, defined per spanwise section:
#      { "rR": 0.19,           # span station y/R
#        "chord": 0.1905,      # local chord [m]
#        "twistDeg": 8.0,      # nose-up twist about the LE [deg]
#        "xLE_c": 0.0,         # LE x-position  / local chord  (sweep)
#        "zLE_c": 0.0,         # LE z-position  / local chord  (dihedral/flap)
#        "airfoil": "naca0012" # 4-digit NACA or path to a Selig .dat file
#      }
#  Between stations chord/twist/xLE/zLE vary linearly; airfoil ordinates are
#  blended linearly at common chord fractions.
#
#  TOPOLOGY (pyHyp-compatible closed sock, no degenerate edges):
#      block 1 : main O-grid   (i = airfoil perimeter, j = span)
#                perimeter: TE_lower -> LE -> TE_upper -> blunt TE seal, the
#                first point repeated to close the O-loop
#      block 2 : tip  cap — Coons-TFI H-grid, its 4 edges are exact subsets of
#                the main end perimeter (LE nose arc / lower / TE seal / upper)
#      block 3 : root cap (same, mirrored)
#      pyHyp's directed-edge normal check rejects ANY collapsed cap edge, so
#      the nose arc (nte+2 points, nte odd) is essential — see README.
#
#  Chordwise spacing is a two-sided tanh (Vinokur): dTE must stay comparable
#  to the blunt-TE seal spacing or the cap corner columns shear during the
#  march; dLE must keep the cap nose arc at nose-radius scale.
#
#  Usage:  blade_surface.py <config.json> <out.fmt>
# -----------------------------------------------------------------------------
import sys, os, json, math
import numpy as np


# ---------------------------------------------------------------------------
# airfoil ordinates
# ---------------------------------------------------------------------------
def naca4(code, x):
    """Upper/lower ordinates of a 4-digit NACA at chord fractions x."""
    m = int(code[0])/100.0; p = int(code[1])/10.0; t = int(code[2:])/100.0
    yt = 5*t*(0.2969*np.sqrt(x) - 0.1260*x - 0.3516*x*x
              + 0.2843*x**3 - 0.1036*x**4)
    if p == 0:
        yc = np.zeros_like(x); dyc = np.zeros_like(x)
    else:
        yc = np.where(x < p, m/p**2*(2*p*x - x*x),
                      m/(1-p)**2*((1-2*p) + 2*p*x - x*x))
        dyc = np.where(x < p, 2*m/p**2*(p - x), 2*m/(1-p)**2*(p - x))
    th = np.arctan(dyc)
    return (x - yt*np.sin(th), yc + yt*np.cos(th),
            x + yt*np.sin(th), yc - yt*np.cos(th))


def load_dat(path, x):
    """Selig-format airfoil (.dat): resample upper/lower at chord fractions x."""
    pts = []
    for line in open(path):
        p = line.split()
        if len(p) == 2:
            try:
                pts.append((float(p[0]), float(p[1])))
            except ValueError:
                continue
    P = np.array(pts)
    P[:, 0] -= P[:, 0].min()
    P[:, 0] /= P[:, 0].max()                     # normalise chord to [0,1]
    ile = int(np.argmin(P[:, 0]))
    up = P[:ile+1][::-1]                         # LE -> TE
    lo = P[ile:]
    yu = np.interp(x, up[:, 0], up[:, 1])
    yl = np.interp(x, lo[:, 0], lo[:, 1])
    return x - 0*x, yu, x - 0*x, yl              # (xu,yu,xl,yl) on given x


def airfoil_ul(name, x, base=None):
    """(xu,yu,xl,yl) at chord fractions x for 'nacaXXXX' or a .dat path."""
    nm = str(name).strip().lower()
    if nm.startswith("naca") and len(nm) == 8 and nm[4:].isdigit():
        return naca4(nm[4:], x)
    path = name if os.path.isabs(name) else os.path.join(base or ".", name)
    return load_dat(path, x)


# ---------------------------------------------------------------------------
# point distributions
# ---------------------------------------------------------------------------
def tanh_dist(n, s0, s1):
    """Vinokur two-sided tanh distribution on [0,1]; s0,s1 = end spacings
    (normalised).  Pass None for either end to get one-sided/uniform."""
    if s0 is None and s1 is None:
        return np.linspace(0.0, 1.0, n)
    if s0 is None:
        return 1.0 - tanh_dist(n, s1, None)[::-1]
    if s1 is None:
        # one-sided: geometric-ish via tanh with free far end
        s1 = 2.0/(n - 1)
    A = math.sqrt(s1/s0)
    B = 1.0/((n - 1)*math.sqrt(s0*s1))
    de = 1.0
    if B > 1.0:
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
        if len(secs) < 2:
            secs = [dict(secs[0]), dict(secs[0])]
            secs[0]["rR"], secs[1]["rR"] = secs[0]["rR"], 1.0
        self.rR = np.array([float(s["rR"]) for s in secs])
        self.chord = np.array([float(s["chord"]) for s in secs])
        self.twist = np.array([float(s.get("twistDeg", 0.0)) for s in secs])
        self.xle = np.array([float(s.get("xLE_c", 0.0)) for s in secs])
        self.zle = np.array([float(s.get("zLE_c", 0.0)) for s in secs])
        self.airfoils = [s.get("airfoil", "naca0012") for s in secs]
        self.base = base

    def at(self, rR, xc):
        """Interpolated (chord, twistDeg, xLE_m, zLE_m, yu, yl) at station."""
        rR = float(np.clip(rR, self.rR[0], self.rR[-1]))
        i = int(np.searchsorted(self.rR, rR, side="right") - 1)
        i = min(max(i, 0), len(self.rR) - 2)
        w = ((rR - self.rR[i])/(self.rR[i+1] - self.rR[i])
             if self.rR[i+1] > self.rR[i] else 0.0)
        chord = (1-w)*self.chord[i] + w*self.chord[i+1]
        twist = (1-w)*self.twist[i] + w*self.twist[i+1]
        xle = ((1-w)*self.xle[i] + w*self.xle[i+1])*chord
        zle = ((1-w)*self.zle[i] + w*self.zle[i+1])*chord
        _, yu0, _, yl0 = airfoil_ul(self.airfoils[i], xc, self.base)
        if self.airfoils[i+1] == self.airfoils[i] or w == 0.0:
            yu, yl = yu0, yl0
        else:
            _, yu1, _, yl1 = airfoil_ul(self.airfoils[i+1], xc, self.base)
            yu = (1-w)*yu0 + w*yu1
            yl = (1-w)*yl0 + w*yl1
        return chord, twist, xle, zle, yu, yl


# ---------------------------------------------------------------------------
# perimeter (closed O-loop) and section placement
# ---------------------------------------------------------------------------
def chord_stations(cfg):
    s = cfg.get("surface", {})
    nchord = int(s.get("nChord", 160))
    nte = int(s.get("nTE", 7))
    if (nte + 2) % 2 == 0:
        raise SystemExit("nTE must be ODD (cap nose-arc symmetry)")
    tecut = float(s.get("teCut", 0.97))
    k = (nte + 1)//2
    dle = s.get("dLE_c", None)
    dte = s.get("dTE_c", 0.003)
    if dle is None:
        dle = 2.5e-3/max(1, k)      # keep the cap nose arc at nose-radius scale
    xc = tecut*tanh_dist(nchord, float(dle)/tecut, float(dte)/tecut)
    return xc, nchord, nte, tecut


def perimeter_2d(xc, yu, yl, nte):
    """Closed O-loop in the section plane (x = chord dir, z = thickness):
    TE_lower -> LE -> TE_upper -> blunt seal (nte interior pts) -> close.
    This winding + span j toward +y makes the pyHyp march normals point OUT."""
    X = np.concatenate([xc[::-1], xc[1:]])
    Z = np.concatenate([yl[::-1], yu[1:]])
    xa = np.linspace(X[-1], X[0], nte + 2)[1:-1]   # seal: upper TE -> lower TE
    za = np.linspace(Z[-1], Z[0], nte + 2)[1:-1]
    X = np.concatenate([X, xa]); Z = np.concatenate([Z, za])
    X = np.concatenate([X, X[:1]]); Z = np.concatenate([Z, Z[:1]])
    return X, Z


def build(cfg, out, base="."):
    pf = Planform(cfg, base)
    s = cfg.get("surface", {})
    xc, nchord, nte, tecut = chord_stations(cfg)
    R = pf.R
    rR0 = float(s.get("rootCut", pf.rR[0]))
    nspan = int(s.get("nSpan", 60))
    span = tanh_dist(nspan, s.get("dRootFrac", None), s.get("dTipFrac", None))
    stations = (rR0 + (pf.rR[-1] - rR0)*span)*R

    nper = 2*nchord - 1 + nte + 1
    Xg = np.zeros((nper, nspan)); Yg = np.zeros((nper, nspan)); Zg = np.zeros((nper, nspan))
    sec2d = []
    for j, yst in enumerate(stations):
        chord, twist, xle, zle, yu, yl = pf.at(yst/R, xc)
        Px, Pz = perimeter_2d(xc, yu, yl, nte)
        # twist about the LE, nose-up positive (rotate x->z by +twist at LE)
        th = math.radians(twist); ct, st = math.cos(th), math.sin(th)
        xr = Px*chord; zr = Pz*chord
        Xs = xle + xr*ct - zr*st*0 - 0*zr  # placeholder, replaced below
        # rotation about LE (origin of section coords): nose-up = LE stays,
        # TE goes DOWN for positive twist:  x' = x cos t + z sin t is nose-down;
        # nose-UP:  x' =  x cos(t) + z sin(t)?  Derive: rotate section by +t
        # about +y axis maps (x,z) -> (x cos t + z sin t, -x sin t + z cos t),
        # which moves TE (+x) toward -z: TE down = LE up = nose-up.  Correct.
        Xs = xle + (xr*ct + zr*st)
        Zs = zle + (-xr*st + zr*ct)
        Xg[:, j] = Xs; Yg[:, j] = yst; Zg[:, j] = Zs
        sec2d.append((Px.copy(), Pz.copy(), chord, th, xle, zle))

    blocks = [(Xg, Yg, Zg)]

    if bool(s.get("closedSock", True)):
        blocks.append(cap_tfi(sec2d[-1], nchord, nte, stations[-1], True))
        blocks.append(cap_tfi(sec2d[0], nchord, nte, stations[0], False))


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
                        " + 2 TFI caps (closed sock)" if len(blocks) > 1 else
                        " (open sock)"))
    return len(blocks)


def cap_tfi(sec, nchord, nte, ystation, isTip):
    """Non-degenerate Coons-TFI cap in the section plane; boundary points are
    exact copies of the main end perimeter (pointReduce-stitchable)."""
    Px, Pz, chord, th, xle, zle = sec
    nth = nte + 2
    k = (nth - 1)//2
    LE = nchord - 1
    upTE = 2*nchord - 2
    arc = list(range(LE - k, LE + k + 1))                  # lower -> upper
    lower = list(range(LE - k, -1, -1))                    # arc end -> TE_lo
    upper = list(range(LE + k, upTE + 1))                  # arc end -> TE_up
    seal = [0] + list(range(upTE + nte, upTE, -1)) + [upTE]  # lo -> up
    il = len(lower)
    x2 = np.zeros((il, nth)); z2 = np.zeros((il, nth))
    for i in range(il):
        u = i/(il - 1)
        for jt in range(nth):
            v = jt/(nth - 1)
            x2[i, jt] = ((1-v)*Px[lower[i]] + v*Px[upper[i]]
                         + (1-u)*Px[arc[jt]] + u*Px[seal[jt]]
                         - ((1-u)*(1-v)*Px[arc[0]] + (1-u)*v*Px[arc[-1]]
                            + u*(1-v)*Px[seal[0]] + u*v*Px[seal[-1]]))
            z2[i, jt] = ((1-v)*Pz[lower[i]] + v*Pz[upper[i]]
                         + (1-u)*Pz[arc[jt]] + u*Pz[seal[jt]]
                         - ((1-u)*(1-v)*Pz[arc[0]] + (1-u)*v*Pz[arc[-1]]
                            + u*(1-v)*Pz[seal[0]] + u*v*Pz[seal[-1]]))
    x2[:, 0] = Px[lower]; z2[:, 0] = Pz[lower]
    x2[:, -1] = Px[upper]; z2[:, -1] = Pz[upper]
    x2[0, :] = Px[arc];   z2[0, :] = Pz[arc]
    x2[-1, :] = Px[seal]; z2[-1, :] = Pz[seal]
    ct, st = math.cos(th), math.sin(th)
    xr = x2*chord; zr = z2*chord
    Xh = xle + (xr*ct + zr*st)
    Zh = zle + (-xr*st + zr*ct)
    Yh = np.full_like(Xh, ystation)
    if isTip:                       # outward normal +y at the tip: reverse i
        Xh = Xh[::-1].copy(); Yh = Yh[::-1].copy(); Zh = Zh[::-1].copy()
    return Xh, Yh, Zh


if __name__ == "__main__":
    cfgFile = sys.argv[1]
    cfg = json.load(open(cfgFile))
    out = sys.argv[2] if len(sys.argv) > 2 else "bladeSurf.fmt"
    build(cfg, out, base=os.path.dirname(os.path.abspath(cfgFile)))
