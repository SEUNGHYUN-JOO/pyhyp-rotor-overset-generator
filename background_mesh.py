#!/usr/bin/env python3
# -----------------------------------------------------------------------------
#  background_mesh.py — structured Cartesian BACKGROUND mesh for overset
#  (chimera) rotor setups, from the same rotor .dat input.
#
#  A single-block PLOT3D brick with a user-defined REFINEMENT BOX around the
#  rotor: uniform spacing inside the box (a multiple of the TIP CHORD,
#  default 0.15 c_tip), geometric growth outward to the domain boundary.
#
#  Axis convention matches the blade tools: x = rotor axis (wake +x),
#  the rotor disk lies in the y-z plane at x ~ 0.
#
#  Keywords (all optional; box extents in units of R):
#      bgSpacing   0.15        # refine-box spacing, multiples of the tip chord
#      bgGrowth    1.12        # geometric growth ratio outside the box
#      bgXmin -4   bgXmax 8    # domain (R units; +x = downstream/wake)
#      bgYmin -4   bgYmax 4
#      bgZmin -4   bgZmax 4
#      refXmin -0.5  refXmax 2.0    # refinement box (R units)
#      refYmin -1.2  refYmax 1.2
#      refZmin -1.2  refZmax 1.2
#
#  Usage:  background_mesh.py <rotor.dat> <out.xyz>
# -----------------------------------------------------------------------------
import sys, os
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from blade_surface import read_input, Planform


def stretch_axis(x0, x1, a, b, h0, ratio):
    """1D distribution on [x0, x1]: uniform spacing h0 inside [a, b],
    geometrically growing spacing (factor 'ratio') outward to both ends.
    The outer segments are rescaled so the end points land exactly."""
    a = max(a, x0); b = min(b, x1)
    # uniform core (exact fit)
    ncore = max(1, int(round((b - a)/h0)))
    core = np.linspace(a, b, ncore + 1)
    hc = (b - a)/ncore

    def outward(from_x, to_x, sign):
        steps = []
        h = hc
        pos = from_x
        while (to_x - pos)*sign > 1e-12:
            h *= ratio
            steps.append(h)
            pos += sign*h
        if not steps:
            return np.array([])
        # rescale all steps so the last point lands exactly on to_x
        steps = np.array(steps)
        steps *= abs(to_x - from_x)/steps.sum()
        return from_x + sign*np.cumsum(steps)

    lo = outward(a, x0, -1.0)[::-1] if a > x0 + 1e-12 else np.array([])
    hi = outward(b, x1, +1.0) if b < x1 - 1e-12 else np.array([])
    return np.concatenate([lo, core, hi])


def build(cfg, out, base="."):
    pf = Planform(cfg, base)
    R = pf.R
    ctip = float(pf.chord[-1])
    bg = cfg.get("background", {})
    h0 = float(bg.get("bgSpacing", 0.15))*ctip
    ratio = float(bg.get("bgGrowth", 1.12))

    dom = [float(bg.get(k, d))*R for k, d in
           (("bgXmin", -4), ("bgXmax", 8), ("bgYmin", -4), ("bgYmax", 4),
            ("bgZmin", -4), ("bgZmax", 4))]
    ref = [float(bg.get(k, d))*R for k, d in
           (("refXmin", -0.5), ("refXmax", 2.0), ("refYmin", -1.2),
            ("refYmax", 1.2), ("refZmin", -1.2), ("refZmax", 1.2))]

    xs = stretch_axis(dom[0], dom[1], ref[0], ref[1], h0, ratio)
    ys = stretch_axis(dom[2], dom[3], ref[2], ref[3], h0, ratio)
    zs = stretch_axis(dom[4], dom[5], ref[4], ref[5], h0, ratio)
    ni, nj, nk = len(xs), len(ys), len(zs)
    sys.stderr.write("[background] tip chord %.4g m -> refine spacing %.4g m "
                     "(%.2f c_tip), growth %.3g\n"
                     % (ctip, h0, h0/ctip, ratio))
    sys.stderr.write("[background] dims %d x %d x %d = %.2fM cells\n"
                     % (ni, nj, nk, (ni-1)*(nj-1)*(nk-1)/1e6))

    X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")   # [i,j,k]
    with open(out, "w") as f:
        f.write("1\n%d %d %d\n" % (ni, nj, nk))
        for A in (X, Y, Z):
            B = np.transpose(A, (2, 1, 0))             # k slowest -> i fastest
            f.write("\n".join("%.10g" % v for v in B.ravel()))
            f.write("\n")
    sys.stderr.write("[background] wrote %s\n" % out)


if __name__ == "__main__":
    inp = sys.argv[1]
    cfg = read_input(inp)
    out = sys.argv[2] if len(sys.argv) > 2 else "backgroundVol.xyz"
    build(cfg, out, base=os.path.dirname(os.path.abspath(inp)))
