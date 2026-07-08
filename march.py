#!/usr/bin/env python3
# -----------------------------------------------------------------------------
#  march.py — extrude a blade_surface.py skin into a wall-resolved structured
#  boundary-layer volume with pyHyp (hyperbolic marching), then reorder the
#  result so that
#
#      i = WALL-NORMAL (marching) direction, i=1 on the blade wall
#      j = airfoil perimeter (TE_lower -> LE -> TE_upper -> seal)
#      k = spanwise (root -> tip, +y)
#
#  Run under a python that has pyhyp (e.g. an MDO-framework conda env):
#      python march.py <input.dat> <surf.fmt> <outVol.xyz>
#
#  march keywords in the .dat input (all lengths in METRES):
#      firstLayer : first cell wall spacing (y+ target)
#      nLayers    : number of marching layers
#      marchDist  : total marching distance
#      splay      : splay factor for open-sock free edges (single block only)
#      volSmoothIter, volBlend, cMax, epsE, epsI, theta, nConstantStart
# -----------------------------------------------------------------------------
import sys, os
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from blade_surface import read_input


def reorder_to_i_normal(vol_in, vol_out):
    """PLOT3D multiblock (i=perimeter, j=span, k=march) ->
    (i=march, j=perimeter, k=span)."""
    v = open(vol_in).read().split()
    nb = int(v[0])
    dims = [(int(v[1+3*b]), int(v[2+3*b]), int(v[3+3*b])) for b in range(nb)]
    off = 1 + 3*nb
    blocks = []
    for (ni, nj, nk) in dims:
        n = ni*nj*nk
        arr = np.array(v[off:off+3*n], dtype=float); off += 3*n
        # PLOT3D formatted ordering: i fastest, then j, then k
        comp = [arr[c*n:(c+1)*n].reshape((nk, nj, ni)) for c in range(3)]
        blocks.append(comp)
    with open(vol_out, "w") as f:
        f.write("%d\n" % nb)
        for (ni, nj, nk) in dims:
            f.write("%d %d %d\n" % (nk, ni, nj))   # (march, perim, span)
        for comp, (ni, nj, nk) in zip(blocks, dims):
            for c in range(3):
                A = comp[c]                        # [k_march, j_span, i_perim]
                B = np.transpose(A, (1, 2, 0))     # [j_span, i_perim, k_march]
                f.write("\n".join("%.10g" % x for x in B.ravel()))
                f.write("\n")
    sys.stderr.write("[march] reordered %s -> %s (i = wall-normal)\n"
                     % (vol_in, vol_out))


def main():
    cfg = read_input(sys.argv[1])
    surf = sys.argv[2]
    outvol = sys.argv[3] if len(sys.argv) > 3 else "bladeVol.xyz"
    m = cfg.get("march", {})

    from pyhyp import pyHyp
    nblocks = int(open(surf).readline())
    bc = ({} if nblocks > 1 else
          {1: {"jLow": "splay", "jHigh": "splay"}})
    options = {
        "inputFile": surf,
        "fileType": "PLOT3D",
        "unattachedEdgesAreSymmetry": False,
        "outerFaceBC": "farfield",
        "autoConnect": True,
        "BC": bc,
        "families": "wall",
        "N": int(m.get("nLayers", 64)),
        "s0": float(m.get("firstLayer", 3e-6)),
        "marchDist": float(m.get("marchDist", 0.1)),
        "nConstantStart": int(m.get("nConstantStart", 1)),
        "ps0": -1.0,
        "pGridRatio": -1.0,
        "cMax": float(m.get("cMax", 3.0)),
        "epsE": float(m.get("epsE", 1.0)),
        "epsI": float(m.get("epsI", 2.0)),
        "theta": float(m.get("theta", 3.0)),
        "volCoef": float(m.get("volCoef", 0.25)),
        "volBlend": float(m.get("volBlend", 0.0005)),
        "volSmoothIter": int(m.get("volSmoothIter", 100)),
        "splay": float(m.get("splay", 0.25)),
        "kspreltol": 1e-4,
    }
    hyp = pyHyp(options=options)
    hyp.run()
    raw = outvol + ".raw"
    hyp.writePlot3D(raw)
    reorder_to_i_normal(raw, outvol)
    sys.stderr.write("[march] wrote %s\n" % outvol)


if __name__ == "__main__":
    main()
