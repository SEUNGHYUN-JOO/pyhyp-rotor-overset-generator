#!/usr/bin/env python3
# -----------------------------------------------------------------------------
#  match_spacing.py — match the OUTERMOST wall-normal (i-direction) cell size
#  of the blade BL mesh to the background refinement-zone spacing (the overset
#  donor-quality rule: outer blade cell equal to or slightly smaller than the
#  background cell).
#
#  pyHyp has no such option: with ps0 = pGridRatio = -1 it distributes N
#  layers over marchDist as a (near-)geometric progression from s0, so the
#  outer spacing is an IMPLICIT result of (s0, N, marchDist).  This script
#  solves the geometric-series design problem instead:
#
#      s0*(r^N - 1)/(r - 1) = marchDist          (total distance)
#      s0*r^(N-1)           = matchFactor * h_bg (target outer spacing)
#
#  ->  r = (marchDist - s0) / (marchDist - h_target),
#      N = 1 + ln(h_target/s0)/ln(r),  then N is rounded and r re-solved for
#  the integer N so the total distance stays exact.
#
#  h_bg = bgSpacing * c_tip (the background refinement-box spacing);
#  matchFactor (keyword, default 0.9) keeps the blade cell slightly smaller.
#
#  Usage:
#      match_spacing.py <rotor.dat>                 # report recommendation
#      match_spacing.py <rotor.dat> --apply         # also patch nLayers in file
#      match_spacing.py <rotor.dat> --check <bladeVol.xyz>
#                                                   # measure the ACTUAL outer
#                                                   # spacing of a marched
#                                                   # volume vs h_bg
# -----------------------------------------------------------------------------
import sys, os, math, re
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from blade_surface import read_input, Planform


def solve_ratio(s0, D, N):
    """Growth ratio r of a geometric series: s0*(r^N-1)/(r-1) = D (r > 1)."""
    r = (D/s0)**(1.0/max(N-1, 1))     # start guess
    for _ in range(200):
        f = s0*(r**N - 1.0)/(r - 1.0) - D
        df = s0*(N*r**(N-1)*(r - 1.0) - (r**N - 1.0))/(r - 1.0)**2
        step = f/df
        r -= step
        if abs(step) < 1e-14:
            break
    return r


def design(cfg, base):
    pf = Planform(cfg, base)
    ctip = float(pf.chord[-1])
    m = cfg.get("march", {})
    bg = cfg.get("background", {})
    s0 = float(m.get("firstLayer", 3e-6))
    D = float(m.get("marchDist", 0.1))
    h_bg = float(bg.get("bgSpacing", 0.15))*ctip
    fac = float(m.get("matchFactor", 0.9))
    h_t = fac*h_bg
    if h_t >= D:
        raise SystemExit("target outer spacing >= marchDist — increase "
                         "marchDist or refine the background")

    r0 = (D - s0)/(D - h_t)
    N0 = 1.0 + math.log(h_t/s0)/math.log(r0)
    best = None
    for N in (int(math.floor(N0)), int(math.ceil(N0)), int(math.ceil(N0)) + 1):
        if N < 4:
            continue
        r = solve_ratio(s0, D, N)
        h_last = s0*r**(N - 1)
        ok = h_last <= h_bg
        cand = (not ok, abs(h_last - h_t), N, r, h_last, ok)
        if best is None or cand < best:
            best = cand
    _, _, N, r, h_last, ok = best
    return {"s0": s0, "marchDist": D, "c_tip": ctip, "h_bg": h_bg,
            "matchFactor": fac, "nLayers": N, "ratio": r, "h_outer": h_last,
            "ok": ok}


def report(d, current_N=None):
    print("background refine spacing  h_bg = %.5g m (bgSpacing x c_tip %.4g)"
          % (d["h_bg"], d["c_tip"]))
    print("first layer s0 = %.4g m,  marchDist = %.4g m" % (d["s0"], d["marchDist"]))
    if current_N:
        rc = solve_ratio(d["s0"], d["marchDist"], current_N)
        hc = d["s0"]*rc**(current_N - 1)
        print("current  nLayers %-4d -> ratio %.4f, outer spacing %.5g m "
              "(%.2f x h_bg)%s" % (current_N, rc, hc, hc/d["h_bg"],
              "" if hc <= d["h_bg"] else "  ** LARGER than background **"))
    print("matched  nLayers %-4d -> ratio %.4f, outer spacing %.5g m "
          "(%.2f x h_bg, target factor %.2f)"
          % (d["nLayers"], d["ratio"], d["h_outer"], d["h_outer"]/d["h_bg"],
             d["matchFactor"]))
    if d["ratio"] > 1.30:
        print("WARNING: growth ratio %.3f > 1.3 — increase marchDist or "
              "nLayers for march robustness" % d["ratio"])


def check_volume(volfile, h_bg):
    v = open(volfile).read().split()
    nb = int(v[0])
    dims = [(int(v[1+3*b]), int(v[2+3*b]), int(v[3+3*b])) for b in range(nb)]
    off = 1 + 3*nb
    worst = 0.0
    for (ni, nj, nk) in dims:          # (march, perim, span)
        n = ni*nj*nk
        P = np.empty((nk, nj, ni, 3))
        for c in range(3):
            P[..., c] = np.array(v[off:off+n], float).reshape(nk, nj, ni)
            off += n
        d_out = np.linalg.norm(P[:, :, -1, :] - P[:, :, -2, :], axis=-1)
        worst = max(worst, float(d_out.max()))
    print("measured max outer i-spacing = %.5g m  (%.2f x h_bg)%s"
          % (worst, worst/h_bg,
             "  OK" if worst <= h_bg else "  ** LARGER than background **"))
    return worst


def main():
    inp = sys.argv[1]
    cfg = read_input(inp)
    base = os.path.dirname(os.path.abspath(inp))
    d = design(cfg, base)
    mm = cfg.get("march", {})
    cur = mm.get("nLayers")
    report(d, current_N=int(cur) if cur else None)
    if int(mm.get("autoMatch", 0)):
        print("autoMatch 1: march.py will FORCE nLayers %d automatically"
              % d["nLayers"])

    if "--check" in sys.argv:
        vol = sys.argv[sys.argv.index("--check") + 1]
        check_volume(vol, d["h_bg"])

    if "--apply" in sys.argv:
        s = open(inp).read()
        if re.search(r"^nLayers\s+\S+", s, re.M):
            s = re.sub(r"^(nLayers\s+)\S+", r"\g<1>%d" % d["nLayers"], s,
                       flags=re.M)
        else:
            s = s.replace("SECTIONS", "nLayers    %d\n\nSECTIONS" % d["nLayers"], 1)
        open(inp, "w").write(s)
        print("applied: nLayers %d written to %s" % (d["nLayers"], inp))


if __name__ == "__main__":
    main()
