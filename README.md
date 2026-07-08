# pyhyp-rotor-overset-generator

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21257648.svg)](https://doi.org/10.5281/zenodo.21257648)

Structured overset (chimera) mesh generation for rotors with
[pyHyp](https://github.com/mdolab/pyhyp): a wall-resolved (y+ ≈ 1)
boundary-layer **blade component mesh** from a plain-text per-section
planform table, full-rotor replication, and a matching structured Cartesian
**background mesh** with a rotor refinement box.

The output is plain formatted **PLOT3D**, so the mesh is solver-agnostic:
convert/import it into OpenFOAM, SU2, CGNS-based codes, or use it as the
near-body component of an overset (chimera) setup.

```
rotor.dat ──▶ blade_surface.py ──▶ skin.fmt ──▶ march.py (pyHyp) ──▶ bladeVol.xyz
               (watertight skin,                     (hyperbolic BL,        │ nBlades > 1
                multiblock surface)               i = wall-normal)      ▼
                                                                   rotorVol.xyz
```

Caradonna-Tung example (ParaView):

| blade skin (single blade) | BL volume — watertight closure |
|---|---|
| ![blade surface](examples/caradonna_tung/image/surface.png) | ![BL volume](examples/caradonna_tung/image/rotor_volume.png) |
| **section slice — O-grid, y+≈1 layers** | **Cartesian background — refinement box** |
| ![section slice](examples/caradonna_tung/image/rotor_volume_slice.png) | ![background](examples/caradonna_tung/image/background_volume.png) |

One-shot (outputs land next to the input file):

```bash
PYHYP_PYTHON=<python-with-pyhyp> ./make_rotor.sh examples/rotor_sc1095/rotor_sc1095.dat
```

## Why

gmsh/snappy-style extruded prisms struggle to produce robust y+≈1 layers on
thin lifting surfaces (self-intersection at the TE, ~60° non-orthogonality
for tetrahedral BLs). A pyHyp hyperbolic march from a structured skin gives
orthogonal hexahedral layers (typically **~2° mean non-orthogonality**) with a
first-cell spacing you set directly in meters. The hard part — a watertight
(watertight) surface topology that pyHyp accepts — is what this package
automates.

## Axis convention (rotor frame)

```
      -x  (suction side / thrust side)
        │
        │      y  (span, tip)
        │    ╱
        │   ╱
        └──╱─────────▶  z   (chordwise: the LE FACES -z, TE toward +z)

  x : rotor axis — rotor wake direction is +x
```

* **twist** is applied about the local **quarter chord (0.25 c)**;
  positive twist = nose-up (LE rotates toward −x)
* **LE_z** moves the leading edge fore/aft along +z (sweep), in **meters**,
  *before* the twist rotation
* volume output index order: **i = wall-normal** (i=1 on the wall),
  j = airfoil perimeter (TE_lower → LE → TE_upper → blunt-TE seal),
  k = span (root → tip, +y)

## Input file

One plain-text `.dat` file: keyword lines followed by a `SECTIONS` table
(`#` starts a comment). Examples are organised one folder per rotor, with
airfoil coordinate files in an `airfoil/` subfolder; generated meshes are
written into the rotor folder (and are git-ignored):

```
examples/
  rotor_sc1095/
    rotor_sc1095.dat        # this input file
    airfoil/sc1095.dat      # airfoil coordinates (paths relative to input)
    bladeSurf.fmt           # generated: single-blade skin
    bladeVol.xyz            # generated: single-blade volume
    rotorVol.xyz            # generated: full rotor (nBlades copies about +x)
    backgroundVol.xyz       # generated: structured overset background
    *_vtk.vtm               # generated: ParaView
```

```
nBlades    2            # full-rotor replication about +x (1 = single blade)
R          1.143        # tip radius [m]

# ---- surface ----
nChord     200          # chordwise points per side                (default 160)
nTE        7            # blunt-TE seal interior points, ODD       (default 7)
teCut      0.96         # TE truncation (x/c)                      (default 0.97)
dTE_c      0.003        # chordwise spacing at the TE (x/c)        (default 0.003)
#dLE_c     auto         # LE chordwise spacing (x/c); default auto = 0.16*r_LE/k
nSpan      70           # spanwise stations                        (default 60)
#dRootFrac auto         # spanwise tanh clustering at the root; default uniform
#dTipFrac  auto         # spanwise tanh clustering at the tip;  default uniform
#rootCut   auto         # span start y/R; default = first SECTIONS r/R
closedSock 1            # 1: watertight tips (TFI caps) / 0: open + splay (default 1)
capDome    0.25         # rounded tip cap, 0 flat .. ~1 half-thickness dome
                        #   (default 0; 0.2-0.4 recommended)
#datSmooth 5            # smoothing passes for tabulated airfoils  (default 5)

# ---- march (pyHyp) ----
firstLayer 2.78e-6      # first wall spacing [m]  (y+ target)
nLayers    76           # marching layers (overridden when autoMatch 1)
marchDist  0.19         # total march distance [m]
autoMatch  1            # 1: force nLayers so the outer wall-normal cell
                        #    matches the background refine spacing / 0: off
#matchFactor 0.9        # target outer cell = matchFactor * h_bg   (default 0.9)
#splay     0.25         # open-end free-edge splay                 (default 0.25)
#volSmoothIter 100      # pyHyp volume smoothing iterations        (default 100)
#volBlend  0.0005       # pyHyp volume blending                    (default 0.0005)
#volCoef   0.25         # pyHyp volume coefficient                 (default 0.25)
#cMax      3.0          # pyHyp marching cMax                      (default 3.0)
#epsE      1.0          # pyHyp explicit smoothing                 (default 1.0)
#epsI      2.0          # pyHyp implicit smoothing                 (default 2.0)
#theta     3.0          # pyHyp implicit blending                  (default 3.0)
#nConstantStart 1       # layers marched at constant s0            (default 1)

# ---- overset background (background_mesh.py) ----
#bgSpacing 0.15         # refine-box spacing [tip chords]          (default 0.15)
#bgGrowth  1.12         # spacing growth ratio outside the box     (default 1.12)
#bgXmin    -4           # domain extents [R]; +x = wake/downstream (defaults:
#bgXmax    8            #   x in [-4, 8], y and z in [-4, 4])
#bgYmin    -4
#bgYmax    4
#bgZmin    -4
#bgZmax    4
#refXmin   -0.5         # refinement box [R] (defaults: 0.5R upstream,
#refXmax   2.0          #   2R of wake downstream, radius 1.2R)
#refYmin   -1.2
#refYmax   1.2
#refZmin   -1.2
#refZmax   1.2

SECTIONS
# r/R    chord[m]   twist[deg]  LE_z[m]  airfoil
0.19     0.1905     8.0         0.0      naca0012
1.00     0.1905     8.0         0.0      naca0012
```

| column | meaning |
|---|---|
| `r/R` | span station y/R |
| `chord` | local chord [m] |
| `twist` | nose-up twist about the local 0.25c [deg] |
| `LE_z` | fore/aft LE position along +z, in METERS (sweep) |
| `airfoil` | `nacaXXXX` (4-digit) or path to a Selig-format `.dat` file |

Airfoil files: both **Selig** (one TE→LE→TE loop) and **Lednicer** (point
counts, then upper/lower LE→TE) formats are auto-detected, so files from the
[UIUC Airfoil Coordinates Database](https://m-selig.ae.illinois.edu/ads/coord_database.html)
work directly. Tabulated ordinates are resampled with a cubic spline in
√x (smooth nose) after `datSmooth` (default 5) endpoint-preserving Laplacian
passes — digitised data carries point-to-point noise that otherwise folds the
hyperbolic march. Paths are resolved relative to the input file.

Examples using airfoil files:
* `examples/rotor_23012/` — NACA 23012 (5-digit, generated analytically)
* `examples/rotor_sc1095/` — Sikorsky SC1095 (UH-60 section, downloaded from
  the UIUC database, Lednicer format)

Between stations chord/twist/LE_z vary linearly and airfoil ordinates are
blended linearly. Commented keywords above show every optional knob with its
default; the same fully annotated block ships in every `examples/*/*.dat`.

## Usage

One-shot (surface -> match advisory -> march -> full rotor -> background -> VTK,
outputs written next to the input file):

```bash
PYHYP_PYTHON=<python-with-pyhyp> ./make_rotor.sh examples/caradonna_tung/caradonna_tung.dat
```

Or step by step:

```bash
# 1. surface (any python3 + numpy)
python3 blade_surface.py examples/caradonna_tung/caradonna_tung.dat skin.fmt

# 2. march (a python with pyHyp installed, e.g. the MDO-lab / DAFoam conda env)
<pyhyp-python> march.py examples/caradonna_tung/caradonna_tung.dat skin.fmt bladeVol.xyz

# 3. optional: ParaView-ready VTK (no VTK library required)
python3 to_vtk.py bladeVol.xyz          # -> bladeVol.vtm + bladeVol/block*.vts
```

`march.py` prints pyHyp's quality table; look for `Normals are consistent!`
and a positive `Min Quality`. The single-blade result is a 3-block PLOT3D
volume (main O-grid + tip cap + root cap) in the i=wall-normal ordering.
With `nBlades > 1` the single-blade volume is additionally replicated by
rotation about the rotor axis (+x) into `rotorVol.xyz` (nBlades x 3 blocks);
the SURFACE stays single-blade — only one blade is ever marched.

## Overset background mesh

`background_mesh.py` (also run by `make_rotor.sh`) builds a structured
single-block Cartesian background for overset assemblies from the same input:
a **refinement box** around the rotor with uniform spacing given in **tip
chords** (default `bgSpacing 0.15`), growing geometrically away from the box
(`bgGrowth`, default 1.12) out to the domain boundary.

```
bgSpacing  0.15         # refine-box spacing [tip chords]
bgGrowth   1.12         # spacing growth ratio outside the box
bgXmin -4   bgXmax 8    # domain extents [R]  (+x = wake/downstream)
bgYmin -4   bgYmax 4
bgZmin -4   bgZmax 4
refXmin -0.5  refXmax 2.0     # refinement box [R]: 0.5R upstream, 2R of
refYmin -1.2  refYmax 1.2     #  wake downstream, radius 1.2R
refZmin -1.2  refZmax 1.2
```

### Matching the blade outer spacing to the background

For overset interpolation quality the OUTERMOST wall-normal cell of the blade
mesh should be equal to or slightly smaller than the background refine
spacing. pyHyp has no such option — the outer spacing is an implicit result
of `(firstLayer, nLayers, marchDist)` — so `match_spacing.py` solves the
geometric-series design problem for you:

```bash
python3 match_spacing.py examples/caradonna_tung/caradonna_tung.dat            # report
python3 match_spacing.py examples/caradonna_tung/caradonna_tung.dat --apply    # write nLayers
python3 match_spacing.py rotor.dat --check bladeVol.xyz   # measure a marched volume
```

It recommends `nLayers` (and reports the implied growth ratio) so the outer
spacing lands at `matchFactor` x `bgSpacing` x c_tip (`matchFactor` keyword,
default 0.9), warns when the ratio exceeds ~1.3, and `--check` measures the
actual outer i-spacing of a generated volume (validated: prediction within
~2% of the marched result). `make_rotor.sh` prints the advisory automatically.

To make this fully automatic, set **`autoMatch 1`** in the rotor `.dat`:
`march.py` then FORCES the matched `nLayers` (any `nLayers` in the file is
overridden), logging e.g.
`autoMatch: nLayers 76 -> 64 (ratio 1.1559, outer 0.0256 m = 0.90 x h_bg)`.

Mind the cell count for small tip chords: at `bgSpacing 0.15` a c_tip = 0.07 m
blade on R = 1 m gives ~15M background cells; the small-chord examples ship
with `bgSpacing 0.4` for compactness.

Note on the march log: the first few levels near the blunt-TE cap corners can
report `Min Quality -1` with tiny negative volumes — the march recovers within
~30 levels **provided the growth ratio stays moderate (~1.15)**. Very
aggressive tests (few layers over a large `marchDist`) collapse instead.

## Topology notes (the parts that are easy to get wrong)

* **Watertight skin**: the skin is a single watertight multiblock quilt — a main O-grid
  (i = perimeter, j = span) plus two **non-degenerate Coons-TFI caps**. pyHyp's
  normal check registers every quad's four *directed* edges and aborts on any
  repeat, so a cap with a collapsed edge (point/camber-line tip) can NEVER
  pass: each cap's four edges are exact point-subsets of the end perimeter
  (LE nose arc / lower surface / blunt-TE seal / upper surface), and the cap
  corners become 4-edge/3-face "LCorner" nodes, which pyHyp's marcher
  supports.
* **`nTE` must be odd** — the LE nose arc has `nTE+2` points placed
  symmetrically about the LE.
* **Chordwise spacing is a two-sided tanh** (Vinokur), not cosine: cosine
  clustering shrinks TE cells to ~1e-4 c, and the cap columns at the 90°
  blunt-TE corner then shear/invert during the march. Keep `dTE_c` comparable
  to the seal spacing (default 0.003) and let `dLE_c` default so the nose arc
  stays at nose-radius scale.
* **Blunt TE**: the TE is truncated at `teCut` and sealed with `nTE` points —
  a sharp TE folds the hyperbolic march.
* The outward orientation of the quilt is verified numerically at write time
  and all blocks are flipped together if needed, so the march direction is
  always out of the body.
* Set `closedSock 0` to fall back to a single open-ended block whose free
  span edges pyHyp splays (quick tests; span ends are then not meshed).
* **Rounded tip cap**: the caps are flat by default; `capDome` (0..~1) bulges
  the cap INTERIOR outboard into a slightly rounded tip (a local
  half-thickness dome, `capDome 1` ~ semicircular cross-section; 0.2-0.4
gives a subtle edge rounding and is the recommended range). Cap
  boundary points are untouched, so the block stitching and pyHyp's normal
  check are unaffected (validated: Caradonna-Tung with `capDome 0.6`, final
  Min Quality 0.41). For a true CAD-defined rounded cap you can also build
  the skin externally (e.g. pyGeo's rounded-tip lifting surfaces, or project
  onto a CAD tip) and feed the resulting multiblock PLOT3D skin straight to
  `march.py` — the watertight-topology rules in this README still apply.
* **Known limitation** — cambered sections with large twist at the blade
  ENDS leave a small number of non-convex corner cells in the end-cap blocks
  near the LE (0 for symmetric sections; ~0.1% of cap cells for SC1095 at
  10° root twist). The march still completes with the usual final quality
  (~0.42); check whether your solver tolerates these cells, reduce the twist
  at the root station, or use `closedSock 0` if the span ends are not needed.

## Requirements

* python3 + numpy (surface generation, VTK export — any OS incl. Windows)
* [pyHyp](https://github.com/mdolab/pyhyp) (volume march only)

### Windows

`blade_surface.py` and `to_vtk.py` are pure python3+numpy and run natively on
Windows. pyHyp itself is Linux/macOS software, so run the march step under
**WSL2** (or a Docker image that ships pyHyp, e.g. the DAFoam/MDO-lab images);
`make_rotor.sh` is a bash script and runs as-is inside WSL. The generated
PLOT3D/VTK files are plain text and portable across OSes.

## Citation

If you use this software in your research, please cite it. GitHub's
*"Cite this repository"* button uses `CITATION.cff`; an archived, DOI-carrying
snapshot of each release is available via Zenodo:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21257648.svg)](https://doi.org/10.5281/zenodo.21257648)

```bibtex
@software{joo_pyhyp_rotor_overset_generator,
  author  = {Joo, Seunghyun},
  title   = {pyhyp-rotor-overset-generator: structured overset mesh
             generation for rotors with pyHyp},
  url     = {https://github.com/SEUNGHYUN-JOO/pyhyp-rotor-overset-generator},
  doi     = {10.5281/zenodo.21257648},
  license = {MIT}
}
```

## License

MIT (this package). pyHyp is a separate project by the MDO Lab, licensed under the Apache License 2.0 (http://www.apache.org/licenses/LICENSE-2.0); it is not bundled here — install it independently and comply with its license.
