# pyhyp-blade-mesh

Structured, wall-resolved (y+ ≈ 1) boundary-layer meshing of rotor/wing
blades with [pyHyp](https://github.com/mdolab/pyhyp) hyperbolic extrusion —
from a **plain-text per-section planform table** to a watertight multiblock
volume.

The output is plain formatted **PLOT3D**, so the mesh is solver-agnostic:
convert/import it into OpenFOAM, SU2, CGNS-based codes, or use it as the
near-body component of an overset (chimera) setup.

```
blade.dat ──▶ blade_surface.py ──▶ skin.fmt ──▶ march.py (pyHyp) ──▶ vol.xyz
               (closed sock,                     (hyperbolic BL,
                multiblock surface)               i = wall-normal)
```

## Why

gmsh/snappy-style extruded prisms struggle to produce robust y+≈1 layers on
thin lifting surfaces (self-intersection at the TE, ~60° non-orthogonality
for tetrahedral BLs). A pyHyp hyperbolic march from a structured skin gives
orthogonal hexahedral layers (typically **~2° mean non-orthogonality**) with a
first-cell spacing you set directly in metres. The hard part — a watertight
("closed sock") surface topology that pyHyp accepts — is what this package
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
* **LE_z** moves the leading edge fore/aft along +z (sweep), in units of the
  **local chord**, *before* the twist rotation
* volume output index order: **i = wall-normal** (i=1 on the wall),
  j = airfoil perimeter (TE_lower → LE → TE_upper → blunt-TE seal),
  k = span (root → tip, +y)

## Input file

One plain-text `.dat` file: keyword lines followed by a `SECTIONS` table
(`#` starts a comment — see `examples/`):

```
R          1.143        # tip radius [m]

nChord     200          # chordwise points per side
nTE        7            # blunt-TE seal interior points (must be ODD)
teCut      0.96         # TE truncation (x/c)
dTE_c      0.003        # chordwise spacing at the TE (x/c)
nSpan      70           # spanwise stations
closedSock 1

firstLayer 2.78e-6      # first wall spacing [m]  (y+ target)
nLayers    76
marchDist  0.19         # total march distance [m]

SECTIONS
# r/R    chord[m]   twist[deg]  LE_z[c]  airfoil
0.19     0.1905     8.0         0.0      naca0012
1.00     0.1905     8.0         0.0      naca0012
```

| column | meaning |
|---|---|
| `r/R` | span station y/R |
| `chord` | local chord [m] |
| `twist` | nose-up twist about the local 0.25c [deg] |
| `LE_z` | fore/aft LE position along +z, in local chords (sweep) |
| `airfoil` | `nacaXXXX` (4-digit) or path to a Selig-format `.dat` file |

See `examples/rotor_23012.dat` for a blade built from a Selig airfoil file
(`examples/naca23012.dat` — a 5-digit section, hence not expressible as
`nacaXXXX`); the path is resolved relative to the input file.

Between stations chord/twist/LE_z vary linearly and airfoil ordinates are
blended linearly. Optional keywords: `dLE_c` (LE chordwise spacing),
`dRootFrac`/`dTipFrac` (spanwise tanh clustering), `rootCut`,
`splay`, `volSmoothIter`, `volBlend`, `cMax`, `epsE`, `epsI`, `theta`,
`nConstantStart`.

## Usage

```bash
# 1. surface (any python3 + numpy)
python3 blade_surface.py examples/caradonna_tung.dat skin.fmt

# 2. march (a python with pyHyp installed, e.g. the MDO-lab / DAFoam conda env)
<pyhyp-python> march.py examples/caradonna_tung.dat skin.fmt bladeVol.xyz

# 3. optional: ParaView-ready VTK (no VTK library required)
python3 to_vtk.py bladeVol.xyz          # -> bladeVol.vtm + bladeVol/block*.vts
```

`march.py` prints pyHyp's quality table; look for `Normals are consistent!`
and a positive `Min Quality`. The result is a 3-block PLOT3D volume (main
O-grid + tip cap + root cap) in the i=wall-normal ordering.

Note on the march log: the first few levels near the blunt-TE cap corners can
report `Min Quality -1` with tiny negative volumes — the march recovers within
~30 levels **provided the growth ratio stays moderate (~1.15)**. Very
aggressive tests (few layers over a large `marchDist`) collapse instead.

## Topology notes (the parts that are easy to get wrong)

* **Closed sock**: the skin is a watertight multiblock quilt — a main O-grid
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

## Requirements

* python3 + numpy (surface generation)
* [pyHyp](https://github.com/mdolab/pyhyp) (volume march only)

## License

MIT (this package). pyHyp is licensed separately (LGPL) by the MDO Lab.
