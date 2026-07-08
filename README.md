# pyhyp-blade-mesh

Structured, wall-resolved (y+ ≈ 1) boundary-layer meshing of rotor/wing
blades with [pyHyp](https://github.com/mdolab/pyhyp) hyperbolic extrusion —
from a **per-section planform definition** to a watertight multiblock volume.

The output is plain formatted **PLOT3D**, so the mesh is solver-agnostic:
convert/import it into OpenFOAM, SU2, CGNS-based codes, or use it as the
near-body component of an overset (chimera) setup.

```
config.json ──▶ blade_surface.py ──▶ skin.fmt ──▶ march.py (pyHyp) ──▶ vol.xyz
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

## Axis convention (blade frame)

```
        z  (thickness, suction side up)
        │
        │     y  (span, tip)
        │    ╱
        │   ╱
        └──╱────────── x  (chord, LE → TE = section wake direction)
```

* positive twist = **leading edge up** (+z), applied **about the LE**
* volume output index order: **i = wall-normal** (i=1 on the wall),
  j = airfoil perimeter (TE_lower → LE → TE_upper → blunt-TE seal),
  k = span (root → tip)

## Planform definition

Arbitrary planforms are built from spanwise sections; properties vary
linearly between stations and airfoil ordinates are blended:

```json
"planform": {
  "R": 1.143,
  "sections": [
    { "rR": 0.19, "chord": 0.1905, "twistDeg": 8.0,
      "xLE_c": 0.0, "zLE_c": 0.0, "airfoil": "naca0012" },
    { "rR": 1.00, "chord": 0.0700, "twistDeg": 0.5,
      "xLE_c": 0.35, "zLE_c": 0.02, "airfoil": "naca0012" }
  ]
}
```

| key | meaning |
|---|---|
| `rR` | span station y/R |
| `chord` | local chord [m] |
| `twistDeg` | nose-up twist about the LE [deg] |
| `xLE_c`, `zLE_c` | LE position in the section plane, in local chords (sweep / dihedral-flap) |
| `airfoil` | `nacaXXXX` (4-digit) or path to a Selig-format `.dat` file |

## Usage

```bash
# 1. surface (any python3 + numpy)
python3 blade_surface.py examples/caradonna_tung.json skin.fmt

# 2. march (a python with pyHyp installed, e.g. the MDO-lab / DAFoam conda env)
<pyhyp-python> march.py examples/caradonna_tung.json skin.fmt bladeVol.xyz
```

`march.py` prints pyHyp's quality table; look for `Normals are consistent!`
and a positive `Min Quality`. The result `bladeVol.xyz` is a 3-block PLOT3D
volume (main O-grid + tip cap + root cap) in the i=wall-normal ordering.

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
* Set `closedSock: false` to fall back to a single open-ended block whose
  free span edges pyHyp splays (useful for quick tests; the span ends are
  then not meshed).

## Requirements

* python3 + numpy (surface generation)
* [pyHyp](https://github.com/mdolab/pyhyp) (volume march only)

## License

MIT (this package). pyHyp is licensed separately (LGPL) by the MDO Lab.
