#!/bin/bash
# -----------------------------------------------------------------------------
#  make_rotor.sh — one-shot rotor mesh generation.  Outputs are written NEXT TO
#  the input file (the per-rotor example folder):
#      bladeSurf.fmt      single-blade skin (PLOT3D surface)
#      bladeVol.xyz       single-blade BL volume (i = wall-normal)
#      rotorVol.xyz       full rotor (nBlades copies rotated about +x)
#      backgroundVol.xyz  structured Cartesian overset background (refine box)
#      *_vtk.vtm          ParaView versions
#
#  Usage:  make_rotor.sh <rotor.dat> [pyhyp-python]
#          PYHYP_PYTHON=<python-with-pyhyp> make_rotor.sh <rotor.dat>
# -----------------------------------------------------------------------------
set -e
IN="$(readlink -f "$1")"
DIR="$(dirname "$IN")"
HERE="$(cd "$(dirname "$0")" && pwd)"
PY="${2:-${PYHYP_PYTHON:-python3}}"

python3 "$HERE/match_spacing.py" "$IN" || true    # overset spacing advisory
python3 "$HERE/blade_surface.py" "$IN" "$DIR/bladeSurf.fmt"
"$PY" "$HERE/march.py" "$IN" "$DIR/bladeSurf.fmt" "$DIR/bladeVol.xyz"
python3 "$HERE/to_vtk.py" "$DIR/bladeSurf.fmt" "$DIR/bladeSurf_vtk"
python3 "$HERE/to_vtk.py" "$DIR/bladeVol.xyz" "$DIR/bladeVol_vtk"
[ -f "$DIR/rotorVol.xyz" ] && python3 "$HERE/to_vtk.py" "$DIR/rotorVol.xyz" "$DIR/rotorVol_vtk"
python3 "$HERE/background_mesh.py" "$IN" "$DIR/backgroundVol.xyz"
python3 "$HERE/to_vtk.py" "$DIR/backgroundVol.xyz" "$DIR/backgroundVol_vtk"
rm -f "$DIR/bladeVol.xyz.raw"
echo "outputs in $DIR"
