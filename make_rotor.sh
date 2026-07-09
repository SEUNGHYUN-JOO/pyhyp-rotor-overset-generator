#!/bin/bash
# -----------------------------------------------------------------------------
#  Author: JOO Seunghyun <chlrh45351@gmail.com>  (2026-07)
#  make_rotor.sh — one-shot rotor mesh generation.  Outputs are written to an
#  output/ folder next to the input file (all PLOT3D volumes use the .x
#  extension):
#      bladeSurf.fmt   single-blade skin (PLOT3D surface)
#      bladeVol.x      single-blade BL volume (i = wall-normal)
#      rotorVol.x      full rotor (nBlades copies rotated about +x)
#      blade1.x ...    per-blade overset component volumes (blade1 = +y span,
#                      then rotating +y -> -z -> -y -> +z about +x)
#      background.x    structured Cartesian overset background (refine box;
#                      bgQuarter 1 = 90-deg sector, diagonal on +y)
#      *_vtk.vtm       ParaView versions
#
#  Usage:  make_rotor.sh <rotor.dat> [pyhyp-python]
#          PYHYP_PYTHON=<python-with-pyhyp> make_rotor.sh <rotor.dat>
# -----------------------------------------------------------------------------
set -e
IN="$(readlink -f "$1")"
DIR="$(dirname "$IN")/output"
mkdir -p "$DIR"
HERE="$(cd "$(dirname "$0")" && pwd)"
PY="${2:-${PYHYP_PYTHON:-python3}}"

python3 "$HERE/match_spacing.py" "$IN" || true    # overset spacing advisory
python3 "$HERE/blade_surface.py" "$IN" "$DIR/bladeSurf.fmt"
"$PY" "$HERE/march.py" "$IN" "$DIR/bladeSurf.fmt" "$DIR/bladeVol.x"
python3 "$HERE/to_vtk.py" "$DIR/bladeSurf.fmt" "$DIR/bladeSurf_vtk"
python3 "$HERE/to_vtk.py" "$DIR/bladeVol.x" "$DIR/bladeVol_vtk"
[ -f "$DIR/rotorVol.x" ] && python3 "$HERE/to_vtk.py" "$DIR/rotorVol.x" "$DIR/rotorVol_vtk"
python3 "$HERE/background_mesh.py" "$IN" "$DIR/background.x"
python3 "$HERE/to_vtk.py" "$DIR/background.x" "$DIR/background_vtk"
rm -f "$DIR/bladeVol.x.raw"
echo "outputs in $DIR"
