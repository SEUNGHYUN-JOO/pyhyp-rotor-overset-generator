#!/usr/bin/env python3
# -----------------------------------------------------------------------------
#  to_vtk.py — convert a (multiblock) formatted PLOT3D volume/surface to VTK
#  XML structured grids (.vts + a .vtm collection) for ParaView.
#
#  Usage:  to_vtk.py <mesh.xyz> [outBase]
#     ->   outBase/block00.vts ... + outBase.vtm
#  Works for both the i=wall-normal volumes from march.py and the surface
#  skins from blade_surface.py (nk=1).  Pure numpy, no VTK library needed.
# -----------------------------------------------------------------------------
import sys, os, base64, struct
import numpy as np


def read_plot3d(fn):
    v = open(fn).read().split()
    nb = int(v[0])
    dims = [(int(v[1+3*b]), int(v[2+3*b]), int(v[3+3*b])) for b in range(nb)]
    off = 1 + 3*nb
    out = []
    for (ni, nj, nk) in dims:
        n = ni*nj*nk
        arr = np.array(v[off:off+3*n], dtype=float); off += 3*n
        # PLOT3D formatted: i fastest, then j, then k
        comp = [arr[c*n:(c+1)*n].reshape((nk, nj, ni)) for c in range(3)]
        out.append((ni, nj, nk, comp))
    return out


def write_vts(fn, ni, nj, nk, comp):
    # VTK structured grid: x fastest -> same as PLOT3D i fastest
    pts = np.empty((nk, nj, ni, 3))
    for c in range(3):
        pts[..., c] = comp[c]
    raw = pts.reshape(-1, 3).astype("<f4").tobytes()
    payload = struct.pack("<I", len(raw)) + raw
    b64 = base64.b64encode(payload).decode()
    with open(fn, "w") as f:
        f.write('<?xml version="1.0"?>\n')
        f.write('<VTKFile type="StructuredGrid" version="0.1" '
                'byte_order="LittleEndian">\n')
        ext = "0 %d 0 %d 0 %d" % (ni-1, nj-1, nk-1)
        f.write('  <StructuredGrid WholeExtent="%s">\n' % ext)
        f.write('    <Piece Extent="%s">\n' % ext)
        f.write('      <Points>\n')
        f.write('        <DataArray type="Float32" NumberOfComponents="3" '
                'format="binary">\n')
        f.write(b64 + "\n")
        f.write('        </DataArray>\n      </Points>\n')
        f.write('    </Piece>\n  </StructuredGrid>\n</VTKFile>\n')


def main():
    fn = sys.argv[1]
    base = sys.argv[2] if len(sys.argv) > 2 else os.path.splitext(fn)[0]
    os.makedirs(base, exist_ok=True)
    blocks = read_plot3d(fn)
    entries = []
    for b, (ni, nj, nk, comp) in enumerate(blocks):
        out = os.path.join(base, "block%02d.vts" % b)
        write_vts(out, ni, nj, nk, comp)
        entries.append(os.path.relpath(out, os.path.dirname(base) or "."))
        sys.stderr.write("[to_vtk] %s (%dx%dx%d)\n" % (out, ni, nj, nk))
    with open(base + ".vtm", "w") as f:
        f.write('<?xml version="1.0"?>\n<VTKFile type="vtkMultiBlockDataSet" '
                'version="1.0">\n  <vtkMultiBlockDataSet>\n')
        for b, e in enumerate(entries):
            f.write('    <DataSet index="%d" file="%s"/>\n' % (b, e))
        f.write('  </vtkMultiBlockDataSet>\n</VTKFile>\n')
    sys.stderr.write("[to_vtk] wrote %s.vtm\n" % base)


if __name__ == "__main__":
    main()
