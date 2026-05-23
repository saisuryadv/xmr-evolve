"""
Minimal mr3_gk stub for the Fortran-only self-contained package.

This file exists ONLY so that evaluate.py (which does
`from mr3_gk import bidiag_svd`) can import successfully.
The actual SVD is performed by the Fortran executable
mr3gk_fortran/mr3gk_run via subprocess + binary I/O.
"""
import os
import struct
import subprocess
import numpy as np

_ROOT = os.path.dirname(os.path.abspath(__file__))
_MR3GK_RUN = os.path.join(_ROOT, "mr3gk_fortran", "mr3gk_run")


def bidiag_svd(d, e):
    """Bidiagonal SVD via the pure-Fortran MR3-GK driver."""
    n = len(d)
    if n == 0:
        return np.zeros(0), np.zeros((0, 0)), np.zeros((0, 0)), 0

    in_path = "/tmp/_mr3gk_in.bin"
    out_path = "/tmp/_mr3gk_out.bin"

    with open(in_path, "wb") as f:
        f.write(struct.pack("i", n))
        f.write(np.asarray(d, dtype=np.float64).tobytes())
        f.write(np.asarray(e, dtype=np.float64).tobytes())

    r = subprocess.run([_MR3GK_RUN, in_path, out_path],
                       capture_output=True, timeout=120)
    if r.returncode != 0:
        raise RuntimeError(
            f"mr3gk_run failed (rc={r.returncode}): "
            + r.stderr.decode(errors="ignore")[:200])

    with open(out_path, "rb") as f:
        info = struct.unpack("i", f.read(4))[0]
        sigma = np.frombuffer(f.read(8 * n), dtype=np.float64).copy()
        U = np.frombuffer(f.read(8 * n * n), dtype=np.float64) \
              .reshape(n, n, order="F").copy()
        V = np.frombuffer(f.read(8 * n * n), dtype=np.float64) \
              .reshape(n, n, order="F").copy()

    return sigma, U, V, int(info)
