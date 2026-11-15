#!/usr/bin/env python3
"""Write a single-band GeoTIFF as raw little-endian Float32, row-major, north row first."""
import sys
import numpy as np
import tifffile

a = tifffile.imread(sys.argv[1]).astype("<f4")
assert a.ndim == 2, a.shape
a.tofile(sys.argv[2])
print(sys.argv[2], a.shape, float(a.min()), float(a.max()))
