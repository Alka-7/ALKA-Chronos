"""Rebuild src/chronos/data/wave_grids.json from the archive.

Each SPHEREx spectral image carries a WCS-WAVE table: wavelength and bandwidth on a 9x9
grid of detector positions. This script reads that table from one exposure per detector
and stores the six grids so the app can estimate wavelengths without downloading files.

Run from the repository root:  python scripts/build_wave_grids.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astroquery.ipac.irsa import Irsa

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from chronos import archive  # noqa: E402

RA, DEC = 270.0, 66.56          # north ecliptic pole deep field: every detector has coverage
OUT = Path(__file__).resolve().parents[1] / "src" / "chronos" / "data" / "wave_grids.json"


def main() -> None:
    coord = SkyCoord(RA, DEC, unit="deg")
    collection = next(c for c in archive.list_collections() if c.endswith("_deep"))
    grids = {}
    for detector, centre_um in archive.DETECTOR_QUERY_UM.items():
        centre = centre_um * 1e-6
        table = Irsa.query_sia(pos=(coord, 5 * u.arcsec), collection=collection,
                               band=(centre * u.m, centre * 1.001 * u.m), maxrec=1)
        url = str(table["access_url"][0])
        path = archive.download_cutout(url, RA, DEC, 0.05)
        with fits.open(path) as hdul:
            row = hdul["WCS-WAVE"].data[0]
            grids[detector] = {
                "x": [int(v) for v in row["X"]],
                "y": [int(v) for v in row["Y"]],
                "wave": row["VALUES"][:, :, 0].round(5).tolist(),         # [y, x]
                "bandwidth": row["VALUES"][:, :, 1].round(5).tolist(),
                "source": Path(url).name,
            }
        print(f"detector {detector}: {Path(url).name}")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(grids, indent=1))
    print("wrote", OUT)


if __name__ == "__main__":
    main()
