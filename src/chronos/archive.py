"""Find and download SPHEREx spectral images from the NASA/IPAC Infrared Science Archive.

Everything here talks to IRSA (https://irsa.ipac.caltech.edu):

* ``search``        - list every SPHEREx exposure of one detector that covers a sky position
* ``find_pairs``    - pick pairs of exposures that are comparable (same detector, nearly the
                      same wavelength at the target) but taken at different times
* ``download_cutout`` - fetch a small cutout around the target instead of the full 66 MB file

SPHEREx detail that drives the design: every detector sits behind a linear variable filter,
so the wavelength an object is seen at depends on where it lands on the detector. Two
exposures are only comparable when the target lands at a similar detector position.
"""
from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests

IMAGE_SIZE = 2040          # SPHEREx detector is 2040 x 2040 pixels
PIXEL_ARCSEC = 6.15        # approximate pixel scale

# A wavelength (micron) that only falls inside one detector, used to query one detector at a time
DETECTOR_QUERY_UM = {1: 0.93, 2: 1.37, 3: 2.04, 4: 3.10, 5: 4.12, 6: 4.70}
DETECTOR_RANGE_UM = {
    1: (0.75, 1.11), 2: (1.10, 1.65), 3: (1.63, 2.44),
    4: (2.40, 3.85), 5: (3.81, 4.44), 6: (4.41, 5.01),
}

# Helium in Earth's upper atmosphere glows at 1.083 um and shows up in band 1 as a bright stripe
# that differs from exposure to exposure. Pairs are not offered at wavelengths close to it.
HELIUM_AIRGLOW_UM = 1.083

_GRID_FILE = Path(__file__).parent / "data" / "wave_grids.json"
_grids: dict | None = None


def detector_label(detector: int) -> str:
    lo, hi = DETECTOR_RANGE_UM[detector]
    return f"Band {detector} ({lo:.2f}-{hi:.2f} um)"


# --------------------------------------------------------------------------------------
# Wavelength lookup
# --------------------------------------------------------------------------------------
def _load_grids() -> dict:
    global _grids
    if _grids is None:
        raw = json.loads(_GRID_FILE.read_text())
        _grids = {
            int(det): {k: np.asarray(grid[k], dtype=float) for k in ("x", "y", "wave", "bandwidth")}
            for det, grid in raw.items()
        }
    return _grids


def wavelength_at(detector: int, x, y):
    """Wavelength and bandwidth (micron) at detector pixel (x, y), 0-based.

    Uses the 9x9 control grid from the WCS-WAVE extension of the SPHEREx files
    (bundled in ``data/wave_grids.json``) with bilinear interpolation.
    """
    grid = _load_grids()[int(detector)]
    gx, gy = grid["x"] - 1.0, grid["y"] - 1.0      # FITS grid is 1-based
    x = np.clip(np.asarray(x, dtype=float), gx[0], gx[-1])
    y = np.clip(np.asarray(y, dtype=float), gy[0], gy[-1])
    ix = np.clip(np.searchsorted(gx, x) - 1, 0, len(gx) - 2)
    iy = np.clip(np.searchsorted(gy, y) - 1, 0, len(gy) - 2)
    fx = (x - gx[ix]) / (gx[ix + 1] - gx[ix])
    fy = (y - gy[iy]) / (gy[iy + 1] - gy[iy])

    def interp(values):
        return (
            values[iy, ix] * (1 - fx) * (1 - fy)
            + values[iy, ix + 1] * fx * (1 - fy)
            + values[iy + 1, ix] * (1 - fx) * fy
            + values[iy + 1, ix + 1] * fx * fy
        )

    return interp(grid["wave"]), interp(grid["bandwidth"])


# --------------------------------------------------------------------------------------
# Where does the target land on each detector image?
# --------------------------------------------------------------------------------------
def _gnomonic(ra, dec, ra0, dec0):
    """Tangent-plane offsets (degrees) of (ra, dec) about (ra0, dec0). All inputs in degrees."""
    ra, dec, ra0, dec0 = (np.radians(np.asarray(v, dtype=float)) for v in (ra, dec, ra0, dec0))
    cos_c = np.sin(dec0) * np.sin(dec) + np.cos(dec0) * np.cos(dec) * np.cos(ra - ra0)
    xi = np.cos(dec) * np.sin(ra - ra0) / cos_c
    eta = (np.cos(dec0) * np.sin(dec) - np.sin(dec0) * np.cos(dec) * np.cos(ra - ra0)) / cos_c
    return np.degrees(xi), np.degrees(eta)


def parse_region(region: str) -> np.ndarray:
    """Turn 'POLYGON ICRS ra1 dec1 ra2 dec2 ...' into an array of the first four corners."""
    numbers = [float(tok) for tok in str(region).split() if _is_number(tok)]
    return np.asarray(numbers[:8], dtype=float).reshape(4, 2)


def _is_number(token: str) -> bool:
    try:
        float(token)
        return True
    except ValueError:
        return False


def estimate_pixel(corners: np.ndarray, center_ra, center_dec, ra: float, dec: float):
    """Approximate detector pixel of (ra, dec) from the image footprint polygons.

    ``corners`` has shape (n, 4, 2): the sky positions of detector pixels
    (0,0), (N,0), (N,N), (0,N) as listed in the archive's ``s_region`` column.
    Accurate to roughly 15 pixels, which is enough to rank exposures by wavelength;
    the exact value is read from each downloaded cutout later.
    """
    corners = np.asarray(corners, dtype=float).reshape(-1, 4, 2)
    center_ra = np.asarray(center_ra, dtype=float)
    center_dec = np.asarray(center_dec, dtype=float)
    px, py = _gnomonic(corners[:, :, 0], corners[:, :, 1], center_ra[:, None], center_dec[:, None])
    tx, ty = _gnomonic(ra, dec, center_ra, center_dec)

    ax, ay = px[:, 1] - px[:, 0], py[:, 1] - py[:, 0]      # detector x direction
    bx, by = px[:, 3] - px[:, 0], py[:, 3] - py[:, 0]      # detector y direction
    dx, dy = tx - px[:, 0], ty - py[:, 0]
    det = ax * by - bx * ay
    u = (dx * by - bx * dy) / det
    v = (ax * dy - dx * ay) / det
    return u * IMAGE_SIZE, v * IMAGE_SIZE


# --------------------------------------------------------------------------------------
# Archive search
# --------------------------------------------------------------------------------------
def list_collections() -> list[str]:
    """Names of the SPHEREx spectral-image collections currently published by IRSA."""
    from astroquery.ipac.irsa import Irsa

    names = [str(n) for n in Irsa.list_collections(filter="spherex")["collection"]]
    return sorted(n for n in names if "cube" not in n and "cal" not in n)


def search(ra: float, dec: float, detector: int, collections: list[str] | None = None) -> pd.DataFrame:
    """All exposures of one detector that cover (ra, dec), oldest first.

    Returns one row per exposure with the observation time, download URL, the estimated
    detector position of the target and the wavelength it is seen at there.
    """
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from astropy.time import Time
    from astroquery.ipac.irsa import Irsa

    if collections is None:
        collections = list_collections()
    coord = SkyCoord(ra, dec, unit="deg")
    centre = DETECTOR_QUERY_UM[int(detector)] * 1e-6
    frames = []
    for name in collections:
        table = Irsa.query_sia(
            pos=(coord, 5 * u.arcsec),
            collection=name,
            band=(centre * u.m, centre * 1.001 * u.m),
        )
        if len(table) == 0:
            continue
        frame = table[["obs_id", "energy_bandpassname", "t_min", "t_exptime", "s_ra", "s_dec",
                       "s_region", "access_url"]].to_pandas()
        frame["collection"] = name
        frames.append(frame)

    columns = ["obs_id", "detector", "mjd", "date", "exptime", "x_est", "y_est", "wave_um",
               "bandwidth_um", "access_url", "collection"]
    if not frames:
        return pd.DataFrame(columns=columns)

    df = pd.concat(frames, ignore_index=True)
    df = df[df["energy_bandpassname"].astype(str).str.endswith(f"D{int(detector)}")]
    df = df.drop_duplicates(subset="access_url").reset_index(drop=True)
    if df.empty:
        return pd.DataFrame(columns=columns)

    corners = np.stack([parse_region(r) for r in df["s_region"]])
    x, y = estimate_pixel(corners, df["s_ra"].to_numpy(), df["s_dec"].to_numpy(), ra, dec)
    wave, bandwidth = wavelength_at(detector, x, y)

    out = pd.DataFrame({
        "obs_id": df["obs_id"].astype(str),
        "detector": int(detector),
        "mjd": df["t_min"].astype(float),
        "date": Time(df["t_min"].to_numpy(dtype=float), format="mjd").strftime("%Y-%m-%d %H:%M"),
        "exptime": df["t_exptime"].astype(float),
        "x_est": x,
        "y_est": y,
        "wave_um": wave,
        "bandwidth_um": bandwidth,
        "access_url": df["access_url"].astype(str),
        "collection": df["collection"],
    })
    return out.sort_values("mjd").reset_index(drop=True)


# --------------------------------------------------------------------------------------
# Choosing comparable pairs
# --------------------------------------------------------------------------------------
def find_pairs(
    obs: pd.DataFrame,
    min_gap_days: float = 30.0,
    max_gap_days: float = float("inf"),
    max_wave_frac: float = 0.25,
    margin_px: float = 60.0,
    top: int = 25,
    avoid_airglow: bool = True,
) -> pd.DataFrame:
    """Pairs of exposures that see the target at nearly the same wavelength.

    * ``min_gap_days``  - the two exposures must be at least this far apart in time
    * ``max_gap_days``  - ... and at most this far apart
    * ``max_wave_frac`` - their wavelengths may differ by at most this fraction of the bandwidth
    * ``margin_px``     - skip exposures where the target is this close to the detector edge
    * ``avoid_airglow`` - skip wavelengths within 1.5 bandwidths of the 1.083 um helium airglow line

    Pairs are ranked longest time baseline first. Each exposure appears in at most two pairs
    so the list offers variety. Columns ``a`` and ``b`` index into ``obs`` (a is the earlier).
    """
    columns = ["a", "b", "gap_days", "dwave_um", "wave_um"]
    if obs.empty:
        return pd.DataFrame(columns=columns)

    inside = (
        obs["x_est"].between(margin_px, IMAGE_SIZE - margin_px)
        & obs["y_est"].between(margin_px, IMAGE_SIZE - margin_px)
    )
    if avoid_airglow:
        inside &= (obs["wave_um"] - HELIUM_AIRGLOW_UM).abs() > 1.5 * obs["bandwidth_um"]
    usable = obs[inside].sort_values("wave_um")
    index = usable.index.to_numpy()
    wave = usable["wave_um"].to_numpy()
    mjd = usable["mjd"].to_numpy()
    tolerance = max_wave_frac * usable["bandwidth_um"].to_numpy()

    rows = []
    for i in range(len(usable)):
        stop = np.searchsorted(wave, wave[i] + tolerance[i], side="right")
        if stop <= i + 1:
            continue
        j = np.arange(i + 1, stop)
        gap = np.abs(mjd[j] - mjd[i])
        wanted = (gap >= min_gap_days) & (gap <= max_gap_days)
        j, gap = j[wanted], gap[wanted]
        if len(j) == 0:
            continue
        best = j[np.argsort(-gap)[:3]]                      # keep the longest baselines for this exposure
        for k in best:
            first, second = (i, k) if mjd[i] <= mjd[k] else (k, i)
            rows.append((index[first], index[second], abs(mjd[k] - mjd[i]),
                         abs(wave[k] - wave[i]), 0.5 * (wave[k] + wave[i])))

    if not rows:
        return pd.DataFrame(columns=columns)

    pairs = pd.DataFrame(rows, columns=columns).sort_values(
        ["gap_days", "dwave_um"], ascending=[False, True])
    used: dict[int, int] = {}
    keep = []
    for row in pairs.itertuples(index=False):
        if used.get(row.a, 0) >= 2 or used.get(row.b, 0) >= 2:
            continue
        used[row.a] = used.get(row.a, 0) + 1
        used[row.b] = used.get(row.b, 0) + 1
        keep.append(row)
        if len(keep) >= top:
            break
    out = pd.DataFrame(keep, columns=columns)
    return out.astype({"a": int, "b": int}).reset_index(drop=True)


# --------------------------------------------------------------------------------------
# Cutout download
# --------------------------------------------------------------------------------------
def cutout_url(access_url: str, ra: float, dec: float, size_deg: float) -> str:
    """IRSA cutout service: append the centre and size to the file URL."""
    return f"{access_url}?center={ra:.6f},{dec:.6f}&size={size_deg:.4f}"


def download_cutout(access_url: str, ra: float, dec: float, size_deg: float,
                    cache_dir: str | Path = "data/cutouts", timeout: float = 120.0,
                    retries: int = 12) -> Path:
    """Download (or reuse) a cutout of one exposure. Returns the local FITS path."""
    url = cutout_url(access_url, ra, dec, size_deg)
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    stem = Path(access_url).name.replace(".fits", "")
    tag = hashlib.sha1(url.encode()).hexdigest()[:10]
    path = cache_dir / f"{stem}_{tag}.fits"
    if path.exists() and path.stat().st_size > 0:
        return path

    # IRSA answers "503 Server busy, retry later" when its cutout service is loaded
    for attempt in range(retries):
        response = requests.get(url, timeout=timeout)
        if response.status_code not in (429, 502, 503, 504):
            break
        if attempt < retries - 1:
            try:
                wait = float(response.headers.get("Retry-After", 5))
            except ValueError:
                wait = 5.0
            time.sleep(min(wait, 30.0))
    response.raise_for_status()
    if not response.content.startswith(b"SIMPLE"):
        raise RuntimeError(f"IRSA did not return a FITS file for {url}")
    partial = path.with_suffix(".part")
    partial.write_bytes(response.content)
    partial.replace(path)
    return path
