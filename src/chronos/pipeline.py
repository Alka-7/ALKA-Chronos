"""High-level steps shared by the command-line tool and the Streamlit app."""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from . import archive, imaging


@dataclass
class BlinkPair:
    """Two epochs of the same sky region on the same grid, ready to display."""

    ra: float
    dec: float
    size_deg: float
    a: np.ndarray                # earlier epoch, aligned, background removed (MJy/sr)
    b: np.ndarray                # later epoch
    meta_a: dict
    meta_b: dict
    name: str = ""

    @property
    def gap_days(self) -> float:
        return abs(self.meta_b["mjd"] - self.meta_a["mjd"])

    def label(self, which: str) -> str:
        meta = self.meta_a if which == "a" else self.meta_b
        return f"{meta['date'][:16].replace('T', ' ')} UTC | {meta['wave_um']:.3f} um"


def _meta(epoch: imaging.Epoch, url: str) -> dict:
    return {
        "date": epoch.date, "mjd": epoch.mjd, "detector": epoch.detector, "obs_id": epoch.obs_id,
        "exptime": epoch.exptime, "wave_um": epoch.wave_um, "bandwidth_um": epoch.bandwidth_um,
        "unit": epoch.unit, "psf_fwhm_arcsec": epoch.extra.get("psf_fwhm_arcsec"),
        "target_pixel": list(epoch.extra.get("target_pixel", (np.nan, np.nan))), "url": url,
    }


def build_pair(ra: float, dec: float, size_deg: float, url_a: str, url_b: str,
               cache_dir: str | Path = "data/cutouts", name: str = "") -> BlinkPair:
    """Download both cutouts, align them on a north-up grid and remove the background.

    The cutouts are requested 1.5x larger than the displayed field so that the rotated
    detector footprint still covers the corners of the north-up frame.
    """
    frame_wcs, shape = imaging.target_frame(ra, dec, size_deg)
    arrays, metas = [], []
    for url in (url_a, url_b):
        path = archive.download_cutout(url, ra, dec, size_deg * 1.5, cache_dir=cache_dir)
        epoch = imaging.load_epoch(path, ra, dec, source=url)
        arrays.append(imaging.align(epoch, frame_wcs, shape))
        metas.append(_meta(epoch, url))
    if metas[0]["mjd"] > metas[1]["mjd"]:
        arrays.reverse()
        metas.reverse()
    return BlinkPair(ra, dec, size_deg, arrays[0], arrays[1], metas[0], metas[1], name=name)


# --------------------------------------------------------------------------------------
# Small example files that can live in the repository
# --------------------------------------------------------------------------------------
def save_pair(pair: BlinkPair, path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    info = {"ra": pair.ra, "dec": pair.dec, "size_deg": pair.size_deg, "name": pair.name,
            "meta_a": pair.meta_a, "meta_b": pair.meta_b}
    np.savez_compressed(path, a=pair.a.astype(np.float32), b=pair.b.astype(np.float32),
                        info=json.dumps(info))
    return path


def load_pair(path: str | Path) -> BlinkPair:
    with np.load(path, allow_pickle=False) as data:
        info = json.loads(str(data["info"]))
        return BlinkPair(info["ra"], info["dec"], info["size_deg"], data["a"].astype(float),
                         data["b"].astype(float), info["meta_a"], info["meta_b"],
                         name=info.get("name", ""))


def render_pair(pair: BlinkPair, low_pct: float = 20.0, high_pct: float = 99.7,
                mode: str = "asinh", palette: str = "Gray", marker: bool = True, width: int = 560):
    """Labelled pictures of both epochs with one shared brightness scale.

    Returns (image_a, image_b, vmin, vmax).
    """
    a01, b01, vmin, vmax = imaging.stretch_pair(pair.a, pair.b, low_pct, high_pct, mode)
    img_a = imaging.annotate(imaging.to_image(a01, palette, width), pair.label("a"),
                             pair.size_deg, marker, tag="A")
    img_b = imaging.annotate(imaging.to_image(b01, palette, width), pair.label("b"),
                             pair.size_deg, marker, tag="B")
    return img_a, img_b, vmin, vmax
