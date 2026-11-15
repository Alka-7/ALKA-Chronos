"""Read SPHEREx cutouts, put two epochs on the same sky grid, and render them for blinking.

Pipeline for one pair:

    load_epoch(path)  ->  Epoch (image, zodiacal model, sky coordinates, date, wavelength)
    target_frame(...) ->  a north-up sky grid centred on the target
    align(epoch, ...) ->  the epoch resampled onto that grid, background removed
    stretch_pair(...) ->  both images scaled with ONE shared brightness range
    render(...)       ->  labelled picture;  make_gif(...) -> blinking animation
"""
from __future__ import annotations

import io
import warnings
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS, FITSFixedWarning
from PIL import Image, ImageDraw, ImageFont

from .archive import PIXEL_ARCSEC

_SIP_PREFIXES = ("A_", "B_", "AP_", "BP_")


@dataclass
class Epoch:
    """One SPHEREx exposure, cut out around the target."""

    image: np.ndarray            # surface brightness, MJy/sr
    zodi: np.ndarray             # zodiacal-light model for the same pixels, MJy/sr
    flags: np.ndarray            # pipeline flag bits per pixel
    wcs: WCS                     # pixel -> sky mapping (with distortion)
    date: str                    # start of exposure, UTC
    mjd: float
    detector: int
    obs_id: str
    exptime: float               # seconds
    wave_um: float               # wavelength at the target position
    bandwidth_um: float
    unit: str
    source: str = ""             # file name or URL the data came from
    extra: dict = field(default_factory=dict)

    @property
    def label(self) -> str:
        return f"{self.date[:16].replace('T', ' ')} UTC | {self.wave_um:.3f} um"


def load_epoch(path: str | Path, ra: float, dec: float, source: str = "") -> Epoch:
    """Read a SPHEREx spectral-image cutout and look up the wavelength at (ra, dec)."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FITSFixedWarning)
        with fits.open(path) as hdul:
            header = hdul["IMAGE"].header
            image = np.asarray(hdul["IMAGE"].data, dtype=float)
            zodi = np.asarray(hdul["ZODI"].data, dtype=float) if "ZODI" in hdul else np.zeros_like(image)
            flags = np.asarray(hdul["FLAGS"].data) if "FLAGS" in hdul else np.zeros(image.shape, int)
            wcs = WCS(header)
            x, y = wcs.all_world2pix(ra, dec, 0)

            # The alternate "W" coordinate system maps pixel -> (wavelength, bandwidth) through
            # the WCS-WAVE lookup table. It must not use the sky distortion terms.
            spectral_header = header.copy()
            for key in [k for k in spectral_header if k.startswith(_SIP_PREFIXES)]:
                del spectral_header[key]
            wave_wcs = WCS(spectral_header, fobj=hdul, key="W")
            wave, bandwidth = wave_wcs.wcs_pix2world([[float(x), float(y)]], 0)[0]

    return Epoch(
        image=image, zodi=zodi, flags=flags, wcs=wcs,
        date=str(header.get("DATE-OBS", "")),
        mjd=float(header.get("MJD-OBS", np.nan)),
        detector=int(header.get("DETECTOR", 0)),
        obs_id=str(header.get("OBSID", "")),
        exptime=float(header.get("XPOSURE", header.get("EXPTIME", np.nan))),
        wave_um=float(wave), bandwidth_um=float(bandwidth),
        unit=str(header.get("BUNIT", "MJy/sr")),
        source=source or Path(path).name,
        extra={"psf_fwhm_arcsec": header.get("PSF_FWHM"), "target_pixel": (float(x), float(y))},
    )


def target_frame(ra: float, dec: float, size_deg: float, pixel_arcsec: float = PIXEL_ARCSEC):
    """A north-up, east-left sky grid centred on the target. Returns (WCS, (ny, nx))."""
    n = max(int(round(size_deg * 3600.0 / pixel_arcsec)), 8)
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    wcs.wcs.crval = [ra, dec]
    wcs.wcs.crpix = [(n + 1) / 2.0, (n + 1) / 2.0]
    wcs.wcs.cdelt = [-pixel_arcsec / 3600.0, pixel_arcsec / 3600.0]
    wcs.wcs.cunit = ["deg", "deg"]
    return wcs, (n, n)


def align(epoch: Epoch, frame_wcs: WCS, shape: tuple[int, int],
          subtract_zodi: bool = True, subtract_median: bool = True) -> np.ndarray:
    """Resample an epoch onto the shared sky grid. Pixels with no data are NaN.

    The zodiacal light changes a lot between dates, so by default the pipeline's zodiacal
    model and any remaining constant background are removed before comparing.
    """
    from reproject import reproject_interp

    data = epoch.image - epoch.zodi if subtract_zodi else epoch.image.copy()
    out, footprint = reproject_interp((data, epoch.wcs), frame_wcs, shape_out=shape)
    out[footprint == 0] = np.nan
    if subtract_median and np.isfinite(out).any():
        _, median, _ = sigma_clipped_stats(out, sigma=3.0, maxiters=5)
        out = out - median
    return out


def stretch_pair(a: np.ndarray, b: np.ndarray, low_pct: float = 20.0, high_pct: float = 99.7,
                 mode: str = "asinh"):
    """Scale both images to 0..1 with one shared brightness range.

    Returns (a01, b01, vmin, vmax). Sharing the range is essential: if each image were
    scaled on its own, an unchanged star could look brighter or fainter between epochs.
    """
    both = np.concatenate([a[np.isfinite(a)], b[np.isfinite(b)]])
    if both.size == 0:
        return np.zeros_like(a), np.zeros_like(b), 0.0, 1.0
    vmin, vmax = np.percentile(both, [low_pct, high_pct])
    if not vmax > vmin:
        vmax = vmin + 1e-6

    def scale(arr):
        x = np.clip((arr - vmin) / (vmax - vmin), 0.0, 1.0)
        if mode == "asinh":
            x = np.arcsinh(10.0 * x) / np.arcsinh(10.0)
        elif mode == "sqrt":
            x = np.sqrt(x)
        return np.where(np.isfinite(arr), x, 0.0)

    return scale(a), scale(b), float(vmin), float(vmax)


def difference(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Later epoch minus earlier epoch (NaN where either has no data)."""
    return b - a


# --------------------------------------------------------------------------------------
# Rendering
# --------------------------------------------------------------------------------------
_PALETTES = {
    "Gray": lambda x: np.stack([x, x, x], axis=-1),
    "Inverted gray": lambda x: np.stack([1 - x, 1 - x, 1 - x], axis=-1),
    "Heat": lambda x: np.stack([np.clip(3 * x, 0, 1), np.clip(3 * x - 1, 0, 1),
                                np.clip(3 * x - 2, 0, 1)], axis=-1),
}
PALETTES = list(_PALETTES)


def _font(size: int):
    try:
        return ImageFont.load_default(size=size)
    except TypeError:                      # older Pillow
        return ImageFont.load_default()


def to_image(arr01: np.ndarray, palette: str = "Gray", width: int = 560) -> Image.Image:
    """Turn a 0..1 array into a picture. Row 0 of the array is south, so flip for display."""
    rgb = (_PALETTES[palette](np.nan_to_num(arr01)) * 255).astype(np.uint8)
    img = Image.fromarray(rgb[::-1])
    scale = max(1, int(round(width / img.width)))
    return img.resize((img.width * scale, img.height * scale), Image.NEAREST)


def diff_image(diff: np.ndarray, limit: float, width: int = 560) -> Image.Image:
    """Blue = fainter in the later epoch, red = brighter, white = unchanged."""
    x = np.clip(np.nan_to_num(diff) / max(limit, 1e-12), -1.0, 1.0)
    rgb = np.stack([1 - np.clip(-x, 0, 1), 1 - np.abs(x), 1 - np.clip(x, 0, 1)], axis=-1)
    rgb[~np.isfinite(diff)] = 0.5
    img = Image.fromarray((rgb[::-1] * 255).astype(np.uint8))
    scale = max(1, int(round(width / img.width)))
    return img.resize((img.width * scale, img.height * scale), Image.NEAREST)


def annotate(img: Image.Image, label: str, size_deg: float, marker: bool = True,
             tag: str = "") -> Image.Image:
    """Add the epoch label, a scale bar, a north/east compass and a target marker."""
    img = img.copy()
    draw = ImageDraw.Draw(img)
    w, h = img.size
    accent, shadow = (255, 214, 10), (0, 0, 0)
    font = _font(max(13, w // 34))

    def text(xy, s, anchor="la"):
        draw.text((xy[0] + 1, xy[1] + 1), s, fill=shadow, font=font, anchor=anchor)
        draw.text(xy, s, fill=accent, font=font, anchor=anchor)

    text((8, 6), f"{tag}  {label}".strip())

    # Scale bar: pick a round length close to a fifth of the field
    field_arcmin = size_deg * 60.0
    bar_arcmin = min((1, 2, 5, 10, 20, 30, 60), key=lambda v: abs(v - field_arcmin / 5.0))
    bar_px = bar_arcmin / field_arcmin * w
    y = h - 14
    draw.line([(10, y), (10 + bar_px, y)], fill=shadow, width=5)
    draw.line([(10, y), (10 + bar_px, y)], fill=accent, width=3)
    text((10, y - 6), f"{bar_arcmin} arcmin", anchor="ld")

    # Compass: north is up, east is left
    cx, cy, arm = w - 22, h - 22, max(26, w // 14)
    draw.line([(cx, cy), (cx, cy - arm)], fill=accent, width=2)
    draw.line([(cx, cy), (cx - arm, cy)], fill=accent, width=2)
    text((cx, cy - arm - 2), "N", anchor="mb")
    text((cx - arm - 4, cy), "E", anchor="rm")

    if marker:                              # open cross-hair so the target itself stays visible
        mx, my, gap, length = w / 2, h / 2, max(8, w // 40), max(14, w // 22)
        for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            draw.line([(mx + dx * gap, my + dy * gap),
                       (mx + dx * (gap + length), my + dy * (gap + length))], fill=accent, width=2)
    return img


def make_gif(frames: list[Image.Image], ms_per_frame: int = 600) -> bytes:
    """An endlessly looping blink animation."""
    buffer = io.BytesIO()
    frames[0].save(buffer, format="GIF", save_all=True, append_images=frames[1:],
                   duration=ms_per_frame, loop=0)
    return buffer.getvalue()


def png_bytes(img: Image.Image) -> bytes:
    buffer = io.BytesIO()
    img.save(buffer, format="PNG")
    return buffer.getvalue()


def side_by_side(left: Image.Image, right: Image.Image, gap: int = 8) -> Image.Image:
    canvas = Image.new("RGB", (left.width + right.width + gap, max(left.height, right.height)), "white")
    canvas.paste(left, (0, 0))
    canvas.paste(right, (left.width + gap, 0))
    return canvas


# --------------------------------------------------------------------------------------
# Alignment check
# --------------------------------------------------------------------------------------
def alignment_offset(a: np.ndarray, b: np.ndarray, max_shift: int = 4):
    """How well do the two aligned images line up?

    Cross-correlates the bright structure of both images over small integer shifts and
    returns (dx, dy, peak) in pixels, where (0, 0) means the star fields coincide.
    A non-zero result on a field of ordinary stars points to an alignment problem.
    """
    def prepare(arr):
        x = np.nan_to_num(arr)
        return np.clip(x - np.percentile(x, 60), 0, np.percentile(x, 99.5))

    pa, pb = prepare(a), prepare(b)
    m = max_shift
    core = pa[m:-m, m:-m]
    best = (0, 0, -np.inf)
    norm = np.sqrt((core ** 2).sum())
    for dy in range(-m, m + 1):
        for dx in range(-m, m + 1):
            window = pb[m + dy: pb.shape[0] - m + dy, m + dx: pb.shape[1] - m + dx]
            score = (core * window).sum() / (norm * np.sqrt((window ** 2).sum()) + 1e-12)
            if score > best[2]:
                best = (dx, dy, float(score))
    return best
