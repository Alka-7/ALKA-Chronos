"""Offline tests: no network needed. Run with  python -m pytest -q"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.wcs import WCS

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from chronos import archive, imaging, pipeline  # noqa: E402


def test_wavelength_grid_covers_each_band():
    for detector, (lo, hi) in archive.DETECTOR_RANGE_UM.items():
        wave, bandwidth = archive.wavelength_at(detector, 1020, np.array([5, 1020, 2035]))
        assert np.all((wave > lo - 0.05) & (wave < hi + 0.05))
        assert np.all(bandwidth > 0)
        assert abs(wave[0] - wave[2]) > 0.3 * (hi - lo)      # wavelength changes along detector y


def test_estimate_pixel_recovers_known_position():
    # Footprint of a real exposure (2025W18_1B_0034_1, band 1). The target's true pixel,
    # measured from the file's own WCS, is about (1230, 2024).
    region = ("POLYGON ICRS 276.21743651268275 63.56993784853953 268.67183205854195 62.842167435941256 "
              "266.6482426906245 66.21664258512845 275.24843557214314 67.0369115771107 "
              "276.21743651268275 63.56993784853953")
    corners = archive.parse_region(region)[None]
    x, y = archive.estimate_pixel(corners, [271.680102], [64.956396], 270.0, 66.56)
    assert abs(x[0] - 1230) < 30 and abs(y[0] - 2024) < 30


def _obs(rows):
    frame = pd.DataFrame(rows, columns=["mjd", "wave_um", "x_est", "y_est"])
    frame["bandwidth_um"] = 0.02
    return frame


def test_find_pairs_requires_similar_wavelength_and_time_gap():
    obs = _obs([
        (60800.0, 0.900, 1000, 1000),
        (60801.0, 0.901, 1000, 1000),    # same wavelength, too close in time
        (61100.0, 0.902, 1000, 1000),    # same wavelength, 300 days later  -> pairs with rows 0 and 1
        (61100.0, 0.990, 1000, 1000),    # different wavelength
        (61200.0, 0.900, 10, 1000),      # too close to the detector edge
    ])
    pairs = archive.find_pairs(obs, min_gap_days=30)
    assert sorted(zip(pairs["a"], pairs["b"])) == [(0, 2), (1, 2)]
    assert (pairs["gap_days"] >= 30).all()
    assert (pairs["dwave_um"] <= 0.25 * 0.02).all()


def test_find_pairs_skips_helium_airglow():
    obs = _obs([(60800.0, 1.083, 1000, 1000), (61100.0, 1.083, 1000, 1000)])
    assert archive.find_pairs(obs).empty
    assert len(archive.find_pairs(obs, avoid_airglow=False)) == 1


def test_target_frame_is_north_up_and_centred():
    wcs, shape = imaging.target_frame(150.0, 20.0, 0.25)
    n = shape[0]
    ra, dec = wcs.all_pix2world((n - 1) / 2, (n - 1) / 2, 0)
    assert abs(ra - 150.0) < 1e-6 and abs(dec - 20.0) < 1e-6
    _, dec_up = wcs.all_pix2world((n - 1) / 2, n - 1, 0)
    ra_right, _ = wcs.all_pix2world(n - 1, (n - 1) / 2, 0)
    assert dec_up > 20.0            # larger row index = north
    assert ra_right < 150.0         # larger column index = west, so east is left


def _fake_epoch(shift=(0.0, 0.0), date="2025-01-01T00:00:00", mjd=60676.0):
    """A synthetic star field on a rotated detector grid."""
    rng = np.random.default_rng(1)
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    wcs.wcs.crval = [150.0, 20.0]
    wcs.wcs.crpix = [100.5, 100.5]
    angle = np.radians(25.0)
    scale = archive.PIXEL_ARCSEC / 3600.0
    wcs.wcs.cd = scale * np.array([[-np.cos(angle), np.sin(angle)], [np.sin(angle), np.cos(angle)]])
    yy, xx = np.mgrid[0:200, 0:200]
    image = np.full((200, 200), 0.3)
    for sx, sy, flux in zip(rng.uniform(20, 180, 40), rng.uniform(20, 180, 40), rng.uniform(1, 8, 40)):
        image += flux * np.exp(-((xx - sx - shift[0]) ** 2 + (yy - sy - shift[1]) ** 2) / (2 * 1.2 ** 2))
    return imaging.Epoch(image=image, zodi=np.full_like(image, 0.3), flags=np.zeros(image.shape, int),
                         wcs=wcs, date=date, mjd=mjd, detector=1, obs_id="fake", exptime=114.0,
                         wave_um=1.0, bandwidth_um=0.03, unit="MJy/sr")


def test_align_puts_two_epochs_on_the_same_grid():
    frame, shape = imaging.target_frame(150.0, 20.0, 0.15)
    a = imaging.align(_fake_epoch(), frame, shape)
    b = imaging.align(_fake_epoch(), frame, shape)
    assert a.shape == shape and np.isfinite(a).all()
    assert abs(np.median(a)) < 0.05                       # background removed
    assert imaging.alignment_offset(a, b)[:2] == (0, 0)


def test_alignment_check_detects_a_shift():
    frame, shape = imaging.target_frame(150.0, 20.0, 0.15)
    a = imaging.align(_fake_epoch(), frame, shape)
    b = np.roll(a, 2, axis=1)
    assert imaging.alignment_offset(a, b)[:2] == (2, 0)


def test_stretch_uses_one_shared_range():
    a = np.array([[0.0, 1.0], [2.0, 3.0]])
    b = a * 2
    a01, b01, vmin, vmax = imaging.stretch_pair(a, b, 0, 100, mode="linear")
    assert (vmin, vmax) == (0.0, 6.0)
    assert a01.max() == 0.5 and b01.max() == 1.0          # same star, twice as bright, looks brighter


def test_bundled_examples_load_and_render(tmp_path):
    for name in ("barnard.npz", "nep_deep_band1.npz"):
        pair = pipeline.load_pair(ROOT / "examples" / name)
        assert pair.a.shape == pair.b.shape and pair.gap_days > 300
        img_a, img_b, vmin, vmax = pipeline.render_pair(pair)
        assert img_a.size == img_b.size and vmax > vmin
        assert len(imaging.make_gif([img_a, img_b])) > 1000
        again = pipeline.load_pair(pipeline.save_pair(pair, tmp_path / name))
        assert np.allclose(np.nan_to_num(again.a), np.nan_to_num(pair.a))
