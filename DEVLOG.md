# Development log

All entries are new hackathon work. No pre-existing research code is used.
Code was written with an AI coding assistant and tested against the live archive.

## 2026-11-14

**Setup.** Repository created, folder structure laid out (`src`, `scripts`, `examples`, `figures`,
`tests`), dependencies listed in `requirements.txt`.

**Data discovery.**
- IRSA's SPHEREx collections are currently `spherex_dr1`, `spherex_dr1_deep`, `spherex_qr3` and
  `spherex_qr3_deep` (the older `spherex_qr2` names in some tutorials no longer return data). The code
  reads the collection list from the archive instead of hard-coding it.
- At the north ecliptic pole (RA 270.0, Dec +66.56) there are 20,898 exposures in DR1 deep alone,
  spanning 2025-04-24 to 2026-05-14; QR3 extends this to 2026-10-04.
- Querying all bands at once took 38 s. Filtering by one detector's wavelength takes about 5 s.

**FITS inspection.** Each spectral image has IMAGE, FLAGS, VARIANCE, ZODI, EPSF and WCS-WAVE
extensions. Units are MJy/sr. The sky WCS is TAN with SIP distortion. Wavelength per pixel comes from
the alternate "W" coordinate system backed by the WCS-WAVE lookup table.

**Decision: use cutouts.** Full files are 66 MB. Appending `?center=RA,DEC&size=DEG` to the file URL
returns a cutout of about 4 MB with all extensions.

**Decision: match wavelength, not just band.** Because of the linear variable filters, two exposures
in the same band can see the target at very different wavelengths. The footprint polygon of each
exposure gives the target's detector position to about 15 pixels, which gives its wavelength. Pairs
must agree to within a quarter of the bandwidth.

**Problem: wavelength lookup disagreed by 0.4%.** `pixel_to_world` on the "W" system applied the sky
distortion terms. Fixed by removing the SIP keywords before building the wavelength WCS.

**Problem: archive returns 503 "Server busy, retry later".** Intermittent on the cutout service. Added
retries that honour the `Retry-After` header, and an on-disk cache.

**Problem: bright stripe in one epoch.** The first pair offered at 1.07 um showed a strong gradient in
one image: the 1.083 um helium airglow line. Pairs within 1.5 bandwidths of that line are now skipped.

**Alignment.** Both epochs are reprojected onto a north-up TAN grid at 6.15 arcsec per pixel. Cutouts
are requested 1.5x larger than the field so rotated detector footprints still fill the frame. An
automatic check cross-correlates the two aligned images and expects the best match at zero shift.

**Validation.**
- North ecliptic pole, band 1, 2025-04-25 vs 2026-10-04: star fields coincide (correlation 0.94).
- Barnard's Star, band 1, 2025-09-01 vs 2026-08-25: measured shift 11.1 arcsec north, 1.8 arcsec west;
  expected from proper motion 10.2 arcsec north, 0.8 arcsec west. The core is saturated, which limits
  the precision of the measurement.

**Interface.** Streamlit app with live archive search, bundled offline examples, in-browser blink
widget, side-by-side view, difference image, observation details and a limitations page.

## To do

- Deploy the app publicly and test the link without university login.
- Confirm the challenge's Resources list and the required submission format.
- Record the demo video or slides.
