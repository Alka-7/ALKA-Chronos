"""Command-line blink comparison: find a pair, align it, save pictures and a small example file.

Examples (run from the repository root):

    python scripts/blink.py                                   # north ecliptic pole, band 1
    python scripts/blink.py --ra 270 --dec 66.56 --band 4 --size 20 --name nep_band4
    python scripts/blink.py --list --band 2                   # only print the observation inventory

Outputs:
    examples/<name>.npz               small aligned pair, loadable by the app without network
    data/inventory/<name>.csv         every exposure found for this target and band (not committed)
    figures/<name>_blink.gif          blink animation
    figures/<name>_side_by_side.png   both epochs next to each other
    figures/<name>_difference.png     later minus earlier
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from chronos import archive, imaging, pipeline  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--ra", type=float, default=270.0, help="right ascension, degrees")
    parser.add_argument("--dec", type=float, default=66.56, help="declination, degrees")
    parser.add_argument("--band", type=int, default=1, choices=range(1, 7), help="SPHEREx detector 1-6")
    parser.add_argument("--size", type=float, default=15.0, help="field of view, arcmin")
    parser.add_argument("--min-gap", type=float, default=30.0, help="minimum days between epochs")
    parser.add_argument("--max-gap", type=float, default=float("inf"), help="maximum days between epochs")
    parser.add_argument("--pair", type=int, default=0, help="which ranked pair to use (0 = longest baseline)")
    parser.add_argument("--name", default="example", help="prefix for output files")
    parser.add_argument("--list", action="store_true", help="print the inventory and pairs, then stop")
    args = parser.parse_args()

    size_deg = args.size / 60.0
    print(f"Searching IRSA for SPHEREx band {args.band} exposures at RA {args.ra}, Dec {args.dec} ...")
    obs = archive.search(args.ra, args.dec, args.band)
    if obs.empty:
        sys.exit("No SPHEREx exposures cover this position in this band.")
    print(f"{len(obs)} exposures from {obs['date'].iloc[0]} to {obs['date'].iloc[-1]} "
          f"({', '.join(sorted(obs['collection'].unique()))})")

    examples = ROOT / "examples"
    examples.mkdir(exist_ok=True)
    inventory = ROOT / "data" / "inventory" / f"{args.name}.csv"
    inventory.parent.mkdir(parents=True, exist_ok=True)
    obs.round({"x_est": 1, "y_est": 1, "wave_um": 4, "bandwidth_um": 4, "mjd": 5}).to_csv(inventory, index=False)
    print("Inventory saved to", inventory.relative_to(ROOT))

    margin = size_deg * 1.5 / 2 * 3600 / archive.PIXEL_ARCSEC + 20
    pairs = None
    for tolerance in (0.25, 0.5, 1.0):
        pairs = archive.find_pairs(obs, args.min_gap, args.max_gap, tolerance, margin)
        if not pairs.empty:
            break
    if pairs is None or pairs.empty:
        sys.exit("No comparable pair found: no two exposures see the target at a similar wavelength.")
    print(f"\nComparable pairs (wavelengths within {tolerance:.2f} of the bandwidth):")
    for i, row in pairs.head(10).iterrows():
        print(f"  [{i}] {obs.loc[row.a, 'date']}  ->  {obs.loc[row.b, 'date']}   "
              f"{row.gap_days:6.1f} days   {row.wave_um:.3f} um (diff {row.dwave_um * 1000:.2f} nm)")
    if args.list:
        return

    chosen = pairs.iloc[args.pair]
    print(f"\nDownloading and aligning pair [{args.pair}] ...")
    pair = pipeline.build_pair(args.ra, args.dec, size_deg, obs.loc[chosen.a, "access_url"],
                               obs.loc[chosen.b, "access_url"], cache_dir=ROOT / "data" / "cutouts",
                               name=args.name)
    for tag, meta in (("A", pair.meta_a), ("B", pair.meta_b)):
        print(f"  {tag}: {meta['date']}  {meta['obs_id']}  {meta['wave_um']:.4f} um "
              f"(bandwidth {meta['bandwidth_um']:.4f})  exposure {meta['exptime']:.0f} s")
    dx, dy, peak = imaging.alignment_offset(pair.a, pair.b)
    print(f"  Alignment check: best match at shift ({dx}, {dy}) pixels, correlation {peak:.2f}"
          f"  -> {'OK' if (dx, dy) == (0, 0) else 'CHECK THIS PAIR'}")

    figures = ROOT / "figures"
    figures.mkdir(exist_ok=True)
    img_a, img_b, vmin, vmax = pipeline.render_pair(pair)
    (figures / f"{args.name}_blink.gif").write_bytes(imaging.make_gif([img_a, img_b]))
    imaging.side_by_side(img_a, img_b).save(figures / f"{args.name}_side_by_side.png")
    imaging.diff_image(imaging.difference(pair.a, pair.b), vmax).save(figures / f"{args.name}_difference.png")
    saved = pipeline.save_pair(pair, examples / f"{args.name}.npz")
    print(f"\nSaved {saved.relative_to(ROOT)} ({saved.stat().st_size / 1024:.0f} kB) and figures/{args.name}_*.")


if __name__ == "__main__":
    main()
