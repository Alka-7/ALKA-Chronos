"""ALKA Chronos - SPHEREx Sky Blink Explorer (Streamlit app).

Run locally:   streamlit run app.py
"""
from __future__ import annotations

import base64
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import streamlit as st
import streamlit.components.v1 as components

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src"))

from chronos import archive, imaging, pipeline  # noqa: E402

st.set_page_config(page_title="ALKA Chronos | SPHEREx Sky Blink Explorer", page_icon="🔭", layout="wide")

PRESETS = {
    "North ecliptic pole deep field": (270.0, 66.56, 15.0, 1),
    "Barnard's Star (fast-moving star)": (269.447, 4.768, 10.0, 1),
    "Custom coordinates": None,
}
EXAMPLES = {
    "Barnard's Star (moves in 1 year)": "barnard.npz",
    "North ecliptic pole (527 days)": "nep_deep_band1.npz",
}


# --------------------------------------------------------------------------------------
# Cached data access
# --------------------------------------------------------------------------------------
@st.cache_data(ttl=6 * 3600, show_spinner=False)
def cached_search(ra: float, dec: float, band: int) -> pd.DataFrame:
    return archive.search(ra, dec, band)


@st.cache_data(ttl=24 * 3600, show_spinner=False, max_entries=40)
def cached_pair(ra: float, dec: float, size_deg: float, url_a: str, url_b: str) -> pipeline.BlinkPair:
    return pipeline.build_pair(ra, dec, size_deg, url_a, url_b, cache_dir=ROOT / "data" / "cutouts")


@st.cache_data(show_spinner=False)
def cached_example(filename: str) -> pipeline.BlinkPair:
    return pipeline.load_pair(ROOT / "examples" / filename)


def pick_pairs(obs: pd.DataFrame, size_deg: float, min_gap: float, max_gap: float):
    """Comparable pairs, relaxing the wavelength tolerance only if nothing matches."""
    margin = size_deg * 1.5 / 2 * 3600 / archive.PIXEL_ARCSEC + 20
    for tolerance in (0.25, 0.5, 1.0):
        pairs = archive.find_pairs(obs, min_gap, max_gap, tolerance, margin)
        if not pairs.empty:
            return pairs, tolerance
    return pairs, tolerance


# --------------------------------------------------------------------------------------
# Blink widget: both pictures are stacked and swapped in the browser, so blinking is instant
# --------------------------------------------------------------------------------------
def blink_widget(img_a, img_b, label_a: str, label_b: str, speed_ms: int) -> None:
    def uri(img) -> str:
        return "data:image/png;base64," + base64.b64encode(imaging.png_bytes(img)).decode()

    size = img_a.width
    html = f"""
    <style>
      body {{ margin:0; font-family: system-ui, sans-serif; color:#31333f; }}
      @media (prefers-color-scheme: dark) {{ body {{ color:#fafafa; }} }}
      #stage {{ position:relative; width:min(100%, {size}px); aspect-ratio:1; background:#000; }}
      #stage img {{ position:absolute; inset:0; width:100%; height:100%; image-rendering:pixelated; }}
      .bar {{ display:flex; flex-wrap:wrap; gap:8px; align-items:center; margin-top:10px; font-size:14px; }}
      button {{ font:inherit; padding:6px 12px; border-radius:8px; border:1px solid #8884; background:#8881;
               color:inherit; cursor:pointer; }}
      button.on {{ background:#ffd60a; color:#000; border-color:#ffd60a; }}
      #now {{ font-weight:600; }}
    </style>
    <div id="stage">
      <img id="imgA" src="{uri(img_a)}" alt="Epoch A">
      <img id="imgB" src="{uri(img_b)}" alt="Epoch B" style="visibility:hidden">
    </div>
    <div class="bar">
      <button id="play" class="on">Pause</button>
      <button id="showA">Show A</button>
      <button id="showB">Show B</button>
      <label>Speed <input id="speed" type="range" min="150" max="2000" step="50" value="{speed_ms}"></label>
      <span id="ms"></span>
    </div>
    <div class="bar"><span>Showing:</span> <span id="now"></span></div>
    <script>
      const A = document.getElementById('imgA'), B = document.getElementById('imgB');
      const play = document.getElementById('play'), speed = document.getElementById('speed');
      const labels = {{A: {json.dumps(label_a)}, B: {json.dumps(label_b)}}};
      let current = 'A', timer = null;
      function show(which) {{
        current = which;
        A.style.visibility = which === 'A' ? 'visible' : 'hidden';
        B.style.visibility = which === 'B' ? 'visible' : 'hidden';
        document.getElementById('now').textContent = which + ' - ' + labels[which];
      }}
      function start() {{
        stop();
        timer = setInterval(() => show(current === 'A' ? 'B' : 'A'), Number(speed.value));
        play.textContent = 'Pause'; play.classList.add('on');
      }}
      function stop() {{
        if (timer) clearInterval(timer);
        timer = null; play.textContent = 'Blink'; play.classList.remove('on');
      }}
      play.onclick = () => timer ? stop() : start();
      document.getElementById('showA').onclick = () => {{ stop(); show('A'); }};
      document.getElementById('showB').onclick = () => {{ stop(); show('B'); }};
      speed.oninput = () => {{
        document.getElementById('ms').textContent = (speed.value / 1000).toFixed(2) + ' s per image';
        if (timer) start();
      }};
      document.addEventListener('keydown', e => {{
        if (e.code === 'Space') {{ e.preventDefault(); stop(); show(current === 'A' ? 'B' : 'A'); }}
      }});
      show('A'); speed.oninput(); start();
    </script>
    """
    components.html(html, height=size + 110)


# --------------------------------------------------------------------------------------
# Sidebar: choose what to look at
# --------------------------------------------------------------------------------------
st.title("ALKA Chronos")
st.caption("SPHEREx Sky Blink Explorer: compare two observations of the same patch of sky, taken at "
           "different times, and blink between them to spot anything that moved or changed brightness.")

with st.sidebar:
    st.header("1. Choose a target")
    source = st.radio("Data source", ["Live NASA/IPAC archive", "Bundled example (no network)"],
                      help="Bundled examples are small aligned pairs stored with the app. "
                           "They load instantly and work offline.")

    if source == "Bundled example (no network)":
        example_name = st.selectbox("Example", list(EXAMPLES))
        searched = None
    else:
        preset = st.selectbox("Target", list(PRESETS))
        default = PRESETS[preset] or (270.0, 66.56, 15.0, 1)
        with st.form("search"):
            ra = st.number_input("Right ascension (degrees)", 0.0, 360.0, float(default[0]), format="%.4f")
            dec = st.number_input("Declination (degrees)", -90.0, 90.0, float(default[1]), format="%.4f")
            size_arcmin = st.slider("Field of view (arcmin)", 6.0, 30.0, float(default[2]), 1.0)
            band = st.selectbox("SPHEREx band", list(archive.DETECTOR_RANGE_UM), index=default[3] - 1,
                                format_func=archive.detector_label)
            gap = st.slider("Days between the two observations", 1, 600, (30, 600))
            submitted = st.form_submit_button("Search SPHEREx", type="primary")
        if submitted:
            st.session_state["query"] = (ra, dec, size_arcmin, band, gap)
        searched = st.session_state.get("query")

    st.header("2. Display")
    mode = st.selectbox("Brightness scaling", ["asinh", "sqrt", "linear"],
                        help="asinh shows faint and bright objects together.")
    high_pct = st.slider("Bright limit (percentile)", 95.0, 100.0, 99.7, 0.1)
    palette = st.selectbox("Colours", imaging.PALETTES)
    marker = st.checkbox("Mark the target position", True)

# --------------------------------------------------------------------------------------
# Get the pair
# --------------------------------------------------------------------------------------
pair = None
tolerance_used = None
if source == "Bundled example (no network)":
    pair = cached_example(EXAMPLES[example_name])
elif searched is None:
    st.info("Pick a target in the sidebar and press **Search SPHEREx**, or switch the data source to "
            "**Bundled example** to try the tool instantly.")
else:
    ra, dec, size_arcmin, band, gap = searched
    size_deg = size_arcmin / 60.0
    try:
        with st.spinner("Searching the NASA/IPAC Infrared Science Archive ..."):
            obs = cached_search(ra, dec, band)
    except Exception as error:                                    # network or archive problem
        st.error(f"Could not reach the archive: {error}")
        st.info("Switch the data source to **Bundled example** to keep exploring offline.")
        st.stop()

    if obs.empty:
        st.warning("SPHEREx has no released exposure covering this position in this band. "
                   "Try another band or position.")
        st.stop()

    pairs, tolerance_used = pick_pairs(obs, size_deg, gap[0], gap[1])
    c1, c2, c3 = st.columns(3)
    c1.metric("Exposures covering the target", f"{len(obs):,}")
    c2.metric("Date range", f"{obs['date'].iloc[0][:7]} to {obs['date'].iloc[-1][:7]}")
    c3.metric("Comparable pairs offered", len(pairs))

    if pairs.empty:
        st.warning("No two exposures in this band see the target at a similar wavelength within the "
                   "chosen time range. Widen the range of days, shrink the field of view, or try another band.")
        st.stop()
    if tolerance_used > 0.25:
        st.warning(f"No closely matched pair exists, so the wavelength tolerance was relaxed to "
                   f"{tolerance_used:.2f} of the bandwidth. Brightness differences are less trustworthy.")

    def describe(i: int) -> str:
        row = pairs.iloc[i]
        return (f"{obs.loc[row.a, 'date']}  →  {obs.loc[row.b, 'date']}   |   {row.gap_days:.0f} days apart"
                f"   |   {row.wave_um:.3f} µm")

    choice = st.selectbox("Pair of observations (longest time gap first)", range(len(pairs)), format_func=describe)
    row = pairs.iloc[choice]
    try:
        with st.spinner("Downloading two cutouts from the archive and aligning them (the archive can be "
                        "busy, this may take up to a minute) ..."):
            pair = cached_pair(ra, dec, size_deg, obs.loc[row.a, "access_url"], obs.loc[row.b, "access_url"])
    except Exception as error:
        st.error(f"Could not download this pair: {error}")
        st.info("Choose another pair, try again in a moment, or use a bundled example.")
        st.stop()

# --------------------------------------------------------------------------------------
# Show it
# --------------------------------------------------------------------------------------
if pair is not None:
    img_a, img_b, vmin, vmax = pipeline.render_pair(pair, high_pct=high_pct, mode=mode,
                                                    palette=palette, marker=marker)
    dx, dy, peak = imaging.alignment_offset(pair.a, pair.b)
    wave_gap = abs(pair.meta_a["wave_um"] - pair.meta_b["wave_um"])

    tab_blink, tab_side, tab_diff, tab_meta, tab_about = st.tabs(
        ["Blink", "Side by side", "Difference", "Observation details", "How to read this"])

    with tab_blink:
        left, right = st.columns([3, 2])
        with left:
            blink_widget(img_a, img_b, pair.label("a"), pair.label("b"), 600)
        with right:
            st.subheader("What to look for")
            st.markdown(
                f"- The two images are **{pair.gap_days:.0f} days apart** and show the same "
                f"{pair.size_deg * 60:.0f} arcmin patch of sky, north up.\n"
                "- Ordinary stars and galaxies **stay put**. Something that **jumps** between frames is "
                "moving across the sky; something that **appears, vanishes or pulses** changed brightness "
                "or is an image artifact.\n"
                "- Press **Space** to flip by hand, or use *Show A* / *Show B*.")
            if (dx, dy) == (0, 0):
                st.success(f"Alignment check passed: the star fields coincide (correlation {peak:.2f}).")
            else:
                st.warning(f"Alignment check: the star fields match best at a shift of ({dx}, {dy}) pixels. "
                           "Treat apparent motion in this pair with caution.")
            st.download_button("Download blink animation (GIF)", imaging.make_gif([img_a, img_b]),
                               file_name="alka_chronos_blink.gif", mime="image/gif")

    with tab_side:
        col_a, col_b = st.columns(2)
        col_a.image(img_a, caption=f"A: {pair.label('a')}")
        col_b.image(img_b, caption=f"B: {pair.label('b')}")

    with tab_diff:
        left, right = st.columns([3, 2])
        diff = imaging.difference(pair.a, pair.b)
        left.image(imaging.diff_image(diff, vmax), caption="Later image minus earlier image")
        right.markdown(
            "**Red** = brighter in the later image, **blue** = fainter, **white** = unchanged.\n\n"
            "A moving object leaves a **blue spot next to a red spot**: it left one place and arrived "
            "at another. Small red/blue specks on ordinary stars are expected, because SPHEREx pixels "
            "are large (6 arcsec) and a star never falls on the pixel grid the same way twice.")

    with tab_meta:
        def link(url: str) -> str:
            return f"[{Path(url).name}]({url})" if url.startswith("http") else url

        rows = {
            "Observation start (UTC)": [m["date"].replace("T", " ")[:19] for m in (pair.meta_a, pair.meta_b)],
            "Observation ID": [m["obs_id"] for m in (pair.meta_a, pair.meta_b)],
            "SPHEREx band (detector)": [str(m["detector"]) for m in (pair.meta_a, pair.meta_b)],
            "Wavelength at target (µm)": [f"{m['wave_um']:.4f}" for m in (pair.meta_a, pair.meta_b)],
            "Bandwidth (µm)": [f"{m['bandwidth_um']:.4f}" for m in (pair.meta_a, pair.meta_b)],
            "Exposure time (s)": [f"{m['exptime']:.0f}" for m in (pair.meta_a, pair.meta_b)],
        }
        st.table(pd.DataFrame(rows, index=["Epoch A", "Epoch B"]).T)
        st.markdown(
            f"- **Target:** RA {pair.ra:.4f}°, Dec {pair.dec:.4f}° | field {pair.size_deg * 60:.0f} arcmin "
            f"| pixel {archive.PIXEL_ARCSEC} arcsec\n"
            f"- **Time between epochs:** {pair.gap_days:.1f} days\n"
            f"- **Wavelength mismatch:** {wave_gap * 1000:.2f} nm "
            f"({wave_gap / pair.meta_a['bandwidth_um']:.0%} of the bandwidth)\n"
            f"- **Shared brightness range:** {vmin:.3f} to {vmax:.3f} MJy/sr after background removal\n"
            f"- **Alignment check:** best shift ({dx}, {dy}) pixels, correlation {peak:.2f}\n"
            f"- **Source files:** {link(pair.meta_a['url'])} and {link(pair.meta_b['url'])}")

    with tab_about:
        st.markdown("""
**What this tool does.** It finds SPHEREx exposures of one sky position in the NASA/IPAC Infrared Science
Archive, picks two taken at different times, puts them on the same north-up sky grid, and lets you blink
between them. This is the classic technique used to discover Pluto: anything in the Solar System shifts
against the fixed background stars.

**How the pair is chosen.** Each SPHEREx detector sits behind a filter whose wavelength changes across
the image, so the same star is seen at a different wavelength depending on where it lands. The tool only
offers pairs from the same detector in which the target lands at nearly the same wavelength.

**How the images are prepared.** The pipeline's zodiacal-light model and any leftover constant background
are subtracted from each image, both are resampled onto one sky grid, and both are shown with one shared
brightness range.

**Limitations - please read before claiming a detection.**
- Away from the centre of the field the two images are at slightly different wavelengths, so small
  brightness differences there are not evidence of real variability.
- SPHEREx pixels are 6 arcsec across and stars are barely resolved, so stars flicker slightly between
  frames even when nothing changed.
- Black pixels inside very bright stars are saturated or masked pixels with no data.
- Band 1 contains a helium glow from Earth's upper atmosphere at 1.083 µm that changes between
  exposures, so the tool does not offer pairs near that wavelength.
- Single bright pixels that appear in only one frame are usually cosmic-ray hits or detector artifacts.
- A real moving object should be confirmed with a third observation and checked against known asteroids.
- A distant planet moves slowly: a body 500 au away shifts by roughly a quarter of a degree over six
  months, almost all of it from Earth's own motion. Use a wide field and epochs months apart.

**Data credit.** SPHEREx spectral images from the NASA/IPAC Infrared Science Archive (IRSA), Caltech/IPAC.
""")
