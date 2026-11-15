
import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from pathlib import Path

st.set_page_config(
    page_title="ALKA Chronos",
    layout="wide"
)

st.title("ALKA Chronos")
st.subheader("SPHEREx Sky Blink Explorer")

st.info(
    "Prototype demonstration using simulated FITS images. "
    "Real SPHEREx observations are not yet integrated."
)

# Locate data relative to this application
BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / "data"

def load_fits(path):
    with fits.open(path) as hdul:
        return (
            np.asarray(hdul[0].data, dtype=float),
            hdul[0].header.copy()
        )

def display_image(image, title, vmin, vmax):
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.imshow(
        image,
        origin="lower",
        cmap="inferno",
        vmin=vmin,
        vmax=vmax
    )
    ax.set_title(title)
    ax.set_xlabel("X pixel")
    ax.set_ylabel("Y pixel")
    st.pyplot(fig)
    plt.close(fig)

file1 = DATA_DIR / "demo_epoch_1.fits"
file2 = DATA_DIR / "demo_epoch_2.fits"

if not file1.exists() or not file2.exists():
    st.error(
        "Demo FITS files not found. "
        "Run the create_demo.py script first."
    )
    st.stop()

image1, header1 = load_fits(file1)
image2, header2 = load_fits(file2)

# Consistent brightness scale for both images
combined = np.concatenate(
    [image1.ravel(), image2.ravel()]
)
vmin, vmax = np.percentile(combined, [5, 99.5])

tab1, tab2, tab3 = st.tabs([
    "Side-by-side",
    "Blink Viewer",
    "Difference Image"
])

with tab1:
    col1, col2 = st.columns(2)

    with col1:
        display_image(image1, "Epoch 1", vmin, vmax)

    with col2:
        display_image(image2, "Epoch 2", vmin, vmax)

with tab2:
    epoch = st.radio(
        "Select observation",
        ["Epoch 1", "Epoch 2"],
        horizontal=True
    )

    selected = image1 if epoch == "Epoch 1" else image2
    display_image(selected, epoch, vmin, vmax)

    st.caption(
        "Switch between epochs to identify changes."
    )

with tab3:
    difference = image2 - image1

    fig, ax = plt.subplots(figsize=(6, 5))

    limit = np.percentile(
        np.abs(difference), 99.5
    )

    ax.imshow(
        difference,
        origin="lower",
        cmap="coolwarm",
        vmin=-limit,
        vmax=limit
    )

    ax.set_title("Epoch 2 − Epoch 1")
    ax.set_xlabel("X pixel")
    ax.set_ylabel("Y pixel")

    st.pyplot(fig)
    plt.close(fig)

with st.expander("Observation Metadata"):
    st.write("Epoch 1:", dict(header1))
    st.write("Epoch 2:", dict(header2))
