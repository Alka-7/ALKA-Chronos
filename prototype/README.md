# ALKA Chronos — SPHEREx Sky Blink Explorer

**NASA International Space Apps Challenge 2026**

**Challenge:** Planet X and SPHEREx

## About the Project

ALKA Chronos is an astronomical image-comparison prototype inspired by NASA's SPHEREx mission.

The application demonstrates how astronomical images taken at different epochs can be compared to identify moving or changing sources.

## Features

- Visualization of astronomical FITS images
- Side-by-side comparison of two observation epochs
- Interactive blink viewer
- Difference-image visualization
- FITS header metadata inspection

## Technologies

Python, NumPy, Astropy, Matplotlib, and Streamlit.

## Running the Application

Install the required packages:

```bash
pip install -r requirements.txt
```

Generate the demonstration images:

```bash
python scripts/create_demo.py
```

Launch the application:

```bash
streamlit run app.py
```

## Data and Scientific Limitations

The current prototype uses simulated FITS images containing stationary background stars and a moving astronomical source.

These demonstration images are not actual SPHEREx observations.

During development, NASA/IPAC IRSA's SPHEREx metadata tables were explored to identify observation metadata and data-access information.

Future versions will integrate real SPHEREx images, WCS alignment, wavelength selection, and automated change detection.

## Project Goal

Our goal is to develop an accessible platform for exploring multi-epoch infrared sky observations and identifying astronomical changes.
