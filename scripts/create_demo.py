
import numpy as np
from astropy.io import fits
from pathlib import Path

# Create the data folder if it does not exist
Path("data").mkdir(exist_ok=True)

# Image dimensions: 160 x 160 pixels
size = 160
y, x = np.mgrid[:size, :size]

# Fixed background stars
stars = [
    (25, 35, 140),
    (55, 120, 190),
    (120, 30, 160),
    (135, 135, 120),
    (85, 85, 80),
]

def create_image(moving_x, seed):
    rng = np.random.default_rng(seed)
    image = rng.normal(10, 1, (size, size))

    # Add stationary stars
    for sx, sy, brightness in stars:
        image += brightness * np.exp(
            -((x - sx)**2 + (y - sy)**2)
            / (2 * 2.2**2)
        )

    # Add a moving source
    image += 220 * np.exp(
        -((x - moving_x)**2 + (y - 75)**2)
        / (2 * 2.2**2)
    )

    return image.astype("float32")

# Two observations of the same simulated sky
for epoch, position in [(1, 65), (2, 73)]:

    image = create_image(position, seed=42)

    header = fits.Header()
    header["OBJECT"] = "SIMULATED SOURCE"
    header["EPOCH"] = epoch
    header["SIMULATE"] = True
    header["BUNIT"] = "arbitrary"

    filename = f"data/demo_epoch_{epoch}.fits"

    fits.writeto(
        filename,
        image,
        header,
        overwrite=True
    )

    print(f"Created: {filename}")
