"""GEMSDOE40 — Euler deconvolution depth-clustering for the DOE GEMS Prize."""

__version__ = "0.1.0"

GRID_HEIGHT = 3730
GRID_WIDTH = 3292
EPSG = 32611
PIXEL_M = 100.0
TRANSFORM = (100.0, 0.0, 243350.0, 0.0, -100.0, 4508550.0)
NODATA_F32 = -3.4028234663852886e38
FOOTPRINT_PX = 5_167_373
TOTAL_PX = GRID_HEIGHT * GRID_WIDTH  # 12_279_160

# Official 1-based GeoTIFF band indices in training_features.tif
# Verified by reading the file's band tags (sha256 4371c82e…).
BANDS = {
    "mag_anom": 1,
    "rtp": 2,
    "tmi_hg": 3,
    "geod_2ndinv": 4,
    "iso_grav_anom_slope": 5,
    "tc": 6,
    "geod_shearrate": 7,
    "geod_dilaterate": 8,
    "tmi_vg": 9,
    "deq_n100a15": 10,
    "iso_grav_anom_vg": 11,
    "det_elev": 12,
    "iso_grav_anom": 13,
    "tmi": 14,
    "depth_to_base_surf": 15,
    "ieq_n100a15": 16,
    "cond_surf": 17,
    "iso_grav_anom_hg": 18,
    "det_elev_slope": 19,
}
