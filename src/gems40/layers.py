"""Official band inventory of the competition feature stack (``training_features.tif``).

Every entry below is the organiser's own band name, band number and description, read directly
from the GeoTIFF band tags of the pinned competition file (sha256
``4371c82e3b8339b807bdffcf4ef59a225520fe2988d521be208ae33743123bc5``, 418,912,844 bytes, 19
float32 bands).  Nothing here is inferred: ``scripts/inspect_data.py`` re-derives the whole table
from the bytes and writes ``data/evidence/band_inventory.json``.

Official source of the file (login-walled):
    https://www.drivendata.org/competitions/306/competition-doe-gems/data/
Byte-pinned public mirror used here (hash-verified):
    data/bridge/gems-geodawn-numerical-features.tif.part-* in github.com/buffedlizard55-lab/GEMSDOE
"""

from __future__ import annotations

BANDS: dict[str, dict] = {
    "mag_anom":          dict(band=1,  category="magnetic_data",   band_name="mag_anom",
                              description="Magnetic anomaly - deviation from expected Earth's magnetic field"),
    "rtp":               dict(band=2,  category="magnetic_data",   band_name="rtp",
                              description="Reduced to pole magnetic data - magnetic anomaly corrected for latitude effects"),
    "tmi_hg":            dict(band=3,  category="magnetic_data",   band_name="tmi_hg",
                              description="Total magnetic intensity horizontal gradient - rate of change in horizontal direction"),
    "geod_2ndinv":       dict(band=4,  category="geodetic_strain", band_name="geod_2ndinv",
                              description="Geodetic second invariant - measure of strain rate tensor magnitude"),
    "iso_grav_anom_slope": dict(band=5, category="gravity_data",   band_name="iso_grav_anom_slope",
                                description="Isostatic gravity anomaly slope - gradient of gravity after isostatic correction"),
    "tc":                dict(band=6,  category="magnetic_data",   band_name="tc",
                              description="Tilt angle or total curvature - magnetic field derivative for edge detection"),
    "geod_shearrate":    dict(band=7,  category="geodetic_strain", band_name="geod_shearrate",
                              description="Geodetic shear rate - rate of angular deformation from GPS/InSAR"),
    "geod_dilaterate":   dict(band=8,  category="geodetic_strain", band_name="geod_dilaterate",
                              description="Geodetic dilatation rate - rate of volumetric strain (expansion/contraction)"),
    "tmi_vg":            dict(band=9,  category="magnetic_data",   band_name="tmi_vg",
                              description="Total magnetic intensity vertical gradient - rate of change in vertical direction"),
    "deq_n100a15":       dict(band=10, category="seismic",         band_name="deq_n100a15",
                              description="Distance to earthquake (n=100km radius, a=15 deg azimuth parameters)"),
    "iso_grav_anom_vg":  dict(band=11, category="gravity_data",    band_name="iso_grav_anom_vg",
                              description="Isostatic gravity anomaly vertical gradient - vertical rate of change"),
    "det_elev":          dict(band=12, category="topographic",     band_name="det_elev",
                              description="Detrended elevation - topography with regional trends removed"),
    "iso_grav_anom":     dict(band=13, category="gravity_data",    band_name="iso_grav_anom",
                              description="Isostatic gravity anomaly - gravity after compensating for topographic mass"),
    "tmi":               dict(band=14, category="magnetic_data",   band_name="tmi",
                              description="Total magnetic intensity - total strength of magnetic field"),
    "depth_to_base_surf": dict(band=15, category="subsurface",    band_name="depth_to_base_surf",
                               description="Depth to basement surface - thickness of sedimentary cover"),
    "ieq_n100a15":       dict(band=16, category="seismic",         band_name="ieq_n100a15",
                              description="Earthquake intensity or density (n=100km radius, a=15 deg parameters)"),
    "cond_surf":         dict(band=17, category="subsurface",      band_name="cond_surf",
                              description="Conductivity surface - electrical conductivity of subsurface"),
    "iso_grav_anom_hg":  dict(band=18, category="gravity_data",    band_name="iso_grav_anom_hg",
                              description="Isostatic gravity anomaly horizontal gradient - horizontal rate of change"),
    "det_elev_slope":    dict(band=19, category="topographic",     band_name="det_elev_slope",
                              description="Detrended elevation slope - gradient of elevation after detrending"),
}

#: Layers used as master potential fields for Euler deconvolution in this repository.
#: The three magnetics are three representations of the same airborne survey; the gravity layer
#: is an independent potential field on the same grid.
EULER_FIELDS: tuple[str, ...] = ("rtp", "tmi", "mag_anom", "iso_grav_anom")

#: Independent potential-field *families* (magnetic vs gravity) used for cross-family conjunction.
FIELD_FAMILIES: dict[str, str] = {
    "rtp": "magnetic", "tmi": "magnetic", "mag_anom": "magnetic",
    "iso_grav_anom": "gravity",
}
