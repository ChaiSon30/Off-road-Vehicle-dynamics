"""
 Pac coefficients for tire on hard-surface.

    import on_road
    Fy = on_road.Fy(alpha_rad, Fz)

ILLUSTRATION ONLY. Read the note at the bottom before using any number
from this file, and before setting it beside off_road.py.
"""
# ---------------------------------------------------------------------
# READ BEFORE USING

# This set shows the SHAPE of tire behaviour: a linear range, a peak, saturation past it, and a friction budget shared with the longitudinal channel. Those features are real and they transfer. The magnitudes do not. It is example data for a road tire on dry pavement.

# It is NOT a comparison partner for off_road.py. Different tire, different surface, different data source. The two files are separate teaching illustrations, and a difference between them cannot be read as a surface effect.
# ---------------------------------------------------------------------

import magic_formula as mf
from magic_formula import TireParams

PARAMS = TireParams(
    name="175/70R13 passenger radial",
    surface="dry pavement, 190 kPa",

    Fz0=3800.0,          # FNOMIN; the file declares itself valid 190-8550 N

    # lateral -- the channel the series runs on
    pCy1=1.4675,         # C > 1, so the curve peaks and then falls away
    pDy1=0.94002,        # friction at 3800 N
    pDy2=-0.17669,       # friction decays with load
    pEy1=0.0040023,
    pEy2=0.00085719,
    pKy1=12.536,         # file gives -12.536: SAE sign convention, flipped here
    pKy2=1.3856,

    # longitudinal -- carried only so the friction ellipse can be drawn
    pCx1=1.5587,
    pDx1=1.09,
    pDx2=-0.079328,
    pEx1=0.27403,
    pEx2=0.10232,
    pEx3=0.074903,
    pKx1=19.733,
    pKx2=0.093405,
    pKx3=0.12433,
)

# Camber, ply-steer and conicity terms in the file are dropped.

Fy = lambda alpha, Fz: mf.Fy(alpha, Fz, PARAMS)
Fx = lambda kappa, Fz: mf.Fx(kappa, Fz, PARAMS)
Fy_peak = lambda Fz: mf.Fy_peak(Fz, PARAMS)
Fx_peak = lambda Fz: mf.Fx_peak(Fz, PARAMS)
cornering_stiffness = lambda Fz: mf.cornering_stiffness(Fz, PARAMS)
Fy_available = lambda Fx_used, Fz: mf.Fy_available(Fx_used, Fz, PARAMS)
combined_slip = lambda alpha, kappa, Fz: mf.combined_slip(alpha, kappa, Fz, PARAMS)
