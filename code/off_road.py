"""
Pac coefficients for tire on soft soil.

    import off_road
    Fy = off_road.Fy(alpha_rad, Fz)
"""
# ---------------------------------------------------------------------
# READ BEFORE USING

# It is not measured. The source curves come from a terramechanics model, so these coefficients inherit every assumption in that model and in the soil parameter tables behind it.

# The two channels are not a matched pair. Lateral is dry sand on an AT 22x7-12; longitudinal is sandy loam on an AT 23x8-12. A friction ellipse built from both is indicative at best.

# The absence of a peak is partly an input, not a finding. The shear law behind the source curves saturates monotonically by construction. Loose sand is believed to behave this way; dense or compacted soil can show a peak and a lower residual strength.

# It is NOT a comparison partner for on_road.py. Different tire, different surface, different data source.
# ---------------------------------------------------------------------

import magic_formula as mf
from magic_formula import TireParams

# Fitted to curves digitised from Figures 5 and 6 of Saunders, White & Compere (2019), ASME IMECE2019-10682. Those curves are Bekker-Wong terramechanics model output, not measured tire data.
PARAMS = TireParams(
    name="AT 22x7-12 / 23x8-12 on soft soil",
    surface="dry sand (lateral) / sandy loam (longitudinal), 96.5 kPa",

    Fz0=1334.5,          # 300 lbf, the reference load of the source figures

    # lateral -- dry sand
    pCy1=1.0000,         # C = 1: PINNED, not fitted. A free fit ran to a bound
                         # because C and E are degenerate on curves with no
                         # peak. C = 1 gives clean saturation.
    pDy1=0.5978,         # friction at 1334 N
    pDy2=-0.0358,        # decays with load, far more weakly than on tarmac
    pEy1=0.5927,
    pEy2=-0.2469,
    pKy1=2.0889,         # cornering stiffness / Fz, in 1/rad
    pKy2=1.0,            # unused: linear_k is True

    # longitudinal -- sandy loam, and it is DRAWBAR PULL, not Fx
    pCx1=1.0000,
    pDx1=0.4441,
    pDx2=-0.0300,
    pEx1=0.5258,
    pEx2=-0.0826,
    pKx1=1.7893,
    SVx=-6.773,          # drawbar pull is net of motion resistance, so the curve does not pass through the origin

    linear_k=True,       # the source loads show no stiffness rolloff, so pK2 is not identifiable; a proportional law is honest
)

Fy = lambda alpha, Fz: mf.Fy(alpha, Fz, PARAMS)
Fx = lambda kappa, Fz: mf.Fx(kappa, Fz, PARAMS)
Fy_peak = lambda Fz: mf.Fy_peak(Fz, PARAMS)
Fx_peak = lambda Fz: mf.Fx_peak(Fz, PARAMS)
cornering_stiffness = lambda Fz: mf.cornering_stiffness(Fz, PARAMS)
Fy_available = lambda Fx_used, Fz: mf.Fy_available(Fx_used, Fz, PARAMS)
combined_slip = lambda alpha, kappa, Fz: mf.combined_slip(alpha, kappa, Fz, PARAMS)


