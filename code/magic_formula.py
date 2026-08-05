"""
The Pacejka Magic Formula.

Holds no coefficients. Import on_road or off_road instead.
"""

from dataclasses import dataclass
import numpy as np

DEG = np.pi / 180.0


@dataclass
class TireParams:
    name: str
    surface: str
    source: str

    Fz0: float                   # nominal load [N], the load dfz is measured against

    # lateral -- the channel this series runs on
    pCy1: float                  # shape factor C
    pDy1: float                  # friction at nominal load
    pDy2: float                  # how friction changes with load (negative)
    pEy1: float                  # curvature into saturation
    pEy2: float
    pKy1: float                  # cornering stiffness magnitude
    pKy2: float                  # load where stiffness stops rising, in units of Fz0

    # longitudinal -- carried so the friction ellipse can be drawn
    pCx1: float
    pDx1: float
    pDx2: float
    pEx1: float
    pEx2: float
    pKx1: float

    pEx3: float = 0.0
    pKx2: float = 0.0
    pKx3: float = 0.0
    SVx: float = 0.0             # vertical offset, non-zero only for drawbar pull
    linear_k: bool = False       # True -> K = pK1*Fz instead of the saturating law


def magic_formula(x, B, C, D, E):
    Bx = B * x      # x is slip angle [rad] or slip ratio [-]
    return D * np.sin(C * np.arctan(Bx - E * (Bx - np.arctan(Bx))))


# ----------------------------------------------------------------------
# Lateral
# ----------------------------------------------------------------------
def lateral_coeffs(Fz, p):
    dfz = (Fz - p.Fz0) / p.Fz0

    mu = p.pDy1 + p.pDy2 * dfz         # friction falls as load rises
    C = p.pCy1
    D = mu * Fz                        # peak force
    E = p.pEy1 + p.pEy2 * dfz

    if p.linear_k:
        K = p.pKy1 * Fz                # stiffness proportional to load
    else: 
        K = p.pKy1 * p.Fz0 * np.sin(2.0 * np.arctan(Fz / (p.pKy2 * p.Fz0)))    # rises with load, flattens, then rolls off

    B = K / (C * D)                    # B is fixed by the initial slope
    return B, C, D, E, K, mu


def Fy(alpha, Fz, p):
    """Lateral force [N]. alpha in RADIANS."""
    B, C, D, E, _, _ = lateral_coeffs(Fz, p)
    return magic_formula(alpha, B, C, D, E)


def cornering_stiffness(Fz, p):
    """dFy/dalpha at alpha = 0, in N/rad. Multiply by DEG for N/deg."""
    return lateral_coeffs(Fz, p)[4]


def Fy_peak(Fz, p):
    """Largest lateral force the curve reaches [N]."""
    return lateral_coeffs(Fz, p)[2]


# ----------------------------------------------------------------------
# Longitudinal
# ----------------------------------------------------------------------
def longitudinal_coeffs(Fz, p):
    dfz = (Fz - p.Fz0) / p.Fz0

    mu = p.pDx1 + p.pDx2 * dfz
    C = p.pCx1
    D = mu * Fz
    E = p.pEx1 + p.pEx2 * dfz + p.pEx3 * dfz**2

    if p.linear_k:
        K = p.pKx1 * Fz
    else:
        K = Fz * (p.pKx1 + p.pKx2 * dfz) * np.exp(p.pKx3 * dfz)

    B = K / (C * D)
    return B, C, D, E, K, mu


def Fx(kappa, Fz, p):
    """Longitudinal force [N]. kappa dimensionless."""
    B, C, D, E, _, _ = longitudinal_coeffs(Fz, p)
    # SVx shifts the curve off the origin; it is zero unless the source data was drawbar pull, which is already net of motion resistance
    return magic_formula(kappa, B, C, D, E) + p.SVx * Fz / p.Fz0


def Fx_peak(Fz, p):
    return longitudinal_coeffs(Fz, p)[2]


# ----------------------------------------------------------------------
# Combined slip
# ----------------------------------------------------------------------
def Fy_available(Fx_used, Fz, p):
    """
    Lateral force left once Fx_used is being asked for, from the friction
    ellipse (Fx/Fx_max)^2 + (Fy/Fy_max)^2 = 1.
    """
    ratio = Fx_used / Fx_peak(Fz, p)
    return Fy_peak(Fz, p) * np.sqrt(1.0 - ratio**2)


def combined_slip(alpha, kappa, Fz, p):
    """Ellipse-scaled combined slip. Returns (Fx, Fy)."""
    fx = Fx(kappa, Fz, p)
    fy = Fy(alpha, Fz, p) * Fy_available(fx, Fz, p) / Fy_peak(Fz, p)
    return np.broadcast_arrays(fx, fy)
