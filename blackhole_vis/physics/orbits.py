"""Planetary dynamics influenced by a Schwarzschild black hole."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .black_hole import SchwarzschildBlackHole
from .constants import CONSTANTS


@dataclass
class Orbit:
    """A sampled orbital trajectory in Cartesian coordinates."""

    xyz: np.ndarray
    semi_major_axis: float
    eccentricity: float
    relativistic_precession: float


def _orbit_rhs(u: float, du_dphi: float, black_hole: SchwarzschildBlackHole, ang_momentum: float) -> float:
    gm = black_hole.gravitational_parameter
    return gm / ang_momentum**2 + 3 * gm / CONSTANTS.speed_of_light**2 * u**2 - u


def integrate_orbit(
    black_hole: SchwarzschildBlackHole,
    semi_major_axis: float,
    eccentricity: float,
    phi_max: float = 12 * math.pi,
    samples: int = 5000,
    inclination: float = math.radians(10.0),
) -> Orbit:
    """Integrate a post-Newtonian orbit around the black hole."""

    if not 0 <= eccentricity < 1:
        raise ValueError("Eccentricity must be between 0 (inclusive) and 1 (exclusive).")
    if semi_major_axis <= 0:
        raise ValueError("Semi-major axis must be positive.")

    gm = black_hole.gravitational_parameter
    ang_momentum = math.sqrt(gm * semi_major_axis * (1 - eccentricity**2))

    phi = np.linspace(0, phi_max, samples)
    dphi = phi[1] - phi[0]

    r_peri = semi_major_axis * (1 - eccentricity)
    u = np.zeros_like(phi)
    du = np.zeros_like(phi)
    u[0] = 1.0 / r_peri
    du[0] = 0.0

    for i in range(1, len(phi)):
        k1 = du[i - 1]
        l1 = _orbit_rhs(u[i - 1], du[i - 1], black_hole, ang_momentum)

        k2 = du[i - 1] + 0.5 * dphi * l1
        l2 = _orbit_rhs(u[i - 1] + 0.5 * dphi * k1, du[i - 1] + 0.5 * dphi * l1, black_hole, ang_momentum)

        k3 = du[i - 1] + 0.5 * dphi * l2
        l3 = _orbit_rhs(u[i - 1] + 0.5 * dphi * k2, du[i - 1] + 0.5 * dphi * l2, black_hole, ang_momentum)

        k4 = du[i - 1] + dphi * l3
        l4 = _orbit_rhs(u[i - 1] + dphi * k3, du[i - 1] + dphi * l3, black_hole, ang_momentum)

        u[i] = u[i - 1] + dphi / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        du[i] = du[i - 1] + dphi / 6.0 * (l1 + 2 * l2 + 2 * l3 + l4)

    r = 1.0 / u
    r[r < black_hole.innermost_stable_circular_orbit] = np.nan

    x = r * np.cos(phi)
    y = r * np.sin(phi)
    z = r * np.sin(inclination) * np.sin(phi)

    xyz = np.vstack([x, y, z]).T

    finite = np.isfinite(r)
    if finite.sum() > 3:
        phi_finite = phi[finite]
        r_finite = r[finite]
        dr_dphi = np.gradient(r_finite, phi_finite)
        d2r_dphi2 = np.gradient(dr_dphi, phi_finite)
        sign_change = (np.hstack([dr_dphi[1:], 0.0]) <= 0) & (np.hstack([dr_dphi[:-1], 0.0]) >= 0)
        maxima = np.where(sign_change & (d2r_dphi2 < 0))[0]
        if len(maxima) >= 2:
            delta_phi = phi_finite[maxima[1]] - phi_finite[maxima[0]]
            relativistic_precession = delta_phi - 2 * math.pi
        else:
            relativistic_precession = float("nan")
    else:
        relativistic_precession = float("nan")

    return Orbit(
        xyz=xyz,
        semi_major_axis=semi_major_axis,
        eccentricity=eccentricity,
        relativistic_precession=relativistic_precession,
    )


def orbital_period(black_hole: SchwarzschildBlackHole, semi_major_axis: float) -> float:
    """Return the relativistic orbital period using Kepler's third law with corrections."""

    gm = black_hole.gravitational_parameter
    newtonian = 2 * math.pi * math.sqrt(semi_major_axis**3 / gm)
    correction = 1 + 3 * gm / (CONSTANTS.speed_of_light**2 * semi_major_axis)
    return newtonian * correction
