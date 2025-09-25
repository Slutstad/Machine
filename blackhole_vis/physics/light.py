"""Light propagation and lensing utilities."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .black_hole import SchwarzschildBlackHole
from .constants import CONSTANTS


@dataclass
class LightRay:
    """A sampled light ray trajectory in Cartesian coordinates."""

    xyz: np.ndarray
    impact_parameter: float
    total_deflection: float


def _closest_approach_radius(black_hole: SchwarzschildBlackHole, impact_parameter: float) -> float:
    """Return the closest approach radius for a photon with ``impact_parameter``."""

    rs = black_hole.schwarzschild_radius
    critical_b = math.sqrt(27) / 2 * rs
    if impact_parameter <= critical_b:
        raise ValueError(
            "Impact parameter is below the critical value; the photon would be captured."
        )

    def f(r: float) -> float:
        term = 1.0 / impact_parameter**2
        return term - (1.0 - rs / r) / r**2

    lower = max(1.0001 * rs, rs * 1.01)
    upper = impact_parameter * 100.0
    if f(upper) <= 0:
        raise ValueError("Upper search bound is insufficient for root finding.")

    for _ in range(200):
        mid = 0.5 * (lower + upper)
        value = f(mid)
        if abs(value) < 1e-12:
            return mid
        if value > 0:
            lower = mid
        else:
            upper = mid
    return mid


def _phi_integrand(r: np.ndarray, rs: float, impact_parameter: float) -> np.ndarray:
    term = 1.0 / impact_parameter**2 - (1.0 - rs / r) / r**2
    term = np.clip(term, 1e-20, None)
    return 1.0 / (r**2 * np.sqrt(term))


def trace_light_ray(
    black_hole: SchwarzschildBlackHole,
    impact_parameter: float,
    r_far_multiplier: float = 80.0,
    samples: int = 2000,
) -> LightRay:
    """Compute the trajectory of a photon deflected by the black hole."""

    rs = black_hole.schwarzschild_radius
    r_min = _closest_approach_radius(black_hole, impact_parameter)
    r_far = impact_parameter * r_far_multiplier
    epsilon = r_min * 1e-6
    r_turn = r_min + epsilon
    r_in = np.linspace(r_far, r_turn, samples // 2)
    r_out = np.linspace(r_turn, r_far, samples // 2)

    integrand_in = _phi_integrand(r_in, rs, impact_parameter)
    integrand_out = _phi_integrand(r_out, rs, impact_parameter)

    phi_in = np.zeros_like(r_in)
    for i in range(1, len(r_in)):
        step = r_in[i] - r_in[i - 1]
        phi_in[i] = phi_in[i - 1] - 0.5 * (integrand_in[i] + integrand_in[i - 1]) * step

    phi_in -= phi_in[-1]

    phi_out = np.zeros_like(r_out)
    for i in range(1, len(r_out)):
        step = r_out[i] - r_out[i - 1]
        phi_out[i] = phi_out[i - 1] + 0.5 * (integrand_out[i] + integrand_out[i - 1]) * step

    phi = np.concatenate([phi_in, phi_out[1:]])
    r = np.concatenate([r_in, r_out[1:]])

    x = r * np.cos(phi)
    y = r * np.sin(phi)
    z = np.zeros_like(x)
    xyz = np.vstack([x, y, z]).T

    total_deflection = 2 * abs(phi[0])
    return LightRay(xyz=xyz, impact_parameter=impact_parameter, total_deflection=total_deflection)


def weak_field_deflection(black_hole: SchwarzschildBlackHole, impact_parameter: float) -> float:
    """Return the weak-field approximation of the deflection angle."""

    return 4.0 * black_hole.gravitational_parameter / (
        CONSTANTS.speed_of_light**2 * impact_parameter
    )


def shapiro_time_delay(
    black_hole: SchwarzschildBlackHole, distance_source: float, distance_observer: float, impact_parameter: float
) -> float:
    """Return the Shapiro time delay for a light signal."""

    rs = black_hole.schwarzschild_radius
    numerator = distance_source + distance_observer + math.sqrt(distance_source**2 + impact_parameter**2)
    numerator *= math.sqrt(distance_observer**2 + impact_parameter**2)
    return (2 * black_hole.gravitational_parameter / CONSTANTS.speed_of_light**3) * math.log(
        numerator / (impact_parameter**2)
    )


def photon_sphere_frequency(black_hole: SchwarzschildBlackHole) -> float:
    """Return the orbital frequency of photons on the photon sphere."""

    r_photon = black_hole.photon_sphere_radius
    return CONSTANTS.speed_of_light / (2 * math.pi * r_photon)
