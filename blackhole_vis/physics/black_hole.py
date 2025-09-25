"""Models describing Schwarzschild black holes."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np

from .constants import CONSTANTS


@dataclass
class SchwarzschildBlackHole:
    """Representation of a non-rotating black hole."""

    mass: float

    def __post_init__(self) -> None:
        if self.mass <= 0:
            raise ValueError("Mass of a black hole must be positive.")

    @property
    def gravitational_parameter(self) -> float:
        """Return the gravitational parameter :math:`GM`."""

        return CONSTANTS.gravitational_constant * self.mass

    @property
    def schwarzschild_radius(self) -> float:
        """Return the Schwarzschild radius (event horizon)."""

        return 2.0 * self.gravitational_parameter / CONSTANTS.speed_of_light**2

    @property
    def photon_sphere_radius(self) -> float:
        """Radius of the unstable photon orbit."""

        return 1.5 * self.schwarzschild_radius

    @property
    def innermost_stable_circular_orbit(self) -> float:
        """Return the radius of the innermost stable circular orbit (ISCO)."""

        return 3.0 * self.schwarzschild_radius

    def gravitational_redshift(self, emitter_radius: float, observer_radius: float | None = None) -> float:
        """Return the gravitational redshift factor between two radii."""

        if emitter_radius <= self.schwarzschild_radius:
            raise ValueError("Emission radius must lie outside the event horizon.")
        if observer_radius is not None and observer_radius <= self.schwarzschild_radius:
            raise ValueError("Observer radius must lie outside the event horizon.")

        r_obs = np.inf if observer_radius is None else observer_radius
        rs = self.schwarzschild_radius
        return np.sqrt((1 - rs / r_obs) / (1 - rs / emitter_radius))

    def tidal_acceleration(self, radius: float, separation: float) -> float:
        """Compute the tidal acceleration between two points."""

        if radius <= self.schwarzschild_radius:
            raise ValueError("Radius must lie outside the event horizon.")
        return 2 * self.gravitational_parameter * separation / radius**3

    def gravitational_time_dilation(self, radius: float | np.ndarray) -> np.ndarray:
        """Return the factor relating local proper time to far-away time."""

        radius_arr = np.asarray(radius)
        if np.any(radius_arr <= self.schwarzschild_radius):
            raise ValueError("Radius must lie outside the event horizon.")
        rs = self.schwarzschild_radius
        return np.sqrt(1 - rs / radius_arr)

    def escape_velocity(self, radius: float) -> float:
        """Return the escape velocity from ``radius`` in m/s."""

        if radius <= 0:
            raise ValueError("Radius must be positive.")
        return np.sqrt(2 * self.gravitational_parameter / radius)

    def light_travel_time(self, path_radii: Iterable[float]) -> float:
        """Estimate the light travel time along a radial path."""

        path_radii = np.asarray(list(path_radii))
        if np.any(path_radii <= self.schwarzschild_radius):
            raise ValueError("All radii must lie outside the event horizon.")
        dilation = self.gravitational_time_dilation(path_radii)
        return np.trapz(1.0 / dilation, path_radii) / CONSTANTS.speed_of_light
