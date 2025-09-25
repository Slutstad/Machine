"""Physical constants used throughout the black hole visualisation package."""

from dataclasses import dataclass


@dataclass(frozen=True)
class PhysicalConstants:
    """Container for the physical constants required by the simulation."""

    gravitational_constant: float = 6.67430e-11  # m^3 kg^-1 s^-2
    speed_of_light: float = 2.99792458e8  # m s^-1
    astronomical_unit: float = 1.495978707e11  # m
    solar_mass: float = 1.98847e30  # kg
    day: float = 86400.0  # s


CONSTANTS = PhysicalConstants()
