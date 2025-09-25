"""Command line entry point for the black hole visualisation project."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import List

from blackhole_vis.physics.black_hole import SchwarzschildBlackHole
from blackhole_vis.physics.constants import CONSTANTS
from blackhole_vis.physics.light import LightRay, trace_light_ray, weak_field_deflection
from blackhole_vis.physics.orbits import Orbit, integrate_orbit, orbital_period
from blackhole_vis.visualization.scene import AccretionDisk, render_scene


def _sample_light_rays(black_hole: SchwarzschildBlackHole) -> List[LightRay]:
    base_impact = 6 * black_hole.schwarzschild_radius
    return [
        trace_light_ray(black_hole, impact_parameter=base_impact * factor, samples=1500)
        for factor in (1.0, 1.5, 2.0)
    ]


def _sample_orbits(black_hole: SchwarzschildBlackHole) -> List[Orbit]:
    a1 = 8 * black_hole.innermost_stable_circular_orbit
    a2 = 16 * black_hole.innermost_stable_circular_orbit
    return [
        integrate_orbit(black_hole, semi_major_axis=a1, eccentricity=0.2, phi_max=8 * math.pi),
        integrate_orbit(black_hole, semi_major_axis=a2, eccentricity=0.4, phi_max=10 * math.pi, inclination=math.radians(25)),
    ]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mass", type=float, default=10.0, help="Black hole mass in solar masses")
    parser.add_argument("--output", type=Path, default=Path("render.png"), help="Path to save the rendered scene")
    parser.add_argument("--show", action="store_true", help="Display the Matplotlib window instead of saving")
    args = parser.parse_args()

    mass_kg = args.mass * CONSTANTS.solar_mass
    black_hole = SchwarzschildBlackHole(mass=mass_kg)

    disk = AccretionDisk(
        inner_radius=black_hole.innermost_stable_circular_orbit * 1.2,
        outer_radius=black_hole.innermost_stable_circular_orbit * 5.0,
        thickness=black_hole.schwarzschild_radius * 0.8,
    )

    rays = _sample_light_rays(black_hole)
    orbits = _sample_orbits(black_hole)

    figure = render_scene(black_hole, disk=disk, light_rays=rays, orbits=orbits, save_path=None if args.show else str(args.output))

    for ray in rays:
        approx = weak_field_deflection(black_hole, ray.impact_parameter)
        print(
            f"Impact parameter: {ray.impact_parameter:.2e} m | deflection (numeric): {ray.total_deflection:.3f} rad | weak-field: {approx:.3f} rad"
        )

    for orbit in orbits:
        period = orbital_period(black_hole, orbit.semi_major_axis)
        print(
            f"Orbit a={orbit.semi_major_axis:.2e} m e={orbit.eccentricity:.2f} | period={period/CONSTANTS.day:.2f} days | precession={orbit.relativistic_precession:.3e} rad"
        )

    if args.show:
        import matplotlib.pyplot as plt

        plt.show()


if __name__ == "__main__":
    main()
