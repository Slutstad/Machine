"""Utilities for visualising black holes and associated physics in 3D."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Iterable, Optional

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 - side effect import required by Matplotlib

from ..physics.black_hole import SchwarzschildBlackHole
from ..physics.light import LightRay
from ..physics.orbits import Orbit


@dataclass
class AccretionDisk:
    """Simple representation of an accretion disk used for visualisation."""

    inner_radius: float
    outer_radius: float
    thickness: float
    particle_count: int = 2000


def generate_accretion_disk(disk: AccretionDisk) -> np.ndarray:
    """Return randomised particle positions representing an accretion disk."""

    radii = np.sqrt(np.random.rand(disk.particle_count))
    radii = disk.inner_radius + (disk.outer_radius - disk.inner_radius) * radii
    angles = 2 * math.pi * np.random.rand(disk.particle_count)
    heights = (np.random.rand(disk.particle_count) - 0.5) * disk.thickness

    x = radii * np.cos(angles)
    y = radii * np.sin(angles)
    z = heights
    return np.vstack([x, y, z]).T


def _draw_event_horizon(ax: Axes3D, black_hole: SchwarzschildBlackHole) -> None:
    radius = black_hole.schwarzschild_radius
    theta = np.linspace(0, math.pi, 64)
    phi = np.linspace(0, 2 * math.pi, 128)
    theta_grid, phi_grid = np.meshgrid(theta, phi)
    x = radius * np.sin(theta_grid) * np.cos(phi_grid)
    y = radius * np.sin(theta_grid) * np.sin(phi_grid)
    z = radius * np.cos(theta_grid)
    ax.plot_surface(x, y, z, color="black", linewidth=0, antialiased=False, alpha=0.9)


def _draw_disk(ax: Axes3D, disk_points: np.ndarray, color: str = "#ff9933") -> None:
    ax.scatter(disk_points[:, 0], disk_points[:, 1], disk_points[:, 2], s=1.0, c=color, alpha=0.5)


def _draw_light_rays(ax: Axes3D, light_rays: Iterable[LightRay]) -> None:
    for ray in light_rays:
        ax.plot(ray.xyz[:, 0], ray.xyz[:, 1], ray.xyz[:, 2], color="#66ccff", linewidth=1.5, alpha=0.8)


def _draw_orbits(ax: Axes3D, orbits: Iterable[Orbit]) -> None:
    for orbit in orbits:
        ax.plot(orbit.xyz[:, 0], orbit.xyz[:, 1], orbit.xyz[:, 2], color="#00ff99", linewidth=1.2)


def render_scene(
    black_hole: SchwarzschildBlackHole,
    disk: Optional[AccretionDisk] = None,
    light_rays: Iterable[LightRay] = (),
    orbits: Iterable[Orbit] = (),
    save_path: Optional[str] = None,
    show_axes: bool = False,
) -> plt.Figure:
    """Render the configured scene and optionally save it to ``save_path``."""

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection="3d")

    _draw_event_horizon(ax, black_hole)

    if disk is not None:
        disk_points = generate_accretion_disk(disk)
        _draw_disk(ax, disk_points)

    _draw_light_rays(ax, light_rays)
    _draw_orbits(ax, orbits)

    horizon = black_hole.schwarzschild_radius
    limit = horizon * 30
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.set_zlim(-limit / 2, limit / 2)

    ax.set_box_aspect((1, 1, 0.5))
    ax.set_facecolor("#000010")
    fig.patch.set_facecolor("black")

    if not show_axes:
        ax.set_axis_off()

    if save_path is not None:
        fig.savefig(save_path, dpi=200, bbox_inches="tight", facecolor=fig.get_facecolor())

    return fig
