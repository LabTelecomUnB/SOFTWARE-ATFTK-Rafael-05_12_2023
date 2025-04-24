#  -*- coding: utf-8 -*-
"""
Author: Rafael R. L. Benevides
"""

# ---------- ---------- ---------- ---------- ---------- ---------- py
import numpy
import stripy
import pyvista

from matplotlib.colors import rgb2hex
from matplotlib import pyplot

# ---------- ---------- ---------- ---------- ---------- ---------- cy
from numpy cimport ndarray


class SphericalMesh:

    def __init__(self, resolution: int = 9, refinement_levels: int = 0):
        self._mesh = stripy.spherical_meshes.uniform_ring_mesh(
            resolution=resolution,
            refinement_levels=refinement_levels)

    def __len__(self):
        return self.theta.shape[0]

    def unstructured_grid(self, ndarray[double, ndim=2] points):
        faces = self._mesh.simplices
        # points = self._mesh.points

        threes = 3 * numpy.ones((faces.shape[0], 1), dtype=numpy.int64)

        cells = numpy.hstack([threes, faces]).flatten()
        cells = numpy.hstack([faces.shape[0], cells])
        cells = numpy.hstack([cells.shape[0], cells])

        POLYHEDRON = pyvista.CellType.POLYHEDRON

        return pyvista.UnstructuredGrid(cells, [POLYHEDRON], points)

    def plot_3D(self,
                point_size: float = 5.0,
                show_faces: bool = True,
                show_edges: bool = True,
                line_width: float = 1.0,
                opacity: float = 0.5):
        plotter = pyvista.Plotter()
        plotter.background_color = 'gray'

        if show_faces:
            # faces = self._mesh.simplices
            # points = self._mesh.points
            #
            # threes = 3 * numpy.ones((faces.shape[0], 1), dtype=numpy.int64)
            #
            # cells = numpy.hstack([threes, faces]).flatten()
            # cells = numpy.hstack([faces.shape[0], cells])
            # cells = numpy.hstack([cells.shape[0], cells])
            #
            # grid = pyvista.UnstructuredGrid(cells,
            #                                 [pyvista.CellType.POLYHEDRON],
            #                                 points)

            plotter.add_mesh(self.unstructured_grid(self._mesh.points),
                             show_edges=show_edges,
                             line_width=line_width,
                             opacity=opacity)

        plotter.add_points(self._mesh.points,
                           render_points_as_spheres=True,
                           point_size=point_size,
                           color=rgb2hex('C0'))

        plotter.add_axes_at_origin(xlabel='x', ylabel='y', zlabel='z', line_width=10)
        plotter.show()

    def plot_2D(self):
        figure, axes = pyplot.subplots()

        axes.plot(numpy.rad2deg(self.phi),
                  numpy.rad2deg(self.theta), marker='o', linestyle='')

        axes.set_xlim(-180, 180)
        axes.set_ylim(0, 180)

        axes.set_xlabel(r'Azimuthal angle - $\phi$ (deg)')
        axes.set_xlabel(r'Polar angle - $\theta$ (deg)')

        axes.set_title(f'Spreading: {self.spreading:.0%}')

        pyplot.show()

    @property
    def theta(self):
        return numpy.pi / 2 - self._mesh.lats

    @property
    def phi(self):
        return self._mesh.lons

    @property
    def points(self):
        return self._mesh.points

    @property
    def x(self):
        return self.points[:, 0]

    @property
    def y(self):
        return self.points[:, 1]

    @property
    def z(self):
        return self.points[:, 2]

    @property
    def spreading(self) -> float:
        areas = self._mesh.areas()
        return areas.min() / areas.max()