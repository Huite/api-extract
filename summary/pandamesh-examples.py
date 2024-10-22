"""
Basic Triangle Example
======================

In this example we'll create some basic geometries and turn them into meshes.
to illustrate some of the mesh generation features that Triangle provides in
combination with polygon, point, and linestring geometries represented by
geopandas.
"""
# %%
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import shapely.geometry as sg

import pandamesh as pm

# %%
# A simple rectangular mesh
# -------------------------
#
# The most simple example is perhaps a rectangle. We'll create a vector
# geometry, store this in a geodataframe, and associate a cell size.

polygon = sg.Polygon(
    [
        [0.0, 0.0],
        [10.0, 0.0],
        [10.0, 10.0],
        [0.0, 10.0],
    ]
)
gdf = gpd.GeoDataFrame(geometry=[polygon])
gdf["cellsize"] = 2.0

# %%
# We'll use this polygon to generate a mesh. We start by initializing a
# TriangleMesher, which is a simple wrapper around the Python bindings to the
# Triangle C-library. This wrapper extracts the coordinates and presents them
# in the appropriate manner for triangle.

mesher = pm.TriangleMesher(gdf)
vertices, triangles = mesher.generate()
pm.plot(vertices, triangles)

# %%
# Defaults
# --------
#
# The TriangleMesher class is initialized with a number of default parameters:

print(mesher)

# %%
# We can change a parameter, and see what effects this has on the mesh:

mesher.conforming_delaunay = False
vertices, triangles = mesher.generate()
pm.plot(vertices, triangles)

# %%
# To generate a mesh with smaller cell sizes, we adjust the geodataframe, and
# recreate the mesher.

gdf["cellsize"] = 1.0
mesher = pm.TriangleMesher(gdf)
vertices, triangles = mesher.generate()
pm.plot(vertices, triangles)
# %%
# Multiple cell size zones
# ------------------------
#
# Multiple zones of cell sizes are supported, as every polygon can be associated
# with a cell size in the geodataframe.

polygon2 = sg.Polygon(
    [
        [10.0, 0.0],
        [20.0, 0.0],
        [20.0, 10.0],
        [10.0, 10.0],
    ]
)
gdf = gpd.GeoDataFrame(geometry=[polygon, polygon2])
gdf["cellsize"] = [2.0, 1.0]

mesher = pm.TriangleMesher(gdf)
vertices, triangles = mesher.generate()
pm.plot(vertices, triangles)
# %%
# Polygons with holes ("donut" geometries)
# ----------------------------------------
#
# Holes in polygons work as expected:

outer = [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]
inner = [(3.0, 3.0), (7.0, 3.0), (7.0, 7.0), (3.0, 7.0)]

donut = sg.Polygon(shell=outer, holes=[inner])
gdf = gpd.GeoDataFrame(geometry=[donut])
gdf["cellsize"] = [2.0]

mesher = pm.TriangleMesher(gdf)
vertices, triangles = mesher.generate()
pm.plot(vertices, triangles)

# %%
# Local refinement
# ----------------
#
# To do local refinement, we need to ensure there is no overlap between the
# polygons. The coordinates of the hole of the outer polygon should match
# exactly with the coordinates of the exterior boundary of the inner polygon.

refined = sg.Polygon(inner)

gdf = gpd.GeoDataFrame(geometry=[donut, refined])
gdf["cellsize"] = [2.0, 0.5]

mesher = pm.TriangleMesher(gdf)
vertices, triangles = mesher.generate()
pm.plot(vertices, triangles)

# %%
# Force points into the triangulation
# -----------------------------------
#
# We may also force points into the triangulation, by adding points to the
# geodataframe. Let's assume we'd like to a series of points at x=1.0, at a
# distance of 0.5.

y = np.arange(0.5, 10.0, 0.5)
x = np.full(y.size, 1.0)
points = gpd.points_from_xy(x, y)

gdf = gpd.GeoDataFrame(geometry=[donut, refined, *points])
gdf["cellsize"] = [2.0, 0.5] + (len(points) * [np.nan])
gdf.plot(facecolor="none")

# %%
# We can now see the points forced in the triangulation, by plotting the
# contents of the geodataframe on top of the generated mesh:

mesher = pm.TriangleMesher(gdf)
vertices, triangles = mesher.generate()

fig, ax = plt.subplots()
pm.plot(vertices, triangles, ax=ax)
gdf.plot(facecolor="none", edgecolor="red", ax=ax)
# %%
# Force linestrings into the triangulation
# ----------------------------------------
#
# We may do the same with linestrings. Here, we will add a vertical line at
# x = 9.0.

line = sg.LineString(
    [
        [9.0, 2.0],
        [9.0, 8.0],
    ]
)
gdf = gpd.GeoDataFrame(geometry=[donut, refined, line, *points])
gdf["cellsize"] = [2.0, 0.5, np.nan] + (len(points) * [np.nan])

mesher = pm.TriangleMesher(gdf)
vertices, triangles = mesher.generate()

fig, ax = plt.subplots()
pm.plot(vertices, triangles, ax=ax)
gdf.plot(facecolor="none", edgecolor="red", ax=ax)

# %%
# Specify cell size along line string
# -----------------------------------
#
# Finally, we may also specify the cell size along the line.

line = sg.LineString([(2.0, 8.0), (8.0, 2.0)])
gdf = gpd.GeoDataFrame(geometry=[polygon, line])
gdf["cellsize"] = [2.0, 0.5]

fig, ax = plt.subplots()

mesher = pm.TriangleMesher(gdf)
vertices, triangles = mesher.generate()
pm.plot(vertices, triangles, ax=ax)
gdf.plot(facecolor="none", edgecolor="red", ax=ax)

# %%
# Conclusion
# ----------
#
# In real use, the vector geometries will be more complex, and not based on
# just a few coordinate pairs. Such cases are presented in the other examples,
# but the same principles apply: we may use polygons with associated cell
# sizes, and linestrings and points to steer the triangulation.


"""
Basic Gmsh Example
==================

In this example we'll create some basic geometries and turn them into meshes.
to illustrate some of the mesh generation features that Gmsh provides in
combination with polygon, point, and linestring geometries represented by
geopandas.

The :py:class:`GmshMesher` supports the geometry show in the basic Triangle
example and has a number of additional features.
"""
# %%
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import shapely.geometry as sg

import pandamesh as pm

# sphinx_gallery_start_ignore
pm.GmshMesher.finalize()
# sphinx_gallery_end_ignore

# %%
# A simple rectangular mesh
# -------------------------
#
# The most simple example is perhaps a rectangle. We'll create a vector
# geometry, store this in a geodataframe, and associate a cell size.

polygon = sg.Polygon(
    [
        [0.0, 0.0],
        [10.0, 0.0],
        [10.0, 10.0],
        [0.0, 10.0],
    ]
)
gdf = gpd.GeoDataFrame(geometry=[polygon])
gdf["cellsize"] = 2.0

# %%
# We'll use this polygon to generate a mesh. We start by initializing a
# TriangleMesher, which is a simple wrapper around the Python bindings to the
# Gmsh C++-library. This wrapper extracts the coordinates and presents them
# in the appropriate manner for Gmsh.

mesher = pm.GmshMesher(gdf)
vertices, triangles = mesher.generate()
pm.plot(vertices, triangles)

# %%
# Before we can instantiate another GmshMesher, we need to ``finalize`` the old
# one.

mesher.finalize()

# %%
# As the name suggests, Triangle only generates triangular meshes. Gmsh is
# capable of generating quadrilateral-dominant meshes, and has a lot more bells
# and whistles for defining cellsizes.

line = sg.LineString([(2.0, 8.0), (8.0, 2.0)])
gdf = gpd.GeoDataFrame(geometry=[polygon, line])
gdf["cellsize"] = [2.0, 0.5]

fig, (ax0, ax1) = plt.subplots(ncols=2)

mesher = pm.TriangleMesher(gdf)
vertices, triangles = mesher.generate()
pm.plot(vertices, triangles, ax=ax0)

mesher = pm.GmshMesher(gdf)
vertices, triangles = mesher.generate()
pm.plot(vertices, triangles, ax=ax1)

# %%
# Gmsh allows for specifying cell sizes in a more flexible way. Triangle (left)
# only supports polygons (regions) with fixed cell sizes and explicitly placed
# vertices. Gmsh is capable of forcing refinement in a larger zone around
# features as is visible around the diagonal (right).
#
# Defaults
# --------
#
# The GmshMesher class is initialized with a number of default parameters:

print(mesher)

mesher.finalize()

# %%
# The parameters of Gmsh differ from Triangle, but they work the same: they can
# be altered after initialization to control the triangulation.
#
# Forcing points, lines, local refinement
# ---------------------------------------
#
# We can force points and lines into the triangulation:

outer = [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]
inner = [(3.0, 3.0), (7.0, 3.0), (7.0, 7.0), (3.0, 7.0)]
donut = sg.Polygon(shell=outer, holes=[inner])
refined = sg.Polygon(inner)

y = np.arange(0.5, 10.0, 0.5)
x = np.full(y.size, 1.0)
points = gpd.points_from_xy(x, y)

line = sg.LineString(
    [
        [9.0, 2.0],
        [9.0, 8.0],
    ]
)

gdf = gpd.GeoDataFrame(geometry=[donut, refined, line, *points])
gdf["cellsize"] = [2.0, 0.5, 2.0] + (len(points) * [2.0])

mesher = pm.GmshMesher(gdf)
vertices, triangles = mesher.generate()
mesher.finalize()

fig, ax = plt.subplots()
pm.plot(vertices, triangles, ax=ax)
gdf.plot(facecolor="none", edgecolor="red", ax=ax)


# Quadrilateral meshes
# --------------------
#
# One of the features of Gmsh is that it is also capable of generating
# quadrilateral (dominant) meshes, by recombining triangles. We can achieve
# this by changing a parameter on the mesher:

gdf = gpd.GeoDataFrame(geometry=[polygon])
gdf["cellsize"] = 2.0
mesher = pm.GmshMesher(gdf)
mesher.recombine_all = True
vertices, faces = mesher.generate()

pm.plot(vertices, faces)

# %%
# Writing to file
# ---------------
# It's also possible to use the Python bindings to write a Gmsh ``.msh`` file.
# This file can be opened using the Gmsh GUI to e.g. inspect the generated
# mesh.

mesher.write("my-mesh.msh")

# %%
# Conclusion
# ----------
#
# In real use, the vector geometries will be more complex, and not based on
# just a few coordinate pairs. Such cases are presented in the other examples,
# but the same principles apply: we may use polygons, linestrings and points
# with associated cell sizes to steer the triangulation; unlike Triangle,
# for Gmsh cell sizes can associated to linestrings and points, not just
# polygons.

# %%


"""
Gmsh Fields Example
===================

Gmsh supports so called "fields" to guide the cell sizes of the generated
meshes. These fields are separate from the geometrical constraints: for
example, a field point does not end up in the generated mesh, but influences
the cell size in its surrounding.

These field geometries can be added via:

* :meth:`pandamesh.GmshMesher.add_threshold_distance_field()`
* :meth:`pandamesh.GmshMesher.add_matheval_distance_field()`
* :meth:`pandamesh.GmshMesher.add_structured_field()`
* :meth:`pandamesh.GmshMesher.add_structured_field_from_dataarray()`,

The examples below demonstrate how to set up these distance fields for meshing.
"""
# %%
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import shapely.geometry as sg

import pandamesh as pm

# sphinx_gallery_start_ignore
pm.GmshMesher.finalize()
# sphinx_gallery_end_ignore

# %%
# Point fields
# ------------
#
# We'll start again with simple rectangular example.

polygon = sg.Polygon(
    [
        [0.0, 0.0],
        [10.0, 0.0],
        [10.0, 10.0],
        [0.0, 10.0],
    ]
)
point = sg.Point([4.0, 4.0])
gdf = gpd.GeoDataFrame(geometry=[polygon])
gdf["cellsize"] = 5.0

mesher = pm.GmshMesher(gdf, shift_origin=False)
mesher.mesh_size_extend_from_boundary = False
mesher.mesh_size_from_curvature = False
mesher.mesh_size_from_points = False

pm.plot(*mesher.generate())

# %%
# Threshold distance fields
# -------------------------
#
# Gmsh supports changing cell sizes gradually, for example as a function of
# distance to a feature. We can add a point, and connect a distance threshold
# field to it:

point = sg.Point([4.0, 4.0])
field = gpd.GeoDataFrame(geometry=[point])
field["dist_min"] = 2.0
field["dist_max"] = 4.0
field["size_min"] = 0.5
field["size_max"] = 2.5
field["spacing"] = np.nan
mesher.add_threshold_distance_field(field)

vertices, faces = mesher.generate()
pm.plot(vertices, faces)

# %%
# Within the ``dist_min`` of the point, all cell sizes have size of at most
# ``size_min``. This changes linearly until ``dist_max`` is reached, at which point
# the cell sizes become ``size_max``.
#
# Fields can be removed via ``.clear_fields()``:

mesher.clear_fields()
vertices, faces = mesher.generate()
pm.plot(vertices, faces)

# %%
# Gmsh only measures distances to point. The ``spacing`` is used to interpolate
# points along lines:

mesher.clear_fields()

line = sg.LineString(
    [
        [3.0, -3.0],
        [3.0, 13.0],
    ]
)
field = gpd.GeoDataFrame(geometry=[line])
field["dist_min"] = 2.0
field["dist_max"] = 4.0
field["size_min"] = 0.5
field["size_max"] = 2.5
field["spacing"] = 2.0
mesher.add_threshold_distance_field(field)

vertices, faces = mesher.generate()
pm.plot(vertices, faces)

# %%
# Note that unlike the mesher input geometries, these geometries may fall
# outside the meshing domain: they only "radiate" a cell size.
#
# Polygons can also be used as field geometries. Distances are measured from
# internal and external boundaries:

mesher.clear_fields()

square = sg.Polygon(
    [
        [3.0, 3.0],
        [7.0, 3.0],
        [7.0, 7.0],
        [3.0, 7.0],
    ]
)
field = gpd.GeoDataFrame(geometry=[square])
field["dist_min"] = 0.5
field["dist_max"] = 1.5
field["size_min"] = 0.3
field["size_max"] = 2.5
field["spacing"] = 1.0
mesher.add_threshold_distance_field(field)

vertices, faces = mesher.generate()
pm.plot(vertices, faces)

# %%
# MathEval distance fields
# ------------------------
#
# Gmsh also supports arbitrary mathematical functions. With Pandamesh, these
# can be easily combined to specify cell size a function to some boundary. For
# example, we can specify cell size as quadratically growing with the distance
# from the left boundary:

mesher.clear_fields()

line = sg.LineString(
    [
        [0.0, 0.0],
        [0.0, 10.0],
    ]
)
field = gpd.GeoDataFrame(geometry=[line])
field["function"] = "distance^2 + 0.3"
field["spacing"] = 1.0
mesher.add_matheval_distance_field(field)

vertices, faces = mesher.generate()
pm.plot(vertices, faces)

# %%
# Note that we should take care to specify a function which is always larger
# than zero in the meshing domain.
#
# Unlike input geometries, fields can be added in a piece by piece manner. The
# distance is always relative to the feature of the geometry in the
# GeoDataFrame row.

second_field = gpd.GeoDataFrame(geometry=[sg.Point([5.0, 5.0])])
second_field["function"] = "max(1/(distance^2), 2.0)"
second_field["spacing"] = np.nan
mesher.add_matheval_distance_field(second_field)

vertices, faces = mesher.generate()
pm.plot(vertices, faces)

# %%
# Structured fields
# -----------------
#
# In some cases, the generated cell size should depend on some physical
# properties of the domain. In geospatial applications, such properties are
# often represented as raster data. These data can be used to guide mesh
# generation as a structured grid. The cell size is prescribed at the grid
# points, and interpolated between.
#
# In the example below, we generate 3 by 3 grid of cell sizes, with small cell
# sizes in the lower left corner, and large cell sizes in the upper right:

mesher.clear_fields()

y, x = np.meshgrid([1.0, 5.0, 9.0], [1.0, 5.0, 9.0], indexing="ij")
distance_from_origin = np.sqrt((x * x + y * y))
cellsize = np.log(distance_from_origin / distance_from_origin.min()) + 0.5
mesher.add_structured_field(
    cellsize=cellsize,
    xmin=x.min(),
    ymin=y.min(),
    dx=1.0,
    dy=1.0,
)
vertices, faces = mesher.generate()

fig, ax = plt.subplots()
pm.plot(vertices, faces, ax=ax)
ax.scatter(x, y)

# %%
# DataArray structured fields
# ---------------------------
#
# These structured fields can also be provided as xarray DataArrays:

mesher.clear_fields()

import xarray as xr

x = np.arange(1.0, 10.0)
y = np.arange(1.0, 10.0)
da = xr.DataArray(np.ones((y.size, x.size)), coords={"y": y, "x": x}, dims=("y", "x"))

mesher.add_structured_field_from_dataarray(da)
vertices, faces = mesher.generate()
pm.plot(vertices, faces)

# %%
# This is arguably the most flexible way of configuring cell sizes, since we
# can easily modify the DataArray values. Note that like the MathEval
# specification, we need to take care to ensure values remain > 0.

mesher.clear_fields()

cos_da = da * np.cos(da["x"]) + 1.1
mesher.add_structured_field_from_dataarray(cos_da)
vertices, faces = mesher.generate()
pm.plot(vertices, faces)
# %%


"""
Geospatial Triangle Example
===========================

In this example we'll illustrate how to generate a mesh from a "real-world"
geospatial vector dataset.
"""
# %%
import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
import shapely.geometry as sg

import pandamesh as pm

# %%
# Overlap
# -------
#
# We will get the data of a GeoJSON file describing the provinces of the
# Netherlands, and select only the name and geometry columns. We'll set the
# coordinate reference system to the Dutch national standard (EPSG:28992).
# Finally we set the name column to be used as index, so we can select
# provinces on name.

provinces = pm.data.provinces_nl().loc[:, ["name", "geometry"]]
provinces = provinces.to_crs("epsg:28992")
provinces.index = provinces["name"]
gdf = provinces.copy()

# %%
# The mesh generation software cannot deal with overlap of polygons. To get rid
# of overlap, we can use the spatial functionality that geopandas provides.
# Let's check the polygons for overlap first.

overlap = gdf.overlay(gdf, how="intersection", keep_geom_type=True)
overlap = overlap.loc[overlap["name_1"] != overlap["name_2"]]

fig, ax = plt.subplots()
gdf.plot(ax=ax)
overlap.plot(edgecolor="red", ax=ax)

# %%
# Clean-up
# --------
#
# There are many small overlaps visible at the province borders.
#
# We can generate a consistent polygon using a unary union.

union = sg.Polygon(gdf.unary_union)
union_gdf = gpd.GeoDataFrame(geometry=[union])
union_gdf["cellsize"] = 10_000.0

# %%
# Unfortunately, the province boundaries of this dataset no do align neatly and
# there are a number of small holes present. Some of these holes are not formed
# by inconsistencies, but by a small number of Belgian exclaves,
# `Baarle-Hertog`_.
#
# Simplify
# --------
#
# We'll ignore the subtleties of international law for now and use geopandas to
# remove all blemishes by:
#
# * squeezing out the holes with ``.buffer``
# * dissolving the buffered polygons into a single polygon with ``.dissolve``
# * simplifying the dissolved polygon to avoid over-refinement with ``.simplify``
#
# This creates a clean, and simpler, geometry.

simplified = gdf.copy()
simplified.geometry = simplified.geometry.buffer(500.0)
simplified["dissolve_column"] = 0
simplified = simplified.dissolve(by="dissolve_column")
simplified.geometry = simplified.geometry.simplify(5_000.0)
simplified["cellsize"] = 10_000.0

simplified.plot()

# %%
# Using this clean geometry, we can generate an unstructured grid with a fairly
# constant cell size.

mesher = pm.TriangleMesher(simplified)
vertices, triangles = mesher.generate()
pm.plot(vertices, triangles)

# %%
# For real work, buffering and simplifying will likely not suffice.
#
# See the preprocessing example to for an overview of common issues and how to
# apply pandamesh's Preprocessor class to resolve them.
#
# Local refinement
# ----------------
#
# To set a zone of refinement, we can define an additional polygon. We need to
# ensure that no overlap occurs in the follwing steps:
#
# * select the geometry of a single province;
# * simplify its geometry to an appropriate level of detail;
# * specify a smaller cell size;
# * remove this province from the enveloping polygon;
# * collect the two polygons in a single geodataframe.

utrecht = gdf.loc[["Utrecht"]]
utrecht.geometry = utrecht.geometry.simplify(2_500.0)
utrecht["cellsize"] = 5000.0

envelope = simplified.overlay(utrecht, how="difference")
refined = pd.concat([envelope, utrecht])
refined.index = [0, 1]
refined.plot(column="name")

# %%
# This results in a mesh with a smaller cell size in the province of Utrecht.

mesher = pm.TriangleMesher(refined)
vertices, triangles = mesher.generate()
pm.plot(vertices, triangles)

# %%
# Conclusion
# ----------
#
# This example provides a taste of how to convert a geospatial vector dataset
# into an unstructured grid with a locally refined part. Real-world data
# generally come with their own idiosyncrasies and inconsistencies. Depending
# on the nature of the necessary fixes, they can be solved with geopandas
# functionality, but sometimes manual editing is required. Fortunately,
# geopandas provides easy input and output for many file formats, which can be
# opened by e.g. QGIS.
#
# .. _Baarle-Hertog: https://en.wikipedia.org/wiki/Baarle-Hertog


"""
Preprocessing
=============

Raw geospatial vector data is often not ready to use directly in mesh
generation:
    
* Polygon data often do not form a valid planar partition: polygons are
  overlapping, or neighboring polygons have small gaps between them.
* Polygon boundaries or linestring segments intersect each other.
* Points may be located on polygon boundaries or lines. Since floating point
  numbers are not exact, points seemingly located on a line are computationally
  just left or just right of the line and form an extremely thin triangle.
* Points may be located extremely close together, thereby generating tiny
  triangles.
  
Such problems either lead to a generated mesh with extremely small elements, or
worse, they lead to a crash of the meshing program. Pandamesh provides a
``Preprocessor`` class to assist with cleaning up some common faults. 

This example will illustrate some common problems and how to resolve them.
"""
# %%
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import shapely
import shapely.geometry as sg

import pandamesh as pm

# sphinx_gallery_start_ignore
pm.GmshMesher.finalize()
# sphinx_gallery_end_ignore

# %%
# Polygons
# --------
#
# When generating a mesh, we often have a general area which may be meshed
# coarsely and an area of interest, which should be meshed more finely.
# Generally, the fine inner zone is located within the coarse outer zone, but
# this requires a hole in the outer zone that exactly matches up with the
# exterior of the inner zone.

outer = sg.Polygon(
    [
        [0.0, 0.0],
        [10.0, 0.0],
        [10.0, 10.0],
        [0.0, 10.0],
    ]
)
inner = sg.Polygon(
    [
        [5.0, 2.0],
        [8.0, 5.0],
        [5.0, 8.0],
        [2.0, 5.0],
    ]
)

gdf = gpd.GeoDataFrame(geometry=[outer, inner])
gdf["cellsize"] = [2.0, 1.0]

fig, (ax0, ax1) = plt.subplots(ncols=2, sharex=True, sharey=True)
gdf.iloc[[0]].plot(ax=ax0)
gdf.iloc[[1]].plot(ax=ax1)

# %%
# In this case, we have two conflicting specified cell sizes in the inner
# square. We can resolve this as follows:

resolved = (
    pm.Preprocessor(geometry=gdf.geometry, values=gdf.cellsize)
    .unify_polygons()
    .to_geodataframe()
).rename(columns={"values": "cellsize"})

# %%
# Note that the Preprocessor supports method chaining, allowing you to flexibly
# execute a set of operations.
#
# The resulting geodataframe's geometries are valid planar partition:

fig, (ax0, ax1) = plt.subplots(ncols=2, sharex=True, sharey=True)
resolved.iloc[[0]].plot(ax=ax0)
resolved.iloc[[1]].plot(ax=ax1)

# %%
# And we can use it directly to generate a mesh:

vertices, faces = pm.TriangleMesher(resolved).generate()
pm.plot(vertices, faces)

# %%
# Alternatively, multiple polygons with the same cell size specification might
# be overlapping

inner0 = shapely.affinity.translate(inner, xoff=-1.0)
inner1 = shapely.affinity.translate(inner, xoff=1.0)
gdf = gpd.GeoDataFrame(geometry=[outer, inner0, inner1])
gdf["cellsize"] = [2.0, 1.0, 1.0]

fig, ax = plt.subplots()
gdf.plot(ax=ax, facecolor="none")
# %%
# These will also be resolved by ``.unify_polygons``.

resolved = (
    pm.Preprocessor(geometry=gdf.geometry, values=gdf.cellsize)
    .unify_polygons()
    .to_geodataframe()
).rename(columns={"values": "cellsize"})

vertices, faces = pm.TriangleMesher(resolved).generate()

fig, ax = plt.subplots()
pm.plot(vertices, faces, ax=ax)
resolved.plot(ax=ax, facecolor="none", edgecolor="red")

# %%
# Note, however, that the internal boundaries of the inner polygons are forced
# into the triangulation. We can rid of these by calling ``.merge_polygons``:

resolved = (
    pm.Preprocessor(geometry=gdf.geometry, values=gdf.cellsize)
    .unify_polygons()
    .merge_polygons()
    .to_geodataframe()
).rename(columns={"values": "cellsize"})

vertices, faces = pm.TriangleMesher(resolved).generate()

fig, ax = plt.subplots()
pm.plot(vertices, faces, ax=ax)
resolved.plot(ax=ax, facecolor="none", edgecolor="red")

# %%
# An alternative problem is when polygons are touching, but do not actually
# share vertices along the boundary.

first = sg.Polygon(
    [
        [0.0, 0.0],
        [10.0, 0.0],
        [10.0, 10.0],
        [0.0, 10.0],
    ]
)
second = sg.Polygon(
    [
        [10.0, 2.0],
        [18.0, 2.0],
        [18.0, 8.0],
        [10.0, 8.0],
    ]
)

gdf = gpd.GeoDataFrame(geometry=[first, second])
gdf["cellsize"] = [4.0, 2.0]

vertices, faces = pm.GmshMesher(gdf, intersecting_edges="warn").generate(finalize=True)
pm.plot(vertices, faces)

# %%
# At x=10.0, the generated triangles are disconnected.
#
# This is caused by the the fact that the polygons do not share an edge:
#
# * The polygon on the left has an edge from (10.0, 0.0) to (10.0, 10.0)
# * The polygon on the right has an edge from (10.0, 2.0) to (10.0, 8.0)
#
# In fact, the vertices of the right polygon are intersecting the (edge) of the
# left polygon. We can identify these intersections with
# :func:`pandamesh.find_edge_intersections`:

intersections = pm.find_edge_intersections(gdf.geometry)

fig, ax = plt.subplots()
pm.plot(vertices, faces, ax=ax)
intersections.plot(ax=ax)

# %%
# Calling ``.unify_polygons()`` ensures that the vertices of touching polygons
# are inserted, such that the polygons share an edge.

resolved = (
    pm.Preprocessor(geometry=gdf.geometry, values=gdf.cellsize)
    .unify_polygons()
    .to_geodataframe()
).rename(columns={"values": "cellsize"})

vertices, faces = pm.TriangleMesher(resolved).generate()
polygon0_coords = shapely.get_coordinates(resolved.geometry[0])

fig, ax = plt.subplots()
pm.plot(vertices, faces, ax=ax)
ax.scatter(*polygon0_coords.T)

# %%
# Lines
# -----
#
# Lines may only be only partially present, or present in holes:

donut = sg.Polygon(
    [
        [0.0, 0.0],
        [10.0, 0.0],
        [10.0, 10.0],
        [0.0, 10.0],
    ],
    holes=[
        [
            [2.0, 5.0],
            [5.0, 8.0],
            [8.0, 5.0],
            [5.0, 2.0],
        ]
    ],
)
line0 = shapely.LineString(
    [
        [-2.0, 0.0],
        [12.0, 10.0],
    ]
)
line1 = shapely.LineString(
    [
        [5.5, 9.0],
        [9.0, 5.5],
    ]
)

gdf = gpd.GeoDataFrame(geometry=[donut, line0, line1])
gdf["cellsize"] = [2.0, 1.0, 1.0]
gdf.plot(edgecolor="k")

# %%
# We can identify these problematic intersections again using
# :func:`pandamesh.find_edge_intersections`:

intersections = pm.find_edge_intersections(gdf.geometry)
fig, ax = plt.subplots()
gdf.plot(ax=ax, facecolor="none")
intersections.plot(ax=ax)

# %%
# A first step is to remove line segments that do not fall in any polygon:

resolved = (
    pm.Preprocessor(geometry=gdf.geometry, values=gdf.cellsize)
    .clip_lines()
    .to_geodataframe()
).rename(columns={"values": "cellsize"})
resolved.plot(edgecolor="k")

# %%
# However, this doesn't create suitable input for meshing. The ``GmshMesher``
# appears to hang on this input, and Triangle generates a grid with very small
# triangles. Pandamesh errors on these intersections by default, but way may
# proceed:

vertices, faces = pm.TriangleMesher(resolved, intersecting_edges="warn").generate()
pm.plot(vertices, faces)

# %%
# A better approach here is to ensure all intersections are present in all
# linework:
#
# * First we clip.
# * Then we call ``unify_lines`` to ensure that the intersection of line0 and
#   line1 at (7.625 6.875) is represented.
# * Next we call ``unify_polygons``. This ensures the intersections of the lines
#   with the poygon exterior is represented as well.
# * The result of ``unify_polygons`` is that the line splits the polygon in two
#   parts. These are merged back together with ``merge_polygons``.
#
# If we plot the vertices of the resolved polygon, we see that the intersection
# vertices have been inserted into the polygon boundaries, and that the tiny
# triangles around the line intersection have disappeared:

resolved = (
    pm.Preprocessor(geometry=gdf.geometry, values=gdf.cellsize)
    .clip_lines()
    .unify_lines()
    .unify_polygons()
    .merge_polygons()
    .to_geodataframe()
).rename(columns={"values": "cellsize"})

vertices, faces = pm.GmshMesher(resolved).generate(finalize=True)
polygon0_coords = shapely.get_coordinates(resolved.geometry[0])

fig, ax = plt.subplots()
pm.plot(vertices, faces, ax=ax)
ax.scatter(*polygon0_coords.T)

# %%
# In some cases, having line segments terminate exactly on polygon boundaries
# still causes trouble. We may also ensure that lines are some distance removed
# from any polygon boundary by providing a distance to ``clip_lines``:

resolved = (
    pm.Preprocessor(geometry=gdf.geometry, values=gdf.cellsize)
    .unify_lines()
    .clip_lines(distance=0.5)
    .to_geodataframe()
).rename(columns={"values": "cellsize"})

vertices, faces = pm.GmshMesher(resolved).generate(finalize=True)
polygon0_coords = shapely.get_coordinates(resolved.geometry[0])

fig, ax = plt.subplots()
pm.plot(vertices, faces, ax=ax)
resolved.plot(facecolor="none", edgecolor="red", ax=ax)
ax.scatter(*polygon0_coords.T)

# %%
# Another pragmatic approach is to convert any line into interpolated points.
# Points cannot intersect each other, which sidesteps a large number of problems.

resolved = (
    pm.Preprocessor(geometry=gdf.geometry, values=gdf.cellsize)
    .interpolate_lines_to_points(distance=0.25)
    .clip_points()
    .to_geodataframe()
).rename(columns={"values": "cellsize"})

vertices, faces = pm.GmshMesher(resolved).generate(finalize=True)

fig, ax = plt.subplots()
pm.plot(vertices, faces, ax=ax)
resolved.plot(facecolor="none", edgecolor="red", ax=ax)

# %%
# Points
# ------
#
# Note that the start and end points of the lines are still on, or very near
# the polygon edges.
#
# We can remove those points by providing a distance to ``clip_points``.

resolved = (
    pm.Preprocessor(geometry=gdf.geometry, values=gdf.cellsize)
    .interpolate_lines_to_points(distance=0.25)
    .clip_points(distance=0.5)
    .to_geodataframe()
).rename(columns={"values": "cellsize"})

vertices, faces = pm.GmshMesher(resolved).generate(finalize=True)

fig, ax = plt.subplots()
pm.plot(vertices, faces, ax=ax)
resolved.plot(facecolor="none", edgecolor="red", ax=ax)

# %%
# A problem with points is that they may be very close together, thereby
# generating very small triangles. Let's generate 200 random points to illustrate:

rng = np.random.default_rng()
points = gpd.points_from_xy(*rng.random((2, 200)) * 10.0)
gdf = gpd.GeoDataFrame(geometry=np.concatenate([[donut], points]))
gdf["cellsize"] = 2.0

resolved = (
    pm.Preprocessor(geometry=gdf.geometry, values=gdf.cellsize)
    .clip_points(distance=0.5)
    .to_geodataframe()
).rename(columns={"values": "cellsize"})

vertices, faces = pm.GmshMesher(resolved).generate(finalize=True)
pm.plot(vertices, faces)

# %%
# We can solve this by snapping points together that are located some distance
# from each other:

resolved = (
    pm.Preprocessor(geometry=gdf.geometry, values=gdf.cellsize)
    .clip_points(distance=0.5)
    .snap_points(distance=0.5)
    .to_geodataframe()
).rename(columns={"values": "cellsize"})

vertices, faces = pm.GmshMesher(resolved).generate(finalize=True)
pm.plot(vertices, faces)

# %%
# Flexibility and composability
# -----------------------------
#
# The Preprocessor class in Pandamesh is designed with flexibility and
# composability in mind through method chaining. By combining various
# preprocessing steps in any order, you can address a wide range of geometric
# issues. For instance, you might start by unifying polygons, then clip lines,
# interpolate them to points, and finally snap those points together.
#
# The steps required depend on the nature of geometrical input, and may require
# experimenting with various methods. The intermediate output can be checked
# and visualized at any moments, by calling ``to_geodataframe``. For example,
# to check the intermediate result after clipping but prior to snapping:

check = (
    pm.Preprocessor(geometry=gdf.geometry, values=gdf.cellsize)
    .clip_points(distance=0.5)
    .to_geodataframe()
)

check.plot(facecolor="none")

# %%
# This also makes it easy to apply the preprocessor in steps. Some steps may be
# relatively costly, such as unifying a large number of detailed polygons. The
# intermediate result can be stored as e.g. a GeoPackage. Then, in a separate
# processing step, the intermediate result can be read again, and other
# processing steps (such as filtering points) can be applied.
