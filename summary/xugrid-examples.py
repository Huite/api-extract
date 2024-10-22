"""
Connectivity
============

A fundamental difference between structured and unstructured grids lies in the
connectivity. This is true for cell to cell connectivity, but also for vertex
(node) connectivity (which set of vertices make up an individual cell). In
structured grids, connectivity is implicit and can be directly derived from row
and column numbers; unstructured grids require explicit connectivity lists.

Xugrid provides a number of methods to derive and extract different kinds of
connectivities, as well as a number of operations which require connectivity
information. These methods and their interrelations are briefly introduced here.

For 2D meshes, the fundamental topological information consists of:

* A list of nodes (vertices): (x, y) coordinate pairs forming points.
* A list of faces (polygons): for every face, a list of index values indicating
  which vertices form its exterior.

Imports
-------

The following imports suffice for the examples.
"""
# %%

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

import xugrid

# %%
# Connectivity arrays
# -------------------
#
# From the fundamental face node connectivity, all other connectivities can be
# derived. These are accessible via the ``grid`` attribute of a XugridDataArray
# or XugridDataset. The are the available (derived) connectivity arrays are
# listed below. Depending on the (ir)regularity of the connectivity, the arrays
# are returned as either (dense) numpy arrays of integers, or as
# ``scipy.sparse.csr_matrix``.
#
# * ``face_node_connectivity``: dense ``(n_face, n_max_nodes_per_face)``
# * ``edge_node_connectivity``: dense ``(n_edge, 2)``
# * ``edge_face_connectivity``: dense ``(n_edge, 2)``
# * ``face_face_connectivity``: sparse
# * ``face_edge_connectivity``: sparse
# * ``node_edge_connectivity``: sparse
# * ``node_face_connectivity``: sparse
#
# Some connectivity arrays are returned in dense form, some in sparse. The
# ``node_edge_connectivity`` is the inverse of the ``edge_node_connectivity``.
# While the edge node connectivity array is very regular -- every edge is
# associated with just two nodes, the node edge connectivity is irregular: a
# node may be associated with just one edge or many and this requires many fill
# values in dense form.
#
# Binary erosion and dilation
# ---------------------------
#
# Binary erosion and dilation are useful operations to e.g. locate boundary
# cells, or to "shrink" some collection of cells. In this example, we start
# with a grid in which all cells are given a value of ``True`` (equal to
# ``1``).
#
# By default, the border value for binary erosion is set to ``False`` (equal to
# ``0``). This means the erosion erodes inwards from the boundaries.

ds = xugrid.data.disk()
uda = xugrid.UgridDataArray(
    xr.full_like(ds.obj["face_z"], True, dtype=bool),
    ds.grids[0],
)
iter2 = uda.ugrid.binary_erosion(iterations=2)
iter5 = uda.ugrid.binary_erosion(iterations=5)

fig, (ax0, ax1) = plt.subplots(ncols=2, figsize=(12, 5))
iter2.ugrid.plot(ax=ax0)
iter5.ugrid.plot(ax=ax1)

# %%
# By default, the border value for binary dilation is **also** set to
# ``False``. This means boundary does not dilate inwards by default.  We'll
# start by setting a single value in the center of the grid to ``True``.

uda = xugrid.UgridDataArray(
    xr.full_like(ds["face_z"].ugrid.obj, False, dtype=bool),
    ds.grids[0],
)
uda[0] = True
uda.ugrid.plot()

# %%
# Now let's run two dilations: one with the default border, and one with the
# alternative border value:

iter1 = uda.ugrid.binary_dilation(iterations=1)
iter1_boundary = uda.ugrid.binary_dilation(iterations=1, border_value=True)

fig, (ax0, ax1) = plt.subplots(ncols=2, figsize=(12, 5))
iter1.ugrid.plot(ax=ax0)
iter1_boundary.ugrid.plot(ax=ax1)

# %%
# Connected Components
# --------------------
#
# Xugrid also wraps :py:func:`scipy.sparse.csgraph.connected_components` to
# analyse connected parts of the mesh.

grid = xugrid.data.xoxo()
uda = xugrid.UgridDataArray(
    xr.DataArray(np.ones(grid.node_face_connectivity.shape[0]), dims=["face"]), grid
)
labeled = uda.ugrid.connected_components()
labeled.ugrid.plot(cmap="RdBu")

# %%
# Centroidal Voronoi Tesselation
# ------------------------------
#
# We can also use connectivity information to derive a centroidal Voronoi
# Tesselation.

voronoi_grid = grid.tesselate_centroidal_voronoi()
xugrid.plot.line(voronoi_grid, color="black")

# %%
# There are two alternative flavors to consider. We can fully ignore the
# exterior and consider only the (interior) centroids. Alternatively, we can
# include intersections of the voronoi edges with the mesh exterior, but
# ignore the original nodes.
#
# Both methods have the benefit of guaranteeing convex Voronoi polygons as
# their output -- provided the input mesh is convex as well! However, neither
# preserves the exterior exactly: the resulting mesh has smaller bounds than
# the original.

centroid_only = grid.tesselate_centroidal_voronoi(add_exterior=False)
convex_exterior = grid.tesselate_centroidal_voronoi(
    add_exterior=True, add_vertices=False
)

fig, (ax0, ax1) = plt.subplots(ncols=2, figsize=(12, 5))
xugrid.plot.line(centroid_only, ax=ax0, color="black")
xugrid.plot.line(convex_exterior, ax=ax1, color="black")

# %%
# Triangulation
# -------------
#
# Triangulation is a commonly required operation: every polygon can be split
# into triangles and triangles are the simplest geometric primitive. This makes
# them very attractive for e.g. visualization.
#
# We can break down one of the Voronoi tesselations from above into triangles:

triangulation = convex_exterior.triangulate()
fig, (ax0, ax1) = plt.subplots(ncols=2, figsize=(12, 5))
xugrid.plot.line(convex_exterior, ax=ax0, color="black")
xugrid.plot.line(triangulation, ax=ax1, color="black")

# %%
# Laplace interpolation
# ---------------------
#
# Laplace interpolation is a simple but powerful method to fill holes in a
# grid. Laplace's equation describes potential flow, such as e.g. steady-state
# heat conduction or steady-state groundwater flow. In this method, we solve
# Laplace's equation for the nodata gaps, with data values functioning as fixed
# potential boundary conditions.
#
# Let's setup a mesh with data exclusively on the left- and rightmost faces of
# the upper and lower parts:

grid = xugrid.data.xoxo()

da = xr.DataArray(
    np.full(283, np.nan),
    dims=[grid.face_dimension],
)
da.data[2] = 0.0
da.data[12] = 0.0
da.data[77] = 10.0
da.data[132] = 10.0

uda = xugrid.UgridDataArray(da, grid)

fig, ax = plt.subplots()
uda.ugrid.plot(ax=ax)
uda.ugrid.plot.line(ax=ax, color="black")

# %%
# We can now use Laplace interpolation to fill the gaps in the grid.

filled = uda.ugrid.laplace_interpolate()
filled.ugrid.plot(cmap="gist_rainbow", vmin=2.5, vmax=7.5)

# %%
# Laplace interpolation can also be used on the nodes of a grid.
# We start by removing 75% of the data. Then we fill it up again using
# interpolation.

disk_nodes = xugrid.data.disk()["node_z"]
disk_emptied = disk_nodes.where(disk_nodes["mesh2d_nNodes"] % 4 == 0)
disk_filled = disk_emptied.ugrid.laplace_interpolate(direct_solve=True)

fig, (ax0, ax1, ax2) = plt.subplots(ncols=3, figsize=(12, 3))
disk_emptied.ugrid.plot.scatter(ax=ax0)
disk_filled.ugrid.plot.scatter(ax=ax1)
disk_filled.ugrid.plot(ax=ax2)

# %%
# Reverse-Cuthill McKee
# ---------------------
#
# For numerical solutions, low "bandwidth" is desirable as this increases
# performance due to more efficient memory access. Xugrid wraps
# :py:func:`scipy.sparse.csgraph.reverse_cuthill_mckee` to reorder
# grids for bandwith reduction.
#
# To illustrate, let's take a look at the connectivity matrix of the Xoxo grid.

grid = xugrid.data.xoxo()
connectivity = grid.face_face_connectivity.toarray() != 0

fig, ax = plt.subplots(figsize=(8, 8))
ax.imshow(connectivity, cmap="Greys")

# %%
# The bandwidth of this matrix is poor. Connections are all over the place: low
# numbered cells are connected to high numbered cells (and vice versa). The
# bandwidth of the reordered grid is much smaller and has much better data
# locality:

renumbered_grid, _ = grid.reverse_cuthill_mckee()
connectivity = renumbered_grid.face_face_connectivity.toarray() != 0

fig, ax = plt.subplots(figsize=(8, 8))
ax.imshow(connectivity, cmap="Greys")

# %%


"""
OverlapRegridder
================

The overlap regridder works in two stages. First, it searches the source grid
for all faces of the target grid, computes the intersections, and stores all
overlaps between source and target faces. This occurs when the regridder is
initialized. Second, the regridder applies the weights: it reduces the
collection of overlapping faces to a single value for the target face.

There are many reductions possible. The best choice generally differs based on
the physical meaning of the variable, or the application. Xugrid provides a
number of reductions, but it's also possible to use a custom reduction
function. This is demonstrated here.

We start with the same example as in the quick overview.
"""
# %%

import matplotlib.pyplot as plt
import numpy as np

import xugrid as xu

# %%
# We'll use a part of a triangular grid with the surface elevation (including
# some bathymetry) of the Netherlands, and a coarser target grid.


def create_grid(bounds, nx, ny):
    """Create a simple grid of triangles covering a rectangle."""
    import numpy as np
    from matplotlib.tri import Triangulation

    xmin, ymin, xmax, ymax = bounds
    dx = (xmax - xmin) / nx
    dy = (ymax - ymin) / ny
    x = np.arange(xmin, xmax + dx, dx)
    y = np.arange(ymin, ymax + dy, dy)
    y, x = [a.ravel() for a in np.meshgrid(y, x, indexing="ij")]
    faces = Triangulation(x, y).triangles
    return xu.Ugrid2d(x, y, -1, faces)


uda = xu.data.elevation_nl().ugrid.sel(
    x=slice(125_000, 225_000), y=slice(440_000, 500_000)
)
grid = create_grid(uda.ugrid.total_bounds, nx=7, ny=6)

# %%

fig, ax = plt.subplots()
uda.ugrid.plot(vmin=-20, vmax=90, cmap="terrain", ax=ax)
grid.plot(ax=ax, color="red")

# %%
# Method comparison
# -----------------
#
# Let's compare the different reduction functions that are available in
# xugrid. We'll create a regridder once for every method, and plot the results
# side by side.
#
# .. note::
#   Sum and results in much higher values. The white in the figures are high
#   values, not no data. In contrast, a geometric mean generally only makes
#   sense for physical quantities with a "true zero": surface elevation is not
#   such quantity, as a datum is an arbitrary level. The xugrid geometric mean
#   returns NaN if reducing over negative values.

functions = [
    "mean",
    "harmonic_mean",
    "geometric_mean",
    "sum",
    "minimum",
    "maximum",
    "mode",
    "median",
    "max_overlap",
]

fig, axes = plt.subplots(nrows=5, ncols=2, figsize=(10, 25), sharey=True, sharex=True)
axes = axes.ravel()

for f, ax in zip(functions, axes):
    regridder = xu.OverlapRegridder(source=uda, target=grid, method=f)
    result = regridder.regrid(uda)
    result.ugrid.plot(vmin=-20, vmax=90, cmap="terrain", ax=ax)
    ax.set_title(f)

# %%
# Relative overlap
# ----------------
#
# For some reductions, the relative degree of overlap with the original source
# cell is required rather than the absolute overlap, e.g. for first-order
# conservative methods, such as conductance:

regridder = xu.RelativeOverlapRegridder(source=uda, target=grid, method="conductance")
result = regridder.regrid(uda)
result.ugrid.plot()

# %%
# Custom reductions
# -----------------
#
# It's also possible to define your own reduction methods. Such a method is
# inserted during the ``.regrid`` call and compiled by `Numba`_ for performance.
#
# A valid reduction method must be compileable by Numba, and takes exactly three
# arguments: ``values``, ``weights``, ``workspace``.
#
# * ``values``: is the array containing the (float) source values.
# * ``weights``: contains the (float) overlap between the target face and the
#   source faces. The size of ``weights`` is equal to the size of ``values``.
# * ``workspace``: used as a temporary workspace of floats. The size of ``work`` is
#   equal to the size of ``values``. (Make sure to zero it beforehand if that's
#   important to your reduction!)
#
# Xugrid regridder reduction functions are implemented in such a way. For a
# example, an area weighted sum could be implemented as follows:


def mean(values, weights, workspace):
    total = 0.0
    weight_sum = 0.0
    for value, weight in zip(values, weights):
        if ~np.isnan(value):
            total += value * weight
            weight_sum += weight
    if weight_sum == 0.0:
        return np.nan
    return total / weight_sum


# %%
# .. note::
#    * Each reduction must return a single float.
#    * Always check for ``np.isnan(value)``: Custom reductions methods must be
#      able to deal with NaN values as these are commonly encountered in datasets
#      as a "no data value".
#    * If Python features are used that are unsupported by Numba, you will get
#      somewhat obscure errors. In such a case, ``numba.njit`` and test your
#      function separately with synthetic values for ``values, weights,
#      workspace``.
#    * The ``workspace`` array is provided to avoid dynamic memory allocations.
#      It is a an array of floats with the same size as ``values`` or
#      ``weights``. You may freely allocate new arrays within the reduction
#      function but it will impact performance. (Methods such as mode or median
#      require a workspace.)
#    * While we could have implemented a weighted mean as:
#      ``np.nansum(values * weights) / np.nansum(weights)``, the function above
#      is efficiently compiled by Numba and does not allocate temporary arrays.
#
# To use our custom method, we provide it at initialization of the
# OverlapRegridder:

regridder = xu.OverlapRegridder(uda, grid, method=mean)
result = regridder.regrid(uda)
result.ugrid.plot(vmin=-20, vmax=90, cmap="terrain")

# %%
# Not every reduction uses the ``weights`` and ``workspace`` arguments. For
# example, a regular sum could only look at the values:


def nansum(values, weights, workspace):
    return np.nansum(values)


# %%
# Custom percentiles
# ------------------
#
# Xugrid provides a number of predefined percentiles (5, 10, 25, 50, 75, 90,
# 95). In case you need a different percentile value, you can use this utility:

p333 = xu.OverlapRegridder.create_percentile_method(33.3)

# %%
# Then, provide it as the regridder method as above:

regridder = xu.OverlapRegridder(uda, grid, method=p333)
result = regridder.regrid(uda)
result.ugrid.plot(vmin=-20, vmax=90, cmap="terrain")

# %%
# .. _Numba: https://numba.pydata.org/


"""
Partitioning
============

Grid partitioning, or domain decomposition, is an important step in setting up
parallellized simulations. Xugrid provides utilities for partitioning a grid
and its associated data, and for merging partitions back into a single whole.
"""
# %%

import matplotlib.pyplot as plt
import numpy as np

import xugrid as xu

# %%
# Create partitions
# -----------------
#
# Xugrid wraps the well known `METIS library`_ via the `pymetis bindings`_.
# METIS is generally used to partition a grid in such a manner that
# communication between parallel processes is minimized.
#
# We'll demonstrate the functionality by diving the elevation example
# into several parts.

uda = xu.data.elevation_nl()
uda.ugrid.plot(vmin=-20, vmax=90, cmap="terrain")

# %%
# Let's start by dividing the grid into four parts:

partitions = uda.ugrid.partition(n_part=4)

fig, axes = plt.subplots(2, 2, figsize=(12.6, 10))
for partition, ax in zip(partitions, axes.ravel()):
    partition.ugrid.plot(ax=ax, vmin=-20, vmax=90, cmap="terrain")

# %%
# Partition the grid
# ------------------
#
# Calling ``.partition`` on a UgridDataArray or UgridDataset will automatically
# partition the grid topology, select all associated data, and create a new
# UgridDataArray or UgridDataset for each partition.
#
# However, in some case, we might prefer to pre-compute the labels, and then
# apply them multiple datasets. To do so, we compute the partition labels from
# the grid. ``label_partitions`` returns a UgridDataArray, with every cell given
# its partition label number.
#
# We can easily plot this data to visualize the partitions:

labels = uda.ugrid.grid.label_partitions(n_part=12)
labels.ugrid.plot()

# %%
# Not quite the twelve provinces of the Netherlands!
#
# However, we may use the labels to partition the data nonetheless:

partitions = uda.ugrid.partition_by_label(labels)

fig, axes = plt.subplots(4, 3, figsize=(15, 15))
for partition, ax in zip(partitions, axes.ravel()):
    partition.ugrid.plot(ax=ax, vmin=-20, vmax=90, cmap="terrain")

# %%
# Since the labels are an ordinary UgridDataArray, we can easily store them in
# a netCDF file and re-use them in another part of a workflow.
#
# Merging partitions
# ------------------
#
# Generally, after partitioning the data we write it as model input and run a
# model in parallel. Many model codes produce output per process. Xugrid can
# merge these partitions back into one whole for post-processing:

merged = xu.merge_partitions(partitions)["elevation"]

merged.ugrid.plot(vmin=-20, vmax=90, cmap="terrain")

# %%
# Partitioning grids without data
# -------------------------------
#
# Of course, we can also partition the grid topology without any associated
# data:

grid = uda.ugrid.grid
grid_parts = grid.partition(n_part=4)

fig, axes = plt.subplots(2, 2, figsize=(12.6, 10))
for part, ax in zip(grid_parts, axes.ravel()):
    part.plot(ax=ax)

# %%
# ... and merge them back into one whole:

merged_grid, _ = xu.Ugrid2d.merge_partitions(grid_parts)
merged_grid.plot()

# %%
# Preserving order
# ----------------
#
# Note that partioning and merging does not preserve order!

uda == merged

# %%
# The topology is equivalent, but the nodes, edges, and faces are in a
# different order. This is because ``merge_partitions`` simply concatenates the
# partitions.
#
# The easiest way to restore the order is by providing an example of the
# original topology. ``reindex_like`` looks at the coordinates of both
# (equivalent!) grids and automatically determines how to reorder:

reordered = merged.ugrid.reindex_like(uda)
uda == reordered

# %%
# Alternatively, we can also assign IDs, carry these along, and use these to
# reorder the data after merging.

uds = xu.UgridDataset(grids=[uda.ugrid.grid])
uds["elevation"] = uda
uds["cell_id"] = ("mesh2d_nFaces", np.arange(len(uda)))

partitions = uds.ugrid.partition(n_part=4)
merged = xu.merge_partitions(partitions)
order = np.argsort(merged["cell_id"].values)
reordered = merged.isel(mesh2d_nFaces=order)

uds["elevation"] == reordered["elevation"]

# %%
# This is required if results are compared with the input, or with results
# stemming from another partitioning, e.g. one with a different number of
# partitions.
#
# .. _METIS library: https://github.com/KarypisLab/METIS
# .. _pymetis bindings: https://github.com/inducer/pymetis


"""
Plot unstructured mesh data
===========================

The labels that are present in xarray's data structures allow for easy creation
of informative plots: think of dates on the x-axis, or geospatial coordinates.
Xarray provides a convenient way of plotting your data provided it is
structured. Xugrid extends these plotting methods to easily make spatial
(x-y) plots of unstructured grids.

Like Xarray's focus for plotting is the DataArray, Xugrid's focus is the
UgridDataArray; like Xarray, if your (extracted) data fits into a pandas
DataFrame, you're better of using pandas tools instead.

As every other method in Xugrid, any logic involving the unstructured topology
is accessed via the ``.ugrid`` accessor on the DataArrays and Datasets;
UgridDatasets and UgridDataArrays behave the same as ordinary Xarray DataArrays
and Datasets otherwise.

Imports
-------

The following imports suffice for the examples.
"""

# %%
import matplotlib.pyplot as plt

import xugrid

# %%
# We'll use a simple synthetic example. This dataset contains data for all
# topological attributes of a two dimensional mesh:
#
# * Nodes: the coordinate pair (x, y) forming a point.
# * Edges: a line or curve bounded by two nodes.
# * Faces: the polygon enclosed by a set of edges.
#
# In this disk example, very similar has been placed on the nodes, edges, and
# faces.

ds = xugrid.data.disk()
ds

# %%
# UgridDataArray
# --------------
#
# Just like Xarray, we can create a plot by selecting a DataArray from the
# Dataset and calling the :py:meth:`UgridDataArray.ugrid.plot()` method.

uda = ds["face_z"]
uda.ugrid.plot()

# %%
# Like Xarray, the axes and the colorbar are labeled automatically using the
# available information.
#
# The convenience function :py:meth:`xugrid.UgridDataArray.ugrid.plot()`
# dispatches on the topological dimension of the variable. In this case, the
# data is associated with the face dimension of the topology. Data located on
# the edges results in a different kind of plot:

ds["edge_z"].ugrid.plot()

# %%
# The method called by default depends on the type of the data:
#
# =============== ===========================
# Dimension       Plotting function
# =============== ===========================
# Face            :py:func:`xugrid.plot.pcolormesh`
# Edge            :py:func:`xugrid.plot.line`
# Node            :py:func:`xugrid.plot.tripcolor`
# =============== ===========================
#
# We can put them side by side to illustrate the differences:

fig, (ax0, ax1, ax2) = plt.subplots(ncols=3, figsize=(11, 3), sharex=True, sharey=True)
ds["face_z"].ugrid.plot(ax=ax0)
ds["edge_z"].ugrid.plot(ax=ax1)
ds["node_z"].ugrid.plot(ax=ax2)

# %%
# We can also exactly control the type of plot we want. For example, to plot
# filled contours for data associated with the face dimension:

ds["face_z"].ugrid.plot.contourf()

# %%
# We can also overlay this data with the edges:

fig, ax = plt.subplots()
ds["face_z"].ugrid.plot.contourf()
ds["face_z"].ugrid.plot.line(color="black")

# %%
# In general, there has to be data associated with the mesh topology before a
# plot can be made. ``plot.line()`` forms an exception to this rule, as the
# location of the edges is meaningful on its own: for this reason
# ``plot.line()`` does not error in the example above.
#
# Other types of plot
# -------------------
#
# The available plotting methods per topology dimension are listed here.
#
# For the **face** dimension:
#
# * :py:func:`xugrid.plot.contour`
# * :py:func:`xugrid.plot.contourf`
# * :py:func:`xugrid.plot.imshow`
# * :py:func:`xugrid.plot.pcolormesh`
# * :py:func:`xugrid.plot.scatter`
# * :py:func:`xugrid.plot.surface`
#
# For the **edge** dimension:
#
# * :py:func:`xugrid.plot.line`
# * :py:func:`xugrid.plot.scatter`
#
# For the **node** dimension:
#
# * :py:func:`xugrid.plot.contour`
# * :py:func:`xugrid.plot.contourf`
# * :py:func:`xugrid.plot.scatter`
# * :py:func:`xugrid.plot.surface`
# * :py:func:`xugrid.plot.tripcolor`
#
# All these (2D) plots are illustrated here for completeness' sake:

fig, axes = plt.subplots(nrows=5, ncols=3, figsize=(10, 15))

ds["face_z"].ugrid.plot.pcolormesh(ax=axes[0, 0])
ds["face_z"].ugrid.plot.contour(ax=axes[1, 0])
ds["face_z"].ugrid.plot.contourf(ax=axes[2, 0])
ds["face_z"].ugrid.plot.imshow(ax=axes[3, 0])
ds["face_z"].ugrid.plot.scatter(ax=axes[4, 0])

ds["edge_z"].ugrid.plot.line(ax=axes[0, 1])
ds["edge_z"].ugrid.plot.scatter(ax=axes[4, 1])

ds["node_z"].ugrid.plot.tripcolor(ax=axes[0, 2])
ds["node_z"].ugrid.plot.contour(ax=axes[1, 2])
ds["node_z"].ugrid.plot.contourf(ax=axes[2, 2])
ds["node_z"].ugrid.plot.scatter(ax=axes[4, 2])

# %%
# The ``surface`` methods generate 3D surface plots:

fig = plt.figure(figsize=plt.figaspect(0.5))
ax0 = fig.add_subplot(1, 2, 1, projection="3d")
ax1 = fig.add_subplot(1, 2, 2, projection="3d")
ds["face_z"].ugrid.plot.surface(ax=ax0)
ds["node_z"].ugrid.plot.surface(ax=ax1)

# %%
# Additional Arguments
# --------------------
#
# Once again like in Xarray, additional arguments are passed to the underlying
# matplotlib function and the additional arguments supported by Xarray can be
# used:

ds["face_z"].ugrid.plot(cmap="RdBu", levels=8, yincrease=False)

# %%
# As a function
# -------------
#
# The plotting methods can also be called as a function, in which case they
# take an xarray DataArray and a xugrid grid as arguments.

grid = ds.ugrid.grids[0]
da = ds.obj["face_z"]

xugrid.plot.pcolormesh(grid, da)

# %%
# Xarray DataArray plots
# ----------------------
#
# As mentioned, apart from the ``.ugrid`` accessor, a UgridDataArray behaves the
# same as an Xarray DataArray. To illustrate, we can select a location
# somewhere in the unstructured topology, and plot the resulting timeseries:

ds = xugrid.data.adh_san_diego()
depth = ds["depth"]
depth.isel(node=1000).plot()

# %%


"""
Quick overview
==============

Here are a number of quick examples of how to get started with xugrid. More
detailed explanation can be found in the rest of the documentation.

We'll start by importing a few essential packages.
"""
# %%

import numpy as np
import xarray as xr

import xugrid as xu

# %%
# Create a UgridDataArray
# -----------------------
#
# There are three ways to create a UgridDataArray:
#
# * From an xarray Dataset containing the grid topology stored according to the
#   UGRID conventions.
# * From a xugrid Ugrid object and an xarray DataArray containing the data.
# * From a UGRID netCDF file, via :py:func:`xugrid.open_dataset`.
#
#
# From xarray Dataset
# ~~~~~~~~~~~~~~~~~~~
#
# xugrid will automatically find the UGRID topological variables, and separate
# them from the main data variables.
#
# Details on the required variables can be found in the `UGRID conventions`_.
# For 1D and 2D UGRID topologies, the required variables are:
#
# * x-coordinates of the nodes
# * y-coordinates of the nodes
# * edge node connectivity (1D) or face node connectivity (2D)
# * a "dummy" variable storing the names of the above variables in its
#   attributes
#
# We'll start by fetching a dataset:

ds = xu.data.adh_san_diego(xarray=True)
ds

# %%
# There are a number of topology coordinates and variables: ``node_x`` and
# ``node_y``, ``mesh2d`` and ``face_node_connectivity``. The dummy variable
# is ``mesh2d`` contains only a 0 for data; its attributes contain a mapping of
# UGRID roles to dataset variables.
#
# We can convert this dataset to a UgridDataset which will automatically
# separate the variables:

uds = xu.UgridDataset(ds)
uds

# %%
# We can then grab one of the data variables as usual for xarray:

elev = uds["elevation"]
elev

# %%
# From Ugrid and DataArray
# ~~~~~~~~~~~~~~~~~~~~~~~~
#
# Alternatively, we can build a Ugrid topology object first from vertices and
# connectivity numpy arrays, for example when using the topology data generated
# by a mesh generator (at which stage there is no data asssociated with the
# nodes, edges, or faces).
#
# There are many ways to construct such arrays, typically via mesh generators
# or Delaunay triangulation, but we will construct two simple triangles and
# some data by hand here:

nodes = np.array([[0, 0], [0, 1.1], [1, 0], [1, 1]])
faces = np.array([[2, 3, 0], [3, 1, 0]])
fill_value = -1

grid = xu.Ugrid2d(nodes[:, 0], nodes[:, 1], fill_value, faces)
da = xr.DataArray(
    data=[1.0, 2.0],
    dims=[grid.face_dimension],
)
uda = xu.UgridDataArray(da, grid)
uda

# %%
# From netCDF file
# ~~~~~~~~~~~~~~~~
#
# :py:func:`xugrid.open_dataset` is demonstrated in the last section of this
# guide. Internally, it opens the netCDF as a regular dataset, then converts it
# as seen in the first example.
#
# Plotting
# --------

elev.ugrid.plot(cmap="viridis")

# %%
# Data selection
# --------------
#
# A UgridDataArray behaves identical to an xarray DataArray:

whole = xu.data.disk()["face_z"]

# %%
# To select based on the topology, use the ``.ugrid`` attribute:

subset = whole.ugrid.sel(y=slice(5.0, None))
subset.ugrid.plot()

# %%
# .. note::
#
#   ``ugrid.sel()`` currently only supports data on the faces for 2D
#   topologies, and data on edges for 1D topologies. More flexibility
#   may be added.
#
# Computation
# -----------
#
# Computation on DataArrays is unchanged from xarray:

uda + 10.0

# %%
# Geopandas
# ---------
#
# Xugrid objects provide a number of conversion functions from and to geopandas
# GeoDataFrames using :py:meth:`xugrid.UgridDataset.from_geodataframe`. Note
# that storing large grids as GeoDataFrames can be very inefficient.

gdf = uda.ugrid.to_geodataframe(name="test")
gdf

# %%
# Conversion from Geopandas is easy too:

xu.UgridDataset.from_geodataframe(gdf)

# %%
# XugridDatasets
# --------------
#
# Like an Xarray Dataset, a UgridDataset is a dict-like container of
# UgridDataArrays. It is required that they share the same grid topology;
# but the individual DataArrays may be located on different aspects of the
# grid (nodes, faces, edges).

xu.data.disk()

# %%
# A UgridDataset may be initialized without data variables, but this requires
# a grid object:

new_uds = xu.UgridDataset(grids=uds.ugrid.grids)
new_uds

# %%
# We can then add variables one-by-one, as we might with an xarray Dataset:

new_uds["elevation"] = elev
new_uds

# %%
# Write netCDF files
# ------------------
#
# Once again like xarray, NetCDF is the recommended file format for xugrid
# objects. Xugrid automatically stores the grid topology according to the UGRID
# conventions and merges it with the main dataset containing the data variables
# before writing.

uds.ugrid.to_netcdf("example-ugrid.nc")
xu.open_dataset("example-ugrid.nc")

# %%
# .. _UGRID Conventions: https://ugrid-conventions.github.io/ugrid-conventions


"""
Regridding overview
===================

`Regridding`_ is the process of converting gridded data from one grid to
another grid. Xugrid provides tools for 2D and 3D regridding of structured
gridded data, represented as xarray objects, as well as (`layered`_)
unstructured gridded data, represented as xugrid objects.

A number of regridding methods are provided, based on area or volume overlap,
as well as interpolation routines. It currently only supports Cartesian
coordinates. See e.g. `xESMF`_ instead for regridding with a spherical Earth
representation (note: EMSF is `not available`_ via conda-forge on Windows).

Here are a number of quick examples of how to get started with regridding.

We'll start by importing a few essential packages.
"""
# %%

import matplotlib.pyplot as plt
import xarray as xr

import xugrid as xu

# %%
# We will take a look at a sample dataset: a triangular grid with the surface
# elevation of the Netherlands.

uda = xu.data.elevation_nl()
uda.ugrid.plot(vmin=-20, vmax=90, cmap="terrain")

# %%
# Xugrid provides several "regridder" classes which can convert gridded data
# from one grid to another grid. Let's generate a very simple coarse mesh that
# covers the entire Netherlands.


def create_grid(bounds, nx, ny):
    """Create a simple grid of triangles covering a rectangle."""
    import numpy as np
    from matplotlib.tri import Triangulation

    xmin, ymin, xmax, ymax = bounds
    dx = (xmax - xmin) / nx
    dy = (ymax - ymin) / ny
    x = np.arange(xmin, xmax + dx, dx)
    y = np.arange(ymin, ymax + dy, dy)
    y, x = [a.ravel() for a in np.meshgrid(y, x, indexing="ij")]
    faces = Triangulation(x, y).triangles
    return xu.Ugrid2d(x, y, -1, faces)


grid = create_grid(uda.ugrid.total_bounds, 7, 7)

# %%
# CentroidLocatorRegridder
# ------------------------
#
# An easy way of regridding is by simply looking in which cell of the original
# the centroids of the new grid fall.

fig, ax = plt.subplots()
uda.ugrid.plot(vmin=-20, vmax=90, cmap="terrain", ax=ax)
grid.plot(ax=ax, color="red")
ax.scatter(*grid.centroids.T, color="red")

# %%
# Xugrid provides the CentroidLocatorRegridder for this:

regridder = xu.CentroidLocatorRegridder(source=uda, target=grid)
result = regridder.regrid(uda)
result.ugrid.plot(vmin=-20, vmax=90, cmap="terrain", edgecolor="red")

# %%
# OverlapRegridder
# ----------------
#
# Such a regridding is not appropriate when the new grid cells are
# so large. Let's try the OverlapOverregridder instead.

regridder = xu.OverlapRegridder(source=uda, target=grid)
mean = regridder.regrid(uda)
mean.ugrid.plot(vmin=-20, vmax=90, cmap="terrain", edgecolor="red")

# %%
# By default, the OverlapRegridder computes an area weighted mean.
# Let's try again, now with the minimum:

regridder = xu.OverlapRegridder(source=uda, target=grid, method="minimum")
minimum = regridder.regrid(uda)
minimum.ugrid.plot(vmin=-20, vmax=90, cmap="terrain", edgecolor="red")

# %%
# Or the maximum:

regridder = xu.OverlapRegridder(source=uda, target=grid, method="maximum")
maximum = regridder.regrid(uda)
maximum.ugrid.plot(vmin=-20, vmax=90, cmap="terrain", edgecolor="red")

# %%
# All regridders also work for multi-dimensional data.
#
# Let's pretend our elevation dataset contains multiple layers, for example to
# denote multiple geological strata. We'll generate five layers, each with a
# thickness of 10.0 meters.

thickness = xr.DataArray(
    data=[10.0, 10.0, 10.0, 10.0, 10.0],
    coords={"layer": [1, 2, 3, 4, 5]},
    dims=["layer"],
)

# %%
# We need to make that the face dimension remains last, so we transpose the
# result.

bottom = (uda - thickness.cumsum("layer")).transpose()
bottom

# %%
# We can feed the result to the regridder, which will automatically regrid over
# all additional dimensions.

mean_bottom = xu.OverlapRegridder(source=bottom, target=grid).regrid(bottom)
mean_bottom

# %%
# Let's take a slice to briefly inspect our original layer bottom elevation,
# and the aggregated mean.

section_y = 475_000.0
section = bottom.ugrid.sel(y=section_y)
section_mean = mean_bottom.ugrid.sel(y=section_y)

fig, (ax0, ax1) = plt.subplots(ncols=2, figsize=(10, 5), sharex=True, sharey=True)
section.plot.line(x="mesh2d_s", hue="layer", ax=ax0)
section_mean.plot.line(x="mesh2d_s", hue="layer", ax=ax1)

# %%
# BarycentricInterpolator
# -----------------------
#
# All examples above show reductions: from a fine grid to a coarse grid.
# However, xugrid also provides interpolation to generate smooth fine
# representations of a coarse grid.
#
# To illustrate, we will zoom in to a part of the Netherlands.

part = uda.ugrid.sel(x=slice(125_000, 225_000), y=slice(440_000, 500_000))
part.ugrid.plot(vmin=-20, vmax=90, cmap="terrain")

# %%
# We can clearly identify the individual triangles that form the grid. To get a
# smooth presentation, we can use the BarycentricInterpolator.
#
# We will generate a fine grid.

grid = create_grid(part.ugrid.total_bounds, nx=100, ny=100)

# %%
# We use the centroids of the fine grid to interpolate between the centroids of
# the triangles.

regridder = xu.BarycentricInterpolator(part, grid)
interpolated = regridder.regrid(part)
interpolated.ugrid.plot(vmin=-20, vmax=90, cmap="terrain")

# %%
# Arbitrary grids
# ---------------
#
# The above examples all feature triangular source and target grids. However,
# the regridders work for any collection of (convex) faces.

grid = create_grid(part.ugrid.total_bounds, nx=20, ny=15)
voronoi_grid = grid.tesselate_centroidal_voronoi()

regridder = xu.CentroidLocatorRegridder(part, voronoi_grid)
result = regridder.regrid(part)

fig, ax = plt.subplots()
result.ugrid.plot(vmin=-20, vmax=90, cmap="terrain")
voronoi_grid.plot(ax=ax, color="red")

# %%
# Re-use
# ------
#
# The most expensive step of the regridding process is finding and computing
# overlaps. A regridder can be used repeatedly, provided the source topology
# is kept the same.

part_other = part - 50.0
result = regridder.regrid(part_other)
result.ugrid.plot(vmin=-20, vmax=90, cmap="terrain")

# %%
# .. _Xarray: https://docs.xarray.dev/en/stable/index.html
# .. _Xugrid: https://deltares.github.io/xugrid/
# .. _Regridding: https://climatedataguide.ucar.edu/climate-tools/regridding-overview
# .. _layered: https://ugrid-conventions.github.io/ugrid-conventions/#3d-layered-mesh-topology
# .. _xESMF: https://xesmf.readthedocs.io/en/latest/index.html
# .. _not available: https://github.com/conda-forge/esmf-feedstock/issues/64

# %%


"""
Select unstructured data
========================

Xarray has flexible tools for label based selection, in the form of ``.sel``
and ``.isel`` for index selection. This works well for structured data since
the orthogonality of the x and y axes is reflected in the axes of the
underlying arrays. This orthogonality does not exist for unstructured grids, as
the data for all faces cannot be stored in a two-dimensional array and is
stored in a one-dimensional array instead.

Xugrid provides tools for convenient spatial selection, primarily via the
``.ugrid.sel`` method; its behavior is comparable to xarray's ``.sel`` method.
The ``.ugrid.sel`` method should only be used for selection in the x or y
dimension. Selections along other dimension (such as time) should be performed
by xarray's ``.sel`` instead (without the ``ugrid`` accessor).

The examples below demonstrate the various ways to select data.

Imports
-------

The following imports suffice for the examples.
"""
# %%
import matplotlib.pyplot as plt
import numpy as np
import shapely

import xugrid as xu

# %%
# We will take a look at a sample dataset: a triangular grid with the surface
# elevation of the Netherlands.

uda = xu.data.elevation_nl()
uda.ugrid.plot(vmin=-20, vmax=90, cmap="terrain")

# %%
# We will start by demonstrating the behavior of ``.ugrid.sel``. This method
# takes several types of arguments, like its xarray equivalent. The return type
# and shape of the selection operation depends on the argument given.
#
# ========== ===========
# Selection  Result type
# ========== ===========
# Subset     xugrid
# Point      xarray
# Line       xarray
# ========== ===========
#
# Grid subset selection
# ---------------------
#
# A subset of the unstructured grid is returned by using slices without a step:

subset = uda.ugrid.sel(x=slice(100_000.0, 200_000.0), y=slice(450_000.0, 550_000.0))
subset.ugrid.plot(vmin=-20, vmax=90, cmap="terrain")

# %%
# The default arguments of ``x`` and ``y`` are: ``slice(None, None)``.
# In such a case the entire grid is returned.

subset = uda.ugrid.sel()
subset.ugrid.plot(vmin=-20, vmax=90, cmap="terrain")

# %%
# .. note::
#
#   ``None`` in a Python slice can be interpreted as "from the start" or "up to
#   and including the end".
#
# This means we can easily select along a single dimension:

subset = uda.ugrid.sel(x=slice(100_000.0, 200_000.0))
subset.ugrid.plot(vmin=-20, vmax=90, cmap="terrain", aspect=1, size=5)

# %%
# Or, using ``None`` if we only care about the start:

subset = uda.ugrid.sel(x=slice(100_000.0, None))
subset.ugrid.plot(vmin=-20, vmax=90, cmap="terrain", aspect=1, size=5)

# %%
# Point selection
# ---------------
#
# Since point data can be represented as an ordinary xarray DataArray with x
# and y coordinates, all point selection result in xarray DataArrays rather
# than UgridDataArrays with an associated unstructured grid topology.
#
# We will use a utility function to show what is selected on the map:


def show_point_selection(uda, da):
    _, ax = plt.subplots()
    uda.ugrid.plot(ax=ax, vmin=-20, vmax=90, cmap="terrain")
    ax.scatter(da["mesh2d_x"], da["mesh2d_y"], color="red")
    ax.set_aspect(1.0)


# %%
# Two values will select a point:

da = uda.ugrid.sel(x=150_000.0, y=463_000.0)
show_point_selection(uda, da)
da

# %%
# Multiple values are broadcasted against each other ("outer indexing").
# If we select by three x values and two y values, the result is a collection
# of six points:

da = uda.ugrid.sel(x=[125_000.0, 150_000.0, 175_000.0], y=[400_000.0, 465_000.0])
show_point_selection(uda, da)
da

# %%
# To select points without broadcasting, use ``.ugrid.sel_points`` instead:

da = uda.ugrid.sel_points(
    x=[125_000.0, 150_000.0, 175_000.0], y=[400_000.0, 430_000.0, 465_000.0]
)
show_point_selection(uda, da)
da

# %%
# We can sample points along a line as well by providing slices **with** a step:

da = uda.ugrid.sel(x=slice(100_000.0, 200_000.0, 10_000.0), y=465_000.0)
show_point_selection(uda, da)
da

# %%
# Two slices with a step results in broadcasting:

da = uda.ugrid.sel(
    x=slice(100_000.0, 200_000.0, 10_000.0), y=slice(400_000.0, 500_000.0, 10_000.0)
)
show_point_selection(uda, da)
da

# %%
# As well as a slice with a step and multiple values:

da = uda.ugrid.sel(x=slice(100_000.0, 200_000.0, 10_000.0), y=[400_000.0, 430_000.0])
show_point_selection(uda, da)
da

# %%
# Line selection
# --------------
#
# Since line data can be represented as an ordinary xarray DataArray with x
# and y coordinates, all line selection result in xarray DataArrays rather
# than UgridDataArrays with an associated unstructured grid topology.
#
# Line selection is performed by finding all faces that are intersected by
# the line.
#
# We start by defining a utility to show the selection again:


def show_line_selection(uda, da, line_x=None, line_y=None):
    _, (ax0, ax1) = plt.subplots(ncols=2, figsize=(10, 5))
    uda.ugrid.plot(ax=ax0, vmin=-20, vmax=90, cmap="terrain")
    da.plot(ax=ax1, x="mesh2d_s")
    if line_x is None:
        ax0.axhline(line_y, color="red")
    elif line_y is None:
        ax0.axvline(line_x, color="red")
    else:
        ax0.plot(line_x, line_y, color="red")
    ax0.set_aspect(1.0)


# %%
# A single value for either x or y in ``.ugrid.sel`` will select values along a
# line:

da = uda.ugrid.sel(y=465_000.0)
show_line_selection(uda, da, line_y=465_000.0)

# %%
# Line segments that are not axis aligned can be selected with
# ``.ugrid.intersect_line``:

da = uda.ugrid.intersect_line(start=(60_000.0, 400_000.0), end=(190_000.0, 475_000.0))
show_line_selection(uda, da, (60_000.0, 190_000.0), (400_000.0, 475_000.0))

# %%
# Linestrings can be selected with ``.ugrid.intersect_linestring``:

linestring = shapely.geometry.LineString(
    [
        (60_000.0, 400_000.0),
        (190_000.0, 400_000.0),
        (120_000.0, 575_000.0),
        (250_000.0, 575_000.0),
    ]
)

da = uda.ugrid.intersect_linestring(linestring)
show_line_selection(uda, da, *shapely.get_coordinates(linestring).T)

# %%
# This will work for any type of shapely line:

ring = shapely.geometry.Point(155_000.0, 463_000).buffer(50_000.0).exterior
da = uda.ugrid.intersect_linestring(ring)
show_line_selection(uda, da, *shapely.get_coordinates(ring).T)

# %%
# Index selection
# ---------------
#
# We may also use ordinary index selection to create a subset. This does not
# require the ``.ugrid`` accessor. For example, to take only the first
# thousands faces:

subset = uda.isel(mesh2d_nFaces=np.arange(1000))
subset.ugrid.plot(vmin=-20, vmax=90, cmap="terrain", aspect=1, size=5)

# %%
# For a 2D topology, selecting faces by an index always results in a valid
# topology. However, selecting by node or edge does not give a guarantee that
# the result forms a valid 2D topology: e.g. if we only select two nodes, or
# only two edges from a face, the result cannot form a valid 2D face.
#
# To avoid generating invalid topologies, xugrid always checks whether the
# result of a selection results in a valid 2D topology and raises an error if
# the result is invalid.
#
# In general, index selection should only be performed on the "core" dimension
# of the UGRID topology. This is the edge dimension for 1D topologies, and the
# face dimension for 2D topologies.

# %%


"""
Vector geometry conversion
==========================

A great deal of geospatial data is available not in gridded form, but in
vector form: as points, lines, and polygons. In the Python data ecosystem,
these geometries and their associated data are generally represented by a
geopandas GeoDataFrame.

Xugrid provides a number of utilities to use such data in combination with
unstructured grids. These are demonstrated below.
"""
# %%

import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd

import xugrid as xu

# %%
# We'll once again use the surface elevation data example.

uda = xu.data.elevation_nl()
uda.ugrid.plot(vmin=-20, vmax=90, cmap="terrain")

# %%
# Conversion to GeoDataFrame
# --------------------------
#
# A UgridDataArray or UgridDataset can be directly converted to a GeoDataFrame,
# provided it only contains a spatial dimension (and not a dimension such as
# time). When calling
# ``.to_geodataframe``, a shapely Polygon is created for every face (cell).

gdf = uda.ugrid.to_geodataframe()
print(gdf)

# %%
# We see that a GeoDataFrame with 5248 rows is created: one row for each face.
#
# Conversion from GeoDataFrame
# ----------------------------
#
# We can also make the opposite conversion: we can create a UgridDataSet from a
# GeoDataFrame.
#
back = xu.UgridDataset.from_geodataframe(gdf)
back

# %%
# .. note::
#   Not every GeoDataFrame can be converted to a ``xugrid`` representation!
#   While an unstructured grid topology is generally always a valid collection
#   of polygon geometries, not every collection of polygon geometries is a
#   valid grid: polygons should be convex and non-overlapping to create a valid
#   unstructured grid.
#
#   Secondly, each polygon fully owns its vertices (nodes), while the face of a
#   UGRID topology shares its nodes with its neighbors. All the vertices of the
#   polygons must therefore be exactly snapped together to form a connected
#   mesh.
#
#   Hence, the ``.from_geodataframe()`` is primarily meant to create ``xugrid``
#   objects from data that were originally created as triangulation or
#   unstructured grid, but that were converted to vector geometry form.
#
# "Rasterizing", or "burning" vector geometries
# ---------------------------------------------
#
# Rasterizing is a common operation when working with raster and vector data.
# While we cannot name the operation "rasterizing" when we're dealing with
# unstructured grids, there is a clearly equivalent operation where we mark
# cells that are covered or touched by a polygon.
#
# In this example, we mark the faces that are covered by a certain province.
#
# We start by re-projecting the provinces dataset to the coordinate reference
# system (CRS), from WGS84 (EPSG:4326) to the Dutch National coordinate system
# (RD New, EPSG: 28992). Then, we give each province a unique id, which we
# burn into the grid.

provinces = xu.data.provinces_nl().to_crs(28992)
provinces["id"] = range(len(provinces))
burned = xu.burn_vector_geometry(provinces, uda, column="id")
burned.ugrid.plot()

# %%
# This makes it very easy to classify and group data. Let's say
# we want to compute the average surface elevation per province:

burned = xu.burn_vector_geometry(provinces, uda, column="id")
uda.groupby(burned).mean()


# %%
# This is a convenient way to create masks for specific regions:

utrecht = provinces[provinces["name"] == "Utrecht"]
burned = xu.burn_vector_geometry(utrecht, uda)
xmin, ymin, xmax, ymax = utrecht.buffer(10_000).total_bounds

fig, ax = plt.subplots()
burned.ugrid.plot(ax=ax)
burned.ugrid.plot.line(ax=ax, edgecolor="black", linewidth=0.5)
utrecht.plot(ax=ax, edgecolor="red", facecolor="none", linewidth=1.5)
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)

# %%
# By default, ``burn_vector_geometry`` will only include grid faces whose
# centroid are located in a polygon. We can also mark all intersected faces
# by setting ``all_touched=True``:

burned = xu.burn_vector_geometry(utrecht, uda, all_touched=True)

fig, ax = plt.subplots()
burned.ugrid.plot(ax=ax)
burned.ugrid.plot.line(ax=ax, edgecolor="black", linewidth=0.5)
utrecht.plot(ax=ax, edgecolor="red", facecolor="none", linewidth=1.5)
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)

# %%
# We can also use such "masks" to e.g. modify specific parts of the grid data:

modified = (uda + 50.0).where(burned == 1, other=uda)
modified.ugrid.plot(vmin=-20, vmax=90, cmap="terrain")

# %%
# Note that ``all_touched=True`` is less suitable when differently valued
# polygons are present that share borders. While the centroid of a face is
# contained by only a single polygon, the area of the polygon may be located
# in more than one polygon. In this case, the results of each polygon will
# overwrite each other.

by_centroid = xu.burn_vector_geometry(provinces, uda, column="id")
by_touch = xu.burn_vector_geometry(provinces, uda, column="id", all_touched=True)

fig, axes = plt.subplots(ncols=2, figsize=(10, 5))
by_centroid.ugrid.plot(ax=axes[0], add_colorbar=False)
by_touch.ugrid.plot(ax=axes[1], add_colorbar=False)

for ax, title in zip(axes, ("centroid", "all touched")):
    burned.ugrid.plot.line(ax=ax, edgecolor="black", linewidth=0.5)
    utrecht.plot(ax=ax, edgecolor="red", facecolor="none", linewidth=1.5)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_title(title)

# %%
# This function can also be used to burn points or lines into the faces of an
# unstructured grid.
#
# The exterior boundaries of the province polygons will provide
# a collection of linestrings that we can burn into the grid:

lines = gpd.GeoDataFrame(geometry=provinces.exterior)
burned = xu.burn_vector_geometry(lines, uda)

fig, ax = plt.subplots()
burned.ugrid.plot(ax=ax)
burned.ugrid.plot.line(ax=ax, edgecolor="black", linewidth=0.5)
provinces.plot(ax=ax, edgecolor="red", facecolor="none", linewidth=1.5)
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)

# %%
# We can also burn points.

province_centroids = gpd.GeoDataFrame(geometry=provinces.centroid)
burned = xu.burn_vector_geometry(province_centroids, uda)

fig, ax = plt.subplots()
burned.ugrid.plot(ax=ax)
provinces.plot(ax=ax, edgecolor="red", facecolor="none")

# %%
# Finally, it's also possible to combine multiple geometry types in a single
# burn operation.

combined = pd.concat([lines, province_centroids])
burned = xu.burn_vector_geometry(combined, uda)
burned.ugrid.plot()

# %%
# Polygonizing
# ------------
#
# We can also do the opposite operation: turn collections of same-valued grid
# faces into vector polygons. Let's classify the elevation data into below and
# above the boundary of 5 m above mean sea level:

classified = uda > 5
polygonized = xu.polygonize(classified)
polygonized.plot(facecolor="none")

# %%
# We see that the results consists of two large polygons, in which the
# triangles of the triangular grid have been merged to form a single polygon,
# and many smaller polygons, some of which correspond one to one to the
# triangles of the grid.
#
# .. note::
#   The produced polygon edges will follow exactly the face boundaries. When
#   the data consists of many unique values (e.g. unbinned elevation data), the
#   result will essentially be one polygon per face. In such cases, it is more
#   efficient to use ``xugrid.UgridDataArray.to_geodataframe``, which directly
#   converts every face to a polygon.
#
# Snap to grid
# ------------
#
# The examples above deal with data on the faces of the grid. However, data can
# also be associated with the nodes or edges of the grid. For example, we can
# snap line data to exactly follow the edges (the face boundaries):

linestrings = gpd.GeoDataFrame(geometry=utrecht.exterior)
snapped_uds, snapped_gdf = xu.snap_to_grid(linestrings, uda, max_snap_distance=1.0)

fig, ax = plt.subplots()
snapped_gdf.plot(ax=ax)
snapped_uds.ugrid.grid.plot(ax=ax, edgecolor="gray", linewidth=0.5)
utrecht.plot(ax=ax, edgecolor="red", facecolor="none", linewidth=0.75)
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)

# %%
# There are (arguably) many ways in which a linestring could be snapped to
# edges. The function above uses the criterion the following criterion: if a
# part of the linestring is located between the centroid of a face and the
# centroid of an edge, it is snapped to that edge.
#
# This sometimes result in less (visually) pleasing results, such as the two
# triangles in the lower left corner which are surrounded by the snapped edges.
# In general, results are best when the line segments are relatively large
# compare to the grid faces.
