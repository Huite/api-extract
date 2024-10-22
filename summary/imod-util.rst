imod.util.empty_2d
==================
Create an empty 2D (x, y) DataArray.

``dx`` and ``dy`` may be provided as:

    * scalar: for equidistant spacing
    * array: for non-equidistant spacing

Note that xarray (and netCDF4) uses midpoint coordinates. ``xmin`` and
``xmax`` are used to generate the appropriate midpoints.

Parameters
----------
dx: float, 1d array of floats
    cell size along x
xmin: float
xmax: float
dy: float, 1d array of floats
    cell size along y
ymin: float
ymax: float

Returns
-------
empty: xr.DataArray
    Filled with NaN.

imod.util.empty_2d_transient
============================
Create an empty transient 2D (time, x, y) DataArray.

``dx`` and ``dy`` may be provided as:

    * scalar: for equidistant spacing
    * array: for non-equidistant spacing

Note that xarray (and netCDF4) uses midpoint coordinates. ``xmin`` and
``xmax`` are used to generate the appropriate midpoints.

Parameters
----------
dx: float, 1d array of floats
    cell size along x
xmin: float
xmax: float
dy: float, 1d array of floats
    cell size along y
ymin: float
ymax: float
time: Any
    One or more of: str, numpy datetime64, pandas Timestamp

Returns
-------
empty: xr.DataArray
    Filled with NaN.

imod.util.empty_3d
==================
Create an empty 2D (x, y) DataArray.

``dx`` and ``dy`` may be provided as:

    * scalar: for equidistant spacing
    * array: for non-equidistant spacing

Note that xarray (and netCDF4) uses midpoint coordinates. ``xmin`` and
``xmax`` are used to generate the appropriate midpoints.

Parameters
----------
dx: float, 1d array of floats
    cell size along x
xmin: float
xmax: float
dy: float, 1d array of floats
    cell size along y
ymin: float
ymax: float
layer: int, sequence of integers, 1d array of integers

Returns
-------
empty: xr.DataArray
    Filled with NaN.

imod.util.empty_3d_transient
============================
Create an empty transient 3D (time, layer, x, y) DataArray.

``dx`` and ``dy`` may be provided as:

    * scalar: for equidistant spacing
    * array: for non-equidistant spacing

Note that xarray (and netCDF4) uses midpoint coordinates. ``xmin`` and
``xmax`` are used to generate the appropriate midpoints.

Parameters
----------
dx: float, 1d array of floats
    cell size along x
xmin: float
xmax: float
dy: float, 1d array of floats
    cell size along y
ymin: float
ymax: float
layer: int, sequence of integers, 1d array of integers
time: Any
    One or more of: str, numpy datetime64, pandas Timestamp

Returns
-------
empty: xr.DataArray
    Filled with NaN.

imod.util.where
===============
Wrapped version of xarray's ``.where``.

This wrapped version does two differently:

Firstly, it prioritizes the dimensions as: ``if_true > if_false > condition``.
``xarray.where(cond, a, b)`` will choose the dimension over ``a`` or ``b``.
This may result in unwanted dimension orders such as ``("y", "x", "layer)``
rather than ``("layer", "y', "x")``.

Secondly, it preserves the NaN values of ``if_true`` by default.  If we
wish to replace all values over 5 by 5, yet keep the NoData parts, this
requires two operations with with xarray's ``where``.

Parameters
----------
condition: DataArray, Dataset
    Locations at which to preserve this object's values. dtype must be `bool`.
if_true : scalar, DataArray or Dataset, optional
    Value to use for locations where ``cond`` is True.
if_false : scalar, DataArray or Dataset, optional
    Value to use for locations where ``cond`` is False.
keep_nan: bool, default: True
    Whether to keep the NaN values in place of ``if_true``.

imod.util.cd
============
Change directory, and change it back after the with block.

Examples
--------
>>> with imod.util.context.cd("docs"):
        do_something_in_docs()

imod.util.ignore_warnings
=========================
Contextmanager to ignore RuntimeWarnings as they are frequently
raised by the Dask delayed scheduler.

Examples
--------
>>> with imod.util.context.ignore_warnings():
        function_that_throws_warnings()

imod.util.to_datetime
=====================
Convert string to datetime. Part of the public API for backwards
compatibility reasons.

Fast performance is important, as this function is used to parse IDF names,
so it being called 100,000 times is a common usecase. Function stored
previously under imod.util.to_datetime.

imod.util.coord_reference
=========================
Extracts dx, xmin, xmax for a coordinate DataArray, where x is any coordinate.

If the DataArray coordinates are nonequidistant, dx will be returned as
1D ndarray instead of float.

Parameters
----------
a : xarray.DataArray of a coordinate

Returns
--------------
tuple
    (dx, xmin, xmax) for a coordinate x

imod.util.spatial_reference
===========================
Extracts spatial reference from DataArray.

If the DataArray coordinates are nonequidistant, dx and dy will be returned
as 1D ndarray instead of float.

Parameters
----------
a : xarray.DataArray

Returns
--------------
tuple
    (dx, xmin, xmax, dy, ymin, ymax)

imod.util.transform
===================
Extract the spatial reference information from the DataArray coordinates,
into an affine.Affine object for writing to e.g. rasterio supported formats.

Parameters
----------
a : xarray.DataArray

Returns
-------
affine.Affine

imod.util.to_ugrid2d
====================
Convert a structured DataArray or Dataset into its UGRID-2D quadrilateral
equivalent.

See:
https://ugrid-conventions.github.io/ugrid-conventions/#2d-flexible-mesh-mixed-triangles-quadrilaterals-etc-topology

Parameters
----------
data: Union[xr.DataArray, xr.Dataset]
    Dataset or DataArray with last two dimensions ("y", "x").
    In case of a Dataset, the 2D topology is defined once and variables are
    added one by one.
    In case of a DataArray, a name is required; a name can be set with:
    ``da.name = "..."``'

Returns
-------
ugrid2d_dataset: xr.Dataset
    The equivalent data, in UGRID-2D quadrilateral form.

imod.util.mdal_compliant_ugrid2d
================================
Ensures the xarray Dataset will be written to a UGRID netCDF that will be
accepted by MDAL.

* Unstacks variables with a layer dimension into separate variables.
* Removes absent entries from the mesh topology attributes.
* Sets encoding to float for datetime variables.

Parameters
----------
dataset: xarray.Dataset
    Dataset to make compliant with MDAL
crs: Any, Optional
    Anything accepted by rasterio.crs.CRS.from_user_input
    Requires ``rioxarray`` installed.

Returns
-------
unstacked: xr.Dataset

imod.util.from_mdal_compliant_ugrid2d
=====================================
Undo some of the changes of ``mdal_compliant_ugrid2d``: re-stack the
layers.

Parameters
----------
dataset: xugrid.UgridDataset

Returns
-------
restacked: xugrid.UgridDataset

