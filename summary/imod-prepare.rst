imod.prepare.Regridder
======================
Object to repeatedly regrid similar objects. Compiles once on first call,
can then be repeatedly called without JIT compilation overhead.

Attributes
----------
method : str, function
    The method to use for regridding. Default available methods are:
    ``{"nearest", "multilinear", mean", "harmonic_mean", "geometric_mean",
    "sum", "minimum", "maximum", "mode", "median", "conductance"}``
ndim_regrid : int, optional
    The number of dimensions over which to regrid. If not provided,
    ``ndim_regrid`` will be inferred. It serves to prevent regridding over an
    unexpected number of dimensions; say you want to regrid over only two
    dimensions. Due to an input error in the coordinates of ``like``, three
    dimensions may be inferred in the first ``.regrid`` call. An error will
    be raised if ndim_regrid not match the number of inferred dimensions.
    Default value is None.
use_relative_weights : bool, optional
    Whether to use relative weights in the regridding method or not.
    Relative weights are defined as: cell_overlap / source_cellsize, for
    every axis.

    This argument should only be used if you are providing your own
    ``method`` as a function, where the function requires relative, rather
    than absolute weights (the provided ``conductance`` method requires
    relative weights, for example). Default value is False.
extra_overlap : integer, optional
    In case of chunked regridding, how many cells of additional overlap is
    necessary. Linear interpolation requires this for example, as it reaches
    beyond cell boundaries to compute values. Default value is 0.

Examples
--------
Initialize the Regridder object:

>>> mean_regridder = imod.prepare.Regridder(method="mean")

Then call the ``regrid`` method to regrid.

>>> result = mean_regridder.regrid(source, like)

The regridder can be re-used if the number of regridding dimensions
match, saving some time by not (re)compiling the regridding method.

>>> second_result = mean_regridder.regrid(second_source, like)

A one-liner is possible for single use:

>>> result = imod.prepare.Regridder(method="mean").regrid(source, like)

It's possible to provide your own methods to the ``Regridder``, provided that
numba can compile them. They need to take the arguments ``values`` and
``weights``. Make sure they deal with ``nan`` values gracefully!

>>> def p30(values, weights):
>>>     return np.nanpercentile(values, 30)

>>> p30_regridder = imod.prepare.Regridder(method=p30)
>>> p30_result = p30_regridder.regrid(source, like)

The Numba developers maintain a list of support Numpy features here:
https://numba.pydata.org/numba-doc/dev/reference/numpysupported.html

In general, however, the provided methods should be adequate for your
regridding needs.

imod.prepare.Regridder Class Members
====================================
   * imod.prepare.Regridder.regrid

imod.prepare.Regridder.regrid
=============================
Regrid ``source`` along dimensions that ``source`` and ``like`` share.
These dimensions will be inferred the first time ``.regrid`` is called
for the Regridder object.

Following xarray conventions, nodata is assumed to ``np.nan``.

Parameters
----------
source : xr.DataArray of floats
like : xr.DataArray of floats
    The like array present what the coordinates should look like.
fill_value : float
    The fill_value. Defaults to np.nan

Returns
-------
result : xr.DataArray
    Regridded result.

imod.prepare.LayerRegridder
===========================
Object to repeatedly regrid layers of similar objects. Compiles
once on first call, can then be repeatedly called without
JIT compilation overhead.

Attributes
----------
method : str, function
    The method to use for regridding. Default available methods are:
    ``{"mean", "harmonic_mean", "geometric_mean", "sum", "minimum",
    "maximum", "mode", "median", "max_overlap"}``

imod.prepare.LayerRegridder Class Members
=========================================
   * imod.prepare.LayerRegridder.regrid

imod.prepare.LayerRegridder.regrid
==================================
Parameters
----------
source : xr.DataArray
    The values of the layered model.
source_top : xr.DataArray
    The vertical location of the layer tops.
destination_top : xr.DataArray
    The vertical location of the layer tops.
source_bottom : xr.DataArray
    The vertical location of the layer bottoms.
destination_bottom : xr.DataArray
    The vertical location of the layer bottoms.

Returns
-------
regridded : xr.DataArray

imod.prepare.Voxelizer
======================
Object to repeatedly voxelize similar objects. Compiles once on first call,
can then be repeatedly called without JIT compilation overhead.

Attributes
----------
method : str, function
    The method to use for regridding. Default available methods are:
    ``{"mean", "harmonic_mean", "geometric_mean", "sum", "minimum",
    "maximum", "mode", "median", "max_overlap"}``

Examples
--------
Usage is similar to the regridding. Initialize the Voxelizer object:

>>> mean_voxelizer = imod.prepare.Voxelizer(method="mean")

Then call the ``voxelize`` method to transform a layered dataset into a
voxel based one. The vertical coordinates of the layers must be provided
by ``top`` and ``bottom``.

>>> mean_voxelizer.voxelize(source, top, bottom, like)

If your data is already voxel based, i.e. the layers have tops and bottoms
that do not differ with x or y, you should use a ``Regridder`` instead.

It's possible to provide your own methods to the ``Regridder``, provided that
numba can compile them. They need to take the arguments ``values`` and
``weights``. Make sure they deal with ``nan`` values gracefully!

>>> def p30(values, weights):
>>>     return np.nanpercentile(values, 30)

>>> p30_voxelizer = imod.prepare.Voxelizer(method=p30)
>>> p30_result = p30_voxelizer.regrid(source, top, bottom, like)

The Numba developers maintain a list of support Numpy features here:
https://numba.pydata.org/numba-doc/dev/reference/numpysupported.html

In general, however, the provided methods should be adequate for your
voxelizing needs.

imod.prepare.Voxelizer Class Members
====================================
   * imod.prepare.Voxelizer.voxelize

imod.prepare.Voxelizer.voxelize
===============================
Parameters
----------
source : xr.DataArray
    The values of the layered model.
top : xr.DataArray
    The vertical location of the layer tops.
bottom : xr.DataArray
    The vertical location of the layer bottoms.
like : xr.DataArray
    An example DataArray providing the coordinates of the voxelized
    results; what it should look like in terms of dimensions, data type,
    and coordinates.

Returns
-------
voxelized : xr.DataArray

imod.prepare.fill
=================
Replace the value of invalid ``da`` cells (indicated by ``invalid``)
using basic nearest neighbour interpolation.

Parameters
----------
da: xr.DataArray with gaps
    array containing missing value
by: str, optional
    dimension over which the array will be filled, one by one.
    See the examples.

invalid: xr.DataArray
    a binary array of same shape as ``da``.
    data value are replaced where invalid is True
    If None (default), uses: `invalid = np.isnan(data)`

Returns
-------
xarray.DataArray
    with the same coordinates as the input.

Examples
--------

A common use case is filling holes in a DataArray, filling it with the
value of its nearest (valid) neighbor:

>>> filled = imod.prepare.fill(da)

In case of a tie (e.g. neighbors in x and y are both one cell removed), the
neighbor in the last dimension is chosen (for rasters, that's generally x).

A typical use case is filling a 3D array (layer, y, x), but only in the
horizontal dimensions. The ``by`` keyword can be used to do this:

>>> filled = imod.prepare.fill(da, by="layer")

In this case, the array is filled by one layer at a time.

imod.prepare.laplace_interpolate
================================
Fills gaps in `source` by interpolating from existing values using Laplace
interpolation.

Parameters
----------
source : xr.DataArray of floats with dims (y, x)
    Data values to interpolate.
ibound : xr.DataArray of bool with dims (y, x)
    Precomputed array which marks where to interpolate.
close : float
    Closure criteration of iterative solver. Should be one to two orders
    of magnitude smaller than desired accuracy.
mxiter : int
    Outer iterations of iterative solver.
iter1 : int
    Inner iterations of iterative solver. Should not exceed 50.
relax : float
    Iterative solver relaxation parameter. Should be between 0 and 1.

Returns
-------
interpolated : xr.DataArray with dims (y, x)
    source, with interpolated values where ibound equals 1

imod.prepare.polygonize
=======================
Polygonize a 2D-DataArray into a GeoDataFrame of polygons.

Parameters
----------
da : xr.DataArray

Returns
-------
polygonized : geopandas.GeoDataFrame

imod.prepare.reproject
======================
Reprojects and/or resamples a 2D xarray DataArray to a
different cellsize or coordinate system.

* To resample to a new cellsize in the same projection: provide only ``like``.
* To only reproject: provide only ``src_crs`` and ``src_crs``.
* To reproject and resample to a specific domain: provide ``src_crs``, ``src_crs``, and ``like``.

Note: when only ``like`` is provided, Cartesian (projected) coordinates are a
ssumed for resampling. In case of non-Cartesian coordinates, specify
``src_crs`` and ``dst_crs`` for correct resampling.

Parameters
----------
source: xarray DataArray
    The DataArray to be resampled and/or reprojected. Must contain dimensions
    ``y`` and ``x``.
like: xarray DataArray
    Example DataArray that shows what the resampled result should look like
    in terms of coordinates. Must contain dimensions ``y`` and ``x``.
src_crs: string, dict, rasterio.crs.CRS
    Coordinate system of ``source``. Options:

    * string: e.g. ``"EPSG:4326"``
    * rasterio.crs.CRS
dst_crs: string, dict, rasterio.crs.CRS
    Coordinate system of result. Options:

    * string: e.g. ``"EPSG:4326"``
    * rasterio.crs.CRS
use_src_attrs: boolean
    If True: Use metadata in ``source.attrs``, as generated by ``xarray.open_rasterio()``, to do
    reprojection.
method: string
    The method to use for resampling/reprojection.
    Defaults to "nearest". GDAL methods are available:

    * nearest
    * bilinear
    * cubic
    * cubic_spline
    * lanczos
    * average
    * mode
    * gauss
    * max
    * min
    * med (50th percentile)
    * q1 (25th percentile)
    * q3 (75th percentile)
reproject_kwargs: dict, optional
    keyword arguments for ``rasterio.warp.reproject()``.

Returns
-------
xarray.DataArray
    Resampled/reprojected DataArray.

Examples
--------
Resample a DataArray ``a`` to a new cellsize, using an existing DataArray ``b``:

>>> c = imod.rasterio.reproject(source=a, like=b)

Resample a DataArray to a new cellsize of 100.0, by creating a ``like`` DataArray first:
(Note that dy must be negative, as is usual for geospatial grids.)

>>> dims = ("y", "x")
>>> coords = {"y": np.arange(200_000.0, 100_000.0, -100.0), "x": np.arange(0.0, 100_000.0, 100.0)}
>>> b = xr.DataArray(data=np.empty((200, 100)), coords=coords, dims=dims)
>>> c = imod.rasterio.reproject(source=a, like=b)

Reproject a DataArray from one coordinate system (WGS84, EPSG:4326) to another (UTM30N, EPSG:32630):

>>> c = imod.rasterio.reproject(source=a, src_crs="EPSG:4326", dst_crs="EPSG:32630")

Get the reprojected DataArray in the desired shape and coordinates by providing ``like``:

>>> c = imod.rasterio.reproject(source=a, like=b, src_crs="EPSG:4326", dst_crs="EPSG:32630")

Open a single band raster, and reproject to RD new coordinate system (EPSG:28992), without explicitly specifying ``src_crs``.
``src_crs`` is taken from ``a.attrs``, so the raster file has to include coordinate system metadata for this to work.

>>> a = rioxarray.open_rasterio("example.tif").squeeze("band")
>>> c = imod.rasterio.reproject(source=a, use_src_attrs=True, dst_crs="EPSG:28992")

In case of a rotated ``source``, provide ``src_transform`` directly or ``use_src_attrs=True`` to rely on generated attributes:

>>> rotated = rioxarray.open_rasterio("rotated_example.tif").squeeze("band")
>>> c = imod.rasterio.reproject(source=rotated, dst_crs="EPSG:28992", reproject_kwargs={"src_transform":affine.Affine(...)})
>>> c = imod.rasterio.reproject(source=rotated, dst_crs="EPSG:28992", use_src_attrs=True)

imod.prepare.rasterize
======================
Rasterize a geopandas GeoDataFrame onto the given
xarray coordinates.

Parameters
----------
geodataframe : geopandas.GeoDataFrame
column : str, int, float
    column name of geodataframe to burn into raster
like : xarray.DataArray
    Example DataArray. The rasterized result will match the shape and
    coordinates of this DataArray.
fill : float, int
    Fill value for nodata areas. Optional, default value is np.nan.
kwargs : additional keyword arguments for rasterio.features.rasterize.
    See: https://rasterio.readthedocs.io/en/stable/api/rasterio.features.html#rasterio.features.rasterize

Returns
-------
rasterized : xarray.DataArray
    Vector data rasterized. Matches shape and coordinates of ``like``.

imod.prepare.gdal_rasterize
===========================
Use GDAL to rasterize a vector file into an xarray.DataArray.

Can be significantly more efficient than rasterize. This doesn't load the
vector data into a GeoDataFrame and loops over the individual shapely
geometries like rasterio.rasterize does, but loops over the features within
GDAL instead.

Parameters
----------
path : str or pathlib.Path
    path to OGR supported vector file (e.g. a shapefile)
column : str
    column name of column to burn into raster
like : xr.DataArray, optional
    example of raster
nodata : int, float; optional
dtype : numpy.dtype, optional
spatial_reference : dict, optional
    Optional dict to avoid allocating the like DataArray. Used if template
    is None. Dict has keys "bounds" and "cellsizes", with:

    * bounds = (xmin, xmax, ymin, ymax)
    * cellsizes = (dx, dy)
all_touched : bool
    If True: all pixels touched by lines or polygons will be updated, not
    just those on the line render path, or whose center point is within the
    polygon. Default value is False.

Returns
-------
rasterized : xr.DataArray

imod.prepare.celltable
======================
Process area of features by rasterizing in a chunkwise manner to limit
memory usage.

Returns a table of cell indices (row, column) with for example feature ID,
and feature area within cell. Essentially returns a COO sparse matrix, but
with duplicate values per cell, since more than one geometry may be present.

The feature area within the cell is approximated by first rasterizing the
feature, and then counting the number of occuring cells. This means the
accuracy of the area depends on the cellsize of the rasterization step.

A celltable is returned, as a ``pandas.DataFrame``. It has the following
columns:

1. ``"row_index"``
2. ``"col_index"``
3. the value of the ``column`` argument
4. ``"area"``

``"row_index"`` and ``"col_index"`` are the indices of the like array in
which the polygon is located. The ``column`` value holds the rasterized
value of the specified column. ``"area"`` contains the area of the
polygon within the cell.

The most convenient way of using this celltable is by specifying a feature
ID as ``column``. After creating a celltable, ``pandas.DataFrame.merge()``
can be used to join additional data on this ID. Refer to the examples.

Parameters
----------
path : str or pathlib.Path
    path to OGR supported vector file (e.g. a shapefile)
column : str
    column name of column to burn into raster
resolution : float
    cellsize at which the rasterization, and determination of area within
    cellsize occurs. Very small values are recommended (e.g. <= 0.5 m).
like : xarray.DataArray
    Example DataArray of where the cells will be located. Used only for the
    coordinates.
dtype: DtypeLike, optional
    datatype of data referred to with "column", defaults to 32-bit integer.
chunksize : int, optional
    The size of the chunksize. Used for both x and y dimension.

Returns
-------
celltable : pandas.DataFrame

Examples
--------
Assume we have a shapefile called ``waterways.shp`` and information on the
model discretization is described by a ``like`` DataArray. The feature ID is
provided by a column in the shapefile called "ID-code". Additionally, this
shapefile also specifies bed hydraulic resistance (c0). For this specific
discretization, we wish to calculate a conductance (area divided by
hydraulic resistance). To do so, we:

1. create a ``celltable``
2. join the additional attributes (such as c0)
3. compute the conductance per feature
4. sum conductances per cell

Import the required packages.

>>> import imod
>>> import geopandas as gpd

Generate the celltable.

>>> celltable = imod.prepare.celltable(
        path="waterways.shp",
        column="ID-code",
        resolution=0.5,
        like=like,
    )

Load the shapefile with geopandas into a ``GeoDataFrame``.

>>> gdf = gpd.read_file("waterways.shp)

Select the relevant columns into a ``pandas.DataFrame`` and merge with the
celltable.

>>> df = gdf[["ID-code", "c0"]]
>>> joined = celltable.merge(gdf, on="ID-code")

We compute the conductance, and sum it per cell using ``pandas`` methods:

>>> joined["conductance"] = joined["area"] / joined["c0"]
>>> summed_conductance = joined.groupby(["row_index", "col_index"], as_index=False)[
        "conductance"
    ].sum()

Finally, turn the result into a DataArray so it can be used as model input:

>>> conductance = imod.prepare.rasterize_celltable(
        table=summed_conductance,
        column="conductance",
        like=like,
    )

imod.prepare.rasterize_celltable
================================
Rasterizes a table, such as produced by ``imod.prepare.spatial.celltable``.
Before rasterization, multiple values should be grouped and aggregated per
cell. Values will be overwritten otherwise.

Parameters
----------
like : xr.DataArray
table : pandas.DataFrame
    with columns: "row_index", "col_index"
column : str, int, float
    column name of values to rasterize

Returns
-------
rasterized : xr.DataArray

imod.prepare.zonal_aggregate_polygons
=====================================
Compute a zonal aggregate of polygon data for (other) polygon geometries,
e.g. a mean, mode, or percentile.

Parameters
----------
path_a : str or pathlib.Path
    path to OGR supported vector file (e.g. a shapefile)
path_b : str or pathlib.Path
    path to OGR supported vector file (e.g. a shapefile)
column_a : str
    column name of path_a. Defines zones of aggregation.
column_b : str
    column name of path_b. Data to aggregate.
like : xarray.DataArray with dims ("y", "x")
    Example DataArray of where the cells will be located. Used only for its
    x and y coordinates.
resolution : float
    cellsize at which the rasterization, sampling, and area measurement occur.
method: Union[str, Callable]
    Aggregation method.
    Anything that is acceptable by a pandas groupby aggregate:
    https://pandas.pydata.org/docs/reference/api/pandas.core.groupby.DataFrameGroupBy.aggregate.html
chunksize : int, optional
    The size of the chunksize. Used for both x and y dimension.

Returns
-------
zonal_aggregates: pandas.DataFrame

imod.prepare.zonal_aggregate_raster
===================================
Compute a zonal aggregate of raster data for polygon geometries, e.g. a mean,
mode, or percentile.

Parameters
----------
path : str or pathlib.Path
    path to OGR supported vector file (e.g. a shapefile). Defines zones
    of aggregation.
column : str
    column name of path, integer IDs defining zones.
raster : xarray.DataArray
    Raster data from which to sample and aggregate data
resolution : float
    cellsize at which the rasterization of polygons and sampling occurs
method : Union[str, Callable]
    Aggregation method.
    Anything that is acceptable by a pandas groupby aggregate:
    https://pandas.pydata.org/docs/reference/api/pandas.core.groupby.DataFrameGroupBy.aggregate.html
chunksize : int, optional
    The size of the chunksize. Used for both x and y dimension.

Returns
-------
zonal_aggregates : pandas.DataFrame

Examples
--------

To compute the mean surface level at polygons of water bodies:

>>> import imod
>>> surface_level = imod.rasterio.open("surface_level.tif")
>>> df = imod.prepare.spatial.zonal_aggregate_raster(
>>>    path="water-bodies.shp",
>>>    column="id",
>>>    raster=surface_level,
>>>    resolution=1.0,
>>>    method="mean",
>>> )

For some functions, like the mode, a function should be passed instead:

>>> import pandas as pd
>>> df = imod.prepare.spatial.zonal_aggregate_raster(
>>>    path="water-bodies.shp",
>>>    column="id",
>>>    raster=surface_level,
>>>    resolution=1.0,
>>>    method=pd.Series.mode,
>>> )

imod.prepare.assign_wells
=========================
Distribute well pumping rate according to filter length when ``k=None``, or
to transmissivity of the sediments surrounding the filter. Minimum
thickness and minimum k should be set to avoid placing wells in clay
layers.

Wells located outside of the grid are removed.

Parameters
----------
wells: pd.DataFrame
    Should contain columns x, y, id, top, bottom, rate.
top: xr.DataArray or xu.UgridDataArray
    Top of the model layers.
bottom: xr.DataArray or xu.UgridDataArray
    Bottom of the model layers.
k: xr.DataArray or xu.UgridDataArray, optional
    Horizontal conductivity of the model layers.
minimum_thickness: float, optional, default: 0.01
minimum_k: float, optional, default: 1.0
    Minimum conductivity
validate: bool
    raise an excpetion if one of the wells is not in the domain
Returns
-------
placed_wells: pd.DataFrame
    Wells with rate subdivided per layer. Contains the original columns of
    ``wells``, as well as layer, overlap, transmissivity.

imod.prepare.get_lower_active_grid_cells
========================================
Returns grid of booleans designating location of the lowermost active grid
cell.

Parameters
----------
active: {DataArray, UgridDataArray}
    Grid of booleans designating active cell.

imod.prepare.get_lower_active_layer_number
==========================================
Returns two-dimensional grid of integers with the layer number of the lower
most active cell.

Parameters
----------
active: {DataArray, UgridDataArray}
    Grid of booleans designating active cell.

imod.prepare.get_upper_active_grid_cells
========================================
Returns grid of booleans designating location of the uppermost active grid
cell.

Parameters
----------
active: {DataArray, UgridDataArray}
    Grid of booleans designating active cell.

imod.prepare.get_upper_active_layer_number
==========================================
Returns planar grid of integers with the layer number of the upper most
active cell.

Parameters
----------
active: {DataArray, UgridDataArray}
    Grid of booleans designating active cell.

imod.prepare.create_layered_top
===============================
Create a top array with a layer dimension, from a top array with no layer
dimension and a bottom array with a layer dimension. The (output) top of
layer n is assigned the bottom of layer n-1.

Parameters
----------
bottom: {DataArray, UgridDataArray}
    Bottoms with layer dimension
top: {DataArray, UgridDataArray}
    Top, without layer dimension

Returns
-------
new_top: {DataArray, UgridDataArray}
    Top with layer dimension.

imod.prepare.ALLOCATION_OPTION
==============================
Enumerator for settings to allocate planar grid with RIV, DRN, GHB, or RCH
cells over the vertical layer dimensions. Numbers match the IDEFLAYER
options in iMOD 5.6.

* ``stage_to_riv_bot``: RIV. Allocate cells spanning from the river stage up
  to the river bottom elevation. This matches the iMOD 5.6 IDEFLAYER = 0
  option.
* ``first_active_to_elevation``: RIV, DRN, GHB. Allocate cells spanning from
  first upper active cell up to the river bottom elevation. This matches the
  iMOD 5.6 IDEFLAYER = -1 option.
* ``stage_to_riv_bot_drn_above``: RIV. Allocate cells spanning from first
  upper active cell up to the river bottom elevation. Method returns both
  allocated cells for a river package as well as a drain package. Cells
  above river stage are allocated as drain cells, cells below are as river
  cells. This matches the iMOD 5.6 IDEFLAYER = 1 option.
* ``at_elevation``: RIV, DRN, GHB. Allocate cells containing the river
  bottom elevation, drain elevation, or head respectively for river, drain
  and general head boundary. This matches the iMOD 5.6 IDEFLAYER = 2
  option.
* ``at_first_active``: RIV, DRN, GHB, RCH. Allocate cells at the upper
  active cells. This has no equivalent option in iMOD 5.6.

Examples
--------

>>> from imod.prepare.topsystem import ALLOCATION_OPTION
>>> setting = ALLOCATION_OPTION.at_first_active
>>> allocated = allocate_rch_cells(setting, active, rate)

imod.prepare.ALLOCATION_OPTION Class Members
============================================
   * imod.prepare.ALLOCATION_OPTION.at_elevation
   * imod.prepare.ALLOCATION_OPTION.at_first_active
   * imod.prepare.ALLOCATION_OPTION.first_active_to_elevation
   * imod.prepare.ALLOCATION_OPTION.name
   * imod.prepare.ALLOCATION_OPTION.stage_to_riv_bot
   * imod.prepare.ALLOCATION_OPTION.stage_to_riv_bot_drn_above
   * imod.prepare.ALLOCATION_OPTION.value

imod.prepare.ALLOCATION_OPTION.at_elevation
===========================================
Enumerator for settings to allocate planar grid with RIV, DRN, GHB, or RCH
cells over the vertical layer dimensions. Numbers match the IDEFLAYER
options in iMOD 5.6.

* ``stage_to_riv_bot``: RIV. Allocate cells spanning from the river stage up
  to the river bottom elevation. This matches the iMOD 5.6 IDEFLAYER = 0
  option.
* ``first_active_to_elevation``: RIV, DRN, GHB. Allocate cells spanning from
  first upper active cell up to the river bottom elevation. This matches the
  iMOD 5.6 IDEFLAYER = -1 option.
* ``stage_to_riv_bot_drn_above``: RIV. Allocate cells spanning from first
  upper active cell up to the river bottom elevation. Method returns both
  allocated cells for a river package as well as a drain package. Cells
  above river stage are allocated as drain cells, cells below are as river
  cells. This matches the iMOD 5.6 IDEFLAYER = 1 option.
* ``at_elevation``: RIV, DRN, GHB. Allocate cells containing the river
  bottom elevation, drain elevation, or head respectively for river, drain
  and general head boundary. This matches the iMOD 5.6 IDEFLAYER = 2
  option.
* ``at_first_active``: RIV, DRN, GHB, RCH. Allocate cells at the upper
  active cells. This has no equivalent option in iMOD 5.6.

Examples
--------

>>> from imod.prepare.topsystem import ALLOCATION_OPTION
>>> setting = ALLOCATION_OPTION.at_first_active
>>> allocated = allocate_rch_cells(setting, active, rate)

imod.prepare.ALLOCATION_OPTION.at_first_active
==============================================
Enumerator for settings to allocate planar grid with RIV, DRN, GHB, or RCH
cells over the vertical layer dimensions. Numbers match the IDEFLAYER
options in iMOD 5.6.

* ``stage_to_riv_bot``: RIV. Allocate cells spanning from the river stage up
  to the river bottom elevation. This matches the iMOD 5.6 IDEFLAYER = 0
  option.
* ``first_active_to_elevation``: RIV, DRN, GHB. Allocate cells spanning from
  first upper active cell up to the river bottom elevation. This matches the
  iMOD 5.6 IDEFLAYER = -1 option.
* ``stage_to_riv_bot_drn_above``: RIV. Allocate cells spanning from first
  upper active cell up to the river bottom elevation. Method returns both
  allocated cells for a river package as well as a drain package. Cells
  above river stage are allocated as drain cells, cells below are as river
  cells. This matches the iMOD 5.6 IDEFLAYER = 1 option.
* ``at_elevation``: RIV, DRN, GHB. Allocate cells containing the river
  bottom elevation, drain elevation, or head respectively for river, drain
  and general head boundary. This matches the iMOD 5.6 IDEFLAYER = 2
  option.
* ``at_first_active``: RIV, DRN, GHB, RCH. Allocate cells at the upper
  active cells. This has no equivalent option in iMOD 5.6.

Examples
--------

>>> from imod.prepare.topsystem import ALLOCATION_OPTION
>>> setting = ALLOCATION_OPTION.at_first_active
>>> allocated = allocate_rch_cells(setting, active, rate)

imod.prepare.ALLOCATION_OPTION.first_active_to_elevation
========================================================
Enumerator for settings to allocate planar grid with RIV, DRN, GHB, or RCH
cells over the vertical layer dimensions. Numbers match the IDEFLAYER
options in iMOD 5.6.

* ``stage_to_riv_bot``: RIV. Allocate cells spanning from the river stage up
  to the river bottom elevation. This matches the iMOD 5.6 IDEFLAYER = 0
  option.
* ``first_active_to_elevation``: RIV, DRN, GHB. Allocate cells spanning from
  first upper active cell up to the river bottom elevation. This matches the
  iMOD 5.6 IDEFLAYER = -1 option.
* ``stage_to_riv_bot_drn_above``: RIV. Allocate cells spanning from first
  upper active cell up to the river bottom elevation. Method returns both
  allocated cells for a river package as well as a drain package. Cells
  above river stage are allocated as drain cells, cells below are as river
  cells. This matches the iMOD 5.6 IDEFLAYER = 1 option.
* ``at_elevation``: RIV, DRN, GHB. Allocate cells containing the river
  bottom elevation, drain elevation, or head respectively for river, drain
  and general head boundary. This matches the iMOD 5.6 IDEFLAYER = 2
  option.
* ``at_first_active``: RIV, DRN, GHB, RCH. Allocate cells at the upper
  active cells. This has no equivalent option in iMOD 5.6.

Examples
--------

>>> from imod.prepare.topsystem import ALLOCATION_OPTION
>>> setting = ALLOCATION_OPTION.at_first_active
>>> allocated = allocate_rch_cells(setting, active, rate)

imod.prepare.ALLOCATION_OPTION.name
===================================
The name of the Enum member.

imod.prepare.ALLOCATION_OPTION.stage_to_riv_bot
===============================================
Enumerator for settings to allocate planar grid with RIV, DRN, GHB, or RCH
cells over the vertical layer dimensions. Numbers match the IDEFLAYER
options in iMOD 5.6.

* ``stage_to_riv_bot``: RIV. Allocate cells spanning from the river stage up
  to the river bottom elevation. This matches the iMOD 5.6 IDEFLAYER = 0
  option.
* ``first_active_to_elevation``: RIV, DRN, GHB. Allocate cells spanning from
  first upper active cell up to the river bottom elevation. This matches the
  iMOD 5.6 IDEFLAYER = -1 option.
* ``stage_to_riv_bot_drn_above``: RIV. Allocate cells spanning from first
  upper active cell up to the river bottom elevation. Method returns both
  allocated cells for a river package as well as a drain package. Cells
  above river stage are allocated as drain cells, cells below are as river
  cells. This matches the iMOD 5.6 IDEFLAYER = 1 option.
* ``at_elevation``: RIV, DRN, GHB. Allocate cells containing the river
  bottom elevation, drain elevation, or head respectively for river, drain
  and general head boundary. This matches the iMOD 5.6 IDEFLAYER = 2
  option.
* ``at_first_active``: RIV, DRN, GHB, RCH. Allocate cells at the upper
  active cells. This has no equivalent option in iMOD 5.6.

Examples
--------

>>> from imod.prepare.topsystem import ALLOCATION_OPTION
>>> setting = ALLOCATION_OPTION.at_first_active
>>> allocated = allocate_rch_cells(setting, active, rate)

imod.prepare.ALLOCATION_OPTION.stage_to_riv_bot_drn_above
=========================================================
Enumerator for settings to allocate planar grid with RIV, DRN, GHB, or RCH
cells over the vertical layer dimensions. Numbers match the IDEFLAYER
options in iMOD 5.6.

* ``stage_to_riv_bot``: RIV. Allocate cells spanning from the river stage up
  to the river bottom elevation. This matches the iMOD 5.6 IDEFLAYER = 0
  option.
* ``first_active_to_elevation``: RIV, DRN, GHB. Allocate cells spanning from
  first upper active cell up to the river bottom elevation. This matches the
  iMOD 5.6 IDEFLAYER = -1 option.
* ``stage_to_riv_bot_drn_above``: RIV. Allocate cells spanning from first
  upper active cell up to the river bottom elevation. Method returns both
  allocated cells for a river package as well as a drain package. Cells
  above river stage are allocated as drain cells, cells below are as river
  cells. This matches the iMOD 5.6 IDEFLAYER = 1 option.
* ``at_elevation``: RIV, DRN, GHB. Allocate cells containing the river
  bottom elevation, drain elevation, or head respectively for river, drain
  and general head boundary. This matches the iMOD 5.6 IDEFLAYER = 2
  option.
* ``at_first_active``: RIV, DRN, GHB, RCH. Allocate cells at the upper
  active cells. This has no equivalent option in iMOD 5.6.

Examples
--------

>>> from imod.prepare.topsystem import ALLOCATION_OPTION
>>> setting = ALLOCATION_OPTION.at_first_active
>>> allocated = allocate_rch_cells(setting, active, rate)

imod.prepare.ALLOCATION_OPTION.value
====================================
The value of the Enum member.

imod.prepare.DISTRIBUTING_OPTION
================================
Enumerator containing settings to distribute 2D conductance grids over
vertical layers for the RIV, DRN or GHB package.

* ``by_corrected_transmissivity``: RIV. Distribute the conductance by
  corrected transmissivities. Crosscut thicknesses are used to compute
  transmissivities. The crosscut thicknesses is computed based on the
  overlap of bottom_elevation over the bottom allocated layer. Same holds
  for the stage and top allocated layer. Furthermore the method corrects
  distribution weights for the mismatch between the midpoints of crosscut
  areas and model layer midpoints. This is the default method in iMOD 5.6,
  thus DISTRCOND = 0.
* ``equally``: RIV, DRN, GHB. Distribute conductances equally over layers.
  This matches iMOD 5.6 DISTRCOND = 1 option.
* ``by_crosscut_thickness``: RIV. Distribute the conductance by crosscut
  thicknesses. The crosscut thicknesses is computed based on the overlap of
  bottom_elevation over the bottom allocated layer. Same holds for the stage
  and top allocated layer. This matches iMOD 5.6 DISTRCOND = 2 option.
* ``by_layer_thickness``: RIV, DRN, GHB. Distribute the conductance by model
  layer thickness. This matches iMOD 5.6 DISTRCOND = 3 option.
* ``by_crosscut_transmissivity``: RIV. Distribute the conductance by
  crosscut transmissivity. Crosscut thicknesses are used to compute
  transmissivities. The crosscut thicknesses is computed based on the
  overlap of bottom_elevation over the bottom allocated layer. Same holds
  for the stage and top allocated layer. This matches iMOD 5.6 DISTRCOND = 4
  option.
* ``by_conductivity``: RIV, DRN, GHB. Distribute the conductance weighted by
  model layer hydraulic conductivities. This matches iMOD 5.6 DISTRCOND = 5
  option.
* ``by_layer_transmissivity``: RIV, DRN, GHB. Distribute the conductance by
  model layer transmissivity. This has no equivalent in iMOD 5.6.

imod.prepare.DISTRIBUTING_OPTION Class Members
==============================================
   * imod.prepare.DISTRIBUTING_OPTION.by_conductivity
   * imod.prepare.DISTRIBUTING_OPTION.by_corrected_transmissivity
   * imod.prepare.DISTRIBUTING_OPTION.by_crosscut_thickness
   * imod.prepare.DISTRIBUTING_OPTION.by_crosscut_transmissivity
   * imod.prepare.DISTRIBUTING_OPTION.by_layer_thickness
   * imod.prepare.DISTRIBUTING_OPTION.by_layer_transmissivity
   * imod.prepare.DISTRIBUTING_OPTION.equally
   * imod.prepare.DISTRIBUTING_OPTION.name
   * imod.prepare.DISTRIBUTING_OPTION.value

imod.prepare.DISTRIBUTING_OPTION.by_conductivity
================================================
Enumerator containing settings to distribute 2D conductance grids over
vertical layers for the RIV, DRN or GHB package.

* ``by_corrected_transmissivity``: RIV. Distribute the conductance by
  corrected transmissivities. Crosscut thicknesses are used to compute
  transmissivities. The crosscut thicknesses is computed based on the
  overlap of bottom_elevation over the bottom allocated layer. Same holds
  for the stage and top allocated layer. Furthermore the method corrects
  distribution weights for the mismatch between the midpoints of crosscut
  areas and model layer midpoints. This is the default method in iMOD 5.6,
  thus DISTRCOND = 0.
* ``equally``: RIV, DRN, GHB. Distribute conductances equally over layers.
  This matches iMOD 5.6 DISTRCOND = 1 option.
* ``by_crosscut_thickness``: RIV. Distribute the conductance by crosscut
  thicknesses. The crosscut thicknesses is computed based on the overlap of
  bottom_elevation over the bottom allocated layer. Same holds for the stage
  and top allocated layer. This matches iMOD 5.6 DISTRCOND = 2 option.
* ``by_layer_thickness``: RIV, DRN, GHB. Distribute the conductance by model
  layer thickness. This matches iMOD 5.6 DISTRCOND = 3 option.
* ``by_crosscut_transmissivity``: RIV. Distribute the conductance by
  crosscut transmissivity. Crosscut thicknesses are used to compute
  transmissivities. The crosscut thicknesses is computed based on the
  overlap of bottom_elevation over the bottom allocated layer. Same holds
  for the stage and top allocated layer. This matches iMOD 5.6 DISTRCOND = 4
  option.
* ``by_conductivity``: RIV, DRN, GHB. Distribute the conductance weighted by
  model layer hydraulic conductivities. This matches iMOD 5.6 DISTRCOND = 5
  option.
* ``by_layer_transmissivity``: RIV, DRN, GHB. Distribute the conductance by
  model layer transmissivity. This has no equivalent in iMOD 5.6.

imod.prepare.DISTRIBUTING_OPTION.by_corrected_transmissivity
============================================================
Enumerator containing settings to distribute 2D conductance grids over
vertical layers for the RIV, DRN or GHB package.

* ``by_corrected_transmissivity``: RIV. Distribute the conductance by
  corrected transmissivities. Crosscut thicknesses are used to compute
  transmissivities. The crosscut thicknesses is computed based on the
  overlap of bottom_elevation over the bottom allocated layer. Same holds
  for the stage and top allocated layer. Furthermore the method corrects
  distribution weights for the mismatch between the midpoints of crosscut
  areas and model layer midpoints. This is the default method in iMOD 5.6,
  thus DISTRCOND = 0.
* ``equally``: RIV, DRN, GHB. Distribute conductances equally over layers.
  This matches iMOD 5.6 DISTRCOND = 1 option.
* ``by_crosscut_thickness``: RIV. Distribute the conductance by crosscut
  thicknesses. The crosscut thicknesses is computed based on the overlap of
  bottom_elevation over the bottom allocated layer. Same holds for the stage
  and top allocated layer. This matches iMOD 5.6 DISTRCOND = 2 option.
* ``by_layer_thickness``: RIV, DRN, GHB. Distribute the conductance by model
  layer thickness. This matches iMOD 5.6 DISTRCOND = 3 option.
* ``by_crosscut_transmissivity``: RIV. Distribute the conductance by
  crosscut transmissivity. Crosscut thicknesses are used to compute
  transmissivities. The crosscut thicknesses is computed based on the
  overlap of bottom_elevation over the bottom allocated layer. Same holds
  for the stage and top allocated layer. This matches iMOD 5.6 DISTRCOND = 4
  option.
* ``by_conductivity``: RIV, DRN, GHB. Distribute the conductance weighted by
  model layer hydraulic conductivities. This matches iMOD 5.6 DISTRCOND = 5
  option.
* ``by_layer_transmissivity``: RIV, DRN, GHB. Distribute the conductance by
  model layer transmissivity. This has no equivalent in iMOD 5.6.

imod.prepare.DISTRIBUTING_OPTION.by_crosscut_thickness
======================================================
Enumerator containing settings to distribute 2D conductance grids over
vertical layers for the RIV, DRN or GHB package.

* ``by_corrected_transmissivity``: RIV. Distribute the conductance by
  corrected transmissivities. Crosscut thicknesses are used to compute
  transmissivities. The crosscut thicknesses is computed based on the
  overlap of bottom_elevation over the bottom allocated layer. Same holds
  for the stage and top allocated layer. Furthermore the method corrects
  distribution weights for the mismatch between the midpoints of crosscut
  areas and model layer midpoints. This is the default method in iMOD 5.6,
  thus DISTRCOND = 0.
* ``equally``: RIV, DRN, GHB. Distribute conductances equally over layers.
  This matches iMOD 5.6 DISTRCOND = 1 option.
* ``by_crosscut_thickness``: RIV. Distribute the conductance by crosscut
  thicknesses. The crosscut thicknesses is computed based on the overlap of
  bottom_elevation over the bottom allocated layer. Same holds for the stage
  and top allocated layer. This matches iMOD 5.6 DISTRCOND = 2 option.
* ``by_layer_thickness``: RIV, DRN, GHB. Distribute the conductance by model
  layer thickness. This matches iMOD 5.6 DISTRCOND = 3 option.
* ``by_crosscut_transmissivity``: RIV. Distribute the conductance by
  crosscut transmissivity. Crosscut thicknesses are used to compute
  transmissivities. The crosscut thicknesses is computed based on the
  overlap of bottom_elevation over the bottom allocated layer. Same holds
  for the stage and top allocated layer. This matches iMOD 5.6 DISTRCOND = 4
  option.
* ``by_conductivity``: RIV, DRN, GHB. Distribute the conductance weighted by
  model layer hydraulic conductivities. This matches iMOD 5.6 DISTRCOND = 5
  option.
* ``by_layer_transmissivity``: RIV, DRN, GHB. Distribute the conductance by
  model layer transmissivity. This has no equivalent in iMOD 5.6.

imod.prepare.DISTRIBUTING_OPTION.by_crosscut_transmissivity
===========================================================
Enumerator containing settings to distribute 2D conductance grids over
vertical layers for the RIV, DRN or GHB package.

* ``by_corrected_transmissivity``: RIV. Distribute the conductance by
  corrected transmissivities. Crosscut thicknesses are used to compute
  transmissivities. The crosscut thicknesses is computed based on the
  overlap of bottom_elevation over the bottom allocated layer. Same holds
  for the stage and top allocated layer. Furthermore the method corrects
  distribution weights for the mismatch between the midpoints of crosscut
  areas and model layer midpoints. This is the default method in iMOD 5.6,
  thus DISTRCOND = 0.
* ``equally``: RIV, DRN, GHB. Distribute conductances equally over layers.
  This matches iMOD 5.6 DISTRCOND = 1 option.
* ``by_crosscut_thickness``: RIV. Distribute the conductance by crosscut
  thicknesses. The crosscut thicknesses is computed based on the overlap of
  bottom_elevation over the bottom allocated layer. Same holds for the stage
  and top allocated layer. This matches iMOD 5.6 DISTRCOND = 2 option.
* ``by_layer_thickness``: RIV, DRN, GHB. Distribute the conductance by model
  layer thickness. This matches iMOD 5.6 DISTRCOND = 3 option.
* ``by_crosscut_transmissivity``: RIV. Distribute the conductance by
  crosscut transmissivity. Crosscut thicknesses are used to compute
  transmissivities. The crosscut thicknesses is computed based on the
  overlap of bottom_elevation over the bottom allocated layer. Same holds
  for the stage and top allocated layer. This matches iMOD 5.6 DISTRCOND = 4
  option.
* ``by_conductivity``: RIV, DRN, GHB. Distribute the conductance weighted by
  model layer hydraulic conductivities. This matches iMOD 5.6 DISTRCOND = 5
  option.
* ``by_layer_transmissivity``: RIV, DRN, GHB. Distribute the conductance by
  model layer transmissivity. This has no equivalent in iMOD 5.6.

imod.prepare.DISTRIBUTING_OPTION.by_layer_thickness
===================================================
Enumerator containing settings to distribute 2D conductance grids over
vertical layers for the RIV, DRN or GHB package.

* ``by_corrected_transmissivity``: RIV. Distribute the conductance by
  corrected transmissivities. Crosscut thicknesses are used to compute
  transmissivities. The crosscut thicknesses is computed based on the
  overlap of bottom_elevation over the bottom allocated layer. Same holds
  for the stage and top allocated layer. Furthermore the method corrects
  distribution weights for the mismatch between the midpoints of crosscut
  areas and model layer midpoints. This is the default method in iMOD 5.6,
  thus DISTRCOND = 0.
* ``equally``: RIV, DRN, GHB. Distribute conductances equally over layers.
  This matches iMOD 5.6 DISTRCOND = 1 option.
* ``by_crosscut_thickness``: RIV. Distribute the conductance by crosscut
  thicknesses. The crosscut thicknesses is computed based on the overlap of
  bottom_elevation over the bottom allocated layer. Same holds for the stage
  and top allocated layer. This matches iMOD 5.6 DISTRCOND = 2 option.
* ``by_layer_thickness``: RIV, DRN, GHB. Distribute the conductance by model
  layer thickness. This matches iMOD 5.6 DISTRCOND = 3 option.
* ``by_crosscut_transmissivity``: RIV. Distribute the conductance by
  crosscut transmissivity. Crosscut thicknesses are used to compute
  transmissivities. The crosscut thicknesses is computed based on the
  overlap of bottom_elevation over the bottom allocated layer. Same holds
  for the stage and top allocated layer. This matches iMOD 5.6 DISTRCOND = 4
  option.
* ``by_conductivity``: RIV, DRN, GHB. Distribute the conductance weighted by
  model layer hydraulic conductivities. This matches iMOD 5.6 DISTRCOND = 5
  option.
* ``by_layer_transmissivity``: RIV, DRN, GHB. Distribute the conductance by
  model layer transmissivity. This has no equivalent in iMOD 5.6.

imod.prepare.DISTRIBUTING_OPTION.by_layer_transmissivity
========================================================
Enumerator containing settings to distribute 2D conductance grids over
vertical layers for the RIV, DRN or GHB package.

* ``by_corrected_transmissivity``: RIV. Distribute the conductance by
  corrected transmissivities. Crosscut thicknesses are used to compute
  transmissivities. The crosscut thicknesses is computed based on the
  overlap of bottom_elevation over the bottom allocated layer. Same holds
  for the stage and top allocated layer. Furthermore the method corrects
  distribution weights for the mismatch between the midpoints of crosscut
  areas and model layer midpoints. This is the default method in iMOD 5.6,
  thus DISTRCOND = 0.
* ``equally``: RIV, DRN, GHB. Distribute conductances equally over layers.
  This matches iMOD 5.6 DISTRCOND = 1 option.
* ``by_crosscut_thickness``: RIV. Distribute the conductance by crosscut
  thicknesses. The crosscut thicknesses is computed based on the overlap of
  bottom_elevation over the bottom allocated layer. Same holds for the stage
  and top allocated layer. This matches iMOD 5.6 DISTRCOND = 2 option.
* ``by_layer_thickness``: RIV, DRN, GHB. Distribute the conductance by model
  layer thickness. This matches iMOD 5.6 DISTRCOND = 3 option.
* ``by_crosscut_transmissivity``: RIV. Distribute the conductance by
  crosscut transmissivity. Crosscut thicknesses are used to compute
  transmissivities. The crosscut thicknesses is computed based on the
  overlap of bottom_elevation over the bottom allocated layer. Same holds
  for the stage and top allocated layer. This matches iMOD 5.6 DISTRCOND = 4
  option.
* ``by_conductivity``: RIV, DRN, GHB. Distribute the conductance weighted by
  model layer hydraulic conductivities. This matches iMOD 5.6 DISTRCOND = 5
  option.
* ``by_layer_transmissivity``: RIV, DRN, GHB. Distribute the conductance by
  model layer transmissivity. This has no equivalent in iMOD 5.6.

imod.prepare.DISTRIBUTING_OPTION.equally
========================================
Enumerator containing settings to distribute 2D conductance grids over
vertical layers for the RIV, DRN or GHB package.

* ``by_corrected_transmissivity``: RIV. Distribute the conductance by
  corrected transmissivities. Crosscut thicknesses are used to compute
  transmissivities. The crosscut thicknesses is computed based on the
  overlap of bottom_elevation over the bottom allocated layer. Same holds
  for the stage and top allocated layer. Furthermore the method corrects
  distribution weights for the mismatch between the midpoints of crosscut
  areas and model layer midpoints. This is the default method in iMOD 5.6,
  thus DISTRCOND = 0.
* ``equally``: RIV, DRN, GHB. Distribute conductances equally over layers.
  This matches iMOD 5.6 DISTRCOND = 1 option.
* ``by_crosscut_thickness``: RIV. Distribute the conductance by crosscut
  thicknesses. The crosscut thicknesses is computed based on the overlap of
  bottom_elevation over the bottom allocated layer. Same holds for the stage
  and top allocated layer. This matches iMOD 5.6 DISTRCOND = 2 option.
* ``by_layer_thickness``: RIV, DRN, GHB. Distribute the conductance by model
  layer thickness. This matches iMOD 5.6 DISTRCOND = 3 option.
* ``by_crosscut_transmissivity``: RIV. Distribute the conductance by
  crosscut transmissivity. Crosscut thicknesses are used to compute
  transmissivities. The crosscut thicknesses is computed based on the
  overlap of bottom_elevation over the bottom allocated layer. Same holds
  for the stage and top allocated layer. This matches iMOD 5.6 DISTRCOND = 4
  option.
* ``by_conductivity``: RIV, DRN, GHB. Distribute the conductance weighted by
  model layer hydraulic conductivities. This matches iMOD 5.6 DISTRCOND = 5
  option.
* ``by_layer_transmissivity``: RIV, DRN, GHB. Distribute the conductance by
  model layer transmissivity. This has no equivalent in iMOD 5.6.

imod.prepare.allocate_drn_cells
===============================
Allocate drain cells from a planar grid across the vertical dimension.
Multiple options are available, which can be selected in ALLOCATION_OPTION.

Parameters
----------
allocation_option: ALLOCATION_OPTION
    Chosen allocation option, can be selected using the ALLOCATION_OPTION
    enumerator.
active: DataArray | UgridDatarray
    Boolean array containing active model cells. For Modflow 6, this is the
    equivalent of ``idomain == 1``.
top: DataArray | UgridDatarray
    Grid containing tops of model layers. If has no layer dimension, is
    assumed as top of upper layer and the other layers are padded with
    bottom values of the overlying model layer.
bottom: DataArray | UgridDatarray
    Grid containing bottoms of model layers.
elevation: DataArray | UgridDatarray
    Planar grid containing drain elevation. Is not allowed to have a layer
    dimension.

Returns
-------
DataArray | UgridDatarray
    Allocated drain cells

Examples
--------

>>> from imod.prepare.topsystem import ALLOCATION_OPTION, allocate_drn_cells
>>> setting = ALLOCATION_OPTION.at_elevation
>>> allocated = allocate_drn_cells(setting, active, top, bottom, stage, drain_elevation)

imod.prepare.allocate_ghb_cells
===============================
Allocate general head boundary (GHB) cells from a planar grid across the
vertical dimension. Multiple options are available, which can be selected in
ALLOCATION_OPTION.

Parameters
----------
allocation_option: ALLOCATION_OPTION
    Chosen allocation option, can be selected using the ALLOCATION_OPTION
    enumerator.
active: DataArray | UgridDatarray
    Boolean array containing active model cells. For Modflow 6, this is the
    equivalent of ``idomain == 1``.
top: DataArray | UgridDatarray
    Grid containing tops of model layers. If has no layer dimension, is
    assumed as top of upper layer and the other layers are padded with
    bottom values of the overlying model layer.
bottom: DataArray | UgridDatarray
    Grid containing bottoms of model layers.
head: DataArray | UgridDatarray
    Planar grid containing general head boundary's head. Is not allowed to
    have a layer dimension.

Returns
-------
DataArray | UgridDatarray
    Allocated general head boundary cells

Examples
--------

>>> from imod.prepare.topsystem import ALLOCATION_OPTION, allocate_ghb_cells
>>> setting = ALLOCATION_OPTION.at_elevation
>>> allocated = allocate_ghb_cells(setting, active, top, bottom, head)

imod.prepare.allocate_rch_cells
===============================
Allocate recharge cells from a planar grid across the vertical dimension.
Multiple options are available, which can be selected in ALLOCATION_OPTION.

Parameters
----------
allocation_option: ALLOCATION_OPTION
    Chosen allocation option, can be selected using the ALLOCATION_OPTION
    enumerator.
active: DataArray | UgridDataArray
    Boolean array containing active model cells. For Modflow 6, this is the
    equivalent of ``idomain == 1``.
rate: DataArray | UgridDataArray
    Array with recharge rates. This will only be used to infer where
    recharge cells are defined.

Returns
-------
DataArray | UgridDataArray
    Allocated recharge cells

Examples
--------

>>> from imod.prepare.topsystem import ALLOCATION_OPTION, allocate_rch_cells
>>> setting = ALLOCATION_OPTION.at_first_active
>>> allocated = allocate_rch_cells(setting, active, rate)

imod.prepare.allocate_riv_cells
===============================
Allocate river cells from a planar grid across the vertical dimension.
Multiple options are available, which can be selected in ALLOCATION_OPTION.

Parameters
----------
allocation_option: ALLOCATION_OPTION
    Chosen allocation option, can be selected using the ALLOCATION_OPTION
    enumerator.
active: DataArray | UgridDatarray
    Boolean array containing active model cells. For Modflow 6, this is the
    equivalent of ``idomain == 1``.
top: DataArray | UgridDatarray
    Grid containing tops of model layers. If has no layer dimension, is
    assumed as top of upper layer and the other layers are padded with
    bottom values of the overlying model layer.
bottom: DataArray | UgridDatarray
    Grid containing bottoms of model layers.
stage: DataArray | UgridDatarray
    Planar grid containing river stages. Is not allowed to have a layer
    dimension.
bottom_elevation: DataArray | UgridDatarray
    Planar grid containing river bottom elevations. Is not allowed to have a
    layer dimension.

Returns
-------
tuple(DataArray | UgridDatarray)
    Allocated river cells

Examples
--------

>>> from imod.prepare.topsystem import ALLOCATION_OPTION, allocate_riv_cells
>>> setting = ALLOCATION_OPTION.stage_to_riv_bot
>>> allocated = allocate_riv_cells(setting, active, top, bottom, stage, bottom_elevation)

imod.prepare.c_leakage
======================
Computes the phreatic leakage resistance.

Parameters
----------
kh : xr.DataArray of floats
    horizontal conductivity of phreatic aquifer
kv : xr.DataArray of floats
    vertical conductivity of phreatic aquifer
c0 : xr.DataArray of floats
    hydraulic bed resistance of water feature
c1 : xr.DataArray of floats
    hydraulic resistance of the first aquitard
D : xr.DataArray of floats
    saturated thickness of the top system
B : xr.DataArray of floats
    water feature wetted perimeter
length : xr.DataArray of floats
    water feature length per cell
dx : xr.DataArray of floats
    cellsize in x
dy : xr.DataArray of floats
    cellsize in y

Returns
-------
c_leakage: xr.DataArray of floats
    Hydraulic resistance of water features corrected for intra-cell
    hydraulic resistance and surface water interaction.

imod.prepare.c_radial
=====================
Ernst's radial resistance term to a drain.

Parameters
----------
L : xr.DataArray of floats
    distance between water features
kh : xr.DataArray of floats
    horizontal conductivity
kv : xr.DataArray of floats
    vertical conductivity
B : xr.DataArray of floats
    water feature wetted perimeter
D : xr.DataArray of floats
    saturated thickness of the top system

Returns
-------
radial_c : xr.DataArray
    Ernst's radial resistance for a drain

imod.prepare.distribute_drn_conductance
=======================================
Function to distribute 2D conductance over vertical layers for the DRN
package. Multiple options are available, which need to be selected in the
DISTRIBUTING_OPTION enumerator.

Parameters
----------
distributing_option : DISTRIBUTING_OPTION
    Distributing option available in the DISTRIBUTING_OPTION enumerator.
allocated: DataArray | UgridDataArray
    3D boolean array with drain cell locations. This can be made with the
    :func:`imod.prepare.allocate_drn_cells` function.
conductance: DataArray | UgridDataArray
    Planar grid with conductances that need to be distributed over layers,
    so grid cannot contain a layer dimension. Can contain a time dimension.
top: DataArray | UgridDataArray
    Model top
bottom: DataArray | UgridDataArray
    Model layer bottoms
k: DataArray | UgridDataArray
    Hydraulic conductivities
elevation: DataArray | UgridDataArray
    Drain elevation

Returns
-------
Conductances distributed over depth.

Examples
--------
>>> from imod.prepare import allocate_drn_cells, distribute_drn_conductance, ALLOCATION_OPTION, DISTRIBUTING_OPTION
>>> allocated = allocate_drn_cells(
    ALLOCATION_OPTION.at_elevation, active, top, bottom, drain_elevation
    )
>>> conductances_distributed = distribute_drn_conductance(
        DISTRIBUTING_OPTION.by_layer_transmissivity, allocated,
        conductance, top, bottom, k, drain_elevation
    )

imod.prepare.distribute_riv_conductance
=======================================
Function to distribute 2D conductance over vertical layers for the RIV
package. Multiple options are available, which need to be selected in the
DISTRIBUTING_OPTION enumerator.

Parameters
----------
distributing_option : DISTRIBUTING_OPTION
    Distributing option available in the DISTRIBUTING_OPTION enumerator.
allocated: DataArray | UgridDataArray
    3D boolean array with river cell locations. This can be made with the
    :func:`imod.prepare.allocate_riv_cells` function.
conductance: DataArray | UgridDataArray
    Planar grid with conductances that need to be distributed over layers,
    so grid cannot contain a layer dimension. Can contain a time dimension.
top: DataArray | UgridDataArray
    Model top
bottom: DataArray | UgridDataArray
    Model layer bottoms
k: DataArray | UgridDataArray
    Hydraulic conductivities
stage: DataArray | UgridDataArray
    Planar grid with river stages, cannot contain a layer dimension. Can
    contain a time dimension.
bottom_elevation: DataArray | UgridDataArray
    Planar grid with river bottom elevations, cannot contain a layer
    dimension. Can contain a time dimension.

Returns
-------
Conductances distributed over depth.

Examples
--------
>>> from imod.prepare import allocate_riv_cells, distribute_riv_conductance, ALLOCATION_OPTION, DISTRIBUTING_OPTION
>>> allocated = allocate_riv_cells(
    ALLOCATION_OPTION.stage_to_riv_bot, active, top, bottom, stage, bottom_elevation
    )
>>> conductances_distributed = distribute_riv_conductance(
        DISTRIBUTING_OPTION.by_corrected_transmissivity, allocated,
        conductance, top, bottom, stage, bottom_elevation, k
    )

