imod.select.cross_section_line
==============================
Obtain an interpolated cross-sectional slice through gridded data.
Utilizing the interpolation functionality in ``xarray``, this function
takes a vertical cross-sectional slice along a line through the given
data on a regular (possibly non-equidistant) grid, which is given as an
`xarray.DataArray` so that we can utilize its coordinate data.

Adapted from Metpy:
https://github.com/Unidata/MetPy/blob/master/metpy/interpolate/slices.py

Parameters
----------
data: `xarray.DataArray` or `xarray.Dataset`
    Three- (or higher) dimensional field(s) to interpolate. The DataArray
    (or each DataArray in the Dataset) must have been parsed by MetPy and
    include both an x and y coordinate dimension and the added ``crs``

    coordinate.
start: (2, ) array_like
    A latitude-longitude pair designating the start point of the cross
    section.
end: (2, ) array_like
    A latitude-longitude pair designating the end point of the cross
    section.

Returns
-------
`xarray.DataArray` or `xarray.Dataset`
    The interpolated cross section, with new dimension "s" along the
    cross-section. The cellsizes along "s" are given in the "ds" coordinate.

imod.select.cross_section_linestring
====================================
Obtain an interpolated cross-sectional slice through gridded data.
Utilizing the interpolation functionality in ``xarray``, this function
takes a vertical cross-sectional slice along a linestring through the given
data on a regular grid, which is given as an `xarray.DataArray` so that
we can utilize its coordinate data.

Adapted from Metpy:
https://github.com/Unidata/MetPy/blob/master/metpy/interpolate/slices.py

Parameters
----------
data: `xarray.DataArray` or `xarray.Dataset`
    Three- (or higher) dimensional field(s) to interpolate. The DataArray
    (or each DataArray in the Dataset) must have been parsed by MetPy and
    include both an x and y coordinate dimension and the added ``crs``

    coordinate.
linestring : shapely.geometry.LineString
    Shapely geometry designating the linestring along which to sample the
    cross section.

    Note that a LineString can easily be taken from a geopandas.GeoDataFrame
    using the .geometry attribute. Please refer to the examples.

Returns
-------
`xarray.DataArray` or `xarray.Dataset`
    The interpolated cross section, with new index dimension along the
    cross-section.

Examples
--------
Load a shapefile (that you might have drawn before using a GIS program),
take a linestring from it, and use it to extract the data for a cross
section.

>>> geodataframe = gpd.read_file("cross_section.shp")
>>> linestring = geodataframe.geometry[0]
>>> section = cross_section_linestring(data, linestring)

Or, construct the linestring directly in Python:

>>> import shapely.geometry as sg
>>> linestring = sg.LineString([(0.0, 1.0), (5.0, 5.0), (7.5, 5.0)])
>>> section = cross_section_linestring(data, linestring)

If you have drawn multiple cross sections within a shapefile, simply loop
over the linestrings:

>>> sections = [cross_section_linestring(data, ls) for ls in geodataframe.geometry]

imod.select.points_in_bounds
============================
Returns whether points specified by keyword arguments fall within the bounds
of ``da``.

Parameters
----------
da : xr.DataArray
points : keyword arguments of coordinate=values
    keyword arguments specifying coordinate and values. Please refer to the
    examples.

Returns
-------
in_bounds : np.array of bools

Examples
--------
Create the DataArray, then use the keyword arguments to define along which
coordinate to check whether the points are within bounds.

>>> nrow, ncol = 3, 4
>>> data = np.arange(12.0).reshape(nrow, ncol)
>>> coords = {"x": [0.5, 1.5, 2.5, 3.5], "y": [2.5, 1.5, 0.5]}
>>> dims = ("y", "x")
>>> da = xr.DataArray(data, coords, dims)
>>> x = [0.4, 2.6]
>>> points_in_bounds(da, x=x)

This works for an arbitrary number of coordinates:

>>> y = [1.3, 2.7]
>>> points_in_bounds(da, x=x, y=y)

imod.select.points_values
=========================
Get values from specified points.

This function will raise a ValueError if the points fall outside of the
bounds of the DataArray to avoid ambiguous behavior. Use the
``imod.select.points_in_bounds`` function to detect these points.

Parameters
----------
da : xr.DataArray
out_of_bounds : {"raise", "warn", "ignore"}, default: "raise"
    What to do if the points are not located in the bounds of the
    DataArray:
    - "raise": raise an exception
    - "warn": raise a warning, and ignore the missing points
    - "ignore": ignore the missing points
points : keyword arguments of coordinate=values
    keyword arguments specifying coordinate and values.
Returns
-------
selection : xr.DataArray

Examples
--------

>>> x = [1.0, 2.2, 3.0]
>>> y = [4.0, 5.6, 7.0]
>>> selection = imod.select.points_values(da, x=x, y=y)

imod.select.points_set_values
=============================
Set values at specified points.

This function will raise a ValueError if the points fall outside of the
bounds of the DataArray to avoid ambiguous behavior. Use the
``imod.select.points_in_bounds`` function to detect these points.

Parameters
----------
da : xr.DataArray
values : (int, float) or array of (int, float)
out_of_bounds : {"raise", "warn", "ignore"}, default: "raise"
    What to do if the points are not located in the bounds of the
    DataArray:
    - "raise": raise an exception
    - "warn": raise a warning, and ignore the missing points
    - "ignore": ignore the missing points
points : keyword arguments of coordinate=values
    keyword arguments specifying coordinate and values.

Returns
-------
da : xr.DataArray
    DataArray with values set at the point locations.

Examples
--------

>>> x = [1.0, 2.2, 3.0]
>>> y = [4.0, 5.6, 7.0]
>>> values = [10.0, 11.0, 12.0]
>>> out = imod.select.points_set_values(da, values, x=x, y=y)

imod.select.points_indices
==========================
Get the indices for points as defined by the arrays x and y.

Not all points may be located in the bounds of the DataArray. By default,
this function raises an error. This behavior can be controlled with the
``out_of_bounds`` argument. If ``out_of_bounds`` is set to "warn" or
"ignore", out of bounds point are removed. Which points have been removed
is visible in the ``index`` coordinate of the resulting selection.

Parameters
----------
da : xr.DataArray
out_of_bounds : {"raise", "warn", "ignore"}, default: "raise"
    What to do if the points are not located in the bounds of the
    DataArray:
    - "raise": raise an exception
    - "warn": raise a warning, and ignore the missing points
    - "ignore": ignore the missing points
points : keyword arguments of coordinates and values

Returns
-------
indices : dict of {coordinate: xr.DataArray with indices}

Examples
--------

To extract values:

>>> x = [1.0, 2.2, 3.0]
>>> y = [4.0, 5.6, 7.0]
>>> indices = imod.select.points_indices(da, x=x, y=y)
>>> ind_y = indices["y"]
>>> ind_x = indices["x"]
>>> selection = da.isel(x=ind_x, y=ind_y)

Or shorter, using dictionary unpacking:

>>> indices = imod.select.points_indices(da, x=x, y=y)
>>> selection = da.isel(**indices)

To set values (in a new array), the following will do the trick:

>>> empty = xr.full_like(da, np.nan)
>>> empty.data[indices["y"].values, indices["x"].values] = values_to_set

Unfortunately, at the time of writing, xarray's .sel method does not
support setting values yet. The method here works for both numpy and dask
arrays, but you'll have to manage dimensions yourself!

The ``imod.select.points_set_values()`` function will take care of the
dimensions.

imod.select.upper_active_layer
==============================
Function to get the upper active layer from ibound xarray.DataArray

Parameters
----------
da : xarray.DataArray
    A 3D DataArray
is_ibound: bool, optional
    If True, ``da`` is interpreted as ibound, with values 0: inactive, 1: active, -1 constant head.
    If False, ``upper_active_layer`` is interpreted as first layer that has data.
    Default is True.
include_constant_head : bool, optional
    If True and ``is_ibound``, also include constant head cells.
    Default is False.

Returns
-------
2d xr.DataArray of layernumber of upper active model layer

imod.select.grid_boundary_xy
============================
Return grid boundary on the xy plane.

Wraps the binary_dilation function.

Parameters
----------
grid : {xarray.DataArray, xugrid.UgridDataArray}
    Grid with either ``x`` and ``y`` dimensions or a face dimesion.

Returns
-------
{xarray.DataArray, xugrid.UgridDataArray}
    2d grid with locations of grid boundaries

imod.select.active_grid_boundary_xy
===================================
Return active boundary cells on the xy plane.

Parameters
----------
active : {xarray.DataArray, xugrid.UgridDataArray}
    Grid with active cells,
    either with ``x`` and ``y`` dimensions or a face dimesion.

Returns
-------
{xarray.DataArray, xugrid.UgridDataArray}
    Locations of active grid boundaries

