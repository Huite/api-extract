xugrid.open_dataarray
=====================
Open an DataArray from a file or file-like object containing a single
data variable.

This is designed to read netCDF files with only one data variable. If
multiple variables are present then a ValueError is raised.

Parameters
----------
filename_or_obj : str, Path, file-like or DataStore
    Strings and Path objects are interpreted as a path to a netCDF file
    or an OpenDAP URL and opened with python-netCDF4, unless the filename
    ends with .gz, in which case the file is gunzipped and opened with
    scipy.io.netcdf (only netCDF3 supported). Byte-strings or file-like
    objects are opened by scipy.io.netcdf (netCDF3) or h5py (netCDF4/HDF).
engine : {"netcdf4", "scipy", "pydap", "h5netcdf", "zarr", None}        , installed backend         or subclass of xarray.backends.BackendEntrypoint, optional
    Engine to use when reading files. If not provided, the default engine
    is chosen based on available dependencies, with a preference for
    "netcdf4".
chunks : int, dict, 'auto' or None, default: None
    If provided, used to load the data into dask arrays.

    - ``chunks='auto'`` will use dask ``auto`` chunking taking into account the
      engine preferred chunks.
    - ``chunks=None`` skips using dask, which is generally faster for
      small arrays.
    - ``chunks=-1`` loads the data with dask using a single chunk for all arrays.
    - ``chunks={}`` loads the data with dask using engine preferred chunks if
      exposed by the backend, otherwise with a single chunk for all arrays.

    See dask chunking for more details.

cache : bool, optional
    If True, cache data loaded from the underlying datastore in memory as
    NumPy arrays when accessed to avoid reading from the underlying data-
    store multiple times. Defaults to True unless you specify the `chunks`
    argument to use dask, in which case it defaults to False. Does not
    change the behavior of coordinates corresponding to dimensions, which
    always load their data from disk into a ``pandas.Index``.
decode_cf : bool, optional
    Whether to decode these variables, assuming they were saved according
    to CF conventions.
mask_and_scale : bool, optional
    If True, replace array values equal to `_FillValue` with NA and scale
    values according to the formula `original_values * scale_factor +
    add_offset`, where `_FillValue`, `scale_factor` and `add_offset` are
    taken from variable attributes (if they exist).  If the `_FillValue` or
    `missing_value` attribute contains multiple values a warning will be
    issued and all array values matching one of the multiple values will
    be replaced by NA. This keyword may not be supported by all the backends.
decode_times : bool, optional
    If True, decode times encoded in the standard NetCDF datetime format
    into datetime objects. Otherwise, leave them encoded as numbers.
    This keyword may not be supported by all the backends.
decode_timedelta : bool, optional
    If True, decode variables and coordinates with time units in
    {"days", "hours", "minutes", "seconds", "milliseconds", "microseconds"}
    into timedelta objects. If False, leave them encoded as numbers.
    If None (default), assume the same value of decode_time.
    This keyword may not be supported by all the backends.
use_cftime: bool, optional
    Only relevant if encoded dates come from a standard calendar
    (e.g. "gregorian", "proleptic_gregorian", "standard", or not
    specified).  If None (default), attempt to decode times to
    ``np.datetime64[ns]`` objects; if this is not possible, decode times to
    ``cftime.datetime`` objects. If True, always decode times to
    ``cftime.datetime`` objects, regardless of whether or not they can be
    represented using ``np.datetime64[ns]`` objects.  If False, always
    decode times to ``np.datetime64[ns]`` objects; if this is not possible
    raise an error. This keyword may not be supported by all the backends.
concat_characters : bool, optional
    If True, concatenate along the last dimension of character arrays to
    form string arrays. Dimensions will only be concatenated over (and
    removed) if they have no corresponding variable and if they are only
    used as the last dimension of character arrays.
    This keyword may not be supported by all the backends.
decode_coords : bool or {"coordinates", "all"}, optional
    Controls which variables are set as coordinate variables:

    - "coordinates" or True: Set variables referred to in the
      ``'coordinates'`` attribute of the datasets or individual variables
      as coordinate variables.
    - "all": Set variables referred to in  ``'grid_mapping'``, ``'bounds'`` and
      other attributes as coordinate variables.

    Only existing variables can be set as coordinates. Missing variables
    will be silently ignored.
drop_variables: str or iterable of str, optional
    A variable or list of variables to exclude from being parsed from the
    dataset. This may be useful to drop variables with problems or
    inconsistent values.
inline_array: bool, default: False
    How to include the array in the dask task graph.
    By default(``inline_array=False``) the array is included in a task by
    itself, and each chunk refers to that task by its key. With
    ``inline_array=True``, Dask will instead inline the array directly
    in the values of the task graph. See :py:func:`dask.array.from_array`.
chunked_array_type: str, optional
    Which chunked array type to coerce the underlying data array to.
    Defaults to 'dask' if installed, else whatever is registered via the `ChunkManagerEnetryPoint` system.
    Experimental API that should not be relied upon.
from_array_kwargs: dict
    Additional keyword arguments passed on to the `ChunkManagerEntrypoint.from_array` method used to create
    chunked arrays, via whichever chunk manager is specified through the `chunked_array_type` kwarg.
    For example if :py:func:`dask.array.Array` objects are used for chunking, additional kwargs will be passed
    to :py:func:`dask.array.from_array`. Experimental API that should not be relied upon.
backend_kwargs: dict
    Additional keyword arguments passed on to the engine open function,
    equivalent to `**kwargs`.
**kwargs: dict
    Additional keyword arguments passed on to the engine open function.
    For example:

    - 'group': path to the netCDF4 group in the given file to open given as
      a str,supported by "netcdf4", "h5netcdf", "zarr".
    - 'lock': resource lock to use when reading data from disk. Only
      relevant when using dask or another form of parallelism. By default,
      appropriate locks are chosen to safely read and write files with the
      currently active dask scheduler. Supported by "netcdf4", "h5netcdf",
      "scipy".

    See engine open function for kwargs accepted by each specific engine.

Notes
-----
This is designed to be fully compatible with `DataArray.to_netcdf`. Saving
using `DataArray.to_netcdf` and then loading with this function will
produce an identical result.

All parameters are passed directly to `xarray.open_dataset`. See that
documentation for further details.

See also
--------
open_dataset

xugrid.open_dataset
===================
Open and decode a dataset from a file or file-like object.

Parameters
----------
filename_or_obj : str, Path, file-like or DataStore
    Strings and Path objects are interpreted as a path to a netCDF file
    or an OpenDAP URL and opened with python-netCDF4, unless the filename
    ends with .gz, in which case the file is gunzipped and opened with
    scipy.io.netcdf (only netCDF3 supported). Byte-strings or file-like
    objects are opened by scipy.io.netcdf (netCDF3) or h5py (netCDF4/HDF).
engine : {"netcdf4", "scipy", "pydap", "h5netcdf", "zarr", None}        , installed backend         or subclass of xarray.backends.BackendEntrypoint, optional
    Engine to use when reading files. If not provided, the default engine
    is chosen based on available dependencies, with a preference for
    "netcdf4". A custom backend class (a subclass of ``BackendEntrypoint``)
    can also be used.
chunks : int, dict, 'auto' or None, default: None
    If provided, used to load the data into dask arrays.

    - ``chunks="auto"`` will use dask ``auto`` chunking taking into account the
      engine preferred chunks.
    - ``chunks=None`` skips using dask, which is generally faster for
      small arrays.
    - ``chunks=-1`` loads the data with dask using a single chunk for all arrays.
    - ``chunks={}`` loads the data with dask using the engine's preferred chunk
      size, generally identical to the format's chunk size. If not available, a
      single chunk for all arrays.

    See dask chunking for more details.
cache : bool, optional
    If True, cache data loaded from the underlying datastore in memory as
    NumPy arrays when accessed to avoid reading from the underlying data-
    store multiple times. Defaults to True unless you specify the `chunks`
    argument to use dask, in which case it defaults to False. Does not
    change the behavior of coordinates corresponding to dimensions, which
    always load their data from disk into a ``pandas.Index``.
decode_cf : bool, optional
    Whether to decode these variables, assuming they were saved according
    to CF conventions.
mask_and_scale : bool or dict-like, optional
    If True, replace array values equal to `_FillValue` with NA and scale
    values according to the formula `original_values * scale_factor +
    add_offset`, where `_FillValue`, `scale_factor` and `add_offset` are
    taken from variable attributes (if they exist).  If the `_FillValue` or
    `missing_value` attribute contains multiple values a warning will be
    issued and all array values matching one of the multiple values will
    be replaced by NA. Pass a mapping, e.g. ``{"my_variable": False}``,
    to toggle this feature per-variable individually.
    This keyword may not be supported by all the backends.
decode_times : bool or dict-like, optional
    If True, decode times encoded in the standard NetCDF datetime format
    into datetime objects. Otherwise, leave them encoded as numbers.
    Pass a mapping, e.g. ``{"my_variable": False}``,
    to toggle this feature per-variable individually.
    This keyword may not be supported by all the backends.
decode_timedelta : bool or dict-like, optional
    If True, decode variables and coordinates with time units in
    {"days", "hours", "minutes", "seconds", "milliseconds", "microseconds"}
    into timedelta objects. If False, leave them encoded as numbers.
    If None (default), assume the same value of decode_time.
    Pass a mapping, e.g. ``{"my_variable": False}``,
    to toggle this feature per-variable individually.
    This keyword may not be supported by all the backends.
use_cftime: bool or dict-like, optional
    Only relevant if encoded dates come from a standard calendar
    (e.g. "gregorian", "proleptic_gregorian", "standard", or not
    specified).  If None (default), attempt to decode times to
    ``np.datetime64[ns]`` objects; if this is not possible, decode times to
    ``cftime.datetime`` objects. If True, always decode times to
    ``cftime.datetime`` objects, regardless of whether or not they can be
    represented using ``np.datetime64[ns]`` objects.  If False, always
    decode times to ``np.datetime64[ns]`` objects; if this is not possible
    raise an error. Pass a mapping, e.g. ``{"my_variable": False}``,
    to toggle this feature per-variable individually.
    This keyword may not be supported by all the backends.
concat_characters : bool or dict-like, optional
    If True, concatenate along the last dimension of character arrays to
    form string arrays. Dimensions will only be concatenated over (and
    removed) if they have no corresponding variable and if they are only
    used as the last dimension of character arrays.
    Pass a mapping, e.g. ``{"my_variable": False}``,
    to toggle this feature per-variable individually.
    This keyword may not be supported by all the backends.
decode_coords : bool or {"coordinates", "all"}, optional
    Controls which variables are set as coordinate variables:

    - "coordinates" or True: Set variables referred to in the
      ``'coordinates'`` attribute of the datasets or individual variables
      as coordinate variables.
    - "all": Set variables referred to in  ``'grid_mapping'``, ``'bounds'`` and
      other attributes as coordinate variables.

    Only existing variables can be set as coordinates. Missing variables
    will be silently ignored.
drop_variables: str or iterable of str, optional
    A variable or list of variables to exclude from being parsed from the
    dataset. This may be useful to drop variables with problems or
    inconsistent values.
inline_array: bool, default: False
    How to include the array in the dask task graph.
    By default(``inline_array=False``) the array is included in a task by
    itself, and each chunk refers to that task by its key. With
    ``inline_array=True``, Dask will instead inline the array directly
    in the values of the task graph. See :py:func:`dask.array.from_array`.
chunked_array_type: str, optional
    Which chunked array type to coerce this datasets' arrays to.
    Defaults to 'dask' if installed, else whatever is registered via the `ChunkManagerEnetryPoint` system.
    Experimental API that should not be relied upon.
from_array_kwargs: dict
    Additional keyword arguments passed on to the `ChunkManagerEntrypoint.from_array` method used to create
    chunked arrays, via whichever chunk manager is specified through the `chunked_array_type` kwarg.
    For example if :py:func:`dask.array.Array` objects are used for chunking, additional kwargs will be passed
    to :py:func:`dask.array.from_array`. Experimental API that should not be relied upon.
backend_kwargs: dict
    Additional keyword arguments passed on to the engine open function,
    equivalent to `**kwargs`.
**kwargs: dict
    Additional keyword arguments passed on to the engine open function.
    For example:

    - 'group': path to the netCDF4 group in the given file to open given as
      a str,supported by "netcdf4", "h5netcdf", "zarr".
    - 'lock': resource lock to use when reading data from disk. Only
      relevant when using dask or another form of parallelism. By default,
      appropriate locks are chosen to safely read and write files with the
      currently active dask scheduler. Supported by "netcdf4", "h5netcdf",
      "scipy".

    See engine open function for kwargs accepted by each specific engine.

Returns
-------
dataset : Dataset
    The newly created dataset.

Notes
-----
``open_dataset`` opens the file with read-only access. When you modify
values of a Dataset, even one linked to files on disk, only the in-memory
copy you are manipulating in xarray is modified: the original file on disk
is never touched.

See Also
--------
open_mfdataset

xugrid.open_mfdataset
=====================
Open multiple files as a single dataset.

If combine='by_coords' then the function ``combine_by_coords`` is used to combine
the datasets into one before returning the result, and if combine='nested' then
``combine_nested`` is used. The filepaths must be structured according to which
combining function is used, the details of which are given in the documentation for
``combine_by_coords`` and ``combine_nested``. By default ``combine='by_coords'``
will be used. Requires dask to be installed. See documentation for
details on dask [1]_. Global attributes from the ``attrs_file`` are used
for the combined dataset.

Parameters
----------
paths : str or nested sequence of paths
    Either a string glob in the form ``"path/to/my/files/*.nc"`` or an explicit list of
    files to open. Paths can be given as strings or as pathlib Paths. If
    concatenation along more than one dimension is desired, then ``paths`` must be a
    nested list-of-lists (see ``combine_nested`` for details). (A string glob will
    be expanded to a 1-dimensional list.)
chunks : int, dict, 'auto' or None, optional
    Dictionary with keys given by dimension names and values given by chunk sizes.
    In general, these should divide the dimensions of each dataset. If int, chunk
    each dimension by ``chunks``. By default, chunks will be chosen to load entire
    input files into memory at once. This has a major impact on performance: please
    see the full documentation for more details [2]_. This argument is evaluated
    on a per-file basis, so chunk sizes that span multiple files will be ignored.
concat_dim : str, DataArray, Index or a Sequence of these or None, optional
    Dimensions to concatenate files along.  You only need to provide this argument
    if ``combine='nested'``, and if any of the dimensions along which you want to
    concatenate is not a dimension in the original datasets, e.g., if you want to
    stack a collection of 2D arrays along a third dimension. Set
    ``concat_dim=[..., None, ...]`` explicitly to disable concatenation along a
    particular dimension. Default is None, which for a 1D list of filepaths is
    equivalent to opening the files separately and then merging them with
    ``xarray.merge``.
combine : {"by_coords", "nested"}, optional
    Whether ``xarray.combine_by_coords`` or ``xarray.combine_nested`` is used to
    combine all the data. Default is to use ``xarray.combine_by_coords``.
compat : {"identical", "equals", "broadcast_equals",               "no_conflicts", "override"}, default: "no_conflicts"
    String indicating how to compare variables of the same name for
    potential conflicts when merging:

     * "broadcast_equals": all values must be equal when variables are
       broadcast against each other to ensure common dimensions.
     * "equals": all values and dimensions must be the same.
     * "identical": all values, dimensions and attributes must be the
       same.
     * "no_conflicts": only values which are not null in both datasets
       must be equal. The returned dataset then contains the combination
       of all non-null values.
     * "override": skip comparing and pick variable from first dataset

preprocess : callable, optional
    If provided, call this function on each dataset prior to concatenation.
    You can find the file-name from which each dataset was loaded in
    ``ds.encoding["source"]``.
engine : {"netcdf4", "scipy", "pydap", "h5netcdf", "zarr", None}        , installed backend         or subclass of xarray.backends.BackendEntrypoint, optional
    Engine to use when reading files. If not provided, the default engine
    is chosen based on available dependencies, with a preference for
    "netcdf4".
data_vars : {"minimal", "different", "all"} or list of str, default: "all"
    These data variables will be concatenated together:
      * "minimal": Only data variables in which the dimension already
        appears are included.
      * "different": Data variables which are not equal (ignoring
        attributes) across all datasets are also concatenated (as well as
        all for which dimension already appears). Beware: this option may
        load the data payload of data variables into memory if they are not
        already loaded.
      * "all": All data variables will be concatenated.
      * list of str: The listed data variables will be concatenated, in
        addition to the "minimal" data variables.
coords : {"minimal", "different", "all"} or list of str, optional
    These coordinate variables will be concatenated together:
     * "minimal": Only coordinates in which the dimension already appears
       are included.
     * "different": Coordinates which are not equal (ignoring attributes)
       across all datasets are also concatenated (as well as all for which
       dimension already appears). Beware: this option may load the data
       payload of coordinate variables into memory if they are not already
       loaded.
     * "all": All coordinate variables will be concatenated, except
       those corresponding to other dimensions.
     * list of str: The listed coordinate variables will be concatenated,
       in addition the "minimal" coordinates.
parallel : bool, default: False
    If True, the open and preprocess steps of this function will be
    performed in parallel using ``dask.delayed``. Default is False.
join : {"outer", "inner", "left", "right", "exact", "override"}, default: "outer"
    String indicating how to combine differing indexes
    (excluding concat_dim) in objects

    - "outer": use the union of object indexes
    - "inner": use the intersection of object indexes
    - "left": use indexes from the first object with each dimension
    - "right": use indexes from the last object with each dimension
    - "exact": instead of aligning, raise `ValueError` when indexes to be
      aligned are not equal
    - "override": if indexes are of same size, rewrite indexes to be
      those of the first object with that dimension. Indexes for the same
      dimension must have the same size in all objects.
attrs_file : str or path-like, optional
    Path of the file used to read global attributes from.
    By default global attributes are read from the first file provided,
    with wildcard matches sorted by filename.
combine_attrs : {"drop", "identical", "no_conflicts", "drop_conflicts",                      "override"} or callable, default: "override"
    A callable or a string indicating how to combine attrs of the objects being
    merged:

    - "drop": empty attrs on returned Dataset.
    - "identical": all attrs must be the same on every object.
    - "no_conflicts": attrs from all objects are combined, any that have
      the same name must also have the same value.
    - "drop_conflicts": attrs from all objects are combined, any that have
      the same name but different values are dropped.
    - "override": skip comparing and copy attrs from the first dataset to
      the result.

    If a callable, it must expect a sequence of ``attrs`` dicts and a context object
    as its only parameters.
**kwargs : optional
    Additional arguments passed on to :py:func:`xarray.open_dataset`. For an
    overview of some of the possible options, see the documentation of
    :py:func:`xarray.open_dataset`

Returns
-------
xarray.Dataset

Notes
-----
``open_mfdataset`` opens files with read-only access. When you modify values
of a Dataset, even one linked to files on disk, only the in-memory copy you
are manipulating in xarray is modified: the original file on disk is never
touched.

See Also
--------
combine_by_coords
combine_nested
open_dataset

Examples
--------
A user might want to pass additional arguments into ``preprocess`` when
applying some operation to many individual files that are being opened. One route
to do this is through the use of ``functools.partial``.

>>> from functools import partial
>>> def _preprocess(x, lon_bnds, lat_bnds):
...     return x.sel(lon=slice(*lon_bnds), lat=slice(*lat_bnds))
...
>>> lon_bnds, lat_bnds = (-110, -105), (40, 45)
>>> partial_func = partial(_preprocess, lon_bnds=lon_bnds, lat_bnds=lat_bnds)
>>> ds = xr.open_mfdataset(
...     "file_*.nc", concat_dim="time", preprocess=partial_func
... )  # doctest: +SKIP

It is also possible to use any argument to ``open_dataset`` together
with ``open_mfdataset``, such as for example ``drop_variables``:

>>> ds = xr.open_mfdataset(
...     "file.nc", drop_variables=["varname_1", "varname_2"]  # any list of vars
... )  # doctest: +SKIP

References
----------

.. [1] https://docs.xarray.dev/en/stable/dask.html
.. [2] https://docs.xarray.dev/en/stable/dask.html#chunking-and-performance

xugrid.open_zarr
================
Load and decode a dataset from a Zarr store.

The `store` object should be a valid store for a Zarr group. `store`
variables must contain dimension metadata encoded in the
`_ARRAY_DIMENSIONS` attribute or must have NCZarr format.

Parameters
----------
store : MutableMapping or str
    A MutableMapping where a Zarr Group has been stored or a path to a
    directory in file system where a Zarr DirectoryStore has been stored.
synchronizer : object, optional
    Array synchronizer provided to zarr
group : str, optional
    Group path. (a.k.a. `path` in zarr terminology.)
chunks : int, dict, 'auto' or None, default: 'auto'
    If provided, used to load the data into dask arrays.

    - ``chunks='auto'`` will use dask ``auto`` chunking taking into account the
      engine preferred chunks.
    - ``chunks=None`` skips using dask, which is generally faster for
      small arrays.
    - ``chunks=-1`` loads the data with dask using a single chunk for all arrays.
    - ``chunks={}`` loads the data with dask using engine preferred chunks if
      exposed by the backend, otherwise with a single chunk for all arrays.

    See dask chunking for more details.
overwrite_encoded_chunks : bool, optional
    Whether to drop the zarr chunks encoded for each variable when a
    dataset is loaded with specified chunk sizes (default: False)
decode_cf : bool, optional
    Whether to decode these variables, assuming they were saved according
    to CF conventions.
mask_and_scale : bool, optional
    If True, replace array values equal to `_FillValue` with NA and scale
    values according to the formula `original_values * scale_factor +
    add_offset`, where `_FillValue`, `scale_factor` and `add_offset` are
    taken from variable attributes (if they exist).  If the `_FillValue` or
    `missing_value` attribute contains multiple values a warning will be
    issued and all array values matching one of the multiple values will
    be replaced by NA.
decode_times : bool, optional
    If True, decode times encoded in the standard NetCDF datetime format
    into datetime objects. Otherwise, leave them encoded as numbers.
concat_characters : bool, optional
    If True, concatenate along the last dimension of character arrays to
    form string arrays. Dimensions will only be concatenated over (and
    removed) if they have no corresponding variable and if they are only
    used as the last dimension of character arrays.
decode_coords : bool, optional
    If True, decode the 'coordinates' attribute to identify coordinates in
    the resulting dataset.
drop_variables : str or iterable, optional
    A variable or list of variables to exclude from being parsed from the
    dataset. This may be useful to drop variables with problems or
    inconsistent values.
consolidated : bool, optional
    Whether to open the store using zarr's consolidated metadata
    capability. Only works for stores that have already been consolidated.
    By default (`consolidate=None`), attempts to read consolidated metadata,
    falling back to read non-consolidated metadata if that fails.

    When the experimental ``zarr_version=3``, ``consolidated`` must be
    either be ``None`` or ``False``.
chunk_store : MutableMapping, optional
    A separate Zarr store only for chunk data.
storage_options : dict, optional
    Any additional parameters for the storage backend (ignored for local
    paths).
decode_timedelta : bool, optional
    If True, decode variables and coordinates with time units in
    {'days', 'hours', 'minutes', 'seconds', 'milliseconds', 'microseconds'}
    into timedelta objects. If False, leave them encoded as numbers.
    If None (default), assume the same value of decode_time.
use_cftime : bool, optional
    Only relevant if encoded dates come from a standard calendar
    (e.g. "gregorian", "proleptic_gregorian", "standard", or not
    specified).  If None (default), attempt to decode times to
    ``np.datetime64[ns]`` objects; if this is not possible, decode times to
    ``cftime.datetime`` objects. If True, always decode times to
    ``cftime.datetime`` objects, regardless of whether or not they can be
    represented using ``np.datetime64[ns]`` objects.  If False, always
    decode times to ``np.datetime64[ns]`` objects; if this is not possible
    raise an error.
zarr_version : int or None, optional
    The desired zarr spec version to target (currently 2 or 3). The default
    of None will attempt to determine the zarr version from ``store`` when
    possible, otherwise defaulting to 2.
chunked_array_type: str, optional
    Which chunked array type to coerce this datasets' arrays to.
    Defaults to 'dask' if installed, else whatever is registered via the `ChunkManagerEntryPoint` system.
    Experimental API that should not be relied upon.
from_array_kwargs: dict, optional
    Additional keyword arguments passed on to the `ChunkManagerEntrypoint.from_array` method used to create
    chunked arrays, via whichever chunk manager is specified through the `chunked_array_type` kwarg.
    Defaults to {'manager': 'dask'}, meaning additional kwargs will be passed eventually to
    :py:func:`dask.array.from_array`. Experimental API that should not be relied upon.

Returns
-------
dataset : Dataset
    The newly created dataset.

See Also
--------
open_dataset
open_mfdataset

References
----------
http://zarr.readthedocs.io/

xugrid.full_like
================
Return a new object with the same shape and type as a given object.

Returned object will be chunked if if the given object is chunked, or if chunks or chunked_array_type are specified.

Parameters
----------
other : DataArray, Dataset or Variable
    The reference object in input
fill_value : scalar or dict-like
    Value to fill the new object with before returning it. If
    other is a Dataset, may also be a dict-like mapping data
    variables to fill values.
dtype : dtype or dict-like of dtype, optional
    dtype of the new array. If a dict-like, maps dtypes to
    variables. If omitted, it defaults to other.dtype.
chunks : int, "auto", tuple of int or mapping of Hashable to int, optional
    Chunk sizes along each dimension, e.g., ``5``, ``"auto"``, ``(5, 5)`` or
    ``{"x": 5, "y": 5}``.
chunked_array_type: str, optional
    Which chunked array type to coerce the underlying data array to.
    Defaults to 'dask' if installed, else whatever is registered via the `ChunkManagerEnetryPoint` system.
    Experimental API that should not be relied upon.
from_array_kwargs: dict, optional
    Additional keyword arguments passed on to the `ChunkManagerEntrypoint.from_array` method used to create
    chunked arrays, via whichever chunk manager is specified through the `chunked_array_type` kwarg.
    For example, with dask as the default chunked array type, this method would pass additional kwargs
    to :py:func:`dask.array.from_array`. Experimental API that should not be relied upon.

Returns
-------
out : same as object
    New object with the same shape and type as other, with the data
    filled with fill_value. Coords will be copied from other.
    If other is based on dask, the new one will be as well, and will be
    split in the same chunks.

Examples
--------
>>> x = xr.DataArray(
...     np.arange(6).reshape(2, 3),
...     dims=["lat", "lon"],
...     coords={"lat": [1, 2], "lon": [0, 1, 2]},
... )
>>> x
<xarray.DataArray (lat: 2, lon: 3)> Size: 48B
array([[0, 1, 2],
       [3, 4, 5]])
Coordinates:
  * lat      (lat) int64 16B 1 2
  * lon      (lon) int64 24B 0 1 2

>>> xr.full_like(x, 1)
<xarray.DataArray (lat: 2, lon: 3)> Size: 48B
array([[1, 1, 1],
       [1, 1, 1]])
Coordinates:
  * lat      (lat) int64 16B 1 2
  * lon      (lon) int64 24B 0 1 2

>>> xr.full_like(x, 0.5)
<xarray.DataArray (lat: 2, lon: 3)> Size: 48B
array([[0, 0, 0],
       [0, 0, 0]])
Coordinates:
  * lat      (lat) int64 16B 1 2
  * lon      (lon) int64 24B 0 1 2

>>> xr.full_like(x, 0.5, dtype=np.double)
<xarray.DataArray (lat: 2, lon: 3)> Size: 48B
array([[0.5, 0.5, 0.5],
       [0.5, 0.5, 0.5]])
Coordinates:
  * lat      (lat) int64 16B 1 2
  * lon      (lon) int64 24B 0 1 2

>>> xr.full_like(x, np.nan, dtype=np.double)
<xarray.DataArray (lat: 2, lon: 3)> Size: 48B
array([[nan, nan, nan],
       [nan, nan, nan]])
Coordinates:
  * lat      (lat) int64 16B 1 2
  * lon      (lon) int64 24B 0 1 2

>>> ds = xr.Dataset(
...     {"a": ("x", [3, 5, 2]), "b": ("x", [9, 1, 0])}, coords={"x": [2, 4, 6]}
... )
>>> ds
<xarray.Dataset> Size: 72B
Dimensions:  (x: 3)
Coordinates:
  * x        (x) int64 24B 2 4 6
Data variables:
    a        (x) int64 24B 3 5 2
    b        (x) int64 24B 9 1 0
>>> xr.full_like(ds, fill_value={"a": 1, "b": 2})
<xarray.Dataset> Size: 72B
Dimensions:  (x: 3)
Coordinates:
  * x        (x) int64 24B 2 4 6
Data variables:
    a        (x) int64 24B 1 1 1
    b        (x) int64 24B 2 2 2
>>> xr.full_like(ds, fill_value={"a": 1, "b": 2}, dtype={"a": bool, "b": float})
<xarray.Dataset> Size: 51B
Dimensions:  (x: 3)
Coordinates:
  * x        (x) int64 24B 2 4 6
Data variables:
    a        (x) bool 3B True True True
    b        (x) float64 24B 2.0 2.0 2.0

See Also
--------
zeros_like
ones_like

xugrid.ones_like
================
Return a new object of ones with the same shape and
type as a given dataarray or dataset.

Parameters
----------
other : DataArray, Dataset, or Variable
    The reference object. The output will have the same dimensions and coordinates as this object.
dtype : dtype, optional
    dtype of the new array. If omitted, it defaults to other.dtype.
chunks : int, "auto", tuple of int or mapping of Hashable to int, optional
    Chunk sizes along each dimension, e.g., ``5``, ``"auto"``, ``(5, 5)`` or
    ``{"x": 5, "y": 5}``.
chunked_array_type: str, optional
    Which chunked array type to coerce the underlying data array to.
    Defaults to 'dask' if installed, else whatever is registered via the `ChunkManagerEnetryPoint` system.
    Experimental API that should not be relied upon.
from_array_kwargs: dict, optional
    Additional keyword arguments passed on to the `ChunkManagerEntrypoint.from_array` method used to create
    chunked arrays, via whichever chunk manager is specified through the `chunked_array_type` kwarg.
    For example, with dask as the default chunked array type, this method would pass additional kwargs
    to :py:func:`dask.array.from_array`. Experimental API that should not be relied upon.

Returns
-------
out : same as object
    New object of ones with the same shape and type as other.

Examples
--------
>>> x = xr.DataArray(
...     np.arange(6).reshape(2, 3),
...     dims=["lat", "lon"],
...     coords={"lat": [1, 2], "lon": [0, 1, 2]},
... )
>>> x
<xarray.DataArray (lat: 2, lon: 3)> Size: 48B
array([[0, 1, 2],
       [3, 4, 5]])
Coordinates:
  * lat      (lat) int64 16B 1 2
  * lon      (lon) int64 24B 0 1 2

>>> xr.ones_like(x)
<xarray.DataArray (lat: 2, lon: 3)> Size: 48B
array([[1, 1, 1],
       [1, 1, 1]])
Coordinates:
  * lat      (lat) int64 16B 1 2
  * lon      (lon) int64 24B 0 1 2

See Also
--------
zeros_like
full_like

xugrid.zeros_like
=================
Return a new object of zeros with the same shape and
type as a given dataarray or dataset.

Parameters
----------
other : DataArray, Dataset or Variable
    The reference object. The output will have the same dimensions and coordinates as this object.
dtype : dtype, optional
    dtype of the new array. If omitted, it defaults to other.dtype.
chunks : int, "auto", tuple of int or mapping of Hashable to int, optional
    Chunk sizes along each dimension, e.g., ``5``, ``"auto"``, ``(5, 5)`` or
    ``{"x": 5, "y": 5}``.
chunked_array_type: str, optional
    Which chunked array type to coerce the underlying data array to.
    Defaults to 'dask' if installed, else whatever is registered via the `ChunkManagerEnetryPoint` system.
    Experimental API that should not be relied upon.
from_array_kwargs: dict, optional
    Additional keyword arguments passed on to the `ChunkManagerEntrypoint.from_array` method used to create
    chunked arrays, via whichever chunk manager is specified through the `chunked_array_type` kwarg.
    For example, with dask as the default chunked array type, this method would pass additional kwargs
    to :py:func:`dask.array.from_array`. Experimental API that should not be relied upon.

Returns
-------
out : DataArray, Dataset or Variable
    New object of zeros with the same shape and type as other.

Examples
--------
>>> x = xr.DataArray(
...     np.arange(6).reshape(2, 3),
...     dims=["lat", "lon"],
...     coords={"lat": [1, 2], "lon": [0, 1, 2]},
... )
>>> x
<xarray.DataArray (lat: 2, lon: 3)> Size: 48B
array([[0, 1, 2],
       [3, 4, 5]])
Coordinates:
  * lat      (lat) int64 16B 1 2
  * lon      (lon) int64 24B 0 1 2

>>> xr.zeros_like(x)
<xarray.DataArray (lat: 2, lon: 3)> Size: 48B
array([[0, 0, 0],
       [0, 0, 0]])
Coordinates:
  * lat      (lat) int64 16B 1 2
  * lon      (lon) int64 24B 0 1 2

>>> xr.zeros_like(x, dtype=float)
<xarray.DataArray (lat: 2, lon: 3)> Size: 48B
array([[0., 0., 0.],
       [0., 0., 0.]])
Coordinates:
  * lat      (lat) int64 16B 1 2
  * lon      (lon) int64 24B 0 1 2

See Also
--------
ones_like
full_like

xugrid.concat
=============
Concatenate xarray objects along a new or existing dimension.

Parameters
----------
objs : sequence of Dataset and DataArray
    xarray objects to concatenate together. Each object is expected to
    consist of variables and coordinates with matching shapes except for
    along the concatenated dimension.
dim : Hashable or Variable or DataArray or pandas.Index
    Name of the dimension to concatenate along. This can either be a new
    dimension name, in which case it is added along axis=0, or an existing
    dimension name, in which case the location of the dimension is
    unchanged. If dimension is provided as a Variable, DataArray or Index, its name
    is used as the dimension to concatenate along and the values are added
    as a coordinate.
data_vars : {"minimal", "different", "all"} or list of Hashable, optional
    These data variables will be concatenated together:
      * "minimal": Only data variables in which the dimension already
        appears are included.
      * "different": Data variables which are not equal (ignoring
        attributes) across all datasets are also concatenated (as well as
        all for which dimension already appears). Beware: this option may
        load the data payload of data variables into memory if they are not
        already loaded.
      * "all": All data variables will be concatenated.
      * list of dims: The listed data variables will be concatenated, in
        addition to the "minimal" data variables.

    If objects are DataArrays, data_vars must be "all".
coords : {"minimal", "different", "all"} or list of Hashable, optional
    These coordinate variables will be concatenated together:
      * "minimal": Only coordinates in which the dimension already appears
        are included.
      * "different": Coordinates which are not equal (ignoring attributes)
        across all datasets are also concatenated (as well as all for which
        dimension already appears). Beware: this option may load the data
        payload of coordinate variables into memory if they are not already
        loaded.
      * "all": All coordinate variables will be concatenated, except
        those corresponding to other dimensions.
      * list of Hashable: The listed coordinate variables will be concatenated,
        in addition to the "minimal" coordinates.
compat : {"identical", "equals", "broadcast_equals", "no_conflicts", "override"}, optional
    String indicating how to compare non-concatenated variables of the same name for
    potential conflicts. This is passed down to merge.

    - "broadcast_equals": all values must be equal when variables are
      broadcast against each other to ensure common dimensions.
    - "equals": all values and dimensions must be the same.
    - "identical": all values, dimensions and attributes must be the
      same.
    - "no_conflicts": only values which are not null in both datasets
      must be equal. The returned dataset then contains the combination
      of all non-null values.
    - "override": skip comparing and pick variable from first dataset
positions : None or list of integer arrays, optional
    List of integer arrays which specifies the integer positions to which
    to assign each dataset along the concatenated dimension. If not
    supplied, objects are concatenated in the provided order.
fill_value : scalar or dict-like, optional
    Value to use for newly missing values. If a dict-like, maps
    variable names to fill values. Use a data array's name to
    refer to its values.
join : {"outer", "inner", "left", "right", "exact"}, optional
    String indicating how to combine differing indexes
    (excluding dim) in objects

    - "outer": use the union of object indexes
    - "inner": use the intersection of object indexes
    - "left": use indexes from the first object with each dimension
    - "right": use indexes from the last object with each dimension
    - "exact": instead of aligning, raise `ValueError` when indexes to be
      aligned are not equal
    - "override": if indexes are of same size, rewrite indexes to be
      those of the first object with that dimension. Indexes for the same
      dimension must have the same size in all objects.
combine_attrs : {"drop", "identical", "no_conflicts", "drop_conflicts",                      "override"} or callable, default: "override"
    A callable or a string indicating how to combine attrs of the objects being
    merged:

    - "drop": empty attrs on returned Dataset.
    - "identical": all attrs must be the same on every object.
    - "no_conflicts": attrs from all objects are combined, any that have
      the same name must also have the same value.
    - "drop_conflicts": attrs from all objects are combined, any that have
      the same name but different values are dropped.
    - "override": skip comparing and copy attrs from the first dataset to
      the result.

    If a callable, it must expect a sequence of ``attrs`` dicts and a context object
    as its only parameters.
create_index_for_new_dim : bool, default: True
    Whether to create a new ``PandasIndex`` object when the objects being concatenated contain scalar variables named ``dim``.

Returns
-------
concatenated : type of objs

See also
--------
merge

Examples
--------
>>> da = xr.DataArray(
...     np.arange(6).reshape(2, 3), [("x", ["a", "b"]), ("y", [10, 20, 30])]
... )
>>> da
<xarray.DataArray (x: 2, y: 3)> Size: 48B
array([[0, 1, 2],
       [3, 4, 5]])
Coordinates:
  * x        (x) <U1 8B 'a' 'b'
  * y        (y) int64 24B 10 20 30

>>> xr.concat([da.isel(y=slice(0, 1)), da.isel(y=slice(1, None))], dim="y")
<xarray.DataArray (x: 2, y: 3)> Size: 48B
array([[0, 1, 2],
       [3, 4, 5]])
Coordinates:
  * x        (x) <U1 8B 'a' 'b'
  * y        (y) int64 24B 10 20 30

>>> xr.concat([da.isel(x=0), da.isel(x=1)], "x")
<xarray.DataArray (x: 2, y: 3)> Size: 48B
array([[0, 1, 2],
       [3, 4, 5]])
Coordinates:
  * x        (x) <U1 8B 'a' 'b'
  * y        (y) int64 24B 10 20 30

>>> xr.concat([da.isel(x=0), da.isel(x=1)], "new_dim")
<xarray.DataArray (new_dim: 2, y: 3)> Size: 48B
array([[0, 1, 2],
       [3, 4, 5]])
Coordinates:
    x        (new_dim) <U1 8B 'a' 'b'
  * y        (y) int64 24B 10 20 30
Dimensions without coordinates: new_dim

>>> xr.concat([da.isel(x=0), da.isel(x=1)], pd.Index([-90, -100], name="new_dim"))
<xarray.DataArray (new_dim: 2, y: 3)> Size: 48B
array([[0, 1, 2],
       [3, 4, 5]])
Coordinates:
    x        (new_dim) <U1 8B 'a' 'b'
  * y        (y) int64 24B 10 20 30
  * new_dim  (new_dim) int64 16B -90 -100

# Concatenate a scalar variable along a new dimension of the same name with and without creating a new index

>>> ds = xr.Dataset(coords={"x": 0})
>>> xr.concat([ds, ds], dim="x")
<xarray.Dataset> Size: 16B
Dimensions:  (x: 2)
Coordinates:
  * x        (x) int64 16B 0 0
Data variables:
    *empty*

>>> xr.concat([ds, ds], dim="x").indexes
Indexes:
    x        Index([0, 0], dtype='int64', name='x')

>>> xr.concat([ds, ds], dim="x", create_index_for_new_dim=False).indexes
Indexes:
    *empty*

xugrid.merge
============
Merge any number of xarray objects into a single Dataset as variables.

Parameters
----------
objects : iterable of Dataset or iterable of DataArray or iterable of dict-like
    Merge together all variables from these objects. If any of them are
    DataArray objects, they must have a name.
compat : {"identical", "equals", "broadcast_equals", "no_conflicts",               "override", "minimal"}, default: "no_conflicts"
    String indicating how to compare variables of the same name for
    potential conflicts:

    - "identical": all values, dimensions and attributes must be the
      same.
    - "equals": all values and dimensions must be the same.
    - "broadcast_equals": all values must be equal when variables are
      broadcast against each other to ensure common dimensions.
    - "no_conflicts": only values which are not null in both datasets
      must be equal. The returned dataset then contains the combination
      of all non-null values.
    - "override": skip comparing and pick variable from first dataset
    - "minimal": drop conflicting coordinates

join : {"outer", "inner", "left", "right", "exact", "override"}, default: "outer"
    String indicating how to combine differing indexes in objects.

    - "outer": use the union of object indexes
    - "inner": use the intersection of object indexes
    - "left": use indexes from the first object with each dimension
    - "right": use indexes from the last object with each dimension
    - "exact": instead of aligning, raise `ValueError` when indexes to be
      aligned are not equal
    - "override": if indexes are of same size, rewrite indexes to be
      those of the first object with that dimension. Indexes for the same
      dimension must have the same size in all objects.

fill_value : scalar or dict-like, optional
    Value to use for newly missing values. If a dict-like, maps
    variable names to fill values. Use a data array's name to
    refer to its values.
combine_attrs : {"drop", "identical", "no_conflicts", "drop_conflicts",                      "override"} or callable, default: "override"
    A callable or a string indicating how to combine attrs of the objects being
    merged:

    - "drop": empty attrs on returned Dataset.
    - "identical": all attrs must be the same on every object.
    - "no_conflicts": attrs from all objects are combined, any that have
      the same name must also have the same value.
    - "drop_conflicts": attrs from all objects are combined, any that have
      the same name but different values are dropped.
    - "override": skip comparing and copy attrs from the first dataset to
      the result.

    If a callable, it must expect a sequence of ``attrs`` dicts and a context object
    as its only parameters.

Returns
-------
Dataset
    Dataset with combined variables from each object.

Examples
--------
>>> x = xr.DataArray(
...     [[1.0, 2.0], [3.0, 5.0]],
...     dims=("lat", "lon"),
...     coords={"lat": [35.0, 40.0], "lon": [100.0, 120.0]},
...     name="var1",
... )
>>> y = xr.DataArray(
...     [[5.0, 6.0], [7.0, 8.0]],
...     dims=("lat", "lon"),
...     coords={"lat": [35.0, 42.0], "lon": [100.0, 150.0]},
...     name="var2",
... )
>>> z = xr.DataArray(
...     [[0.0, 3.0], [4.0, 9.0]],
...     dims=("time", "lon"),
...     coords={"time": [30.0, 60.0], "lon": [100.0, 150.0]},
...     name="var3",
... )

>>> x
<xarray.DataArray 'var1' (lat: 2, lon: 2)> Size: 32B
array([[1., 2.],
       [3., 5.]])
Coordinates:
  * lat      (lat) float64 16B 35.0 40.0
  * lon      (lon) float64 16B 100.0 120.0

>>> y
<xarray.DataArray 'var2' (lat: 2, lon: 2)> Size: 32B
array([[5., 6.],
       [7., 8.]])
Coordinates:
  * lat      (lat) float64 16B 35.0 42.0
  * lon      (lon) float64 16B 100.0 150.0

>>> z
<xarray.DataArray 'var3' (time: 2, lon: 2)> Size: 32B
array([[0., 3.],
       [4., 9.]])
Coordinates:
  * time     (time) float64 16B 30.0 60.0
  * lon      (lon) float64 16B 100.0 150.0

>>> xr.merge([x, y, z])
<xarray.Dataset> Size: 256B
Dimensions:  (lat: 3, lon: 3, time: 2)
Coordinates:
  * lat      (lat) float64 24B 35.0 40.0 42.0
  * lon      (lon) float64 24B 100.0 120.0 150.0
  * time     (time) float64 16B 30.0 60.0
Data variables:
    var1     (lat, lon) float64 72B 1.0 2.0 nan 3.0 5.0 nan nan nan nan
    var2     (lat, lon) float64 72B 5.0 nan 6.0 nan nan nan 7.0 nan 8.0
    var3     (time, lon) float64 48B 0.0 nan 3.0 4.0 nan 9.0

>>> xr.merge([x, y, z], compat="identical")
<xarray.Dataset> Size: 256B
Dimensions:  (lat: 3, lon: 3, time: 2)
Coordinates:
  * lat      (lat) float64 24B 35.0 40.0 42.0
  * lon      (lon) float64 24B 100.0 120.0 150.0
  * time     (time) float64 16B 30.0 60.0
Data variables:
    var1     (lat, lon) float64 72B 1.0 2.0 nan 3.0 5.0 nan nan nan nan
    var2     (lat, lon) float64 72B 5.0 nan 6.0 nan nan nan 7.0 nan 8.0
    var3     (time, lon) float64 48B 0.0 nan 3.0 4.0 nan 9.0

>>> xr.merge([x, y, z], compat="equals")
<xarray.Dataset> Size: 256B
Dimensions:  (lat: 3, lon: 3, time: 2)
Coordinates:
  * lat      (lat) float64 24B 35.0 40.0 42.0
  * lon      (lon) float64 24B 100.0 120.0 150.0
  * time     (time) float64 16B 30.0 60.0
Data variables:
    var1     (lat, lon) float64 72B 1.0 2.0 nan 3.0 5.0 nan nan nan nan
    var2     (lat, lon) float64 72B 5.0 nan 6.0 nan nan nan 7.0 nan 8.0
    var3     (time, lon) float64 48B 0.0 nan 3.0 4.0 nan 9.0

>>> xr.merge([x, y, z], compat="equals", fill_value=-999.0)
<xarray.Dataset> Size: 256B
Dimensions:  (lat: 3, lon: 3, time: 2)
Coordinates:
  * lat      (lat) float64 24B 35.0 40.0 42.0
  * lon      (lon) float64 24B 100.0 120.0 150.0
  * time     (time) float64 16B 30.0 60.0
Data variables:
    var1     (lat, lon) float64 72B 1.0 2.0 -999.0 3.0 ... -999.0 -999.0 -999.0
    var2     (lat, lon) float64 72B 5.0 -999.0 6.0 -999.0 ... 7.0 -999.0 8.0
    var3     (time, lon) float64 48B 0.0 -999.0 3.0 4.0 -999.0 9.0

>>> xr.merge([x, y, z], join="override")
<xarray.Dataset> Size: 144B
Dimensions:  (lat: 2, lon: 2, time: 2)
Coordinates:
  * lat      (lat) float64 16B 35.0 40.0
  * lon      (lon) float64 16B 100.0 120.0
  * time     (time) float64 16B 30.0 60.0
Data variables:
    var1     (lat, lon) float64 32B 1.0 2.0 3.0 5.0
    var2     (lat, lon) float64 32B 5.0 6.0 7.0 8.0
    var3     (time, lon) float64 32B 0.0 3.0 4.0 9.0

>>> xr.merge([x, y, z], join="inner")
<xarray.Dataset> Size: 64B
Dimensions:  (lat: 1, lon: 1, time: 2)
Coordinates:
  * lat      (lat) float64 8B 35.0
  * lon      (lon) float64 8B 100.0
  * time     (time) float64 16B 30.0 60.0
Data variables:
    var1     (lat, lon) float64 8B 1.0
    var2     (lat, lon) float64 8B 5.0
    var3     (time, lon) float64 16B 0.0 4.0

>>> xr.merge([x, y, z], compat="identical", join="inner")
<xarray.Dataset> Size: 64B
Dimensions:  (lat: 1, lon: 1, time: 2)
Coordinates:
  * lat      (lat) float64 8B 35.0
  * lon      (lon) float64 8B 100.0
  * time     (time) float64 16B 30.0 60.0
Data variables:
    var1     (lat, lon) float64 8B 1.0
    var2     (lat, lon) float64 8B 5.0
    var3     (time, lon) float64 16B 0.0 4.0

>>> xr.merge([x, y, z], compat="broadcast_equals", join="outer")
<xarray.Dataset> Size: 256B
Dimensions:  (lat: 3, lon: 3, time: 2)
Coordinates:
  * lat      (lat) float64 24B 35.0 40.0 42.0
  * lon      (lon) float64 24B 100.0 120.0 150.0
  * time     (time) float64 16B 30.0 60.0
Data variables:
    var1     (lat, lon) float64 72B 1.0 2.0 nan 3.0 5.0 nan nan nan nan
    var2     (lat, lon) float64 72B 5.0 nan 6.0 nan nan nan 7.0 nan 8.0
    var3     (time, lon) float64 48B 0.0 nan 3.0 4.0 nan 9.0

>>> xr.merge([x, y, z], join="exact")
Traceback (most recent call last):
...
ValueError: cannot align objects with join='exact' where ...

Raises
------
xarray.MergeError
    If any variables with the same name have conflicting values.

See also
--------
concat
combine_nested
combine_by_coords

xugrid.merge_partitions
=======================
Merge topology and data, partitioned along UGRID dimensions, into a single
UgridDataset.

UGRID topologies and variables are merged if they share a name. Topologies
and variables must be present in *all* partitions. Dimension names must
match.

Variables are omitted from the merged result if non-UGRID dimensions do not
match in size.

Parameters
----------
partitions : sequence of UgridDataset or UgridDataArray
merge_ugrid_chunks: bool, default is True.
    Whether to merge chunks along the UGRID topology dimensions.

Returns
-------
merged : UgridDataset

xugrid.burn_vector_geometry
===========================
Burn vector geometries (points, lines, polygons) into a Ugrid2d mesh.

If no ``column`` argument is provided, a value of 1.0 will be burned in to
the mesh.

Parameters
----------
gdf: geopandas.GeoDataFrame
    Polygons, points, and/or lines to be burned into the grid.
like: UgridDataArray, UgridDataset, or Ugrid2d
    Grid to burn the vector data into.
column: str, optional
    Name of the geodataframe column of which to the values to burn into
    grid.
fill: int, float, optional, default value ``np.nan``.
    Fill value for nodata areas.
all_touched: bool, optional, default value ``False``.
    All mesh faces (cells) touched by polygons will be updated, not just
    those whose center point is within the polygon.

Returns
-------
burned: UgridDataArray

xugrid.earcut_triangulate_polygons
==================================
Break down polygons using mapbox_earcut, and create a mesh from the
resulting triangles.

If no ``column`` argument is provided, the polygon index will be assigned
to the grid faces.

Parameters
----------
polygons: geopandas.GeoDataFrame
    Polygons to convert to triangles.
column: str, optional
    Name of the geodataframe column of which to the values to assign
    to the grid faces.

Returns
-------
triangulated: UgridDataArray

xugrid.polygonize
=================
Polygonize a UgridDataArray.

This function creates vector polygons for all connected regions of cells
(faces) in the Ugrid2d topology sharing a common value.

The produced polygon edges will follow exactly the cell boundaries. When
the data consists of many unique values (e.g. unbinned elevation data), the
result will essentially be one polygon per face. In such cases, it is much
more efficient to use ``xugrid.UgridDataArray.to_geodataframe``, which
directly converts every cell to a polygon. This function is meant for data
with relatively few unique values such as classification results.

Parameters
----------
uda: UgridDataArray
    The DataArray should only contain the face dimension. Additional
    dimensions, such as time, are not allowed.

Returns
-------
polygonized: GeoDataFrame

xugrid.UgridDataArray.ugrid
===========================
UGRID Accessor. This "accessor" makes operations using the UGRID
topology available.

xugrid.UgridDataArray.from_structured
=====================================
Create a UgridDataArray from a (structured) xarray DataArray.

The spatial dimensions are flattened into a single UGRID face dimension.

By default, this method looks for the "x" and "y" coordinates and assumes
they are one-dimensional. To convert rotated or curvilinear coordinates,
provide the names of the x and y coordinates.

Parameters
----------
da: xr.DataArray
    Last two dimensions must be ``("y", "x")``.
x: str, default: None
    Which coordinate to use as the UGRID x-coordinate.
y: str, default: None
    Which coordinate to use as the UGRID y-coordinate.

Returns
-------
unstructured: UgridDataArray

xugrid.UgridDataArray.from_data
===============================
Create a UgridDataArray from a grid and a 1D array of values.

Parameters
----------
data: array like
    Values for this array. Must be a ``numpy.ndarray`` or castable to
    it.
grid: Ugrid1d, Ugrid2d
facet: str
    With which facet to associate the data. Options for Ugrid1d are,
    ``"node"`` or ``"edge"``. Options for Ugrid2d are ``"node"``,
    ``"edge"``, or ``"face"``.

Returns
-------
uda: UgridDataArray

xugrid.UgridDataset.ugrid
=========================
UGRID Accessor. This "accessor" makes operations using the UGRID
topology available.

xugrid.UgridDataset.from_geodataframe
=====================================
Convert a geodataframe into the appropriate Ugrid topology and dataset.

Parameters
----------
geodataframe: gpd.GeoDataFrame

Returns
-------
dataset: UGridDataset

xugrid.UgridDatasetAccessor
===========================
Helper class that provides a standard way to create an ABC using
inheritance.

xugrid.UgridDatasetAccessor Class Members
=========================================
   * xugrid.UgridDatasetAccessor.assign_edge_coords
   * xugrid.UgridDatasetAccessor.assign_face_coords
   * xugrid.UgridDatasetAccessor.assign_node_coords
   * xugrid.UgridDatasetAccessor.bounds
   * xugrid.UgridDatasetAccessor.clip_box
   * xugrid.UgridDatasetAccessor.crs
   * xugrid.UgridDatasetAccessor.grid
   * xugrid.UgridDatasetAccessor.intersect_line
   * xugrid.UgridDatasetAccessor.intersect_linestring
   * xugrid.UgridDatasetAccessor.name
   * xugrid.UgridDatasetAccessor.names
   * xugrid.UgridDatasetAccessor.partition
   * xugrid.UgridDatasetAccessor.partition_by_label
   * xugrid.UgridDatasetAccessor.rasterize
   * xugrid.UgridDatasetAccessor.rasterize_like
   * xugrid.UgridDatasetAccessor.reindex_like
   * xugrid.UgridDatasetAccessor.rename
   * xugrid.UgridDatasetAccessor.sel_points
   * xugrid.UgridDatasetAccessor.set_crs
   * xugrid.UgridDatasetAccessor.set_node_coords
   * xugrid.UgridDatasetAccessor.to_crs
   * xugrid.UgridDatasetAccessor.to_dataset
   * xugrid.UgridDatasetAccessor.to_geodataframe
   * xugrid.UgridDatasetAccessor.to_netcdf
   * xugrid.UgridDatasetAccessor.to_nonperiodic
   * xugrid.UgridDatasetAccessor.to_periodic
   * xugrid.UgridDatasetAccessor.to_zarr
   * xugrid.UgridDatasetAccessor.topology
   * xugrid.UgridDatasetAccessor.total_bounds

xugrid.UgridDatasetAccessor.assign_edge_coords
==============================================
Assign edge coordinates from the grid to the object.

Returns a new object with all the original data in addition to the new
node coordinates of the grid.

Returns
-------
assigned: UgridDataset

xugrid.UgridDatasetAccessor.assign_face_coords
==============================================
Assign face coordinates from the grid to the object.

Returns a new object with all the original data in addition to the new
node coordinates of the grid.

Returns
-------
assigned: UgridDataset

xugrid.UgridDatasetAccessor.assign_node_coords
==============================================
Assign node coordinates from the grid to the object.

Returns a new object with all the original data in addition to the new
node coordinates of the grid.

Returns
-------
assigned: UgridDataset

xugrid.UgridDatasetAccessor.bounds
==================================
Mapping from grid name to tuple containing ``minx, miny, maxx, maxy``
values of the grid's node coordinates for every grid in the dataset.

xugrid.UgridDatasetAccessor.clip_box
====================================
Clip the DataArray or Dataset by a bounding box.

Parameters
----------
xmin: float
ymin: float
xmax: float
ymax: float

-------
clipped:
    xugrid.UgridDataArray or xugrid.UgridDataset

xugrid.UgridDatasetAccessor.crs
===============================
The Coordinate Reference System (CRS) represented as a ``pyproj.CRS`` object.

Returns None if the CRS is not set.

Returns
-------
crs: dict
    A dictionary containing the names of the grids and their CRS.

xugrid.UgridDatasetAccessor.grid
================================
Returns the single UGRID topology in this dataset. Raises a TypeError if
the dataset contains more than one topology.

xugrid.UgridDatasetAccessor.intersect_line
==========================================
Intersect a line with the grid of this data, and fetch the values of
the intersected faces.

Parameters
----------
obj: xr.DataArray or xr.Dataset
start: sequence of two floats
    coordinate pair (x, y), designating the start point of the line.
end: sequence of two floats
    coordinate pair (x, y), designating the end point of the line.

Returns
-------
intersection: xr.Dataset
    The name of the topology is prefixed in the x, y and s
    (spatium=distance) coordinates.

xugrid.UgridDatasetAccessor.intersect_linestring
================================================
Intersect the grid along a collection of linestrings. Returns a new Dataset
with the values for each intersected segment.

Parameters
----------
linestring: shapely.LineString

Returns
-------
intersection: xr.Dataset
    The name of the topology is prefixed in the x, y and s
    (spatium=distance) coordinates.

xugrid.UgridDatasetAccessor.name
================================
Returns name of the single UGRID topology in this dataset. Raises a
TypeError if the dataset contains more than one topology.

xugrid.UgridDatasetAccessor.names
=================================
Names of all the UGRID topologies in the dataset.

xugrid.UgridDatasetAccessor.partition
=====================================
Partition a grid into a given number of parts.

Parameters
----------
n_part: integer
    The number of parts to partition the mesh.

Returns
-------
partitioned: list of partitions

xugrid.UgridDatasetAccessor.partition_by_label
==============================================
Partition a grid by labels.

Parameters
----------
labels: np.ndarray of integers labeling each face.

Returns
-------
partitioned: list of partitions

xugrid.UgridDatasetAccessor.rasterize
=====================================
Rasterize all face data on 2D unstructured grids by sampling.

Parameters
----------
resolution: float
    Spacing in x and y.

Returns
-------
rasterized: xr.Dataset

xugrid.UgridDatasetAccessor.rasterize_like
==========================================
Rasterize unstructured all face data on 2D unstructured grids by
sampling on the x and y coordinates of ``other``.

Parameters
----------
resolution: float
    Spacing in x and y.
other: Union[xr.DataArray, xr.Dataset]
    Object to take x and y coordinates from.

Returns
-------
rasterized: xr.Dataset

xugrid.UgridDatasetAccessor.reindex_like
========================================
Conform this object to match the topology of another object. The
topologies must be exactly equivalent: only the order of the nodes,
edges, and faces may differ.

Topologies are matched by name, and dimension names must match for
equivalent topologies.

Parameters
----------
other: Ugrid1d, Ugrid2d, UgridDataArray, UgridDataset
obj: DataArray or Dataset
tolerance: float, default value 0.0.
    Maximum distance between inexact coordinate matches.

Returns
-------
reindexed: UgridDataset

xugrid.UgridDatasetAccessor.rename
==================================
Give a new name to the UGRID topology and update the associated
coordinate and dimension names in the Dataset.

Parameters
----------
new_name_or_name_dict: str or dict
    If the argument is a string, the new name of the topology. This
    only works if the dataset contains a single UGRID topology. If the
    argument is a dict, it used as a mapping from old names to new
    names.

xugrid.UgridDatasetAccessor.sel_points
======================================
Select points in the unstructured grid.

Out-of-bounds points are ignored. They may be identified via the
``index`` coordinate of the returned selection.

Parameters
----------
x: ndarray of floats with shape ``(n_points,)``
y: ndarray of floats with shape ``(n_points,)``
out_of_bounds: str, default ``"warn"``
    What to do when points are located outside of any feature:

    * raise: raise a ValueError.
    * ignore: return ``fill_value`` for the out of bounds points.
    * warn: give a warning and return NaN for the out of bounds points.
    * drop: drop the out of bounds points. They may be identified
      via the ``index`` coordinate of the returned selection.
fill_value: scalar, DataArray, Dataset, or callable, optional, default: np.nan
    Value to assign to out-of-bounds points if out_of_bounds is warn
    or ignore. Forwarded to xarray's ``.where()`` method.

Returns
-------
points: Union[xr.DataArray, xr.Dataset]
    The name of the topology is prefixed in the x, y coordinates.

xugrid.UgridDatasetAccessor.set_crs
===================================
Set the Coordinate Reference System (CRS) of a UGRID topology.

NOTE: The underlying geometries are not transformed to this CRS. To
transform the geometries to a new CRS, use the ``to_crs`` method.

Parameters
----------
crs : pyproj.CRS, optional if `epsg` is specified
    The value can be anything accepted
    by :meth:`pyproj.CRS.from_user_input() <pyproj.crs.CRS.from_user_input>`,
    such as an authority string (eg "EPSG:4326") or a WKT string.
epsg : int, optional if `crs` is specified
    EPSG code specifying the projection.
allow_override : bool, default False
    If the the UGRID topology already has a CRS, allow to replace the
    existing CRS, even when both are not equal.
topology: str, optional
    Name of the grid topology in which to set the CRS.
    Sets the CRS for all grids if left unspecified.

xugrid.UgridDatasetAccessor.set_node_coords
===========================================
Given names of x and y coordinates of the nodes of an object, set them
as the coordinates in the grid.

Parameters
----------
node_x: str
    Name of the x coordinate of the nodes in the object.
node_y: str
    Name of the y coordinate of the nodes in the object.
topology: str, optional
    Name of the grid topology in which to set the node_x and node_y
    coordinates. Can be omitted if the UgridDataset contains only a
    single grid.

xugrid.UgridDatasetAccessor.to_crs
==================================
Transform geometries to a new coordinate reference system.
Transform all geometries in an active geometry column to a different coordinate
reference system. The ``crs`` attribute on the current Ugrid must
be set. Either ``crs`` or ``epsg`` may be specified for output.

This method will transform all points in all objects. It has no notion
of projecting the cells. All segments joining points are assumed to be
lines in the current projection, not geodesics. Objects crossing the
dateline (or other projection boundary) will have undesirable behavior.

Parameters
----------
crs : pyproj.CRS, optional if `epsg` is specified
    The value can be anything accepted by
    :meth:`pyproj.CRS.from_user_input() <pyproj.crs.CRS.from_user_input>`,
    such as an authority string (eg "EPSG:4326") or a WKT string.
epsg : int, optional if `crs` is specified
    EPSG code specifying output projection.
topology: str, optional
    Name of the grid topology to reproject.
    Reprojects all grids if left unspecified.

xugrid.UgridDatasetAccessor.to_dataset
======================================
Convert this UgridDataset into a standard
xarray.Dataset.

The UGRID topology information is added as standard data variables.

Parameters
----------
optional_attributes: bool, default: False.
    Whether to generate the UGRID optional attributes.

Returns
-------
dataset: UgridDataset

xugrid.UgridDatasetAccessor.to_geodataframe
===========================================
Convert data and topology of one facet (node, edge, face) of the grid
to a geopandas GeoDataFrame. This also determines the geometry type of
the geodataframe:

* node: point
* edge: line
* face: polygon

Parameters
----------
name: str
    Name to give to the array (required if unnamed).
dim_order:
    Hierarchical dimension order for the resulting dataframe. Array content is
    transposed to this order and then written out as flat vectors in contiguous
    order, so the last dimension in this list will be contiguous in the resulting
    DataFrame. This has a major influence on which operations are efficient on the
    resulting dataframe.

    If provided, must include all dimensions of this DataArray. By default,
    dimensions are sorted according to the DataArray dimensions order.

Returns
-------
geodataframe: gpd.GeoDataFrame

xugrid.UgridDatasetAccessor.to_netcdf
=====================================
Write dataset contents to a UGRID compliant netCDF file.

This function wraps :py:meth:`xr.Dataset.to_netcdf`; it adds the UGRID
variables and coordinates to a standard xarray Dataset, then writes the
result to a netCDF.

All arguments are forwarded to :py:meth:`xr.Dataset.to_netcdf`.

xugrid.UgridDatasetAccessor.to_nonperiodic
==========================================
Convert the grid from a periodic grid (where the rightmost boundary shares its
nodes with the leftmost boundary) to an aperiodic grid, where the leftmost nodes
are separate from the rightmost nodes.

Parameters
----------
xmax: float
    The x-value of the newly created rightmost boundary nodes.

Returns
-------
nonperiodic: UgridDataset

xugrid.UgridDatasetAccessor.to_periodic
=======================================
Convert every grid to a periodic grid, where the rightmost boundary
shares its nodes with the leftmost boundary.

Returns
-------
periodic: UgridDataset

xugrid.UgridDatasetAccessor.to_zarr
===================================
Write dataset contents to a UGRID compliant Zarr file.

This function wraps :py:meth:`xr.Dataset.to_zarr`; it adds the UGRID
variables and coordinates to a standard xarray Dataset, then writes the
result to a Zarr file.

All arguments are forwarded to :py:meth:`xr.Dataset.to_zarr`.

xugrid.UgridDatasetAccessor.topology
====================================
Mapping from names to UGRID topologies.

xugrid.UgridDatasetAccessor.total_bounds
========================================
Returns a tuple containing ``minx, miny, maxx, maxy`` values for the
bounds of the dataset as a whole. Currently does not check whether the
coordinate reference systems (CRS) of the grids in the dataset match.

xugrid.UgridDataArrayAccessor
=============================
This "accessor" makes operations using the UGRID topology available via the
``.ugrid`` attribute for UgridDataArrays and UgridDatasets.

xugrid.UgridDataArrayAccessor Class Members
===========================================
   * xugrid.UgridDataArrayAccessor.assign_edge_coords
   * xugrid.UgridDataArrayAccessor.assign_face_coords
   * xugrid.UgridDataArrayAccessor.assign_node_coords
   * xugrid.UgridDataArrayAccessor.binary_dilation
   * xugrid.UgridDataArrayAccessor.binary_erosion
   * xugrid.UgridDataArrayAccessor.bounds
   * xugrid.UgridDataArrayAccessor.clip_box
   * xugrid.UgridDataArrayAccessor.connected_components
   * xugrid.UgridDataArrayAccessor.crs
   * xugrid.UgridDataArrayAccessor.grids
   * xugrid.UgridDataArrayAccessor.interpolate_na
   * xugrid.UgridDataArrayAccessor.intersect_line
   * xugrid.UgridDataArrayAccessor.intersect_linestring
   * xugrid.UgridDataArrayAccessor.laplace_interpolate
   * xugrid.UgridDataArrayAccessor.name
   * xugrid.UgridDataArrayAccessor.names
   * xugrid.UgridDataArrayAccessor.partition
   * xugrid.UgridDataArrayAccessor.partition_by_label
   * xugrid.UgridDataArrayAccessor.plot
   * xugrid.UgridDataArrayAccessor.rasterize
   * xugrid.UgridDataArrayAccessor.rasterize_like
   * xugrid.UgridDataArrayAccessor.reindex_like
   * xugrid.UgridDataArrayAccessor.rename
   * xugrid.UgridDataArrayAccessor.reverse_cuthill_mckee
   * xugrid.UgridDataArrayAccessor.sel
   * xugrid.UgridDataArrayAccessor.sel_points
   * xugrid.UgridDataArrayAccessor.set_crs
   * xugrid.UgridDataArrayAccessor.set_node_coords
   * xugrid.UgridDataArrayAccessor.to_crs
   * xugrid.UgridDataArrayAccessor.to_dataset
   * xugrid.UgridDataArrayAccessor.to_geodataframe
   * xugrid.UgridDataArrayAccessor.to_netcdf
   * xugrid.UgridDataArrayAccessor.to_nonperiodic
   * xugrid.UgridDataArrayAccessor.to_periodic
   * xugrid.UgridDataArrayAccessor.to_zarr
   * xugrid.UgridDataArrayAccessor.topology
   * xugrid.UgridDataArrayAccessor.total_bounds

xugrid.UgridDataArrayAccessor.assign_edge_coords
================================================
Assign edge coordinates from the grid to the object.

Returns a new object with all the original data in addition to the new
node coordinates of the grid.

Returns
-------
assigned: UgridDataset

xugrid.UgridDataArrayAccessor.assign_face_coords
================================================
Assign face coordinates from the grid to the object.

Returns a new object with all the original data in addition to the new
node coordinates of the grid.

Returns
-------
assigned: UgridDataset

xugrid.UgridDataArrayAccessor.assign_node_coords
================================================
Assign node coordinates from the grid to the object.

Returns a new object with all the original data in addition to the new
node coordinates of the grid.

Returns
-------
assigned: UgridDataset

xugrid.UgridDataArrayAccessor.binary_dilation
=============================================
Binary dilation can be used on a boolean array to expand the "shape" of
features.

Compare with :py:func:`scipy.ndimage.binary_dilation`.

Parameters
----------
iterations: int, default: 1
mask: 1d array of bool, optional
border_value: bool, default value: False

Returns
-------
dilated: UgridDataArray

xugrid.UgridDataArrayAccessor.binary_erosion
============================================
Binary erosion can be used on a boolean array to shrink the "shape" of
features.

Compare with :py:func:`scipy.ndimage.binary_erosion`.

Parameters
----------
iterations: int, default: 1
mask: 1d array of bool, optional
border_value: bool, default value: False

Returns
-------
eroded: UgridDataArray

xugrid.UgridDataArrayAccessor.bounds
====================================
Mapping from grid name to tuple containing ``minx, miny, maxx, maxy``
values of the grid's node coordinates.

xugrid.UgridDataArrayAccessor.connected_components
==================================================
Every edge or face is given a component number. If all are connected,
all will have the same number.

Wraps :py:func:`scipy.sparse.csgraph.connected_components``.

Returns
-------
labelled: UgridDataArray

xugrid.UgridDataArrayAccessor.crs
=================================
The Coordinate Reference System (CRS) represented as a ``pyproj.CRS`` object.

Returns None if the CRS is not set.

Returns
-------
crs: dict
    A dictionary containing the name of the grid and its CRS.

xugrid.UgridDataArrayAccessor.grids
===================================
The UGRID topology of this DataArry, as a list. Included for
consistency with UgridDataset.

xugrid.UgridDataArrayAccessor.interpolate_na
============================================
Fill in NaNs by interpolating.

This function automatically finds the UGRID dimension and broadcasts
over the other dimensions.

Parameters
----------
method: str, default is "nearest"
    Currently the only supported method.
max_distance: nonnegative float, optional.
    Use ``None`` for no maximum distance.

Returns
-------
filled: UgridDataArray of floats

xugrid.UgridDataArrayAccessor.intersect_line
============================================
Intersect a line with the grid of this data, and fetch the values of
the intersected faces.

Parameters
----------
obj: xr.DataArray or xr.Dataset
start: sequence of two floats
    coordinate pair (x, y), designating the start point of the line.
end: sequence of two floats
    coordinate pair (x, y), designating the end point of the line.

Returns
-------
intersection: xr.DataArray
    The length along the line is returned as the "s" coordinate.

xugrid.UgridDataArrayAccessor.intersect_linestring
==================================================
Intersect the grid along a collection of linestrings. Returns a new DataArray
with the values for each intersected segment.

Parameters
----------
linestring: shapely.LineString

Returns
-------
intersection: xr.DataArray
    The length along the linestring is returned as the "s" coordinate.

xugrid.UgridDataArrayAccessor.laplace_interpolate
=================================================
Fill in NaNs by using Laplace interpolation.

This function automatically finds the UGRID dimension and broadcasts
over the other dimensions.

This solves Laplace's equation where where there is no data, with data
values functioning as fixed potential boundary conditions.

Note that an iterative solver method will be required for large grids.
In this case, some experimentation with the solver settings may be
required to find a converging solution of sufficient accuracy. Refer to
the documentation of :py:func:`scipy.sparse.linalg.spilu` and
:py:func:`scipy.sparse.linalg.cg`.

Data can be interpolated from nodes or faces. Direct interpolation of edge
associated data is not allowed. Instead, create node associated data first,
then translate that data to the edges.

Parameters
----------
xy_weights: bool, default False.
    Wether to use the inverse of the centroid to centroid distance in
    the coefficient matrix. If ``False``, defaults to uniform
    coefficients of 1 so that each face connection has equal weight.
direct_solve: bool, optional, default ``False``
    Whether to use a direct or an iterative solver or a conjugate gradient
    solver. Direct method provides an exact answer, but are unsuitable
    for large problems.
delta: float, default 0.0.
    ILU0 preconditioner non-diagonally dominant correction.
relax: float, default 0.0.
    Modified ILU0 preconditioner relaxation factor.
rtol: float, optional, default 1.0e-5.
    Convergence tolerance for ``scipy.sparse.linalg.cg``.
atol: float, optional, default 0.0.
    Convergence tolerance for ``scipy.sparse.linalg.cg``.
maxiter: int, default 500.
    Maximum number of iterations for ``scipy.sparse.linalg.cg``.

Returns
-------
filled: UgridDataArray of floats

xugrid.UgridDataArrayAccessor.name
==================================
Name of the UGRID topology of this DataArray.

xugrid.UgridDataArrayAccessor.names
===================================
Name of the UGRID topology, as a list. Included for consistency with
UgridDataset.

xugrid.UgridDataArrayAccessor.plot
==================================
Enables use of plot functions as attributes.
For example UgridDataArray.ugrid.plot.pcolormesh()

xugrid.UgridDataArrayAccessor.rasterize
=======================================
Rasterize unstructured grid by sampling.

Parameters
----------
resolution: float
    Spacing in x and y.

Returns
-------
rasterized: xr.DataArray

xugrid.UgridDataArrayAccessor.rasterize_like
============================================
Rasterize unstructured grid by sampling on the x and y coordinates
of ``other``.

Parameters
----------
resolution: float
    Spacing in x and y.
other: Union[xr.DataArray, xr.Dataset]
    Object to take x and y coordinates from.

Returns
-------
rasterized: xr.DataArray

xugrid.UgridDataArrayAccessor.reindex_like
==========================================
Conform this object to match the topology of another object. The
topologies must be exactly equivalent: only the order of the nodes,
edges, and faces may differ.

Dimension names must match for equivalent topologies.

Parameters
----------
other: Ugrid1d, Ugrid2d, UgridDataArray, UgridDataset
obj: DataArray or Dataset
tolerance: float, default value 0.0.
    Maximum distance between inexact coordinate matches.

Returns
-------
reindexed: UgridDataArray

xugrid.UgridDataArrayAccessor.rename
====================================
Give a new name to the UGRID topology and update the associated
coordinate and dimension names in the DataArray.

Parameters
----------
name: str
    The new name of the topology.

xugrid.UgridDataArrayAccessor.reverse_cuthill_mckee
===================================================
Reduces bandwith of the connectivity matrix.

Wraps :py:func:`scipy.sparse.csgraph.reverse_cuthill_mckee`.

Returns
-------
reordered: Union[UgridDataArray, UgridDataset]

xugrid.UgridDataArrayAccessor.sel
=================================
Return a new object, a subselection in the UGRID x and y coordinates.

The indexing for x and y always occurs orthogonally, i.e.:
``.sel(x=[0.0, 5.0], y=[10.0, 15.0])`` results in a four points. For
vectorized indexing (equal to ``zip``ing through x and y), see
``.sel_points``.

Depending on the nature of the x and y indexers, a xugrid or xarray
object is returned:

* slice without step: ``x=slice(-100, 100)``: returns xugrid object, a
  part of the unstructured grid.
* slice with step: ``x=slice(-100, 100, 10)``: returns xarray object, a
  series of points (x=[-100, -90, -80, ..., 90, 100]).
* a scalar: ``x=5.0``: returns xarray object, a point.
* an array: ``x=[1.0, 15.0, 17.0]``: returns xarray object, a series of
  points.

Parameters
----------
x: float, 1d array, slice
y: float, 1d array, slice

Returns
-------
selection: Union[UgridDataArray, UgridDataset, xr.DataArray, xr.Dataset]

xugrid.UgridDataArrayAccessor.sel_points
========================================
Select points in the unstructured grid.

Out-of-bounds points are ignored. They may be identified via the
``index`` coordinate of the returned selection.

Parameters
----------
x: ndarray of floats with shape ``(n_points,)``
y: ndarray of floats with shape ``(n_points,)``
out_of_bounds: str, default: "warn"
    What to do when points are located outside of any feature:

    * raise: raise a ValueError.
    * ignore: return ``fill_value`` for the out of bounds points.
    * warn: give a warning and return NaN for the out of bounds points.
    * drop: drop the out of bounds points. They may be identified
      via the ``index`` coordinate of the returned selection.
fill_value: scalar, DataArray, Dataset, or callable, optional, default: np.nan
    Value to assign to out-of-bounds points if out_of_bounds is warn
    or ignore. Forwarded to xarray's ``.where()`` method.

Returns
-------
points: Union[xr.DataArray, xr.Dataset]

xugrid.UgridDataArrayAccessor.set_crs
=====================================
Set the Coordinate Reference System (CRS) of a UGRID topology.

NOTE: The underlying geometries are not transformed to this CRS. To
transform the geometries to a new CRS, use the ``to_crs`` method.

Parameters
----------
crs : pyproj.CRS, optional if `epsg` is specified
    The value can be anything accepted
    by :meth:`pyproj.CRS.from_user_input() <pyproj.crs.CRS.from_user_input>`,
    such as an authority string (eg "EPSG:4326") or a WKT string.
epsg : int, optional if `crs` is specified
    EPSG code specifying the projection.
allow_override : bool, default False
    If the the UGRID topology already has a CRS, allow to replace the
    existing CRS, even when both are not equal.

xugrid.UgridDataArrayAccessor.set_node_coords
=============================================
Given names of x and y coordinates of the nodes of an object, set them
as the coordinates in the grid.

Parameters
----------
node_x: str
    Name of the x coordinate of the nodes in the object.
node_y: str
    Name of the y coordinate of the nodes in the object.

xugrid.UgridDataArrayAccessor.to_crs
====================================
Transform geometries to a new coordinate reference system.
Transform all geometries in an active geometry column to a different coordinate
reference system. The ``crs`` attribute on the current Ugrid must
be set. Either ``crs`` or ``epsg`` may be specified for output.

This method will transform all points in all objects. It has no notion
of projecting the cells. All segments joining points are assumed to be
lines in the current projection, not geodesics. Objects crossing the
dateline (or other projection boundary) will have undesirable behavior.

Parameters
----------
crs : pyproj.CRS, optional if `epsg` is specified
    The value can be anything accepted by
    :meth:`pyproj.CRS.from_user_input() <pyproj.crs.CRS.from_user_input>`,
    such as an authority string (eg "EPSG:4326") or a WKT string.
epsg : int, optional if `crs` is specified
    EPSG code specifying output projection.

xugrid.UgridDataArrayAccessor.to_dataset
========================================
Convert this UgridDataArray or UgridDataset into a standard
xarray.Dataset.

The UGRID topology information is added as standard data variables.

Parameters
----------
optional_attributes: bool, default: False.
    Whether to generate the UGRID optional attributes.

Returns
-------
dataset: UgridDataset

xugrid.UgridDataArrayAccessor.to_geodataframe
=============================================
Convert data and topology of one facet (node, edge, face) of the grid
to a geopandas GeoDataFrame. This also determines the geometry type of
the geodataframe:

* node: point
* edge: line
* face: polygon

Parameters
----------
dim: str
    node, edge, or face dimension. Inferred for DataArray.
name: str
    Name to give to the array (required if unnamed).
dim_order:
    Hierarchical dimension order for the resulting dataframe. Array content is
    transposed to this order and then written out as flat vectors in contiguous
    order, so the last dimension in this list will be contiguous in the resulting
    DataFrame. This has a major influence on which operations are efficient on the
    resulting dataframe.

    If provided, must include all dimensions of this DataArray. By default,
    dimensions are sorted according to the DataArray dimensions order.

Returns
-------
geodataframe: gpd.GeoDataFrame

xugrid.UgridDataArrayAccessor.to_nonperiodic
============================================
Convert this grid from a periodic grid (where the rightmost boundary shares its
nodes with the leftmost boundary) to an aperiodic grid, where the leftmost nodes
are separate from the rightmost nodes.

Parameters
----------
xmax: float
    The x-value of the newly created rightmost boundary nodes.

Returns
-------
nonperiodic: UgridDataArray

xugrid.UgridDataArrayAccessor.to_periodic
=========================================
Convert this grid to a periodic grid, where the rightmost boundary
shares its nodes with the leftmost boundary.

Returns
-------
periodic: UgridDataArray

xugrid.UgridDataArrayAccessor.topology
======================================
Mapping from name to UGRID topology.

xugrid.UgridDataArrayAccessor.total_bounds
==========================================
Returns a tuple containing ``minx, miny, maxx, maxy`` values of the grid's
node coordinates.

xugrid.snap_nodes
=================
Snap neigbhoring vertices together that are located within a maximum
snapping distance from each other.

If vertices are located within a maximum distance, some of them are snapped
to their neighbors ("targets"), thereby guaranteeing a minimum distance
between nodes in the result. The determination of whether a point becomes a
target itself or gets snapped to another point is primarily based on the
order in which points are processed and their spatial relationships.

This function also return an inverse index array. In case of a connectivity
array, ``inverse`` can be used to index into, yielding the updated
numbers. E.g.:

``updated_face_nodes = inverse[face_nodes]``

Parameters
----------
x: 1D nd array of floats of size N
y: 1D nd array of floats of size N
max_snap_distance: float

Returns
-------
inverse: 1D nd array of ints of size N
    Inverse index array: the new vertex number for every old vertex. Is
    None when no vertices within max_distance of each other.
x_snapped: 1D nd array of floats of size M
    Returns a copy of ``x`` when no vertices within max_distance of each
    other.
y_snapped: 1D nd array of floats of size M
    Returns a copy of ``y`` when no vertices within max_distance of each
    other.

xugrid.snap_to_grid
===================
Snap a collection of lines to a grid.

A line is included and snapped to a grid edge when the line separates
the centroid of the cell with the centroid of the edge.

Parameters
----------
lines: gpd.GeoDataFrame
    Line data. Geometry colum should contain exclusively LineStrings.
grid: xr.DataArray or xu.UgridDataArray of integers
    Grid of cells to snap lines to.
max_snap_distance: float

Returns
-------
uds: UgridDataset
    Snapped line geometries as edges in a Ugrid2d topology. Contains a
    ``line_index`` variable identifying the original geodataframe line.
gdf: gpd.GeoDataFrame
    Snapped line geometries.

xugrid.create_snap_to_grid_dataframe
====================================
Create a dataframe required to snap line geometries to a Ugrid2d topology.

A line is included and snapped to a grid edge when the line separates
the centroid of the cell with the centroid of the edge.

Parameters
----------
lines: gpd.GeoDataFrame
    Line data. Geometry colum should contain exclusively LineStrings.
grid: xugrid.Ugrid2d
    Grid of cells to snap lines to.
max_snap_distance: float

Returns
-------
result: pd.DataFrame
    DataFrame with columns:

    * line_index: the index of the geodataframe geometry.
    * edge_index: the index of the
    * x0: start x-coordinate of edge segment.
    * y0: start y-coordinate of edge segment.
    * x1: end x-coordinate of edge segment.
    * y1: end y-coordinate of edge segment.
    * length: length of the edge.

Examples
--------
First create data frame:

>>> snapping_df = create_snap_to_grid_dataframe(lines, grid2d, max_snap_distance=0.5)

Use the ``line_index`` column to assign values from ``lines`` to this new dataframe:

>>> snapping_df["my_variable"] = lines["my_variable"].iloc[snapping_df["line_index"]].to_numpy()

Run some reduction on the variable, to create an aggregated value per grid edge:

>>> aggregated = snapping_df.groupby("edge_index").sum()

Assign the aggregated values to a Ugrid2d topology:

>>> new = xu.full_like(edge_data, np.nan)
>>> new.data[aggregated.index] = aggregated["my_variable"]

xugrid.BarycentricInterpolator
==============================
The BaryCentricInterpolator searches the centroid of every face of the
target grid in the source grid. It finds by which source faces the centroid
is surrounded (via its centroidal voronoi tesselation), and computes
barycentric weights which can be used for to interpolate smoothly between
the values associated with the source faces.

Parameters
----------
source: Ugrid2d, UgridDataArray
target: Ugrid2d, UgridDataArray

xugrid.BarycentricInterpolator Class Members
============================================
   * xugrid.BarycentricInterpolator.from_dataset
   * xugrid.BarycentricInterpolator.regrid
   * xugrid.BarycentricInterpolator.to_dataset
   * xugrid.BarycentricInterpolator.weights_as_dataframe

xugrid.BarycentricInterpolator.from_dataset
===========================================
Reconstruct the regridder from a dataset with source, target indices
and weights.

xugrid.BarycentricInterpolator.regrid
=====================================
Regrid the data from a DataArray from its old grid topology to the new
target topology.

Automatically regrids over additional dimensions (e.g. time).

Supports lazy evaluation for dask arrays inside the DataArray.

Parameters
----------
object: UgridDataArray or xarray.DataArray

Returns
-------
regridded: UgridDataArray or xarray.DataArray

xugrid.BarycentricInterpolator.to_dataset
=========================================
Store the computed weights and target in a dataset for re-use.

xugrid.BarycentricInterpolator.weights_as_dataframe
===================================================
Return the weights as a three column dataframe:

* source index
* target index
* weight

Returns
-------
weights: pd.DataFrame

xugrid.CentroidLocatorRegridder
===============================
The CentroidLocatorRegridded regrids by searching the source grid for the
centroids of the target grid.

If a centroid is exactly located on an edge between two faces, the value of
either face may be used.

Parameters
----------
source: Ugrid2d, UgridDataArray
target: Ugrid2d, UgridDataArray
weights: Optional[MatrixCOO]

xugrid.CentroidLocatorRegridder Class Members
=============================================
   * xugrid.CentroidLocatorRegridder.from_dataset
   * xugrid.CentroidLocatorRegridder.regrid
   * xugrid.CentroidLocatorRegridder.to_dataset
   * xugrid.CentroidLocatorRegridder.weights_as_dataframe

xugrid.CentroidLocatorRegridder.from_dataset
============================================
Reconstruct the regridder from a dataset with source, target indices
and weights.

xugrid.OverlapRegridder
=======================
The OverlapRegridder regrids by computing which target faces overlap with
which source faces. It stores the area of overlap, which can be used in
multiple ways to aggregate the values associated with the source faces.

Currently supported aggregation methods are:

* ``"mean"``
* ``"harmonic_mean"``
* ``"geometric_mean"``
* ``"sum"``
* ``"minimum"``
* ``"maximum"``
* ``"mode"``
* ``"median"``
* ``"max_overlap"``
* percentiles 5, 10, 25, 50, 75, 90, 95: as ``"p5"``, ``"p10"``, etc.

Custom aggregation functions are also supported, if they can be compiled by
Numba. See the User Guide.

Any percentile method can be created via:
``method = OverlapRegridder.create_percentile_methode(percentile)``
See the examples.

Parameters
----------
source: Ugrid2d, UgridDataArray
target: Ugrid2d, UgridDataArray
method: str, function, optional
    Default value is ``"mean"``.

Examples
--------
Create an OverlapRegridder to regrid with mean:

>>> regridder = OverlapRegridder(source_grid, target_grid, method="mean")
>>> regridder.regrid(source_data)

Setup a custom percentile method and apply it:

>>> p33_3 = OverlapRegridder.create_percentile_method(33.3)
>>> regridder = OverlapRegridder(source_grid, target_grid, method=p33_3)
>>> regridder.regrid(source_data)

xugrid.OverlapRegridder Class Members
=====================================
   * xugrid.OverlapRegridder.from_dataset
   * xugrid.OverlapRegridder.regrid
   * xugrid.OverlapRegridder.to_dataset
   * xugrid.OverlapRegridder.weights_as_dataframe

xugrid.OverlapRegridder.from_dataset
====================================
Reconstruct the regridder from a dataset with source, target indices
and weights.

xugrid.RelativeOverlapRegridder
===============================
The RelativeOverlapRegridder regrids by computing which target faces
overlap with which source faces. It stores the area of overlap, which can
be used in multiple ways to aggregate the values associated with the source
faces. Unlike the OverlapRegridder, the intersection area is divided by the
total area of the source face. This is required for e.g. first-order
conserative regridding.

Currently supported aggregation methods are:

* ``"max_overlap"``

Custom aggregation functions are also supported, if they can be compiled by
Numba. See the User Guide.

Parameters
----------
source: Ugrid2d, UgridDataArray
target: Ugrid2d, UgridDataArray
method: str, function, optional
    Default value is "first_order_conservative".

xugrid.RelativeOverlapRegridder Class Members
=============================================
   * xugrid.RelativeOverlapRegridder.from_dataset
   * xugrid.RelativeOverlapRegridder.regrid
   * xugrid.RelativeOverlapRegridder.to_dataset
   * xugrid.RelativeOverlapRegridder.weights_as_dataframe

xugrid.RelativeOverlapRegridder.from_dataset
============================================
Reconstruct the regridder from a dataset with source, target indices
and weights.

xugrid.plot.contour
===================
Create a contour plot of a 2D UgridDataArray.

Wraps :py:func:`matplotlib:matplotlib.pyplot.tricontour`.


Parameters
----------
topology : Union[Ugrid1d, Ugrid2d, Tuple[FloatArray, FloatArray], Triangulation]
    Mesh topology.
darray : DataArray
    Must be two-dimensional, unless creating faceted plots.
figsize : tuple, optional
    A tuple (width, height) of the figure in inches.
    Mutually exclusive with ``size`` and ``ax``.
aspect : scalar, optional
    Aspect ratio of plot, so that ``aspect * size`` gives the *width* in
    inches. Only used if a ``size`` is provided.
size : scalar, optional
    If provided, create a new figure for the plot with the given size:
    *height* (in inches) of each plot. See also: ``aspect``.
ax : matplotlib axes object, optional
    Axes on which to plot. By default, use the current axes.
    Mutually exclusive with ``size`` and ``figsize``.
row : string, optional
    If passed, make row faceted plots on this dimension name.
col : string, optional
    If passed, make column faceted plots on this dimension name.
col_wrap : int, optional
    Use together with ``col`` to wrap faceted plots.
xticks, yticks : array-like, optional
    Specify tick locations for *x*- and *y*-axis.
xlim, ylim : array-like, optional
    Specify *x*- and *y*-axis limits.
xincrease : None, True, or False, optional
    Should the values on the *x* axis be increasing from left to right?
    If ``None``, use the default for the Matplotlib function.
yincrease : None, True, or False, optional
    Should the values on the *y* axis be increasing from top to bottom?
    If ``None``, use the default for the Matplotlib function.
add_colorbar : bool, optional
    Add colorbar to axes.
add_labels : bool, optional
    Use xarray metadata to label axes.
norm : matplotlib.colors.Normalize, optional
    If ``norm`` has ``vmin`` or ``vmax`` specified, the corresponding
    kwarg must be ``None``.
vmin, vmax : float, optional
    Values to anchor the colormap, otherwise they are inferred from the
    data and other keyword arguments. When a diverging dataset is inferred,
    setting one of these values will fix the other by symmetry around
    ``center``. Setting both values prevents use of a diverging colormap.
    If discrete levels are provided as an explicit list, both of these
    values are ignored.
cmap : matplotlib colormap name or colormap, optional
    The mapping from data values to color space. If not provided, this
    will be either be ``'viridis'`` (if the function infers a sequential
    dataset) or ``'RdBu_r'`` (if the function infers a diverging dataset).
    See :doc:`Choosing Colormaps in Matplotlib <matplotlib:tutorials/colors/colormaps>`
    for more information.
    If *seaborn* is installed, ``cmap`` may also be a
    `seaborn color palette <https://seaborn.pydata.org/tutorial/color_palettes.html>`_.
    Note: if ``cmap`` is a seaborn color palette and the plot type
    is not ``'contour'`` or ``'contourf'``, ``levels`` must also be specified.
colors : str or array-like of color-like, optional
    A single color or a sequence of colors. If the plot type is not ``'contour'``
    or ``'contourf'``, the ``levels`` argument is required.
center : float, optional
    The value at which to center the colormap. Passing this value implies
    use of a diverging colormap. Setting it to ``False`` prevents use of a
    diverging colormap.
robust : bool, optional
    If ``True`` and ``vmin`` or ``vmax`` are absent, the colormap range is
    computed with 2nd and 98th percentiles instead of the extreme values.
extend : {'neither', 'both', 'min', 'max'}, optional
    How to draw arrows extending the colorbar beyond its limits. If not
    provided, ``extend`` is inferred from ``vmin``, ``vmax`` and the data limits.
levels : int or array-like, optional
    Split the colormap (``cmap``) into discrete color intervals. If an integer
    is provided, "nice" levels are chosen based on the data range: this can
    imply that the final number of levels is not exactly the expected one.
    Setting ``vmin`` and/or ``vmax`` with ``levels=N`` is equivalent to
    setting ``levels=np.linspace(vmin, vmax, N)``.
infer_intervals : bool, optional
    Only applies to pcolormesh. If ``True``, the coordinate intervals are
    passed to pcolormesh. If ``False``, the original coordinates are used
    (this can be useful for certain map projections). The default is to
    always infer intervals, unless the mesh is irregular and plotted on
    a map projection.
subplot_kws : dict, optional
    Dictionary of keyword arguments for Matplotlib subplots. Only used
    for 2D and faceted plots.
    (see :py:meth:`matplotlib:matplotlib.figure.Figure.add_subplot`).
cbar_ax : matplotlib axes object, optional
    Axes in which to draw the colorbar.
cbar_kwargs : dict, optional
    Dictionary of keyword arguments to pass to the colorbar
    (see :meth:`matplotlib:matplotlib.figure.Figure.colorbar`).
**kwargs : optional
    Additional keyword arguments to wrapped Matplotlib function.
Returns
-------
artist :
    The same type of primitive artist that the wrapped Matplotlib
    function returns.

xugrid.plot.contourf
====================
Create a filled contour plot of a 2D UgridDataArray.

Wraps :py:func:`matplotlib:matplotlib.pyplot.tricontourf`.


Parameters
----------
topology : Union[Ugrid1d, Ugrid2d, Tuple[FloatArray, FloatArray], Triangulation]
    Mesh topology.
darray : DataArray
    Must be two-dimensional, unless creating faceted plots.
figsize : tuple, optional
    A tuple (width, height) of the figure in inches.
    Mutually exclusive with ``size`` and ``ax``.
aspect : scalar, optional
    Aspect ratio of plot, so that ``aspect * size`` gives the *width* in
    inches. Only used if a ``size`` is provided.
size : scalar, optional
    If provided, create a new figure for the plot with the given size:
    *height* (in inches) of each plot. See also: ``aspect``.
ax : matplotlib axes object, optional
    Axes on which to plot. By default, use the current axes.
    Mutually exclusive with ``size`` and ``figsize``.
row : string, optional
    If passed, make row faceted plots on this dimension name.
col : string, optional
    If passed, make column faceted plots on this dimension name.
col_wrap : int, optional
    Use together with ``col`` to wrap faceted plots.
xticks, yticks : array-like, optional
    Specify tick locations for *x*- and *y*-axis.
xlim, ylim : array-like, optional
    Specify *x*- and *y*-axis limits.
xincrease : None, True, or False, optional
    Should the values on the *x* axis be increasing from left to right?
    If ``None``, use the default for the Matplotlib function.
yincrease : None, True, or False, optional
    Should the values on the *y* axis be increasing from top to bottom?
    If ``None``, use the default for the Matplotlib function.
add_colorbar : bool, optional
    Add colorbar to axes.
add_labels : bool, optional
    Use xarray metadata to label axes.
norm : matplotlib.colors.Normalize, optional
    If ``norm`` has ``vmin`` or ``vmax`` specified, the corresponding
    kwarg must be ``None``.
vmin, vmax : float, optional
    Values to anchor the colormap, otherwise they are inferred from the
    data and other keyword arguments. When a diverging dataset is inferred,
    setting one of these values will fix the other by symmetry around
    ``center``. Setting both values prevents use of a diverging colormap.
    If discrete levels are provided as an explicit list, both of these
    values are ignored.
cmap : matplotlib colormap name or colormap, optional
    The mapping from data values to color space. If not provided, this
    will be either be ``'viridis'`` (if the function infers a sequential
    dataset) or ``'RdBu_r'`` (if the function infers a diverging dataset).
    See :doc:`Choosing Colormaps in Matplotlib <matplotlib:tutorials/colors/colormaps>`
    for more information.
    If *seaborn* is installed, ``cmap`` may also be a
    `seaborn color palette <https://seaborn.pydata.org/tutorial/color_palettes.html>`_.
    Note: if ``cmap`` is a seaborn color palette and the plot type
    is not ``'contour'`` or ``'contourf'``, ``levels`` must also be specified.
colors : str or array-like of color-like, optional
    A single color or a sequence of colors. If the plot type is not ``'contour'``
    or ``'contourf'``, the ``levels`` argument is required.
center : float, optional
    The value at which to center the colormap. Passing this value implies
    use of a diverging colormap. Setting it to ``False`` prevents use of a
    diverging colormap.
robust : bool, optional
    If ``True`` and ``vmin`` or ``vmax`` are absent, the colormap range is
    computed with 2nd and 98th percentiles instead of the extreme values.
extend : {'neither', 'both', 'min', 'max'}, optional
    How to draw arrows extending the colorbar beyond its limits. If not
    provided, ``extend`` is inferred from ``vmin``, ``vmax`` and the data limits.
levels : int or array-like, optional
    Split the colormap (``cmap``) into discrete color intervals. If an integer
    is provided, "nice" levels are chosen based on the data range: this can
    imply that the final number of levels is not exactly the expected one.
    Setting ``vmin`` and/or ``vmax`` with ``levels=N`` is equivalent to
    setting ``levels=np.linspace(vmin, vmax, N)``.
infer_intervals : bool, optional
    Only applies to pcolormesh. If ``True``, the coordinate intervals are
    passed to pcolormesh. If ``False``, the original coordinates are used
    (this can be useful for certain map projections). The default is to
    always infer intervals, unless the mesh is irregular and plotted on
    a map projection.
subplot_kws : dict, optional
    Dictionary of keyword arguments for Matplotlib subplots. Only used
    for 2D and faceted plots.
    (see :py:meth:`matplotlib:matplotlib.figure.Figure.add_subplot`).
cbar_ax : matplotlib axes object, optional
    Axes in which to draw the colorbar.
cbar_kwargs : dict, optional
    Dictionary of keyword arguments to pass to the colorbar
    (see :meth:`matplotlib:matplotlib.figure.Figure.colorbar`).
**kwargs : optional
    Additional keyword arguments to wrapped Matplotlib function.
Returns
-------
artist :
    The same type of primitive artist that the wrapped Matplotlib
    function returns.

xugrid.plot.imshow
==================
Image plot of 2D DataArray.
Wraps :py:func:`matplotlib:matplotlib.pyplot.imshow`.

This rasterizes the grid before plotting. Pass a ``resolution`` keyword to
control the rasterization resolution.


Parameters
----------
topology : Union[Ugrid1d, Ugrid2d, Tuple[FloatArray, FloatArray], Triangulation]
    Mesh topology.
darray : DataArray
    Must be two-dimensional, unless creating faceted plots.
figsize : tuple, optional
    A tuple (width, height) of the figure in inches.
    Mutually exclusive with ``size`` and ``ax``.
aspect : scalar, optional
    Aspect ratio of plot, so that ``aspect * size`` gives the *width* in
    inches. Only used if a ``size`` is provided.
size : scalar, optional
    If provided, create a new figure for the plot with the given size:
    *height* (in inches) of each plot. See also: ``aspect``.
ax : matplotlib axes object, optional
    Axes on which to plot. By default, use the current axes.
    Mutually exclusive with ``size`` and ``figsize``.
row : string, optional
    If passed, make row faceted plots on this dimension name.
col : string, optional
    If passed, make column faceted plots on this dimension name.
col_wrap : int, optional
    Use together with ``col`` to wrap faceted plots.
xticks, yticks : array-like, optional
    Specify tick locations for *x*- and *y*-axis.
xlim, ylim : array-like, optional
    Specify *x*- and *y*-axis limits.
xincrease : None, True, or False, optional
    Should the values on the *x* axis be increasing from left to right?
    If ``None``, use the default for the Matplotlib function.
yincrease : None, True, or False, optional
    Should the values on the *y* axis be increasing from top to bottom?
    If ``None``, use the default for the Matplotlib function.
add_colorbar : bool, optional
    Add colorbar to axes.
add_labels : bool, optional
    Use xarray metadata to label axes.
norm : matplotlib.colors.Normalize, optional
    If ``norm`` has ``vmin`` or ``vmax`` specified, the corresponding
    kwarg must be ``None``.
vmin, vmax : float, optional
    Values to anchor the colormap, otherwise they are inferred from the
    data and other keyword arguments. When a diverging dataset is inferred,
    setting one of these values will fix the other by symmetry around
    ``center``. Setting both values prevents use of a diverging colormap.
    If discrete levels are provided as an explicit list, both of these
    values are ignored.
cmap : matplotlib colormap name or colormap, optional
    The mapping from data values to color space. If not provided, this
    will be either be ``'viridis'`` (if the function infers a sequential
    dataset) or ``'RdBu_r'`` (if the function infers a diverging dataset).
    See :doc:`Choosing Colormaps in Matplotlib <matplotlib:tutorials/colors/colormaps>`
    for more information.
    If *seaborn* is installed, ``cmap`` may also be a
    `seaborn color palette <https://seaborn.pydata.org/tutorial/color_palettes.html>`_.
    Note: if ``cmap`` is a seaborn color palette and the plot type
    is not ``'contour'`` or ``'contourf'``, ``levels`` must also be specified.
colors : str or array-like of color-like, optional
    A single color or a sequence of colors. If the plot type is not ``'contour'``
    or ``'contourf'``, the ``levels`` argument is required.
center : float, optional
    The value at which to center the colormap. Passing this value implies
    use of a diverging colormap. Setting it to ``False`` prevents use of a
    diverging colormap.
robust : bool, optional
    If ``True`` and ``vmin`` or ``vmax`` are absent, the colormap range is
    computed with 2nd and 98th percentiles instead of the extreme values.
extend : {'neither', 'both', 'min', 'max'}, optional
    How to draw arrows extending the colorbar beyond its limits. If not
    provided, ``extend`` is inferred from ``vmin``, ``vmax`` and the data limits.
levels : int or array-like, optional
    Split the colormap (``cmap``) into discrete color intervals. If an integer
    is provided, "nice" levels are chosen based on the data range: this can
    imply that the final number of levels is not exactly the expected one.
    Setting ``vmin`` and/or ``vmax`` with ``levels=N`` is equivalent to
    setting ``levels=np.linspace(vmin, vmax, N)``.
infer_intervals : bool, optional
    Only applies to pcolormesh. If ``True``, the coordinate intervals are
    passed to pcolormesh. If ``False``, the original coordinates are used
    (this can be useful for certain map projections). The default is to
    always infer intervals, unless the mesh is irregular and plotted on
    a map projection.
subplot_kws : dict, optional
    Dictionary of keyword arguments for Matplotlib subplots. Only used
    for 2D and faceted plots.
    (see :py:meth:`matplotlib:matplotlib.figure.Figure.add_subplot`).
cbar_ax : matplotlib axes object, optional
    Axes in which to draw the colorbar.
cbar_kwargs : dict, optional
    Dictionary of keyword arguments to pass to the colorbar
    (see :meth:`matplotlib:matplotlib.figure.Figure.colorbar`).
**kwargs : optional
    Additional keyword arguments to wrapped Matplotlib function.
Returns
-------
artist :
    The same type of primitive artist that the wrapped Matplotlib
    function returns.

xugrid.plot.line
================
None

Parameters
----------
topology : Union[Ugrid1d, Ugrid2d, Tuple[FloatArray, FloatArray], Triangulation]
    Mesh topology.
darray : DataArray
    Must be two-dimensional, unless creating faceted plots.
figsize : tuple, optional
    A tuple (width, height) of the figure in inches.
    Mutually exclusive with ``size`` and ``ax``.
aspect : scalar, optional
    Aspect ratio of plot, so that ``aspect * size`` gives the *width* in
    inches. Only used if a ``size`` is provided.
size : scalar, optional
    If provided, create a new figure for the plot with the given size:
    *height* (in inches) of each plot. See also: ``aspect``.
ax : matplotlib axes object, optional
    Axes on which to plot. By default, use the current axes.
    Mutually exclusive with ``size`` and ``figsize``.
row : string, optional
    If passed, make row faceted plots on this dimension name.
col : string, optional
    If passed, make column faceted plots on this dimension name.
col_wrap : int, optional
    Use together with ``col`` to wrap faceted plots.
xticks, yticks : array-like, optional
    Specify tick locations for *x*- and *y*-axis.
xlim, ylim : array-like, optional
    Specify *x*- and *y*-axis limits.
xincrease : None, True, or False, optional
    Should the values on the *x* axis be increasing from left to right?
    If ``None``, use the default for the Matplotlib function.
yincrease : None, True, or False, optional
    Should the values on the *y* axis be increasing from top to bottom?
    If ``None``, use the default for the Matplotlib function.
add_colorbar : bool, optional
    Add colorbar to axes.
add_labels : bool, optional
    Use xarray metadata to label axes.
norm : matplotlib.colors.Normalize, optional
    If ``norm`` has ``vmin`` or ``vmax`` specified, the corresponding
    kwarg must be ``None``.
vmin, vmax : float, optional
    Values to anchor the colormap, otherwise they are inferred from the
    data and other keyword arguments. When a diverging dataset is inferred,
    setting one of these values will fix the other by symmetry around
    ``center``. Setting both values prevents use of a diverging colormap.
    If discrete levels are provided as an explicit list, both of these
    values are ignored.
cmap : matplotlib colormap name or colormap, optional
    The mapping from data values to color space. If not provided, this
    will be either be ``'viridis'`` (if the function infers a sequential
    dataset) or ``'RdBu_r'`` (if the function infers a diverging dataset).
    See :doc:`Choosing Colormaps in Matplotlib <matplotlib:tutorials/colors/colormaps>`
    for more information.
    If *seaborn* is installed, ``cmap`` may also be a
    `seaborn color palette <https://seaborn.pydata.org/tutorial/color_palettes.html>`_.
    Note: if ``cmap`` is a seaborn color palette and the plot type
    is not ``'contour'`` or ``'contourf'``, ``levels`` must also be specified.
colors : str or array-like of color-like, optional
    A single color or a sequence of colors. If the plot type is not ``'contour'``
    or ``'contourf'``, the ``levels`` argument is required.
center : float, optional
    The value at which to center the colormap. Passing this value implies
    use of a diverging colormap. Setting it to ``False`` prevents use of a
    diverging colormap.
robust : bool, optional
    If ``True`` and ``vmin`` or ``vmax`` are absent, the colormap range is
    computed with 2nd and 98th percentiles instead of the extreme values.
extend : {'neither', 'both', 'min', 'max'}, optional
    How to draw arrows extending the colorbar beyond its limits. If not
    provided, ``extend`` is inferred from ``vmin``, ``vmax`` and the data limits.
levels : int or array-like, optional
    Split the colormap (``cmap``) into discrete color intervals. If an integer
    is provided, "nice" levels are chosen based on the data range: this can
    imply that the final number of levels is not exactly the expected one.
    Setting ``vmin`` and/or ``vmax`` with ``levels=N`` is equivalent to
    setting ``levels=np.linspace(vmin, vmax, N)``.
infer_intervals : bool, optional
    Only applies to pcolormesh. If ``True``, the coordinate intervals are
    passed to pcolormesh. If ``False``, the original coordinates are used
    (this can be useful for certain map projections). The default is to
    always infer intervals, unless the mesh is irregular and plotted on
    a map projection.
subplot_kws : dict, optional
    Dictionary of keyword arguments for Matplotlib subplots. Only used
    for 2D and faceted plots.
    (see :py:meth:`matplotlib:matplotlib.figure.Figure.add_subplot`).
cbar_ax : matplotlib axes object, optional
    Axes in which to draw the colorbar.
cbar_kwargs : dict, optional
    Dictionary of keyword arguments to pass to the colorbar
    (see :meth:`matplotlib:matplotlib.figure.Figure.colorbar`).
**kwargs : optional
    Additional keyword arguments to wrapped Matplotlib function.
Returns
-------
artist :
    The same type of primitive artist that the wrapped Matplotlib
    function returns.

xugrid.plot.pcolormesh
======================
Create a pseudocolor mesh plot of a 2D UgridDataArray.

Wraps matplotlib PolyCollection.


Parameters
----------
topology : Union[Ugrid1d, Ugrid2d, Tuple[FloatArray, FloatArray], Triangulation]
    Mesh topology.
darray : DataArray
    Must be two-dimensional, unless creating faceted plots.
figsize : tuple, optional
    A tuple (width, height) of the figure in inches.
    Mutually exclusive with ``size`` and ``ax``.
aspect : scalar, optional
    Aspect ratio of plot, so that ``aspect * size`` gives the *width* in
    inches. Only used if a ``size`` is provided.
size : scalar, optional
    If provided, create a new figure for the plot with the given size:
    *height* (in inches) of each plot. See also: ``aspect``.
ax : matplotlib axes object, optional
    Axes on which to plot. By default, use the current axes.
    Mutually exclusive with ``size`` and ``figsize``.
row : string, optional
    If passed, make row faceted plots on this dimension name.
col : string, optional
    If passed, make column faceted plots on this dimension name.
col_wrap : int, optional
    Use together with ``col`` to wrap faceted plots.
xticks, yticks : array-like, optional
    Specify tick locations for *x*- and *y*-axis.
xlim, ylim : array-like, optional
    Specify *x*- and *y*-axis limits.
xincrease : None, True, or False, optional
    Should the values on the *x* axis be increasing from left to right?
    If ``None``, use the default for the Matplotlib function.
yincrease : None, True, or False, optional
    Should the values on the *y* axis be increasing from top to bottom?
    If ``None``, use the default for the Matplotlib function.
add_colorbar : bool, optional
    Add colorbar to axes.
add_labels : bool, optional
    Use xarray metadata to label axes.
norm : matplotlib.colors.Normalize, optional
    If ``norm`` has ``vmin`` or ``vmax`` specified, the corresponding
    kwarg must be ``None``.
vmin, vmax : float, optional
    Values to anchor the colormap, otherwise they are inferred from the
    data and other keyword arguments. When a diverging dataset is inferred,
    setting one of these values will fix the other by symmetry around
    ``center``. Setting both values prevents use of a diverging colormap.
    If discrete levels are provided as an explicit list, both of these
    values are ignored.
cmap : matplotlib colormap name or colormap, optional
    The mapping from data values to color space. If not provided, this
    will be either be ``'viridis'`` (if the function infers a sequential
    dataset) or ``'RdBu_r'`` (if the function infers a diverging dataset).
    See :doc:`Choosing Colormaps in Matplotlib <matplotlib:tutorials/colors/colormaps>`
    for more information.
    If *seaborn* is installed, ``cmap`` may also be a
    `seaborn color palette <https://seaborn.pydata.org/tutorial/color_palettes.html>`_.
    Note: if ``cmap`` is a seaborn color palette and the plot type
    is not ``'contour'`` or ``'contourf'``, ``levels`` must also be specified.
colors : str or array-like of color-like, optional
    A single color or a sequence of colors. If the plot type is not ``'contour'``
    or ``'contourf'``, the ``levels`` argument is required.
center : float, optional
    The value at which to center the colormap. Passing this value implies
    use of a diverging colormap. Setting it to ``False`` prevents use of a
    diverging colormap.
robust : bool, optional
    If ``True`` and ``vmin`` or ``vmax`` are absent, the colormap range is
    computed with 2nd and 98th percentiles instead of the extreme values.
extend : {'neither', 'both', 'min', 'max'}, optional
    How to draw arrows extending the colorbar beyond its limits. If not
    provided, ``extend`` is inferred from ``vmin``, ``vmax`` and the data limits.
levels : int or array-like, optional
    Split the colormap (``cmap``) into discrete color intervals. If an integer
    is provided, "nice" levels are chosen based on the data range: this can
    imply that the final number of levels is not exactly the expected one.
    Setting ``vmin`` and/or ``vmax`` with ``levels=N`` is equivalent to
    setting ``levels=np.linspace(vmin, vmax, N)``.
infer_intervals : bool, optional
    Only applies to pcolormesh. If ``True``, the coordinate intervals are
    passed to pcolormesh. If ``False``, the original coordinates are used
    (this can be useful for certain map projections). The default is to
    always infer intervals, unless the mesh is irregular and plotted on
    a map projection.
subplot_kws : dict, optional
    Dictionary of keyword arguments for Matplotlib subplots. Only used
    for 2D and faceted plots.
    (see :py:meth:`matplotlib:matplotlib.figure.Figure.add_subplot`).
cbar_ax : matplotlib axes object, optional
    Axes in which to draw the colorbar.
cbar_kwargs : dict, optional
    Dictionary of keyword arguments to pass to the colorbar
    (see :meth:`matplotlib:matplotlib.figure.Figure.colorbar`).
**kwargs : optional
    Additional keyword arguments to wrapped Matplotlib function.
Returns
-------
artist :
    The same type of primitive artist that the wrapped Matplotlib
    function returns.

xugrid.plot.scatter
===================
None

Parameters
----------
topology : Union[Ugrid1d, Ugrid2d, Tuple[FloatArray, FloatArray], Triangulation]
    Mesh topology.
darray : DataArray
    Must be two-dimensional, unless creating faceted plots.
figsize : tuple, optional
    A tuple (width, height) of the figure in inches.
    Mutually exclusive with ``size`` and ``ax``.
aspect : scalar, optional
    Aspect ratio of plot, so that ``aspect * size`` gives the *width* in
    inches. Only used if a ``size`` is provided.
size : scalar, optional
    If provided, create a new figure for the plot with the given size:
    *height* (in inches) of each plot. See also: ``aspect``.
ax : matplotlib axes object, optional
    Axes on which to plot. By default, use the current axes.
    Mutually exclusive with ``size`` and ``figsize``.
row : string, optional
    If passed, make row faceted plots on this dimension name.
col : string, optional
    If passed, make column faceted plots on this dimension name.
col_wrap : int, optional
    Use together with ``col`` to wrap faceted plots.
xticks, yticks : array-like, optional
    Specify tick locations for *x*- and *y*-axis.
xlim, ylim : array-like, optional
    Specify *x*- and *y*-axis limits.
xincrease : None, True, or False, optional
    Should the values on the *x* axis be increasing from left to right?
    If ``None``, use the default for the Matplotlib function.
yincrease : None, True, or False, optional
    Should the values on the *y* axis be increasing from top to bottom?
    If ``None``, use the default for the Matplotlib function.
add_colorbar : bool, optional
    Add colorbar to axes.
add_labels : bool, optional
    Use xarray metadata to label axes.
norm : matplotlib.colors.Normalize, optional
    If ``norm`` has ``vmin`` or ``vmax`` specified, the corresponding
    kwarg must be ``None``.
vmin, vmax : float, optional
    Values to anchor the colormap, otherwise they are inferred from the
    data and other keyword arguments. When a diverging dataset is inferred,
    setting one of these values will fix the other by symmetry around
    ``center``. Setting both values prevents use of a diverging colormap.
    If discrete levels are provided as an explicit list, both of these
    values are ignored.
cmap : matplotlib colormap name or colormap, optional
    The mapping from data values to color space. If not provided, this
    will be either be ``'viridis'`` (if the function infers a sequential
    dataset) or ``'RdBu_r'`` (if the function infers a diverging dataset).
    See :doc:`Choosing Colormaps in Matplotlib <matplotlib:tutorials/colors/colormaps>`
    for more information.
    If *seaborn* is installed, ``cmap`` may also be a
    `seaborn color palette <https://seaborn.pydata.org/tutorial/color_palettes.html>`_.
    Note: if ``cmap`` is a seaborn color palette and the plot type
    is not ``'contour'`` or ``'contourf'``, ``levels`` must also be specified.
colors : str or array-like of color-like, optional
    A single color or a sequence of colors. If the plot type is not ``'contour'``
    or ``'contourf'``, the ``levels`` argument is required.
center : float, optional
    The value at which to center the colormap. Passing this value implies
    use of a diverging colormap. Setting it to ``False`` prevents use of a
    diverging colormap.
robust : bool, optional
    If ``True`` and ``vmin`` or ``vmax`` are absent, the colormap range is
    computed with 2nd and 98th percentiles instead of the extreme values.
extend : {'neither', 'both', 'min', 'max'}, optional
    How to draw arrows extending the colorbar beyond its limits. If not
    provided, ``extend`` is inferred from ``vmin``, ``vmax`` and the data limits.
levels : int or array-like, optional
    Split the colormap (``cmap``) into discrete color intervals. If an integer
    is provided, "nice" levels are chosen based on the data range: this can
    imply that the final number of levels is not exactly the expected one.
    Setting ``vmin`` and/or ``vmax`` with ``levels=N`` is equivalent to
    setting ``levels=np.linspace(vmin, vmax, N)``.
infer_intervals : bool, optional
    Only applies to pcolormesh. If ``True``, the coordinate intervals are
    passed to pcolormesh. If ``False``, the original coordinates are used
    (this can be useful for certain map projections). The default is to
    always infer intervals, unless the mesh is irregular and plotted on
    a map projection.
subplot_kws : dict, optional
    Dictionary of keyword arguments for Matplotlib subplots. Only used
    for 2D and faceted plots.
    (see :py:meth:`matplotlib:matplotlib.figure.Figure.add_subplot`).
cbar_ax : matplotlib axes object, optional
    Axes in which to draw the colorbar.
cbar_kwargs : dict, optional
    Dictionary of keyword arguments to pass to the colorbar
    (see :meth:`matplotlib:matplotlib.figure.Figure.colorbar`).
**kwargs : optional
    Additional keyword arguments to wrapped Matplotlib function.
Returns
-------
artist :
    The same type of primitive artist that the wrapped Matplotlib
    function returns.

xugrid.plot.surface
===================
Create a surface plot of a 2D UgridDataArray.

Wraps :py:func:`matplotlib:mplot3d:plot_trisurf`.


Parameters
----------
topology : Union[Ugrid1d, Ugrid2d, Tuple[FloatArray, FloatArray], Triangulation]
    Mesh topology.
darray : DataArray
    Must be two-dimensional, unless creating faceted plots.
figsize : tuple, optional
    A tuple (width, height) of the figure in inches.
    Mutually exclusive with ``size`` and ``ax``.
aspect : scalar, optional
    Aspect ratio of plot, so that ``aspect * size`` gives the *width* in
    inches. Only used if a ``size`` is provided.
size : scalar, optional
    If provided, create a new figure for the plot with the given size:
    *height* (in inches) of each plot. See also: ``aspect``.
ax : matplotlib axes object, optional
    Axes on which to plot. By default, use the current axes.
    Mutually exclusive with ``size`` and ``figsize``.
row : string, optional
    If passed, make row faceted plots on this dimension name.
col : string, optional
    If passed, make column faceted plots on this dimension name.
col_wrap : int, optional
    Use together with ``col`` to wrap faceted plots.
xticks, yticks : array-like, optional
    Specify tick locations for *x*- and *y*-axis.
xlim, ylim : array-like, optional
    Specify *x*- and *y*-axis limits.
xincrease : None, True, or False, optional
    Should the values on the *x* axis be increasing from left to right?
    If ``None``, use the default for the Matplotlib function.
yincrease : None, True, or False, optional
    Should the values on the *y* axis be increasing from top to bottom?
    If ``None``, use the default for the Matplotlib function.
add_colorbar : bool, optional
    Add colorbar to axes.
add_labels : bool, optional
    Use xarray metadata to label axes.
norm : matplotlib.colors.Normalize, optional
    If ``norm`` has ``vmin`` or ``vmax`` specified, the corresponding
    kwarg must be ``None``.
vmin, vmax : float, optional
    Values to anchor the colormap, otherwise they are inferred from the
    data and other keyword arguments. When a diverging dataset is inferred,
    setting one of these values will fix the other by symmetry around
    ``center``. Setting both values prevents use of a diverging colormap.
    If discrete levels are provided as an explicit list, both of these
    values are ignored.
cmap : matplotlib colormap name or colormap, optional
    The mapping from data values to color space. If not provided, this
    will be either be ``'viridis'`` (if the function infers a sequential
    dataset) or ``'RdBu_r'`` (if the function infers a diverging dataset).
    See :doc:`Choosing Colormaps in Matplotlib <matplotlib:tutorials/colors/colormaps>`
    for more information.
    If *seaborn* is installed, ``cmap`` may also be a
    `seaborn color palette <https://seaborn.pydata.org/tutorial/color_palettes.html>`_.
    Note: if ``cmap`` is a seaborn color palette and the plot type
    is not ``'contour'`` or ``'contourf'``, ``levels`` must also be specified.
colors : str or array-like of color-like, optional
    A single color or a sequence of colors. If the plot type is not ``'contour'``
    or ``'contourf'``, the ``levels`` argument is required.
center : float, optional
    The value at which to center the colormap. Passing this value implies
    use of a diverging colormap. Setting it to ``False`` prevents use of a
    diverging colormap.
robust : bool, optional
    If ``True`` and ``vmin`` or ``vmax`` are absent, the colormap range is
    computed with 2nd and 98th percentiles instead of the extreme values.
extend : {'neither', 'both', 'min', 'max'}, optional
    How to draw arrows extending the colorbar beyond its limits. If not
    provided, ``extend`` is inferred from ``vmin``, ``vmax`` and the data limits.
levels : int or array-like, optional
    Split the colormap (``cmap``) into discrete color intervals. If an integer
    is provided, "nice" levels are chosen based on the data range: this can
    imply that the final number of levels is not exactly the expected one.
    Setting ``vmin`` and/or ``vmax`` with ``levels=N`` is equivalent to
    setting ``levels=np.linspace(vmin, vmax, N)``.
infer_intervals : bool, optional
    Only applies to pcolormesh. If ``True``, the coordinate intervals are
    passed to pcolormesh. If ``False``, the original coordinates are used
    (this can be useful for certain map projections). The default is to
    always infer intervals, unless the mesh is irregular and plotted on
    a map projection.
subplot_kws : dict, optional
    Dictionary of keyword arguments for Matplotlib subplots. Only used
    for 2D and faceted plots.
    (see :py:meth:`matplotlib:matplotlib.figure.Figure.add_subplot`).
cbar_ax : matplotlib axes object, optional
    Axes in which to draw the colorbar.
cbar_kwargs : dict, optional
    Dictionary of keyword arguments to pass to the colorbar
    (see :meth:`matplotlib:matplotlib.figure.Figure.colorbar`).
**kwargs : optional
    Additional keyword arguments to wrapped Matplotlib function.
Returns
-------
artist :
    The same type of primitive artist that the wrapped Matplotlib
    function returns.

xugrid.plot.tripcolor
=====================
None

Parameters
----------
topology : Union[Ugrid1d, Ugrid2d, Tuple[FloatArray, FloatArray], Triangulation]
    Mesh topology.
darray : DataArray
    Must be two-dimensional, unless creating faceted plots.
figsize : tuple, optional
    A tuple (width, height) of the figure in inches.
    Mutually exclusive with ``size`` and ``ax``.
aspect : scalar, optional
    Aspect ratio of plot, so that ``aspect * size`` gives the *width* in
    inches. Only used if a ``size`` is provided.
size : scalar, optional
    If provided, create a new figure for the plot with the given size:
    *height* (in inches) of each plot. See also: ``aspect``.
ax : matplotlib axes object, optional
    Axes on which to plot. By default, use the current axes.
    Mutually exclusive with ``size`` and ``figsize``.
row : string, optional
    If passed, make row faceted plots on this dimension name.
col : string, optional
    If passed, make column faceted plots on this dimension name.
col_wrap : int, optional
    Use together with ``col`` to wrap faceted plots.
xticks, yticks : array-like, optional
    Specify tick locations for *x*- and *y*-axis.
xlim, ylim : array-like, optional
    Specify *x*- and *y*-axis limits.
xincrease : None, True, or False, optional
    Should the values on the *x* axis be increasing from left to right?
    If ``None``, use the default for the Matplotlib function.
yincrease : None, True, or False, optional
    Should the values on the *y* axis be increasing from top to bottom?
    If ``None``, use the default for the Matplotlib function.
add_colorbar : bool, optional
    Add colorbar to axes.
add_labels : bool, optional
    Use xarray metadata to label axes.
norm : matplotlib.colors.Normalize, optional
    If ``norm`` has ``vmin`` or ``vmax`` specified, the corresponding
    kwarg must be ``None``.
vmin, vmax : float, optional
    Values to anchor the colormap, otherwise they are inferred from the
    data and other keyword arguments. When a diverging dataset is inferred,
    setting one of these values will fix the other by symmetry around
    ``center``. Setting both values prevents use of a diverging colormap.
    If discrete levels are provided as an explicit list, both of these
    values are ignored.
cmap : matplotlib colormap name or colormap, optional
    The mapping from data values to color space. If not provided, this
    will be either be ``'viridis'`` (if the function infers a sequential
    dataset) or ``'RdBu_r'`` (if the function infers a diverging dataset).
    See :doc:`Choosing Colormaps in Matplotlib <matplotlib:tutorials/colors/colormaps>`
    for more information.
    If *seaborn* is installed, ``cmap`` may also be a
    `seaborn color palette <https://seaborn.pydata.org/tutorial/color_palettes.html>`_.
    Note: if ``cmap`` is a seaborn color palette and the plot type
    is not ``'contour'`` or ``'contourf'``, ``levels`` must also be specified.
colors : str or array-like of color-like, optional
    A single color or a sequence of colors. If the plot type is not ``'contour'``
    or ``'contourf'``, the ``levels`` argument is required.
center : float, optional
    The value at which to center the colormap. Passing this value implies
    use of a diverging colormap. Setting it to ``False`` prevents use of a
    diverging colormap.
robust : bool, optional
    If ``True`` and ``vmin`` or ``vmax`` are absent, the colormap range is
    computed with 2nd and 98th percentiles instead of the extreme values.
extend : {'neither', 'both', 'min', 'max'}, optional
    How to draw arrows extending the colorbar beyond its limits. If not
    provided, ``extend`` is inferred from ``vmin``, ``vmax`` and the data limits.
levels : int or array-like, optional
    Split the colormap (``cmap``) into discrete color intervals. If an integer
    is provided, "nice" levels are chosen based on the data range: this can
    imply that the final number of levels is not exactly the expected one.
    Setting ``vmin`` and/or ``vmax`` with ``levels=N`` is equivalent to
    setting ``levels=np.linspace(vmin, vmax, N)``.
infer_intervals : bool, optional
    Only applies to pcolormesh. If ``True``, the coordinate intervals are
    passed to pcolormesh. If ``False``, the original coordinates are used
    (this can be useful for certain map projections). The default is to
    always infer intervals, unless the mesh is irregular and plotted on
    a map projection.
subplot_kws : dict, optional
    Dictionary of keyword arguments for Matplotlib subplots. Only used
    for 2D and faceted plots.
    (see :py:meth:`matplotlib:matplotlib.figure.Figure.add_subplot`).
cbar_ax : matplotlib axes object, optional
    Axes in which to draw the colorbar.
cbar_kwargs : dict, optional
    Dictionary of keyword arguments to pass to the colorbar
    (see :meth:`matplotlib:matplotlib.figure.Figure.colorbar`).
**kwargs : optional
    Additional keyword arguments to wrapped Matplotlib function.
Returns
-------
artist :
    The same type of primitive artist that the wrapped Matplotlib
    function returns.

xugrid.Ugrid1d
==============
This class stores the topological data of a "1-D unstructured grid": a
collection of connected line elements, such as a river network.

Parameters
----------
node_x: ndarray of floats
node_y: ndarray of floats
fill_value: int
edge_node_connectivity: ndarray of integers
name: string, optional
    Network name. Defaults to "network1d".
dataset: xr.Dataset, optional
indexes: Dict[str, str], optional
    When a dataset is provided, a mapping from the UGRID role to the dataset
    variable name. E.g. {"face_x": "mesh2d_face_lon"}.
projected: bool, optional
    Whether node_x and node_y are longitude and latitude or projected x and
    y coordinates. Used to write the appropriate standard_name in the
    coordinate attributes.
crs: Any, optional
    Coordinate Reference System of the geometry objects. Can be anything accepted by
    :meth:`pyproj.CRS.from_user_input() <pyproj.crs.CRS.from_user_input>`,
    such as an authority string (eg "EPSG:4326") or a WKT string.
attrs: Dict[str, str], optional
    UGRID topology attributes. Should not be provided together with
    dataset: if other names are required, update the dataset instead.
    A name entry is ignored, as name is given explicitly.
start_index: int, 0 or 1, default is 0.
    Start index of the connectivity arrays. Must match the start index
    of the provided face_node_connectivity and edge_node_connectivity.

xugrid.Ugrid1d Class Members
============================
   * xugrid.Ugrid1d.assign_edge_coords
   * xugrid.Ugrid1d.assign_node_coords
   * xugrid.Ugrid1d.bounds
   * xugrid.Ugrid1d.contract_vertices
   * xugrid.Ugrid1d.coords
   * xugrid.Ugrid1d.copy
   * xugrid.Ugrid1d.create_data_array
   * xugrid.Ugrid1d.dimensions
   * xugrid.Ugrid1d.dims
   * xugrid.Ugrid1d.directed_node_node_connectivity
   * xugrid.Ugrid1d.edge_bounds
   * xugrid.Ugrid1d.edge_coordinates
   * xugrid.Ugrid1d.edge_dimension
   * xugrid.Ugrid1d.edge_node_coordinates
   * xugrid.Ugrid1d.edge_x
   * xugrid.Ugrid1d.edge_y
   * xugrid.Ugrid1d.fill_value
   * xugrid.Ugrid1d.find_ugrid_dim
   * xugrid.Ugrid1d.from_dataset
   * xugrid.Ugrid1d.from_geodataframe
   * xugrid.Ugrid1d.from_meshkernel
   * xugrid.Ugrid1d.from_shapely
   * xugrid.Ugrid1d.get_connectivity_matrix
   * xugrid.Ugrid1d.get_coordinates
   * xugrid.Ugrid1d.isel
   * xugrid.Ugrid1d.merge_partitions
   * xugrid.Ugrid1d.mesh
   * xugrid.Ugrid1d.meshkernel
   * xugrid.Ugrid1d.n_edge
   * xugrid.Ugrid1d.n_node
   * xugrid.Ugrid1d.node_coordinates
   * xugrid.Ugrid1d.node_dimension
   * xugrid.Ugrid1d.node_edge_connectivity
   * xugrid.Ugrid1d.node_node_connectivity
   * xugrid.Ugrid1d.plot
   * xugrid.Ugrid1d.reindex_like
   * xugrid.Ugrid1d.rename
   * xugrid.Ugrid1d.sel
   * xugrid.Ugrid1d.set_crs
   * xugrid.Ugrid1d.set_node_coords
   * xugrid.Ugrid1d.start_index
   * xugrid.Ugrid1d.to_crs
   * xugrid.Ugrid1d.to_shapely
   * xugrid.Ugrid1d.topological_sort_by_dfs
   * xugrid.Ugrid1d.topology_dimension
   * xugrid.Ugrid1d.topology_subset

xugrid.Ugrid1d.assign_edge_coords
=================================
Assign node coordinates from the grid to the object.

Returns a new object with all the original data in addition to the new
node coordinates of the grid.

Parameters
----------
obj: xr.DataArray or xr.Dataset

Returns
-------
assigned (same type as obj)

xugrid.Ugrid1d.assign_node_coords
=================================
Assign node coordinates from the grid to the object.

Returns a new object with all the original data in addition to the new
node coordinates of the grid.

Parameters
----------
obj: xr.DataArray or xr.Dataset

Returns
-------
assigned (same type as obj)

xugrid.Ugrid1d.bounds
=====================
Returns a tuple with the node bounds: xmin, ymin, xmax, ymax

xugrid.Ugrid1d.contract_vertices
================================
Return a simplified network topology by removing all nodes that are
not listed in ``indices``.

Parameters
----------
indices: np.ndarray of integers

Returns
-------
contracted: Ugrid1d

xugrid.Ugrid1d.coords
=====================
Dictionary for grid coordinates.

xugrid.Ugrid1d.copy
===================
Create a deepcopy.

xugrid.Ugrid1d.create_data_array
================================
Create a UgridDataArray from this grid and a 1D array of values.

Parameters
----------
data: array like
    Values for this array. Must be a ``numpy.ndarray`` or castable to
    it.
grid: Ugrid1d, Ugrid2d
facet: str
    With which facet to associate the data. Options for Ugrid1d are,
    ``"node"`` or ``"edge"``. Options for Ugrid2d are ``"node"``,
    ``"edge"``, or ``"face"``.

Returns
-------
uda: UgridDataArray

xugrid.Ugrid1d.dimensions
=========================
Mapping from UGRID dimension names to lengths.

This property will be changed to return a type more consistent with
DataArray.dims in the future, i.e. a set of dimension names.

xugrid.Ugrid1d.dims
===================
Set of UGRID dimension names: node dimension, edge dimension.

xugrid.Ugrid1d.directed_node_node_connectivity
==============================================
Directed node to node connectivity.

The connectivity is represented as an adjacency matrix in CSR format,
with the row and column indices as a (0-based) node index. The data of
the matrix contains the edge index as every connection is formed by an
edge.

Returns
-------
connectivity: csr_matrix

xugrid.Ugrid1d.edge_bounds
==========================
Returns a numpy array with columns ``minx, miny, maxx, maxy``,
describing the bounds of every edge in the grid.

Returns
-------
edge_bounds: np.ndarray of shape (n_edge, 4)

xugrid.Ugrid1d.edge_coordinates
===============================
Centroid (x,y) coordinates of every edge in the UGRID topology

xugrid.Ugrid1d.edge_dimension
=============================
Name of edge dimension

xugrid.Ugrid1d.edge_node_coordinates
====================================
Node coordinates for every edge, shape: ``n_edge, 2, 2``.

xugrid.Ugrid1d.edge_x
=====================
x-coordinate of every edge in the UGRID topology

xugrid.Ugrid1d.edge_y
=====================
y-coordinate of every edge in the UGRID topology

xugrid.Ugrid1d.fill_value
=========================
Fill value for UGRID connectivity arrays.

xugrid.Ugrid1d.find_ugrid_dim
=============================
Find the UGRID dimension that is present in the object.

xugrid.Ugrid1d.from_dataset
===========================
Extract the 1D UGRID topology information from an xarray Dataset.

Parameters
----------
dataset: xr.Dataset
    Dataset containing topology information stored according to UGRID conventions.

Returns
-------
grid: Ugrid1dAdapter

xugrid.Ugrid1d.from_geodataframe
================================
Convert geodataframe of linestrings into a UGRID1D topology.

Parameters
----------
geodataframe: geopandas GeoDataFrame

Returns
-------
topology: Ugrid1d

xugrid.Ugrid1d.from_meshkernel
==============================
Create a 1D UGRID topology from a MeshKernel Mesh1d object.

Parameters
----------
mesh: MeshKernel.Mesh2d
name: str
    Mesh name. Defaults to "network1d".
projected: bool
    Whether node_x and node_y are longitude and latitude or projected x and
    y coordinates. Used to write the appropriate standard_name in the
    coordinate attributes.
crs: Any, optional
    Coordinate Reference System of the geometry objects. Can be anything accepted by
    :meth:`pyproj.CRS.from_user_input() <pyproj.crs.CRS.from_user_input>`,
    such as an authority string (eg "EPSG:4326") or a WKT string.

Returns
-------
grid: Ugrid1d

xugrid.Ugrid1d.from_shapely
===========================
Convert an array of shapely linestrings to UGRID1D topology.

Parameters
----------
geometry: np.ndarray of shapely linestrings
crs: Any, optional
    Coordinate Reference System of the geometry objects. Can be anything accepted by
    :meth:`pyproj.CRS.from_user_input() <pyproj.crs.CRS.from_user_input>`,
    such as an authority string (eg "EPSG:4326") or a WKT string.

xugrid.Ugrid1d.get_connectivity_matrix
======================================
Return the connectivity matrix for the specified UGRID dimension.

xugrid.Ugrid1d.get_coordinates
==============================
Return the coordinates for the specified UGRID dimension.

xugrid.Ugrid1d.isel
===================
Select based on node or edge.

Edge selection always results in a valid UGRID topology. Node selection
may result in invalid topologies (incomplete edges), and will error in
such a case.

Parameters
----------
indexers: dict of str to np.ndarray of integers or bools
return_index: bool, optional
    Whether to return node_index, edge_index.

Returns
-------
obj: xr.Dataset or xr.DataArray
grid: Ugrid2d
indexes: dict
    Dictionary with keys node dimension, edge dimension and values
    their respective index. Only returned if return_index is True.

xugrid.Ugrid1d.merge_partitions
===============================
Merge grid partitions into a single whole.

Duplicate edges are included only once, and removed from subsequent
partitions before merging.

Parameters
----------
grids: sequence of Ugrid1d

Returns
-------
merged: Ugrid1d

xugrid.Ugrid1d.mesh
===================
Create if needed, and return meshkernel Mesh1d object.

Returns
-------
mesh: meshkernel.Mesh1d

xugrid.Ugrid1d.meshkernel
=========================
Create if needed, and return meshkernel MeshKernel instance.

Returns
-------
meshkernel: meshkernel.MeshKernel

xugrid.Ugrid1d.n_edge
=====================
Number of edges in the UGRID topology

xugrid.Ugrid1d.n_node
=====================
Number of nodes (vertices) in the UGRID topology

xugrid.Ugrid1d.node_coordinates
===============================
Coordinates (x, y) of the nodes (vertices)

xugrid.Ugrid1d.node_dimension
=============================
Name of node dimension

xugrid.Ugrid1d.node_edge_connectivity
=====================================
Node to edge connectivity.

Returns
-------
connectivity: csr_matrix

xugrid.Ugrid1d.node_node_connectivity
=====================================
Node to node connectivity.

The connectivity is represented as an adjacency matrix in CSR format,
with the row and column indices as a (0-based) node index. The data of
the matrix contains the edge index as every connection is formed by an
edge.

Returns
-------
connectivity: csr_matrix

xugrid.Ugrid1d.plot
===================
Plot the edges of the mesh.

Parameters
----------
**kwargs : optional
    Additional keyword arguments to ``matplotlib.pyplot.line``.

xugrid.Ugrid1d.reindex_like
===========================
Conform a DataArray or Dataset to match the topology of another Ugrid1D
topology. The topologies must be exactly equivalent: only the order of
the nodes and edges may differ.

Parameters
----------
other: Ugrid1d
obj: DataArray or Dataset
tolerance: float, default value 0.0.
    Maximum distance between inexact coordinate matches.

Returns
-------
reindexed: DataArray or Dataset

xugrid.Ugrid1d.rename
=====================
Create a new grid with all variables named according to the default
naming conventions.

xugrid.Ugrid1d.sel
==================
Select a selection of edges, based on edge centroids.

Parameters
----------
x: slice
y: slice

Returns
-------
dimension: str
as_ugrid: bool
index: 1d array of integers
coords: dict

xugrid.Ugrid1d.set_crs
======================
Set the Coordinate Reference System (CRS) of a UGRID topology.

NOTE: The underlying geometries are not transformed to this CRS. To
transform the geometries to a new CRS, use the ``to_crs`` method.

Parameters
----------
crs : pyproj.CRS, optional if `epsg` is specified
    The value can be anything accepted
    by :meth:`pyproj.CRS.from_user_input() <pyproj.crs.CRS.from_user_input>`,
    such as an authority string (eg "EPSG:4326") or a WKT string.
epsg : int, optional if `crs` is specified
    EPSG code specifying the projection.
allow_override : bool, default False
    If the the UGRID topology already has a CRS, allow to replace the
    existing CRS, even when both are not equal.

xugrid.Ugrid1d.set_node_coords
==============================
Given names of x and y coordinates of the nodes of an object, set them
as the coordinates in the grid.

Parameters
----------
node_x: str
    Name of the x coordinate of the nodes in the object.
node_y: str
    Name of the y coordinate of the nodes in the object.

xugrid.Ugrid1d.start_index
==========================
Start index for UGRID connectivity arrays.

xugrid.Ugrid1d.to_crs
=====================
Transform geometries to a new coordinate reference system.
Transform all geometries in an active geometry column to a different coordinate
reference system. The ``crs`` attribute on the current Ugrid must
be set. Either ``crs`` or ``epsg`` may be specified for output.

This method will transform all points in all objects. It has no notion
of projecting the cells. All segments joining points are assumed to be
lines in the current projection, not geodesics. Objects crossing the
dateline (or other projection boundary) will have undesirable behavior.

Parameters
----------
crs : pyproj.CRS, optional if `epsg` is specified
    The value can be anything accepted by
    :meth:`pyproj.CRS.from_user_input() <pyproj.crs.CRS.from_user_input>`,
    such as an authority string (eg "EPSG:4326") or a WKT string.
epsg : int, optional if `crs` is specified
    EPSG code specifying output projection.

xugrid.Ugrid1d.to_shapely
=========================
Convert UGRID topology to shapely objects.

* nodes: points
* edges: linestrings

Parameters
----------
dim: str
    Node or edge dimension.

Returns
-------
geometry: ndarray of shapely.Geometry

xugrid.Ugrid1d.topological_sort_by_dfs
======================================
Return an array of vertices in topological order.

Returns
-------
sorted_vertices: np.ndarray of integer

xugrid.Ugrid1d.topology_dimension
=================================
Highest dimensionality of the geometric elements: 1

xugrid.Ugrid1d.topology_subset
==============================
Create a new UGRID1D topology for a subset of this topology.

Parameters
----------
edge_index: 1d array of integers or bool
    Edges of the subset.
return_index: bool, optional
    Whether to return node_index, edge_index.

Returns
-------
subset: Ugrid1d
indexes: dict
    Dictionary with keys node dimension and edge dimension and values
    their respective index. Only returned if return_index is True.

xugrid.Ugrid2d
==============
This class stores the topological data of a 2-D unstructured grid.

Parameters
----------
node_x: ndarray of floats
node_y: ndarray of floats
fill_value: int
face_node_connectivity: ndarray of integers
name: string, optional
    Mesh name. Defaults to "mesh2d".
edge_node_connectivity: ndarray of integers, optional
dataset: xr.Dataset, optional
indexes: Dict[str, str], optional
    When a dataset is provided, a mapping from the UGRID role to the dataset
    variable name. E.g. {"face_x": "mesh2d_face_lon"}.
projected: bool, optional
    Whether node_x and node_y are longitude and latitude or projected x and
    y coordinates. Used to write the appropriate standard_name in the
    coordinate attributes.
crs: Any, optional
    Coordinate Reference System of the geometry objects. Can be anything accepted by
    :meth:`pyproj.CRS.from_user_input() <pyproj.crs.CRS.from_user_input>`,
    such as an authority string (eg "EPSG:4326") or a WKT string.
attrs: Dict[str, str], optional
    UGRID topology attributes. Should not be provided together with
    dataset: if other names are required, update the dataset instead.
    A name entry is ignored, as name is given explicitly.
start_index: int, 0 or 1, default is 0.
    Start index of the connectivity arrays. Must match the start index
    of the provided face_node_connectivity and edge_node_connectivity.

xugrid.Ugrid2d Class Members
============================
   * xugrid.Ugrid2d.area
   * xugrid.Ugrid2d.assign_edge_coords
   * xugrid.Ugrid2d.assign_face_coords
   * xugrid.Ugrid2d.assign_node_coords
   * xugrid.Ugrid2d.boundary_node_connectivity
   * xugrid.Ugrid2d.bounding_polygon
   * xugrid.Ugrid2d.bounds
   * xugrid.Ugrid2d.celltree
   * xugrid.Ugrid2d.centroid_triangulation
   * xugrid.Ugrid2d.centroids
   * xugrid.Ugrid2d.circumcenters
   * xugrid.Ugrid2d.compute_barycentric_weights
   * xugrid.Ugrid2d.coords
   * xugrid.Ugrid2d.copy
   * xugrid.Ugrid2d.create_data_array
   * xugrid.Ugrid2d.dimensions
   * xugrid.Ugrid2d.dims
   * xugrid.Ugrid2d.directed_node_node_connectivity
   * xugrid.Ugrid2d.earcut_triangulate_polygons
   * xugrid.Ugrid2d.edge_bounds
   * xugrid.Ugrid2d.edge_coordinates
   * xugrid.Ugrid2d.edge_dimension
   * xugrid.Ugrid2d.edge_face_connectivity
   * xugrid.Ugrid2d.edge_node_connectivity
   * xugrid.Ugrid2d.edge_node_coordinates
   * xugrid.Ugrid2d.edge_x
   * xugrid.Ugrid2d.edge_y
   * xugrid.Ugrid2d.exterior_edges
   * xugrid.Ugrid2d.exterior_faces
   * xugrid.Ugrid2d.face_bounds
   * xugrid.Ugrid2d.face_coordinates
   * xugrid.Ugrid2d.face_dimension
   * xugrid.Ugrid2d.face_edge_connectivity
   * xugrid.Ugrid2d.face_face_connectivity
   * xugrid.Ugrid2d.face_node_coordinates
   * xugrid.Ugrid2d.face_x
   * xugrid.Ugrid2d.face_y
   * xugrid.Ugrid2d.fill_value
   * xugrid.Ugrid2d.find_ugrid_dim
   * xugrid.Ugrid2d.from_dataset
   * xugrid.Ugrid2d.from_geodataframe
   * xugrid.Ugrid2d.from_meshkernel
   * xugrid.Ugrid2d.from_shapely
   * xugrid.Ugrid2d.from_structured
   * xugrid.Ugrid2d.from_structured_bounds
   * xugrid.Ugrid2d.from_structured_intervals1d
   * xugrid.Ugrid2d.from_structured_intervals2d
   * xugrid.Ugrid2d.from_structured_multicoord
   * xugrid.Ugrid2d.get_connectivity_matrix
   * xugrid.Ugrid2d.get_coordinates
   * xugrid.Ugrid2d.intersect_edges
   * xugrid.Ugrid2d.intersect_line
   * xugrid.Ugrid2d.intersect_linestring
   * xugrid.Ugrid2d.isel
   * xugrid.Ugrid2d.label_partitions
   * xugrid.Ugrid2d.locate_bounding_box
   * xugrid.Ugrid2d.locate_points
   * xugrid.Ugrid2d.merge_partitions
   * xugrid.Ugrid2d.mesh
   * xugrid.Ugrid2d.meshkernel
   * xugrid.Ugrid2d.n_edge
   * xugrid.Ugrid2d.n_face
   * xugrid.Ugrid2d.n_max_node_per_face
   * xugrid.Ugrid2d.n_node
   * xugrid.Ugrid2d.node_coordinates
   * xugrid.Ugrid2d.node_dimension
   * xugrid.Ugrid2d.node_edge_connectivity
   * xugrid.Ugrid2d.node_face_connectivity
   * xugrid.Ugrid2d.node_node_connectivity
   * xugrid.Ugrid2d.partition
   * xugrid.Ugrid2d.perimeter
   * xugrid.Ugrid2d.plot
   * xugrid.Ugrid2d.rasterize
   * xugrid.Ugrid2d.rasterize_like
   * xugrid.Ugrid2d.reindex_like
   * xugrid.Ugrid2d.rename
   * xugrid.Ugrid2d.reverse_cuthill_mckee
   * xugrid.Ugrid2d.sel
   * xugrid.Ugrid2d.sel_points
   * xugrid.Ugrid2d.set_crs
   * xugrid.Ugrid2d.set_node_coords
   * xugrid.Ugrid2d.start_index
   * xugrid.Ugrid2d.tesselate_centroidal_voronoi
   * xugrid.Ugrid2d.tesselate_circumcenter_voronoi
   * xugrid.Ugrid2d.to_crs
   * xugrid.Ugrid2d.to_nonperiodic
   * xugrid.Ugrid2d.to_periodic
   * xugrid.Ugrid2d.to_shapely
   * xugrid.Ugrid2d.topology_dimension
   * xugrid.Ugrid2d.topology_subset
   * xugrid.Ugrid2d.triangulate
   * xugrid.Ugrid2d.triangulation
   * xugrid.Ugrid2d.validate_edge_node_connectivity
   * xugrid.Ugrid2d.voronoi_topology

xugrid.Ugrid2d.area
===================
Area of every face.

xugrid.Ugrid2d.assign_face_coords
=================================
Assign face coordinates from the grid to the object.

Returns a new object with all the original data in addition to the new
node coordinates of the grid.

Parameters
----------
obj: xr.DataArray or xr.Dataset

Returns
-------
assigned (same type as obj)

xugrid.Ugrid2d.boundary_node_connectivity
=========================================
Boundary node connectivity

Returns
-------
connectivity: ndarray of integers with shape ``(n_boundary_edge, 2)``

xugrid.Ugrid2d.bounding_polygon
===============================
Construct the bounding polygon of the grid. This polygon may include
holes if the grid also contains holes.

xugrid.Ugrid2d.celltree
=======================
Initializes the celltree if needed, and returns celltree.

A celltree is a search structure for spatial lookups in unstructured grids.

xugrid.Ugrid2d.centroid_triangulation
=====================================
Triangulation of centroidal voronoi tesselation.

Required for e.g. contouring face data, which takes triangles and
associated values at the triangle vertices.

Returns
-------
vertices: ndarray of floats with shape ``(n_centroids, 2)``
face_node_connectivity: ndarray of integers with shape ``(n_triangle, 3)``
    Describes face node connectivity of triangle topology.
face_index: 1d array of integers

xugrid.Ugrid2d.centroids
========================
Centroid (x, y) of every face.

Returns
-------
centroids: ndarray of floats with shape ``(n_face, 2)``

xugrid.Ugrid2d.circumcenters
============================
Circumenter (x, y) of every face; only works for fully triangular
grids.

xugrid.Ugrid2d.compute_barycentric_weights
==========================================
Find in which face the points are located, and compute the barycentric
weight for every vertex of the face.

Parameters
----------
points: ndarray of floats with shape ``(n_point, 2)``

Returns
-------
face_index: ndarray of integers with shape ``(n_points,)``
weights: ndarray of floats with shape ```(n_points, n_max_node)``

xugrid.Ugrid2d.coords
=====================
Dictionary for grid coordinates.

xugrid.Ugrid2d.create_data_array
================================
Create a UgridDataArray from this grid and a 1D array of values.

Parameters
----------
data: array like
    Values for this array. Must be a ``numpy.ndarray`` or castable to
    it.
grid: Ugrid1d, Ugrid2d
facet: str
    With which facet to associate the data. Options for Ugrid1d are,
    ``"node"`` or ``"edge"``. Options for Ugrid2d are ``"node"``,
    ``"edge"``, or ``"face"``.

Returns
-------
uda: UgridDataArray

xugrid.Ugrid2d.dims
===================
Set of UGRID dimension names: node dimension, edge dimension, face_dimension.

xugrid.Ugrid2d.earcut_triangulate_polygons
==========================================
Break down polygons using mapbox_earcut, and create a mesh from the
resulting triangles.

Parameters
----------
polygons: ndarray of shapely polygons
return_index: bool, default is False.

Returns
-------
grid: xugrid.Ugrid2d
index: ndarray of integer, optional
    The polygon index for each triangle. Only provided if ``return_index``
    is True.

xugrid.Ugrid2d.edge_face_connectivity
=====================================
Edge to face connectivity. An edge may belong to a single face
(exterior edge), or it may be shared by two faces (interior edge).

An exterior edge will contain a FILL_VALUE of -1 for the second column.

Returns
-------
connectivity: ndarray of integers with shape ``(n_edge, 2)``.

xugrid.Ugrid2d.edge_node_connectivity
=====================================
Edge to node connectivity. Every edge consists of a connection between
two nodes.

Returns
-------
connectivity: ndarray of integers with shape ``(n_edge, 2)``.

xugrid.Ugrid2d.exterior_edges
=============================
Get all exterior edges, i.e. edges with no other face.

Returns
-------
edge_index: 1d array of integers

xugrid.Ugrid2d.exterior_faces
=============================
Get all exterior faces, i.e. faces with an unshared edge.

Returns
-------
face_index: 1d array of integers

xugrid.Ugrid2d.face_bounds
==========================
Returns a numpy array with columns ``minx, miny, maxx, maxy``,
describing the bounds of every face in the grid.

Returns
-------
face_bounds: np.ndarray of shape (n_face, 4)

xugrid.Ugrid2d.face_coordinates
===============================
Centroid (x, y) of every face.

Returns
-------
centroids: ndarray of floats with shape ``(n_face, 2)``

xugrid.Ugrid2d.face_dimension
=============================
Return the name of the face dimension.

xugrid.Ugrid2d.face_edge_connectivity
=====================================
Face to edge connectivity.

Returns
-------
connectivity: csr_matrix

xugrid.Ugrid2d.face_face_connectivity
=====================================
Face to face connectivity. Derived from shared edges.

The connectivity is represented as an adjacency matrix in CSR format,
with the row and column indices as a (0-based) face index. The data of
the matrix contains the edge index as every connection is formed by a
shared edge.

Returns
-------
connectivity: csr_matrix

xugrid.Ugrid2d.face_node_coordinates
====================================
Node coordinates of every face.

"Fill node" coordinates are set as NaN.

Returns
-------
face_node_coordinates: ndarray of floats with shape ``(n_face, n_max_node_per_face, 2)``

xugrid.Ugrid2d.face_x
=====================
x-coordinate of centroid of every face

xugrid.Ugrid2d.face_y
=====================
y-coordinate of centroid of every face

xugrid.Ugrid2d.from_dataset
===========================
Extract the 2D UGRID topology information from an xarray Dataset.

Parameters
----------
dataset: xr.Dataset
    Dataset containing topology information stored according to UGRID conventions.

Returns
-------
grid: Ugrid1dAdapter

xugrid.Ugrid2d.from_geodataframe
================================
Convert a geodataframe of polygons to UGRID2D topology.

Parameters
----------
geodataframe: geopandas GeoDataFrame

Returns
-------
topology: Ugrid2d

xugrid.Ugrid2d.from_meshkernel
==============================
Create a 2D UGRID topology from a MeshKernel Mesh2d object.

Parameters
----------
mesh: MeshKernel.Mesh2d
name: str
    Mesh name. Defaults to "mesh2d".
projected: bool
    Whether node_x and node_y are longitude and latitude or projected x and
    y coordinates. Used to write the appropriate standard_name in the
    coordinate attributes.
crs: Any, optional
    Coordinate Reference System of the geometry objects. Can be anything accepted by
    :meth:`pyproj.CRS.from_user_input() <pyproj.crs.CRS.from_user_input>`,
    such as an authority string (eg "EPSG:4326") or a WKT string.

Returns
-------
grid: Ugrid2d

xugrid.Ugrid2d.from_shapely
===========================
Convert an array of shapely polygons to UGRID2D topology.

Parameters
----------
geometry: np.ndarray of shapely polygons
crs: Any, optional
    Coordinate Reference System of the geometry objects. Can be anything accepted by
    :meth:`pyproj.CRS.from_user_input() <pyproj.crs.CRS.from_user_input>`,
    such as an authority string (eg "EPSG:4326") or a WKT string.

Returns
-------
topology: Ugrid2d

xugrid.Ugrid2d.from_structured
==============================
Create a Ugrid2d topology from an axis-aligned rectilinear structured topology.

This method assumes the coordinates are 1D.

Use ``from_structured_multicoord`` for 2D x and y coordinates, e.g. for
(approximated) curvilinear and rotated structured topologies.

Parameters
----------
data: xr.DataArray or xr.Dataset
x: str, optional
    Name of the 1D coordinate to use as the UGRID x-coordinate.
y: str, optional
    Name of the 1D coordinate to use as the UGRID y-coordinate.

Returns
-------
grid: Ugrid2d

xugrid.Ugrid2d.from_structured_bounds
=====================================
Create a Ugrid2d topology from a structured topology based on 1D bounds.

The bounds contain the lower and upper cell boundary for each cell.

Parameters
----------
x_bounds: np.ndarray of shape (M, 2)
    x-coordinate bounds for N row and M columns.
y_bounds: np.ndarray of shape (N, 2)
    y-coordinate bounds for N row and M columns.

Returns
-------
grid: Ugrid2d

xugrid.Ugrid2d.from_structured_intervals1d
==========================================
Create a Ugrid2d topology from a structured topology based on 1D intervals.

Parameters
----------
x_intervals: np.ndarray of shape (M + 1,)
    x-coordinate interval values for N row and M columns.
y_intervals: np.ndarray of shape (N + 1,)
    y-coordinate interval values for N row and M columns.

xugrid.Ugrid2d.from_structured_intervals2d
==========================================
Create a Ugrid2d topology from a structured topology based on 2D intervals.

Parameters
----------
x_intervals: np.ndarray of shape shape (N + 1, M + 1)
    x-coordinate interval values for N row and M columns.
y_intervals: np.ndarray of shape shape (N + 1, M + 1)
    y-coordinate interval values for N row and M columns.

xugrid.Ugrid2d.from_structured_multicoord
=========================================
Create a Ugrid2d topology from a structured topology, including rotated
and (approximated) curvilinear topologies.

This method assumes the coordinates are 2D.

Use ``from_structured`` for 1D x and y coordinates, which is generally
the case for axis-aligned rectilinear topologies (most rasters).

Parameters
----------
data: xr.DataArray or xr.Dataset
x: str
    Name of the 2D coordinate to use as the UGRID x-coordinate.
y: str
    Name of the 2D coordinate to use as the UGRID y-coordinate.

Returns
-------
grid: Ugrid2d

xugrid.Ugrid2d.get_connectivity_matrix
======================================
Return the connectivity matrix for the specified UGRID dimension.

xugrid.Ugrid2d.get_coordinates
==============================
Return the coordinates for the specified UGRID dimension.

xugrid.Ugrid2d.intersect_edges
==============================
Find in which face edges are located and compute the intersection with
the face edges.

Parameters
----------
edges: ndarray of floats with shape ``(n_edge, 2, 2)``
    The first dimensions represents the different edges.
    The second dimensions represents the start and end of every edge.
    The third dimensions reresent the x and y coordinate of every vertex.

Returns
-------
edge_index: ndarray of integers with shape ``(n_intersection,)``
face_index: ndarray of integers with shape ``(n_intersection,)``
intersections: ndarray of float with shape ``(n_intersection, 2, 2)``

xugrid.Ugrid2d.intersect_line
=============================
Intersect a line with this grid, and fetch the values of the
intersected faces.

Parameters
----------
obj: xr.DataArray or xr.Dataset
start: sequence of two floats
    coordinate pair (x, y), designating the start point of the line.
end: sequence of two floats
    coordinate pair (x, y), designating the end point of the line.

Returns
-------
selection: xr.DataArray or xr.Dataset
    The name of the topology is prefixed in the x, y and s
    (spatium=distance) coordinates.

xugrid.Ugrid2d.intersect_linestring
===================================
Intersect linestrings with this grid, and fetch the values of the
intersected faces.

Parameters
----------
obj: xr.DataArray or xr.Dataset
linestring: shapely.geometry.lineString

Returns
-------
selection: xr.DataArray or xr.Dataset
    The name of the topology is prefixed in the x, y and s
    (spatium=distance) coordinates.

xugrid.Ugrid2d.isel
===================
Select based on node, edge, or face.

Face selection always results in a valid UGRID topology.
Node or edge selection may result in invalid topologies (incomplete
faces), and will error in such a case.

Parameters
----------
indexers: dict of str to np.ndarray of integers or bools
return_index: bool, optional
    Whether to return node_index, edge_index, face_index.

Returns
-------
obj: xr.Dataset or xr.DataArray
grid: Ugrid2d
indexes: dict
    Dictionary with keys node dimension, edge dimension, face dimension
    and values their respective index. Only returned if return_index is
    True.

xugrid.Ugrid2d.label_partitions
===============================
Generate partition labesl for this grid topology using METIS:
https://github.com/KarypisLab/METIS

This method utilizes the pymetis Python bindings:
https://github.com/inducer/pymetis

Parameters
----------
n_part: integer
    The number of parts to partition the mesh.

Returns
-------
partition_labels: UgridDataArray of integers

xugrid.Ugrid2d.locate_bounding_box
==================================
Find which faces are located in the bounding box. The centroids of the
faces are used.

Parameters
----------
xmin: float,
ymin: float,
xmax: float,
ymax: float

Returns
-------
face_index: ndarray of bools with shape ``(n_face,)``

xugrid.Ugrid2d.locate_points
============================
Find in which face points are located.

Parameters
----------
points: ndarray of floats with shape ``(n_point, 2)``

Returns
-------
face_index: ndarray of integers with shape ``(n_points,)``

xugrid.Ugrid2d.merge_partitions
===============================
Merge grid partitions into a single whole.

Duplicate faces are included only once, and removed from subsequent
partitions before merging.

Parameters
----------
grids: sequence of Ugrid2d

Returns
-------
merged: Ugrid2d

xugrid.Ugrid2d.mesh
===================
Create if needed, and return meshkernel Mesh2d object.

Returns
-------
mesh: meshkernel.Mesh2d

xugrid.Ugrid2d.meshkernel
=========================
Create if needed, and return meshkernel MeshKernel instance.

Returns
-------
meshkernel: meshkernel.MeshKernel

xugrid.Ugrid2d.n_face
=====================
Return the number of faces in the UGRID2D topology.

xugrid.Ugrid2d.n_max_node_per_face
==================================
Return the maximum number of nodes that a face can contain in the
UGRID2D topology.

xugrid.Ugrid2d.node_face_connectivity
=====================================
Node to face connectivity. Inverted from face node connectivity.

Returns
-------
connectivity: csr_matrix

xugrid.Ugrid2d.partition
========================
Partition this grid topology using METIS:
https://github.com/KarypisLab/METIS

This method utilizes the pymetis Python bindings:
https://github.com/inducer/pymetis

Parameters
----------
n_part: integer
    The number of parts to partition the mesh.

Returns
-------
partitions

xugrid.Ugrid2d.perimeter
========================
Perimeter length of every face.

xugrid.Ugrid2d.rasterize
========================
Rasterize unstructured grid by sampling.

x and y coordinates are generated from the bounds of the UGRID2D
topology and the provided resolution.

Parameters
----------
resolution: float
    Spacing in x and y.
bounds: tuple of four floats, optional
    xmin, ymin, xmax, ymax

Returns
-------
x: 1d array of floats with shape ``(ncol,)``
y: 1d array of floats with shape ``(nrow,)``
face_index: 1d array of integers with shape ``(nrow * ncol,)``

xugrid.Ugrid2d.rasterize_like
=============================
Rasterize unstructured grid by sampling on the x and y coordinates.

Parameters
----------
x: 1d array of floats with shape ``(ncol,)``
y: 1d array of floats with shape ``(nrow,)``

Returns
-------
x: 1d array of floats with shape ``(ncol,)``
y: 1d array of floats with shape ``(nrow,)``
face_index: 1d array of integers with shape ``(nrow * ncol,)``

xugrid.Ugrid2d.reindex_like
===========================
Conform a DataArray or Dataset to match the topology of another Ugrid2D
topology. The topologies must be exactly equivalent: only the order of
the nodes, edges, and faces may differ.

Parameters
----------
other: Ugrid2d
obj: DataArray or Dataset
tolerance: float, default value 0.0.
    Maximum distance between inexact coordinate matches.

Returns
-------
reindexed: DataArray or Dataset

xugrid.Ugrid2d.reverse_cuthill_mckee
====================================
Reduces bandwith of the connectivity matrix.

Wraps :py:func:`scipy.sparse.csgraph.reverse_cuthill_mckee`.

Returns
-------
reordered: Ugrid2d

xugrid.Ugrid2d.sel
==================
Find selection in the UGRID x and y coordinates.

The indexing for x and y always occurs orthogonally, i.e.:
``.sel(x=[0.0, 5.0], y=[10.0, 15.0])`` results in a four points. For
vectorized indexing (equal to ``zip``ing through x and y), see
``.sel_points``.

Parameters
----------
obj: xr.DataArray or xr.Dataset
x: float, 1d array, slice
y: float, 1d array, slice

Returns
-------
dimension: str
as_ugrid: bool
index: 1d array of integers
coords: dict

xugrid.Ugrid2d.sel_points
=========================
Select points in the unstructured grid.


Parameters
----------
x: 1d array of floats with shape ``(n_points,)``
y: 1d array of floats with shape ``(n_points,)``
obj: xr.DataArray or xr.Dataset
out_of_bounds: str, default ``"warn"``
    What to do when points are located outside of any feature:

    * raise: raise a ValueError.
    * ignore: return ``fill_value`` for the out of bounds points.
    * warn: give a warning and return NaN for the out of bounds points.
    * drop: drop the out of bounds points. They may be identified
      via the ``index`` coordinate of the returned selection.
fill_value: scalar, DataArray, Dataset, or callable, optional, default: np.nan
    Value to assign to out-of-bounds points if out_of_bounds is warn
    or ignore. Forwarded to xarray's ``.where()`` method.

Returns
-------
selection: xr.DataArray or xr.Dataset
    The name of the topology is prefixed in the x, y coordinates.

xugrid.Ugrid2d.tesselate_centroidal_voronoi
===========================================
Create a centroidal Voronoi tesselation of this UGRID2D topology.

Such a tesselation is not guaranteed to produce convex cells. To ensure
convexity, set ``add_vertices=False`` -- this will result in a
different exterior, however.

Parameters
----------
add_exterior: bool, default: True
add_vertices: bool, default: True
skip_concave: bool, default: False

Returns
-------
tesselation: Ugrid2d

xugrid.Ugrid2d.tesselate_circumcenter_voronoi
=============================================
Create a circumcenter Voronoi tesselation of this UGRID2D topology.

Such a tesselation is not guaranteed to produce convex cells. To ensure
convexity, set ``add_vertices=False`` -- this will result in a
different exterior, however.

Parameters
----------
add_exterior: bool, default: True
add_vertices: bool, default: True
skip_concave: bool, default: False

Returns
-------
tesselation: Ugrid2d

xugrid.Ugrid2d.to_nonperiodic
=============================
Convert this grid from a periodic grid (where the rightmost boundary shares its
nodes with the leftmost boundary) to an aperiodic grid, where the leftmost nodes
are separate from the rightmost nodes.

Parameters
----------
xmax: float
    The x-value of the newly created rightmost boundary nodes.
obj: xr.DataArray or xr.Dataset

Returns
-------
nonperiodic_grid: Ugrid2d
aligned: xr.DataArray or xr.Dataset

xugrid.Ugrid2d.to_periodic
==========================
Convert this grid to a periodic grid, where the rightmost nodes are
equal to the leftmost nodes. Note: for this to work, the y-coordinates
on the left boundary must match those on the right boundary exactly.

Returns
-------
periodic_grid: Ugrid2d
aligned: xr.DataArray or xr.Dataset

xugrid.Ugrid2d.to_shapely
=========================
Convert UGRID topology to shapely objects.

* nodes: points
* edges: linestrings
* faces: polygons

Parameters
----------
dim: str
    Node, edge, or face dimension.

Returns
-------
geometry: ndarray of shapely.Geometry

xugrid.Ugrid2d.topology_dimension
=================================
Highest dimensionality of the geometric elements: 2

xugrid.Ugrid2d.topology_subset
==============================
Create a new UGRID1D topology for a subset of this topology.

Parameters
----------
face_index: 1d array of integers or bool
    Edges of the subset.
return_index: bool, optional
    Whether to return node_index, edge_index, face_index.

Returns
-------
subset: Ugrid2d
indexes: dict
    Dictionary with keys node dimension, edge dimension, face dimension
    and values their respective index. Only returned if return_index is
    True.

xugrid.Ugrid2d.triangulate
==========================
Triangulate this UGRID2D topology, breaks more complex polygons down
into triangles.

Returns
-------
triangles: Ugrid2d

xugrid.Ugrid2d.triangulation
============================
Triangulation of the UGRID2D topology.

Returns
-------
triangulation: tuple
    Contains node_x, node_y, triangle face_node_connectivity.
triangle_face_connectivity: 1d array of integers
    Identifies the original face for every triangle.

xugrid.Ugrid2d.validate_edge_node_connectivity
==============================================
Mark valid edges, by comparing face_node_connectivity and
edge_node_connectivity. Edges that are not part of a face, as well as
duplicate edges are marked ``False``.

An error is raised if the face_node_connectivity defines more unique
edges than the edge_node_connectivity.

Returns
-------
valid: np.ndarray of bool
    Marks for every edge whether it is valid.

Examples
--------
To purge invalid edges and associated data from a dataset that contains
un-associated or duplicate edges:

>>> uds = xugrid.open_dataset("example.nc")
>>> valid = uds.ugrid.grid.validate_edge_node_connectivity()
>>> purged = uds.isel({grid.edge_dimension: valid})

xugrid.Ugrid2d.voronoi_topology
===============================
Centroidal Voronoi tesselation of this UGRID2D topology.

Returns
-------
vertices: ndarray of floats with shape ``(n_centroids, 2)``
face_node_connectivity: csr_matrix
    Describes face node connectivity of voronoi topology.
face_index: 1d array of integers

xugrid.UgridRolesAccessor
=========================
Xarray Dataset "accessor" to retrieve the names of UGRID variables.

Examples
--------
To get a list of the UGRID dummy variables in the dataset:

>>> dataset.ugrid_roles.topology

To get the names of the connectivity variables in the dataset:

>>> dataset.ugrid_roles.connectivity

Names can also be accessed directly through the topology:

>>> dataset.ugrid_roles["mesh2d"]["node_dimension"]

xugrid.UgridRolesAccessor Class Members
=======================================
   * xugrid.UgridRolesAccessor.connectivity
   * xugrid.UgridRolesAccessor.coordinates
   * xugrid.UgridRolesAccessor.dimensions
   * xugrid.UgridRolesAccessor.topology

xugrid.UgridRolesAccessor.connectivity
======================================
Get the names of the variables containing the UGRID connectivity data.

    * face_node_connectivity
    * edge_node_connectivity
    * face_edge_connectivity
    * edge_face_connectivity

Returns
-------
connectivity: Dict[str, Dict[str, str]]

xugrid.UgridRolesAccessor.coordinates
=====================================
Get the names of the coordinate variables from the topology attributes.

Returns a dictionary with the coordinates for the UGRID coordinates:

    * node coordinates
    * edge coordinates
    * face coordinates

Multiple coordinates may be defined. The coordinates are grouped by
their role (x or y).

Returns
-------
coordinates: dict[str, dict[str, Tuple[List[str]]]]

xugrid.UgridRolesAccessor.dimensions
====================================
Get the dimension names from the topology attributes and infer them
from connectivity arrays or coordinates.

Returns a dictionary with the UGRID dimensions per topology:

    * node dimension
    * edge dimension
    * face dimension

Returns
-------
dimensions: dict[str, dict[str, str]]

xugrid.UgridRolesAccessor.topology
==================================
Get the names of the topology dummy variables, marked by a CF-role of
``mesh_topology``.

Returns
-------
topology: List[str]

