imod.mf6.open_hds
=================
Open modflow6 heads (.hds) file.

The data is lazily read per timestep and automatically converted into
(dense) xr.DataArrays or xu.UgridDataArrays, for DIS and DISV respectively.
The conversion is done via the information stored in the Binary Grid file
(GRB).


Parameters
----------
hds_path: Union[str, pathlib.Path]
grb_path: Union[str, pathlib.Path]
dry_nan: bool, default value: False.
    Whether to convert dry values to NaN.
simulation_start_time : Optional datetime
    The time and date correpsonding to the beginning of the simulation.
    Use this to convert the time coordinates of the output array to
    calendar time/dates. time_unit must also be present if this argument is present.
time_unit: Optional str
    The time unit MF6 is working in, in string representation.
    Only used if simulation_start_time was provided.
    Admissible values are:
    ns -> nanosecond
    ms -> microsecond
    s -> second
    m -> minute
    h -> hour
    d -> day
    w -> week
    Units "month" or "year" are not supported, as they do not represent unambiguous timedelta values durations.

Returns
-------
head: Union[xr.DataArray, xu.UgridDataArray]

imod.mf6.open_cbc
=================
Open modflow6 cell-by-cell (.cbc) file.

The data is lazily read per timestep and automatically converted into
(dense) xr.DataArrays or xu.UgridDataArrays, for DIS and DISV respectively.
The conversion is done via the information stored in the Binary Grid file
(GRB).

The ``flowja`` argument controls whether the flow-ja-face array (if present)
is returned in grid form as "as is". By default ``flowja=False`` and the
array is returned in "grid form", meaning:

    * DIS: in right, front, and lower face flow. All flows are placed in
      the cell.
    * DISV: in horizontal and lower face flow.the horizontal flows are
      placed on the edges and the lower face flow is placed on the faces.

When ``flowja=True``, the flow-ja-face array is returned as it is found in
the CBC file, with a flow for every cell to cell connection. Additionally,
a ``connectivity`` DataArray is returned describing for every cell (n) its
connected cells (m).

Parameters
----------
cbc_path: str, pathlib.Path
    Path to the cell-by-cell flows file
grb_path: str, pathlib.Path
    Path to the binary grid file
flowja: bool, default value: False
    Whether to return the flow-ja-face values "as is" (``True``) or in a
    grid form (``False``).
simulation_start_time : Optional datetime
    The time and date correpsonding to the beginning of the simulation.
    Use this to convert the time coordinates of the output array to
    calendar time/dates. time_unit must also be present if this argument is present.
time_unit: Optional str
    The time unit MF6 is working in, in string representation.
    Only used if simulation_start_time was provided.
    Admissible values are:
    ns -> nanosecond
    ms -> microsecond
    s -> second
    m -> minute
    h -> hour
    d -> day
    w -> week
    Units "month" or "year" are not supported, as they do not represent unambiguous timedelta values durations.
merge_to_dataset: bool, default value: False
    Merge output to dataset.

Returns
-------
cbc_content: xr.Dataset | Dict[str, xr.DataArray]
    DataArray contains float64 data of the budgets, with dimensions ("time",
    "layer", "y", "x").

Examples
--------

Open a cbc file:

>>> import imod
>>> cbc_content = imod.mf6.open_cbc("budgets.cbc", "my-model.grb")

Check the contents:

>>> print(cbc_content.keys())

Get the drainage budget, compute a time mean for the first layer:

>>> drn_budget = cbc_content["drn]
>>> mean = drn_budget.sel(layer=1).mean("time")

imod.mf6.read_cbc_headers
=========================
Read all the header data from a cell-by-cell (.cbc) budget file.

All budget data for a MODFLOW6 model is stored in a single file. This
function collects all header data, as well as the starting byte position of
the actual budget data.

This function groups the headers per TEXT record (e.g. "flow-ja-face",
"drn", etc.). The headers are stored as a list of named tuples.
flow-ja-face, storage-ss, and storage-sy are written using IMETH=1, all
others with IMETH=6.

Parameters
----------
cbc_path: str, pathlib.Path
    Path to the budget file.

Returns
-------
headers: Dict[List[UnionImeth1Header, Imeth6Header]]
    Dictionary containing a list of headers per TEXT record in the budget
    file.

imod.mf6.Modflow6Simulation
===========================
A MutableMapping is a generic container for associating
key/value pairs.

This class provides concrete generic implementations of all
methods except for __getitem__, __setitem__, __delitem__,
__iter__, and __len__.

imod.mf6.Modflow6Simulation Class Members
=========================================
   * imod.mf6.Modflow6Simulation.clear
   * imod.mf6.Modflow6Simulation.clip_box
   * imod.mf6.Modflow6Simulation.create_time_discretization
   * imod.mf6.Modflow6Simulation.dump
   * imod.mf6.Modflow6Simulation.get
   * imod.mf6.Modflow6Simulation.items
   * imod.mf6.Modflow6Simulation.keys
   * imod.mf6.Modflow6Simulation.mask_all_models
   * imod.mf6.Modflow6Simulation.open_concentration
   * imod.mf6.Modflow6Simulation.open_flow_budget
   * imod.mf6.Modflow6Simulation.open_head
   * imod.mf6.Modflow6Simulation.open_transport_budget
   * imod.mf6.Modflow6Simulation.pop
   * imod.mf6.Modflow6Simulation.popitem
   * imod.mf6.Modflow6Simulation.regrid_like
   * imod.mf6.Modflow6Simulation.render
   * imod.mf6.Modflow6Simulation.run
   * imod.mf6.Modflow6Simulation.setdefault
   * imod.mf6.Modflow6Simulation.split
   * imod.mf6.Modflow6Simulation.update
   * imod.mf6.Modflow6Simulation.values
   * imod.mf6.Modflow6Simulation.write

imod.mf6.Modflow6Simulation.clear
=================================
D.clear() -> None.  Remove all items from D.

imod.mf6.Modflow6Simulation.clip_box
====================================
Clip a simulation by a bounding box (time, layer, y, x).

Slicing intervals may be half-bounded, by providing None:

* To select 500.0 <= x <= 1000.0:
  ``clip_box(x_min=500.0, x_max=1000.0)``.
* To select x <= 1000.0: ``clip_box(x_min=None, x_max=1000.0)``
  or ``clip_box(x_max=1000.0)``.
* To select x >= 500.0: ``clip_box(x_min = 500.0, x_max=None.0)``
  or ``clip_box(x_min=1000.0)``.

Parameters
----------
time_min: optional
time_max: optional
layer_min: optional, int
layer_max: optional, int
x_min: optional, float
x_max: optional, float
y_min: optional, float
y_max: optional, float
states_for_boundary : optional, Dict[pkg_name:str, boundary_values:Union[xr.DataArray, xu.UgridDataArray]]

Returns
-------
clipped : Simulation

imod.mf6.Modflow6Simulation.create_time_discretization
======================================================
Collect all unique times from model packages and additional given
`times`. These unique times are used as stress periods in the model. All
stress packages must have the same starting time. Function creates
TimeDiscretization object which is set to self["time_discretization"]

The time discretization in imod-python works as follows:

- The datetimes of all packages you send in are always respected
- Subsequently, the input data you use is always included fully as well
- All times are treated as starting times for the stress: a stress is
  always applied until the next specified date
- For this reason, a final time is required to determine the length of
  the last stress period
- Additional times can be provided to force shorter stress periods &
  more detailed output
- Every stress has to be defined on the first stress period (this is a
  modflow requirement)

Or visually (every letter a date in the time axes):

>>> recharge a - b - c - d - e - f
>>> river    g - - - - h - - - - j
>>> times    - - - - - - - - - - - i
>>> model    a - b - c h d - e - f i

with the stress periods defined between these dates. I.e. the model
times are the set of all times you include in the model.

Parameters
----------
additional_times : str, datetime; or iterable of str, datetimes.
    Times to add to the time discretization. At least one single time
    should be given, which will be used as the ending time of the
    simulation.

Note
----
To set the other parameters of the TimeDiscretization object, you have
to set these to the object after calling this function.

Example
-------
>>> simulation = imod.mf6.Modflow6Simulation("example")
>>> simulation.create_time_discretization(times=["2000-01-01", "2000-01-02"])
>>> # Set number of timesteps
>>> simulation["time_discretization"]["n_timesteps"] = 5

imod.mf6.Modflow6Simulation.dump
================================
Dump simulation to files. Writes a model definition as .TOML file, which
points to data for each package. Each package is stored as a separate
NetCDF. Structured grids are saved as regular NetCDFs, unstructured
grids are saved as UGRID NetCDF. Structured grids are always made GDAL
compliant, unstructured grids can be made MDAL compliant optionally.

Parameters
----------
directory: str or Path, optional
    directory to dump simulation into. Defaults to current working directory.
validate: bool, optional
    Whether to validate simulation data. Defaults to True.
mdal_compliant: bool, optional
    Convert data with
    :func:`imod.prepare.spatial.mdal_compliant_ugrid2d` to MDAL
    compliant unstructured grids. Defaults to False.
crs: Any, optional
    Anything accepted by rasterio.crs.CRS.from_user_input
    Requires ``rioxarray`` installed.

imod.mf6.Modflow6Simulation.get
===============================
D.get(k[,d]) -> D[k] if k in D, else d.  d defaults to None.

imod.mf6.Modflow6Simulation.items
=================================
D.items() -> a set-like object providing a view on D's items

imod.mf6.Modflow6Simulation.keys
================================
D.keys() -> a set-like object providing a view on D's keys

imod.mf6.Modflow6Simulation.mask_all_models
===========================================
This function applies a mask to all models in a simulation, provided they use
the same discretization. The  method parameter "mask" is an idomain-like array.
Masking will overwrite idomain with the mask where the mask is 0 or -1.
Where the mask is 1, the original value of idomain will be kept.
Masking will update the packages accordingly, blanking their input where needed,
and is therefore not a reversible operation.

Parameters
----------
mask: xr.DataArray, xu.UgridDataArray of ints
    idomain-like integer array. 1 sets cells to active, 0 sets cells to inactive,
    -1 sets cells to vertical passthrough

imod.mf6.Modflow6Simulation.open_concentration
==============================================
Open concentration of finished simulation, requires that the ``run``
method has been called.

The data is lazily read per timestep and automatically converted into
(dense) xr.DataArrays or xu.UgridDataArrays, for DIS and DISV
respectively. The conversion is done via the information stored in the
Binary Grid file (GRB).

Parameters
----------
species_ls: list of strings, default value: None.
    List of species names, which will be used to concatenate the
    concentrations along the ``"species"`` dimension, in case the
    simulation has multiple species and thus multiple transport models.
    If None, transport model names will be used as species names.
dry_nan: bool, default value: False.
    Whether to convert dry values to NaN.

Returns
-------
concentration: Union[xr.DataArray, xu.UgridDataArray]

Examples
--------
Make sure you write and run your model first

>>> simulation.write(path/to/model)
>>> simulation.run()

Then open concentrations:

>>> concentration = simulation.open_concentration()

imod.mf6.Modflow6Simulation.open_flow_budget
============================================
Open flow budgets of finished simulation, requires that the ``run``
method has been called.

The data is lazily read per timestep and automatically converted into
(dense) xr.DataArrays or xu.UgridDataArrays, for DIS and DISV
respectively. The conversion is done via the information stored in the
Binary Grid file (GRB).

The ``flowja`` argument controls whether the flow-ja-face array (if
present) is returned in grid form as "as is". By default
``flowja=False`` and the array is returned in "grid form", meaning:

    * DIS: in right, front, and lower face flow. All flows are placed in
      the cell.
    * DISV: in horizontal and lower face flow.the horizontal flows are
      placed on the edges and the lower face flow is placed on the faces.

When ``flowja=True``, the flow-ja-face array is returned as it is found in
the CBC file, with a flow for every cell to cell connection. Additionally,
a ``connectivity`` DataArray is returned describing for every cell (n) its
connected cells (m).

Parameters
----------
flowja: bool, default value: False
    Whether to return the flow-ja-face values "as is" (``True``) or in a
    grid form (``False``).

Returns
-------
budget: Dict[str, xr.DataArray|xu.UgridDataArray]
    DataArray contains float64 data of the budgets, with dimensions ("time",
    "layer", "y", "x").

Examples
--------
Make sure you write and run your model first

>>> simulation.write(path/to/model)
>>> simulation.run()

Then open budgets:

>>> budget = simulation.open_flow_budget()

Check the contents:

>>> print(budget.keys())

Get the drainage budget, compute a time mean for the first layer:

>>> drn_budget = budget["drn]
>>> mean = drn_budget.sel(layer=1).mean("time")

imod.mf6.Modflow6Simulation.open_head
=====================================
Open heads of finished simulation, requires that the ``run`` method has
been called.

The data is lazily read per timestep and automatically converted into
(dense) xr.DataArrays or xu.UgridDataArrays, for DIS and DISV
respectively. The conversion is done via the information stored in the
Binary Grid file (GRB).

Parameters
----------
dry_nan: bool, default value: False.
    Whether to convert dry values to NaN.
simulation_start_time : Optional datetime
    The time and date correpsonding to the beginning of the simulation.
    Use this to convert the time coordinates of the output array to
    calendar time/dates. time_unit must also be present if this argument is present.
time_unit: Optional str
    The time unit MF6 is working in, in string representation.
    Only used if simulation_start_time was provided.
    Admissible values are:
    ns -> nanosecond
    ms -> microsecond
    s -> second
    m -> minute
    h -> hour
    d -> day
    w -> week
    Units "month" or "year" are not supported, as they do not represent unambiguous timedelta values durations.

Returns
-------
head: Union[xr.DataArray, xu.UgridDataArray]

Examples
--------
Make sure you write and run your model first

>>> simulation.write(path/to/model)
>>> simulation.run()

Then open heads:

>>> head = simulation.open_head()

imod.mf6.Modflow6Simulation.open_transport_budget
=================================================
Open transport budgets of finished simulation, requires that the ``run``
method has been called.

The data is lazily read per timestep and automatically converted into
(dense) xr.DataArrays or xu.UgridDataArrays, for DIS and DISV
respectively. The conversion is done via the information stored in the
Binary Grid file (GRB).

Parameters
----------
species_ls: list of strings, default value: None.
    List of species names, which will be used to concatenate the
    concentrations along the ``"species"`` dimension, in case the
    simulation has multiple species and thus multiple transport models.
    If None, transport model names will be used as species names.

Returns
-------
budget: Dict[str, xr.DataArray|xu.UgridDataArray]
    DataArray contains float64 data of the budgets, with dimensions ("time",
    "layer", "y", "x").

imod.mf6.Modflow6Simulation.pop
===============================
D.pop(k[,d]) -> v, remove specified key and return the corresponding value.
If key is not found, d is returned if given, otherwise KeyError is raised.

imod.mf6.Modflow6Simulation.popitem
===================================
D.popitem() -> (k, v), remove and return some (key, value) pair
as a 2-tuple; but raise KeyError if D is empty.

imod.mf6.Modflow6Simulation.regrid_like
=======================================
This method creates a new simulation object. The models contained in the new simulation are regridded versions
of the models in the input object (this).
Time discretization and solver settings are copied.

Parameters
----------
regridded_simulation_name: str
    name given to the output simulation
target_grid: xr.DataArray or  xu.UgridDataArray
    discretization onto which the models  in this simulation will be regridded
validate: bool
    set to true to validate the regridded packages

Returns
-------
a new simulation object with regridded models

imod.mf6.Modflow6Simulation.render
==================================
Renders simulation namefile

imod.mf6.Modflow6Simulation.run
===============================
Run Modflow 6 simulation. This method runs a subprocess calling
``mf6path``. This argument is set to ``mf6``, which means the Modflow 6
executable is expected to be added to your PATH environment variable.
:doc:`See this writeup how to add Modflow 6 to your PATH on Windows </examples/mf6/index>`

Note that the ``write`` method needs to be called before this method is
called.

Parameters
----------
mf6path: Union[str, Path]
    Path to the Modflow 6 executable. Defaults to calling ``mf6``.

Examples
--------
Make sure you write your model first

>>> simulation.write(path/to/model)
>>> simulation.run()

imod.mf6.Modflow6Simulation.setdefault
======================================
D.setdefault(k[,d]) -> D.get(k,d), also set D[k]=d if k not in D

imod.mf6.Modflow6Simulation.split
=================================
Split a simulation in different partitions using a submodel_labels array.

The submodel_labels array defines how a simulation will be split. The array should have the same topology as
the domain being split i.e. similar shape as a layer in the domain. The values in the array indicate to
which partition a cell belongs. The values should be zero or greater.

The method return a new simulation containing all the split models and packages

imod.mf6.Modflow6Simulation.update
==================================
D.update([E, ]**F) -> None.  Update D from mapping/iterable E and F.
If E present and has a .keys() method, does:     for k in E: D[k] = E[k]
If E present and lacks .keys() method, does:     for (k, v) in E: D[k] = v
In either case, this is followed by: for k, v in F.items(): D[k] = v

imod.mf6.Modflow6Simulation.values
==================================
D.values() -> an object providing a view on D's values

imod.mf6.Modflow6Simulation.write
=================================
Write Modflow6 simulation, including assigned groundwater flow and
transport models.

Parameters
----------
directory: str, pathlib.Path
    Directory to write Modflow 6 simulation to.
binary: ({True, False}, optional)
    Whether to write time-dependent input for stress packages as binary
    files, which are smaller in size, or more human-readable text files.
validate: ({True, False}, optional)
    Whether to validate the Modflow6 simulation, including models, at
    write. If True, erronous model input will throw a
    ``ValidationError``.
absolute_paths: ({True, False}, optional)
    True if all paths written to the mf6 inputfiles should be absolute.

imod.mf6.GroundwaterFlowModel
=============================
The GroundwaterFlowModel (GWF) simulates flow of (liquid) groundwater.
More information can be found here:
https://water.usgs.gov/water-resources/software/MODFLOW-6/mf6io_6.4.2.pdf#page=27

Parameters
----------

listing_file: Optional[str] = None
    name of the listing file to create for this GWF model. If not specified,
    then the name of the list file will be the basename of the GWF model
    name file and the 'lst' extension.
print_input: bool = False
    keyword to indicate that the list of all model stress package
    information will be written to the listing file immediately after it is
    read.
print_flows: bool = False
    keyword to indicate that the list of all model package flow rates will
    be printed to the listing file for every stress period time step in
    which "BUDGET PRINT" is specified in Output Control.
save_flows: bool = False
    indicate that all model package flow terms will be written to the file
    specified with "BUDGET FILEOUT" in Output Control.
newton: bool = False
    activates the Newton-Raphson formulation for groundwater flow between
    connected, convertible groundwater cells and stress packages that
    support calculation of Newton-Raphson terms for groundwater exchanges.
under_relaxation: bool = False,
    indicates whether the groundwater head in a cell will be under-relaxed when
    water levels fall below the bottom of the model below any given cell. By
    default, Newton-Raphson UNDER_RELAXATION is not applied.

imod.mf6.GroundwaterFlowModel Class Members
===========================================
   * imod.mf6.GroundwaterFlowModel.clear
   * imod.mf6.GroundwaterFlowModel.clip_box
   * imod.mf6.GroundwaterFlowModel.dump
   * imod.mf6.GroundwaterFlowModel.get
   * imod.mf6.GroundwaterFlowModel.is_clipping_supported
   * imod.mf6.GroundwaterFlowModel.is_regridding_supported
   * imod.mf6.GroundwaterFlowModel.is_splitting_supported
   * imod.mf6.GroundwaterFlowModel.items
   * imod.mf6.GroundwaterFlowModel.keys
   * imod.mf6.GroundwaterFlowModel.mask_all_packages
   * imod.mf6.GroundwaterFlowModel.pop
   * imod.mf6.GroundwaterFlowModel.popitem
   * imod.mf6.GroundwaterFlowModel.purge_empty_packages
   * imod.mf6.GroundwaterFlowModel.regrid_like
   * imod.mf6.GroundwaterFlowModel.setdefault
   * imod.mf6.GroundwaterFlowModel.update
   * imod.mf6.GroundwaterFlowModel.update_buoyancy_package
   * imod.mf6.GroundwaterFlowModel.values
   * imod.mf6.GroundwaterFlowModel.write

imod.mf6.GroundwaterFlowModel.clip_box
======================================
Clip a model by a bounding box (time, layer, y, x).

Slicing intervals may be half-bounded, by providing None:

* To select 500.0 <= x <= 1000.0:
  ``clip_box(x_min=500.0, x_max=1000.0)``.
* To select x <= 1000.0: ``clip_box(x_min=None, x_max=1000.0)``
  or ``clip_box(x_max=1000.0)``.
* To select x >= 500.0: ``clip_box(x_min = 500.0, x_max=None.0)``
  or ``clip_box(x_min=1000.0)``.

Parameters
----------
time_min: optional
time_max: optional
layer_min: optional, int
layer_max: optional, int
x_min: optional, float
x_max: optional, float
y_min: optional, float
y_max: optional, float
state_for_boundary: optional, float

imod.mf6.GroundwaterFlowModel.dump
==================================
Dump simulation to files. Writes a model definition as .TOML file, which
points to data for each package. Each package is stored as a separate
NetCDF. Structured grids are saved as regular NetCDFs, unstructured
grids are saved as UGRID NetCDF. Structured grids are always made GDAL
compliant, unstructured grids can be made MDAL compliant optionally.

Parameters
----------
directory: str or Path
    directory to dump simulation into.
modelname: str
    modelname, will be used to create a subdirectory.
validate: bool, optional
    Whether to validate simulation data. Defaults to True.
mdal_compliant: bool, optional
    Convert data with
    :func:`imod.prepare.spatial.mdal_compliant_ugrid2d` to MDAL
    compliant unstructured grids. Defaults to False.
crs: Any, optional
    Anything accepted by rasterio.crs.CRS.from_user_input
    Requires ``rioxarray`` installed.

imod.mf6.GroundwaterFlowModel.is_clipping_supported
===================================================
Returns True if all the packages in the model supports clipping. If one
of the packages in the model does not support clipping, it returns the
name of the first one.

imod.mf6.GroundwaterFlowModel.is_regridding_supported
=====================================================
Returns True if all the packages in the model supports regridding. If one
of the packages in the model does not support regridding, it returns the
name of the first one.

imod.mf6.GroundwaterFlowModel.is_splitting_supported
====================================================
Returns True if all the packages in the model supports splitting. If one
of the packages in the model does not support splitting, it returns the
name of the first one.

imod.mf6.GroundwaterFlowModel.mask_all_packages
===============================================
This function applies a mask to all packages in a model. The mask must
be presented as an idomain-like integer array that has 0 (inactive) or
-1 (vertical passthrough) values in filtered cells and 1 in active
cells.
Masking will overwrite idomain with the mask where the mask is 0 or -1.
Where the mask is 1, the original value of idomain will be kept. Masking
will update the packages accordingly, blanking their input where needed,
and is therefore not a reversible operation.

Parameters
----------
mask: xr.DataArray, xu.UgridDataArray of ints
    idomain-like integer array. 1 sets cells to active, 0 sets cells to inactive,
    -1 sets cells to vertical passthrough

imod.mf6.GroundwaterFlowModel.purge_empty_packages
==================================================
This function removes empty packages from the model.

imod.mf6.GroundwaterFlowModel.regrid_like
=========================================
Creates a model by regridding the packages of this model to another discretization.
It regrids all the arrays in the package using the default regridding methods.
At the moment only regridding to a different planar grid is supported, meaning
``target_grid`` has different ``"x"`` and ``"y"`` or different ``cell2d`` coords.

Parameters
----------
target_grid: xr.DataArray or xu.UgridDataArray
    a grid defined over the same discretization as the one we want to regrid the package to
validate: bool
    set to true to validate the regridded packages
regrid_context: Optional RegridderWeightsCache
    stores regridder weights for different regridders. Can be used to speed up regridding,
    if the same regridders are used several times for regridding different arrays.

Returns
-------
a model with similar packages to the input model, and with all the data-arrays regridded to another discretization,
similar to the one used in input argument "target_grid"

imod.mf6.GroundwaterFlowModel.update
====================================
D.update([E, ]**F) -> None.  Update D from mapping/iterable E and F.
If E present and has a .keys() method, does:     for k in E: D[k] = E[k]
If E present and lacks .keys() method, does:     for (k, v) in E: D[k] = v
In either case, this is followed by: for k, v in F.items(): D[k] = v

imod.mf6.GroundwaterFlowModel.update_buoyancy_package
=====================================================
If the simulation is partitioned, then the buoyancy package, if present,
must be updated for the renamed transport models.

imod.mf6.GroundwaterFlowModel.write
===================================
Write model namefile
Write packages

imod.mf6.GroundwaterTransportModel
==================================
The GroundwaterTransportModel (GWT) simulates transport of a single solute
species flowing in groundwater.
More information can be found here:
https://water.usgs.gov/water-resources/software/MODFLOW-6/mf6io_6.4.2.pdf#page=172

Parameters
----------

listing_file: Optional[str] = None
    name of the listing file to create for this GWT model. If not specified,
    then the name of the list file will be the basename of the GWT model
    name file and the 'lst' extension.
print_input: bool = False
    if True, indicates that the list of exchange entries will be echoed to
    the listing file immediately after it is read.
print_flows: bool = False
    if True, indicates that the list of exchange flow rates will be printed
    to the listing file for every stress period in which "SAVE BUDGET" is
    specified in Output Control
save_flows: bool = False,
    if True, indicates that all model package flow terms will be written to
    the file specified with "BUDGET FILEOUT" in Output Control.

imod.mf6.GroundwaterTransportModel Class Members
================================================
   * imod.mf6.GroundwaterTransportModel.clear
   * imod.mf6.GroundwaterTransportModel.clip_box
   * imod.mf6.GroundwaterTransportModel.dump
   * imod.mf6.GroundwaterTransportModel.get
   * imod.mf6.GroundwaterTransportModel.is_clipping_supported
   * imod.mf6.GroundwaterTransportModel.is_regridding_supported
   * imod.mf6.GroundwaterTransportModel.is_splitting_supported
   * imod.mf6.GroundwaterTransportModel.items
   * imod.mf6.GroundwaterTransportModel.keys
   * imod.mf6.GroundwaterTransportModel.mask_all_packages
   * imod.mf6.GroundwaterTransportModel.pop
   * imod.mf6.GroundwaterTransportModel.popitem
   * imod.mf6.GroundwaterTransportModel.purge_empty_packages
   * imod.mf6.GroundwaterTransportModel.regrid_like
   * imod.mf6.GroundwaterTransportModel.setdefault
   * imod.mf6.GroundwaterTransportModel.update
   * imod.mf6.GroundwaterTransportModel.values
   * imod.mf6.GroundwaterTransportModel.write

imod.mf6.GroundwaterTransportModel.clip_box
===========================================
Clip a model by a bounding box (time, layer, y, x).

Slicing intervals may be half-bounded, by providing None:

* To select 500.0 <= x <= 1000.0:
  ``clip_box(x_min=500.0, x_max=1000.0)``.
* To select x <= 1000.0: ``clip_box(x_min=None, x_max=1000.0)``
  or ``clip_box(x_max=1000.0)``.
* To select x >= 500.0: ``clip_box(x_min = 500.0, x_max=None.0)``
  or ``clip_box(x_min=1000.0)``.

Parameters
----------
time_min: optional
time_max: optional
layer_min: optional, int
layer_max: optional, int
x_min: optional, float
x_max: optional, float
y_min: optional, float
y_max: optional, float
state_for_boundary: optional, float

imod.mf6.StructuredDiscretization
=================================
Discretization information for structered grids is specified using the file.
(DIS6) Only one discretization input file (DISU6, DISV6 or DIS6) can be
specified for a model.
https://water.usgs.gov/water-resources/software/MODFLOW-6/mf6io_6.0.4.pdf#page=35

Parameters
----------
top: array of floats (xr.DataArray)
    is the top elevation for each cell in the top model layer.
bottom: array of floats (xr.DataArray)
    is the bottom elevation for each cell.
idomain: array of integers (xr.DataArray)
    Indicates the existence status of a cell. Horizontal discretization
    information will be derived from the x and y coordinates of the
    DataArray. If the idomain value for a cell is 0, the cell does not exist
    in the simulation. Input and output values will be read and written for
    the cell, but internal to the program, the cell is excluded from the
    solution. If the idomain value for a cell is 1, the cell exists in the
    simulation. if the idomain value for a cell is -1, the cell does not
    exist in the simulation. Furthermore, the first existing cell above will
    be connected to the first existing cell below. This type of cell is
    referred to as a "vertical pass through" cell.
validate: {True, False}
    Flag to indicate whether the package should be validated upon
    initialization. This raises a ValidationError if package input is
    provided in the wrong manner. Defaults to True.

imod.mf6.StructuredDiscretization Class Members
===============================================
   * imod.mf6.StructuredDiscretization.clip_box
   * imod.mf6.StructuredDiscretization.from_file
   * imod.mf6.StructuredDiscretization.get_non_grid_data
   * imod.mf6.StructuredDiscretization.is_empty
   * imod.mf6.StructuredDiscretization.mask
   * imod.mf6.StructuredDiscretization.regrid_like
   * imod.mf6.StructuredDiscretization.to_netcdf

imod.mf6.StructuredDiscretization.clip_box
==========================================
Clip a package by a bounding box (time, layer, y, x).

Slicing intervals may be half-bounded, by providing None:

* To select 500.0 <= x <= 1000.0:
  ``clip_box(x_min=500.0, x_max=1000.0)``.
* To select x <= 1000.0: ``clip_box(x_min=None, x_max=1000.0)``
  or ``clip_box(x_max=1000.0)``.
* To select x >= 500.0: ``clip_box(x_min = 500.0, x_max=None.0)``
  or ``clip_box(x_min=1000.0)``.

Parameters
----------
time_min: optional
time_max: optional
layer_min: optional, int
layer_max: optional, int
x_min: optional, float
x_max: optional, float
y_min: optional, float
y_max: optional, float
top: optional, GridDataArray
bottom: optional, GridDataArray
state_for_boundary: optional, GridDataArray


Returns
-------
clipped: Package

imod.mf6.StructuredDiscretization.from_file
===========================================
Loads an imod mf6 package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.mf6.Package.dataset.to_netcdf(),
as the checks upon package initialization are not done again!

Parameters
----------
path : str, pathlib.Path
    Path to the file.
**kwargs : keyword arguments
    Arbitrary keyword arguments forwarded to ``xarray.open_dataset()``, or
    ``xarray.open_zarr()``.
Refer to the examples.

Returns
-------
package : imod.mf6.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.mf6.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.mf6.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.mf6.StructuredDiscretization.get_non_grid_data
===================================================
This function copies the attributes of a dataset that are scalars, such as options.

parameters
----------
grid_names: list of str
    the names of the attribbutes of a dataset that are grids.

imod.mf6.StructuredDiscretization.is_empty
==========================================
Returns True if the package is empty- for example if it contains only no-data values.

imod.mf6.StructuredDiscretization.mask
======================================
Mask values outside of domain.

Floating values outside of the condition are set to NaN (nodata).
Integer values outside of the condition are set to 0 (inactive in
MODFLOW terms).

Parameters
----------
mask: xr.DataArray, xu.UgridDataArray of ints
    idomain-like integer array. 1 sets cells to active, 0 sets cells to inactive,
    -1 sets cells to vertical passthrough

Returns
-------
masked: Package
    The package with part masked.

imod.mf6.StructuredDiscretization.regrid_like
=============================================
Creates a package of the same type as this package, based on another
discretization. It regrids all the arrays in this package to the desired
discretization, and leaves the options unmodified. At the moment only
regridding to a different planar grid is supported, meaning
``target_grid`` has different ``"x"`` and ``"y"`` or different
``cell2d`` coords.

The default regridding methods are specified in the ``_regrid_method``
attribute of the package. These defaults can be overridden using the
input parameters of this function.

Examples
--------
To regrid the npf package with a non-default method for the k-field, call regrid_like with these arguments:

>>> regridder_types = imod.mf6.regrid.NodePropertyFlowRegridMethod(k=(imod.RegridderType.OVERLAP, "mean"))
>>> new_npf = npf.regrid_like(like,  RegridderWeightsCache, regridder_types)


Parameters
----------
target_grid: xr.DataArray or xu.UgridDataArray
    a grid defined over the same discretization as the one we want to regrid the package to.
regrid_context: RegridderWeightsCache, optional
    stores regridder weights for different regridders. Can be used to speed up regridding,
    if the same regridders are used several times for regridding different arrays.
regridder_types: RegridMethodType, optional
    dictionary mapping arraynames (str) to a tuple of regrid type (a specialization class of BaseRegridder) and function name (str)
    this dictionary can be used to override the default mapping method.

Returns
-------
a package with the same options as this package, and with all the data-arrays regridded to another discretization,
similar to the one used in input argument "target_grid"

imod.mf6.StructuredDiscretization.to_netcdf
===========================================
Write dataset contents to a netCDF file.
Custom encoding rules can be provided on package level by overriding the _netcdf_encoding in the package

imod.mf6.VerticesDiscretization
===============================
Discretization by Vertices (DISV).

Parameters
----------
top: array of floats (xu.UgridDataArray)
bottom: array of floats (xu.UgridDataArray)
idomain: array of integers (xu.UgridDataArray)
validate: {True, False}
    Flag to indicate whether the package should be validated upon
    initialization. This raises a ValidationError if package input is
    provided in the wrong manner. Defaults to True.

imod.mf6.VerticesDiscretization Class Members
=============================================
   * imod.mf6.VerticesDiscretization.clip_box
   * imod.mf6.VerticesDiscretization.from_file
   * imod.mf6.VerticesDiscretization.get_non_grid_data
   * imod.mf6.VerticesDiscretization.is_empty
   * imod.mf6.VerticesDiscretization.mask
   * imod.mf6.VerticesDiscretization.regrid_like
   * imod.mf6.VerticesDiscretization.to_netcdf

imod.mf6.VerticesDiscretization.from_file
=========================================
Loads an imod mf6 package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.mf6.Package.dataset.to_netcdf(),
as the checks upon package initialization are not done again!

Parameters
----------
path : str, pathlib.Path
    Path to the file.
**kwargs : keyword arguments
    Arbitrary keyword arguments forwarded to ``xarray.open_dataset()``, or
    ``xarray.open_zarr()``.
Refer to the examples.

Returns
-------
package : imod.mf6.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.mf6.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.mf6.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.mf6.TimeDiscretization
===========================
Timing for all models of the simulation is controlled by the Temporal
Discretization (TDIS) Package.
https://water.usgs.gov/water-resources/software/MODFLOW-6/mf6io_6.0.4.pdf#page=17

Parameters
----------
timestep_duration: float
    is the length of a stress period. (PERLEN)
n_timesteps: int, optional
    is the number of time steps in a stress period (nstp).
    Default value: 1
timestep_multiplier: float, optional
    is the multiplier for the length of successive time steps. The length of
    a time step is calculated by multiplying the length of the previous time
    step by timestep_multiplier (TSMULT).
    Default value: 1.0
validate: {True, False}
    Flag to indicate whether the package should be validated upon
    initialization. This raises a ValidationError if package input is
    provided in the wrong manner. Defaults to True.

imod.mf6.TimeDiscretization Class Members
=========================================
   * imod.mf6.TimeDiscretization.clip_box
   * imod.mf6.TimeDiscretization.from_file
   * imod.mf6.TimeDiscretization.get_non_grid_data
   * imod.mf6.TimeDiscretization.is_empty
   * imod.mf6.TimeDiscretization.mask
   * imod.mf6.TimeDiscretization.regrid_like
   * imod.mf6.TimeDiscretization.to_netcdf

imod.mf6.TimeDiscretization.from_file
=====================================
Loads an imod mf6 package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.mf6.Package.dataset.to_netcdf(),
as the checks upon package initialization are not done again!

Parameters
----------
path : str, pathlib.Path
    Path to the file.
**kwargs : keyword arguments
    Arbitrary keyword arguments forwarded to ``xarray.open_dataset()``, or
    ``xarray.open_zarr()``.
Refer to the examples.

Returns
-------
package : imod.mf6.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.mf6.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.mf6.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.mf6.OutputControl
======================
The Output Control Option determines how and when heads, budgets and/or
concentrations are printed to the listing file and/or written to a separate
binary output file.
https://water.usgs.gov/water-resources/software/MODFLOW-6/mf6io_6.4.2.pdf#page=53

Currently the settings "first", "last", "all", and "frequency"
are supported, the "steps" setting is not supported, because of
its ragged nature. Furthermore, only one setting per stress period
can be specified in imod-python.

Parameters
----------
save_head : {string, integer}, or xr.DataArray of {string, integer}, optional
    String or integer indicating output control for head file (.hds)
    If string, should be one of ["first", "last", "all"].
    If integer, interpreted as frequency.
save_budget : {string, integer}, or xr.DataArray of {string, integer}, optional
    String or integer indicating output control for cell budgets (.cbc)
    If string, should be one of ["first", "last", "all"].
    If integer, interpreted as frequency.
save_concentration : {string, integer}, or xr.DataArray of {string, integer}, optional
    String or integer indicating output control for concentration file (.ucn)
    If string, should be one of ["first", "last", "all"].
    If integer, interpreted as frequency.
validate: {True, False}
    Flag to indicate whether the package should be validated upon
    initialization. This raises a ValidationError if package input is
    provided in the wrong manner. Defaults to True.

Examples
--------
To specify a mix of both 'frequency' and 'first' setting,
we need to specify an array with both integers and strings.
For this we need to create a numpy object array first,
otherwise xarray converts all to strings automatically.

>>> time = [np.datetime64("2000-01-01"), np.datetime64("2000-01-02")]
>>> data = np.array(["last", 5], dtype="object")
>>> save_head = xr.DataArray(data, coords={"time": time}, dims=("time"))
>>> oc = imod.mf6.OutputControl(save_head=save_head, save_budget=None, save_concentration=None)

imod.mf6.OutputControl Class Members
====================================
   * imod.mf6.OutputControl.clip_box
   * imod.mf6.OutputControl.from_file
   * imod.mf6.OutputControl.get_non_grid_data
   * imod.mf6.OutputControl.is_empty
   * imod.mf6.OutputControl.mask
   * imod.mf6.OutputControl.regrid_like
   * imod.mf6.OutputControl.to_netcdf

imod.mf6.OutputControl.from_file
================================
Loads an imod mf6 package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.mf6.Package.dataset.to_netcdf(),
as the checks upon package initialization are not done again!

Parameters
----------
path : str, pathlib.Path
    Path to the file.
**kwargs : keyword arguments
    Arbitrary keyword arguments forwarded to ``xarray.open_dataset()``, or
    ``xarray.open_zarr()``.
Refer to the examples.

Returns
-------
package : imod.mf6.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.mf6.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.mf6.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.mf6.Solution
=================
Iterative Model Solution.
The model solution will solve all of the models that are added to it, as
specified in the simulation name file, and will include Numerical Exchanges,
if they are present. The iterative model solution requires specification of
both nonlinear and linear settings.
https://water.usgs.gov/water-resources/software/MODFLOW-6/mf6io_6.0.4.pdf#page=147

Three predifined solutions settings are available: SolutionPresetSimple,
SolutionPresetModerate and SolutionPresetComplex. When using one of the
predefined solutions only the print_option, csv_output, and no_ptc have to
be defined. The default values for each are described below.

Parameters
----------
modelnames: list of str
    Which models to solve in this solution. Only models of the same type
    (GWF or GWT) should be added to the same solution.
outer_dvclose: float
    real value defining the head change criterion for convergence of the
    outer (nonlinear) iterations, in units of length. When the maximum
    absolute value of the head change at all nodes during an iteration is
    less than or equal to outer_dvclose, iteration stops. Commonly,
    outer_dvclose equals 0.01.
    SolutionPresetSimple: 0.001
    SolutionPresetModerate: 0.01
    SolutionPresetComplex: 0.1
outer_maximum: int
    integer value defining the maximum number of outer (nonlinear)
    iterations - that is, calls to the solution routine. For a linear
    problem outer_maximum should be 1.
    SolutionPresetSimple: 25
    SolutionPresetModerate: 50
    SolutionPresetComplex: 100
inner_maximum: int
    integer value defining the maximum number of inner (linear) iterations.
    The number typically depends on the characteristics of the matrix
    solution scheme being used. For nonlinear problems, inner_maximum
    usually ranges from 60 to 600; a value of 100 will be sufficient for
    most linear problems.
    SolutionPresetSimple: 50
    SolutionPresetModerate: 100
    SolutionPresetComplex: 500
inner_dvclose: float
    real value defining the head change criterion for convergence of the
    inner (linear) iterations, in units of length. When the maximum absolute
    value of the head change at all nodes during an iteration is less than
    or equal to inner_dvclose, the matrix solver assumes convergence.
    Commonly, inner_dvclose is set an order of magnitude less than the
    outer_dvclose value.
    SolutionPresetSimple: 0.001
    SolutionPresetModerate: 0.01
    SolutionPresetComplex: 0.1
inner_rclose: float
    real value that defines the flow residual tolerance for convergence of
    the IMS linear solver and specific flow residual criteria used. This
    value represents the maximum allowable residual at any single node.
    Value is in units of length cubed per time, and must be consistent with
    MODFLOW 6 length and time units. Usually a value of 1.0 × 10−1 is
    sufficient for the flow-residual criteria when meters and seconds are
    the defined MODFLOW 6 length and time.
    SolutionPresetSimple: 0.1
    SolutionPresetModerate: 0.1
    SolutionPresetComplex: 0.1
linear_acceleration: str
    options: {"cg", "bicgstab"}
    a keyword that defines the linear acceleration method used by the
    default IMS linear solvers. CG - preconditioned conjugate gradient
    method. BICGSTAB - preconditioned bi-conjugate gradient stabilized
    method.
    SolutionPresetSimple: "cg"
    SolutionPresetModerate: "bicgstab"
    SolutionPresetComplex: "bicgstab"
under_relaxation: str, optional
    options: {None, "simple", "cooley", "dbd"}
    is an optional keyword that defines the nonlinear relative_rclose
    schemes used. By default under_relaxation is not used.
    None - relative_rclose is not used.
    simple - Simple relative_rclose scheme with a fixed relaxation factor is
    used.
    cooley - Cooley relative_rclose scheme is used.
    dbd - delta-bar-delta relative_rclose is used.
    Note that the relative_rclose schemes are used in conjunction with
    problems that use the Newton-Raphson formulation, however, experience
    has indicated that the Cooley relative_rclose and damping work well also
    for the Picard scheme with the wet/dry options of MODFLOW 6.
    Default value: None
    SolutionPresetSimple: None
    SolutionPresetModerate: "dbd"
    SolutionPresetComplex: "dbd"
under_relaxation_theta: float, optional
    real value defining the reduction factor for the learning rate
    (underrelaxation term) of the delta-bar-delta algorithm. The value of
    under relaxation theta is between zero and one. If the change in the
    variable (head) is of opposite sign to that of the previous iteration,
    the relative_rclose term is reduced by a factor of under relaxation
    theta. The value usually ranges from 0.3 to 0.9; a value of 0.7 works
    well for most problems. under relaxation theta only needs to be
    specified if under relaxation is dbd.
    Default value: None
    SolutionPresetSimple: 0.0
    SolutionPresetModerate: 0.9
    SolutionPresetComplex: 0.8
under_relaxation_kappa: float, optional
    real value defining the increment for the learning rate (relative_rclose
    term) of the delta-bar-delta algorithm. The value of under relaxation
    kappa is between zero and one. If the change in the variable (head) is
    of the same sign to that of the previous iteration, the relative_rclose
    term is increased by an increment of under_relaxation_kappa. The value
    usually ranges from 0.03 to 0.3; a value of 0.1 works well for most
    problems. under relaxation kappa only needs to be specified if under
    relaxation is dbd.
    Default value: None
    SolutionPresetSimple: 0.0
    SolutionPresetModerate: 0.0001
    SolutionPresetComplex: 0.0001
under_relaxation_gamma: float, optional
    real value defining the history or memory term factor of the
    delta-bardelta algorithm. under relaxation gamma is between zero and 1
    but cannot be equal to one. When under relaxation gamma is zero, only
    the most recent history (previous iteration value) is maintained. As
    under relaxation gamma is increased, past history of iteration changes
    has greater influence on the memory term. The memory term is maintained
    as an exponential average of past changes. Retaining some past history
    can overcome granular behavior in the calculated function surface and
    therefore helps to overcome cyclic patterns of nonconvergence. The value
    usually ranges from 0.1 to 0.3; a value of 0.2 works well for most
    problems. under relaxation gamma only needs to be specified if under
    relaxation is not none.
    Default value: None
    SolutionPresetSimple: 0.0
    SolutionPresetModerate: 0.0
    SolutionPresetComplex: 0.0
under_relaxation_momentum: float, optional
    real value defining the fraction of past history changes that is added
    as a momentum term to the step change for a nonlinear iteration. The
    value of under relaxation momentum is between zero and one. A large
    momentum term should only be used when small learning rates are
    expected. Small amounts of the momentum term help convergence. The value
    usually ranges from 0.0001 to 0.1; a value of 0.001 works well for most
    problems. under relaxation momentum only needs to be specified if under
    relaxation is dbd.
    Default value: None
    SolutionPresetSimple: 0.0
    SolutionPresetModerate: 0.0
    SolutionPresetComplex: 0.0
backtracking_number: int, optional
    integer value defining the maximum number of backtracking iterations
    allowed for residual reduction computations. If backtracking number = 0
    then the backtracking iterations are omitted. The value usually ranges
    from 2 to 20; a value of 10 works well for most problems.
    Default value: None
    SolutionPresetSimple: 0
    SolutionPresetModerate: 0
    SolutionPresetComplex: 20
backtracking_tolerance: float, optional
    real value defining the tolerance for residual change that is allowed
    for residual reduction computations. backtracking tolerance should not
    be less than one to avoid getting stuck in local minima. A large value
    serves to check for extreme residual increases, while a low value serves
    to control step size more severely. The value usually ranges from 1.0 to
    106; a value of 104 works well for most problems but lower values like
    1.1 may be required for harder problems. backtracking tolerance only
    needs to be specified if backtracking_number is greater than zero.
    Default value: None
    SolutionPresetSimple: 0.0
    SolutionPresetModerate: 0.0
    SolutionPresetComplex: 1.05
backtracking_reduction_factor: float, optional
    real value defining the reduction in step size used for residual
    reduction computations. The value of backtracking reduction factor is
    between 142 MODFLOW 6 - Description of Input and Output zero and one.
    The value usually ranges from 0.1 to 0.3; a value of 0.2 works well for
    most problems. backtracking_reduction_factor only needs to be specified
    if backtracking number is greater than zero.
    Default value: None
    SolutionPresetSimple: 0.0
    SolutionPresetModerate: 0.0
    SolutionPresetComplex: 0.1
backtracking_residual_limit: float, optional
    real value defining the limit to which the residual is reduced with
    backtracking. If the residual is smaller than
    backtracking_residual_limit, then further backtracking is not performed.
    A value of 100 is suitable for large problems and residual reduction to
    smaller values may only slow down computations. backtracking residual
    limit only needs to be specified if backtracking_number is greater than
    zero.
    Default value: None
    SolutionPresetSimple: 0.0
    SolutionPresetModerate: 0.0
    SolutionPresetComplex: 0.002
rclose_option: str, optional
    options: {"strict", "l2norm_rclose", "relative_rclose"}
    an optional keyword that defines the specific flow residual criterion
    used.
    strict: an optional keyword that is used to specify that inner rclose
    represents a infinity-norm (absolute convergence criteria) and that the
    head and flow convergence criteria must be met on the first inner
    iteration (this criteria is equivalent to the criteria used by the
    MODFLOW-2005 PCG package (Hill, 1990)).
    l2norm_rclose: an optionalkeyword that is used to specify that inner
    rclose represents a l-2 norm closure criteria instead of a infinity-norm
    (absolute convergence criteria). When l2norm_rclose is specified, a
    reasonable initial inner rclose value is 0.1 times the number of active
    cells when meters and seconds are the defined MODFLOW 6 length and time.
    relative_rclose:  an optional keyword that is used to specify that
    inner_rclose represents a relative L-2 Norm reduction closure criteria
    instead of a infinity-Norm (absolute convergence criteria). When
    relative_rclose is specified, a reasonable initial inner_rclose value is
    1.0 * 10-4 and convergence is achieved for a given inner (linear)
    iteration when ∆h ≤ inner_dvclose and the current L-2 Norm is ≤ the
    product of the relativ_rclose and the initial L-2 Norm for the current
    inner (linear) iteration. If rclose_option is not specified, an absolute
    residual (infinity-norm) criterion is used.
    Default value: None
    SolutionPresetSimple: "strict"
    SolutionPresetModerate: "strict"
    SolutionPresetComplex: "strict"
relaxation_factor: float, optional
    optional real value that defines the relaxation factor used by the
    incomplete LU factorization preconditioners (MILU(0) and MILUT).
    relaxation_factor is unitless and should be greater than or equal to 0.0
    and less than or equal to 1.0. relaxation_factor Iterative Model
    Solution 143 values of about 1.0 are commonly used, and experience
    suggests that convergence can be optimized in some cases with relax
    values of 0.97. A relaxation_factor value of 0.0 will result in either
    ILU(0) or ILUT preconditioning (depending on the value specified for
    preconditioner_levels and/or preconditioner_drop_tolerance). By default,
    relaxation_factor is zero.
    Default value: None
    SolutionPresetSimple: 0.0
    SolutionPresetModerate: 0
    SolutionPresetComplex: 0.0
preconditioner_levels: int, optional
    optional integer value defining the level of fill for ILU decomposition
    used in the ILUT and MILUT preconditioners. Higher levels of fill
    provide more robustness but also require more memory. For optimal
    performance, it is suggested that a large level of fill be applied (7 or
    8) with use of a drop tolerance. Specification of a
    preconditioner_levels value greater than zero results in use of the ILUT
    preconditioner. By default, preconditioner_levels is zero and the
    zero-fill incomplete LU factorization preconditioners (ILU(0) and
    MILU(0)) are used.
    Default value: None
    SolutionPresetSimple: 0
    SolutionPresetModerate: 0
    SolutionPresetComplex: 5
preconditioner_drop_tolerance: float, optional
    optional real value that defines the drop tolerance used to drop
    preconditioner terms based on the magnitude of matrix entries in the
    ILUT and MILUT preconditioners. A value of 10−4 works well for most
    problems. By default, preconditioner_drop_tolerance is zero and the
    zero-fill incomplete LU factorization preconditioners (ILU(0) and
    MILU(0)) are used.
    Default value: None
    SolutionPresetSimple: 0
    SolutionPresetModerate: 0.0
    SolutionPresetComplex: 0.0001
number_orthogonalizations: int, optional
    optional integer value defining the interval used to explicitly
    recalculate the residual of the flow equation using the solver
    coefficient matrix, the latest head estimates, and the right hand side.
    For problems that benefit from explicit recalculation of the residual, a
    number between 4 and 10 is appropriate. By default,
    number_orthogonalizations is zero.
    Default value: None
    SolutionPresetSimple: 0
    SolutionPresetModerate: 0
    SolutionPresetComplex: 2
scaling_method: str
    options: {None, "diagonal", "l2norm"}
    an optional keyword that defines the matrix scaling approach used. By
    default, matrix scaling is not applied.
    None - no matrix scaling applied.
    diagonal - symmetric matrix scaling using the POLCG preconditioner
    scaling method in Hill (1992).
    l2norm - symmetric matrix scaling using the L2 norm.
    Default value: None
reordering_method: str
    options: {None, "rcm", "md"}
    an optional keyword that defines the matrix reordering approach used. By
    default, matrix reordering is not applied.
    None - original ordering.
    rcm - reverse Cuthill McKee ordering.
    md - minimum degree ordering
    Default value: None
print_option: str
    options: {"none", "summary", "all"}
    is a flag that controls printing of convergence information from the
    solver.
    None - means print nothing.
    summary - means print only the total
    number of iterations and nonlinear residual reduction summaries.
    all - means print linear matrix solver convergence information to the
    solution listing file and model specific linear matrix solver
    convergence information to each model listing file in addition to
    SUMMARY information.
    Default value: "summary"
outer_csvfile: str, optional
    None if no csv is to be written for the output, str of filename
    if csv is to be written.
    Default value: None
inner_csvfile: str, optional
    None if no csv is to be written for the output, str of filename
    if csv is to be written.
    Default value: None
no_ptc: str, optional, either None, "all", or "first".
    is a flag that is used to disable pseudo-transient continuation (PTC).
    Option only applies to steady-state stress periods for models using the
    Newton-Raphson formulation. For many problems, PTC can significantly
    improve convergence behavior for steady-state simulations, and for this
    reason it is active by default. In some cases, however, PTC can worsen
    the convergence behavior, especially when the initial conditions are
    similar to the solution. When the initial conditions are similar to, or
    exactly the same as, the solution and convergence is slow, then this NO
    PTC option should be used to deactivate PTC. This NO PTC option should
    also be used in order to compare convergence behavior with other MODFLOW
    versions, as PTC is only available in MODFLOW 6.
ats_outer_maximum_fraction: float, optional.
    real value defining the fraction of the maximum allowable outer iterations
    used with the Adaptive Time Step (ATS) capability if it is active. If this
    value is set to zero by the user, then this solution will have no effect
    on ATS behavior. This value must be greater than or equal to zero and less
    than or equal to 0.5 or the program will terminate with an error. If it is
    not specified by the user, then it is assigned a default value of one
    third. When the number of outer iterations for this solution is less than
    the product of this value and the maximum allowable outer iterations, then
    ATS will increase the time step length by a factor of DTADJ in the ATS
    input file. When the number of outer iterations for this solution is
    greater than the maximum allowable outer iterations minus the product of
    this value and the maximum allowable outer iterations, then the ATS (if
    active) will decrease the time step length by a factor of 1 / DTADJ.

    Default value is None.
validate: {True, False}
    Flag to indicate whether the package should be validated upon
    initialization. This raises a ValidationError if package input is
    provided in the wrong manner. Defaults to True.

imod.mf6.Solution Class Members
===============================
   * imod.mf6.Solution.clip_box
   * imod.mf6.Solution.from_file
   * imod.mf6.Solution.get_non_grid_data
   * imod.mf6.Solution.is_empty
   * imod.mf6.Solution.mask
   * imod.mf6.Solution.regrid_like
   * imod.mf6.Solution.to_netcdf

imod.mf6.Solution.from_file
===========================
Loads an imod mf6 package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.mf6.Package.dataset.to_netcdf(),
as the checks upon package initialization are not done again!

Parameters
----------
path : str, pathlib.Path
    Path to the file.
**kwargs : keyword arguments
    Arbitrary keyword arguments forwarded to ``xarray.open_dataset()``, or
    ``xarray.open_zarr()``.
Refer to the examples.

Returns
-------
package : imod.mf6.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.mf6.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.mf6.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.mf6.ApiPackage
===================
The API package can be used to interact with the dll-version of Modflow (libMF6.dll).
libMF6 implement the XMI api, which can be used among other things to get
or set internal MF6 arrays mid-simulation.
For more information, see:

https://doi.org/10.1016/j.envsoft.2021.105257
https://modflow6.readthedocs.io/en/stable/_mf6io/gwf-api.html
https://modflow6.readthedocs.io/en/stable/_mf6io/gwt-api.html

Parameters
----------
maxbound: int
    The number of cells for which information will be queried or set with api calls.
print_input: ({True, False}, optional)
    keyword to indicate that the list of constant head information will
    be written to the listing file immediately after it is read. Default is
    False.
print_flows: ({True, False}, optional)
    Indicates that the list of constant head flow rates will be printed to
    the listing file for every stress period time step in which "BUDGET
    PRINT" is specified in Output Control. If there is no Output Control
    option and PRINT FLOWS is specified, then flow rates are printed for the
    last time step of each stress period.
    Default is False.
save_flows: ({True, False}, optional)
    Indicates that constant head flow terms will be written to the file
    specified with "BUDGET FILEOUT" in Output Control. Default is False.

.. note::
   This package can be added to both flow and transport models.

imod.mf6.ApiPackage Class Members
=================================
   * imod.mf6.ApiPackage.clip_box
   * imod.mf6.ApiPackage.from_file
   * imod.mf6.ApiPackage.get_non_grid_data
   * imod.mf6.ApiPackage.is_empty
   * imod.mf6.ApiPackage.mask
   * imod.mf6.ApiPackage.regrid_like
   * imod.mf6.ApiPackage.to_netcdf

imod.mf6.ApiPackage.from_file
=============================
Loads an imod mf6 package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.mf6.Package.dataset.to_netcdf(),
as the checks upon package initialization are not done again!

Parameters
----------
path : str, pathlib.Path
    Path to the file.
**kwargs : keyword arguments
    Arbitrary keyword arguments forwarded to ``xarray.open_dataset()``, or
    ``xarray.open_zarr()``.
Refer to the examples.

Returns
-------
package : imod.mf6.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.mf6.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.mf6.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.mf6.Buoyancy
=================
Buoyancy package. This package must be included when performing variable
density simulation.

Note that ``reference_density`` is a single value, but
``density_concentration_slope``, ``reference_concentration`` and
``modelname`` require an entry for every active species. Please refer to the
examples.

Parameters
----------
reference_density: float,
    fluid reference density used in the equation of state.
density_concentration_slope: sequence of floats
    Slope of the (linear) density concentration line used in the density
    equation of state.
reference_concentration: sequence of floats
    Reference concentration used in the density equation of state.
modelname: sequence of strings,
    Name of the GroundwaterTransport (GWT) model used for the
    concentrations.
species: sequence of str,
    Name of the species used to calculate a density value.
hhformulation_rhs: bool, optional.
    use the variable-density hydraulic head formulation and add off-diagonal
    terms to the right-hand. This option will prevent the BUY Package from
    adding asymmetric terms to the flow matrix. Default value is ``False``.
densityfile:
    name of the binary output file to write density information. The density
    file has the same format as the head file. Density values will be
    written to the density file whenever heads are written to the binary
    head file. The settings for controlling head output are contained in the
    Output Control option.
validate: {True, False}
    Flag to indicate whether the package should be validated upon
    initialization. This raises a ValidationError if package input is
    provided in the wrong manner. Defaults to True.

Examples
--------

The buoyancy input for a single species called "salinity", which is
simulated by a GWT model called "gwt-1" are specified as follows:

>>> buy = imod.mf6.Buoyance(
...     reference_density=1000.0,
...     density_concentration_slope=[0.025],
...     reference_concentration=[0.0],
...     modelname=["gwt-1"],
...     species=["salinity"],
... )

Multiple species can be specified by presenting multiple values with an
associated species coordinate. Two species, called "c1" and "c2", simulated
by the GWT models "gwt-1" and "gwt-2" are specified as:

>>> coords = {"species": ["c1", "c2"]}
>>> buy = imod.mf6.Buoyance(
...     reference_density=1000.0,
...     density_concentration_slope=[0.025, 0.01],
...     reference_concentration=[0.0, 0.0],
...     modelname=["gwt-1", "gwt-2"],
...     species=["c1", "c2"],
... )

imod.mf6.Buoyancy Class Members
===============================
   * imod.mf6.Buoyancy.clip_box
   * imod.mf6.Buoyancy.from_file
   * imod.mf6.Buoyancy.get_non_grid_data
   * imod.mf6.Buoyancy.get_transport_model_names
   * imod.mf6.Buoyancy.is_empty
   * imod.mf6.Buoyancy.mask
   * imod.mf6.Buoyancy.regrid_like
   * imod.mf6.Buoyancy.to_netcdf
   * imod.mf6.Buoyancy.update_transport_models

imod.mf6.Buoyancy.from_file
===========================
Loads an imod mf6 package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.mf6.Package.dataset.to_netcdf(),
as the checks upon package initialization are not done again!

Parameters
----------
path : str, pathlib.Path
    Path to the file.
**kwargs : keyword arguments
    Arbitrary keyword arguments forwarded to ``xarray.open_dataset()``, or
    ``xarray.open_zarr()``.
Refer to the examples.

Returns
-------
package : imod.mf6.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.mf6.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.mf6.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.mf6.Buoyancy.get_transport_model_names
===========================================
Returns the names of the transport  models used by this buoyancy package.

imod.mf6.Buoyancy.update_transport_models
=========================================
The names of the transport models can change in some cases, for example
when partitioning. Use this function to update the names of the
transport models.

imod.mf6.ConstantHead
=====================
Constant-Head package. Any number of CHD Packages can be specified for a
single groundwater flow model; however, an error will occur if a CHD Package
attempts to make a GWF cell a constant-head cell when that cell has already
been designated as a constant-head cell either within the present CHD
Package or within another CHD Package. In previous MODFLOW versions, it was
not possible to convert a constant-head cell to an active cell. Once a cell
was designated as a constant-head cell, it remained a constant-head cell
until the end of the end of the simulation. In MODFLOW 6 a constant-head
cell will become active again if it is not included as a constant-head cell
in subsequent stress periods. Previous MODFLOW versions allowed
specification of SHEAD and EHEAD, which were the starting and ending
prescribed heads for a stress period. Linear interpolation was used to
calculate a value for each time step. In MODFLOW 6 only a single head value
can be specified for any constant-head cell in any stress period. The
time-series functionality must be used in order to interpolate values to
individual time steps.

Parameters
----------
head: array of floats (xr.DataArray)
    Is the head at the boundary.
print_input: ({True, False}, optional)
    keyword to indicate that the list of constant head information will
    be written to the listing file immediately after it is read. Default is
    False.
concentration: array of floats (xr.DataArray, optional)
    if this flow package is used in simulations also involving transport, then this array is used
    as the  concentration for inflow over this boundary.
concentration_boundary_type: ({"AUX", "AUXMIXED"}, optional)
    if this flow package is used in simulations also involving transport, then this keyword specifies
    how outflow over this boundary is computed.
print_flows: ({True, False}, optional)
    Indicates that the list of constant head flow rates will be printed to
    the listing file for every stress period time step in which "BUDGET
    PRINT" is specified in Output Control. If there is no Output Control
    option and PRINT FLOWS is specified, then flow rates are printed for the
    last time step of each stress period.
    Default is False.
save_flows: ({True, False}, optional)
    Indicates that constant head flow terms will be written to the file
    specified with "BUDGET FILEOUT" in Output Control. Default is False.
observations: [Not yet supported.]
    Default is None.
validate: {True, False}
    Flag to indicate whether the package should be validated upon
    initialization. This raises a ValidationError if package input is
    provided in the wrong manner. Defaults to True.
repeat_stress: Optional[xr.DataArray] of datetimes
    Used to repeat data for e.g. repeating stress periods such as
    seasonality without duplicating the values. The DataArray should have
    dimensions ``("repeat", "repeat_items")``. The ``repeat_items``
    dimension should have size 2: the first value is the "key", the second
    value is the "value". For the "key" datetime, the data of the "value"
    datetime will be used. Can also be set with a dictionary using the
    ``set_repeat_stress`` method.

imod.mf6.ConstantHead Class Members
===================================
   * imod.mf6.ConstantHead.clip_box
   * imod.mf6.ConstantHead.from_file
   * imod.mf6.ConstantHead.get_non_grid_data
   * imod.mf6.ConstantHead.is_empty
   * imod.mf6.ConstantHead.mask
   * imod.mf6.ConstantHead.regrid_like
   * imod.mf6.ConstantHead.render
   * imod.mf6.ConstantHead.set_repeat_stress
   * imod.mf6.ConstantHead.to_netcdf
   * imod.mf6.ConstantHead.write

imod.mf6.ConstantHead.from_file
===============================
Loads an imod mf6 package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.mf6.Package.dataset.to_netcdf(),
as the checks upon package initialization are not done again!

Parameters
----------
path : str, pathlib.Path
    Path to the file.
**kwargs : keyword arguments
    Arbitrary keyword arguments forwarded to ``xarray.open_dataset()``, or
    ``xarray.open_zarr()``.
Refer to the examples.

Returns
-------
package : imod.mf6.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.mf6.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.mf6.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.mf6.ConstantHead.render
============================
Render fills in the template only, doesn't write binary data

imod.mf6.ConstantHead.set_repeat_stress
=======================================
Set repeat stresses: re-use data of earlier periods.

Parameters
----------
times: Dict of datetime-like to datetime-like.
    The data of the value datetime is used for the key datetime.

imod.mf6.ConstantHead.write
===========================
writes the blockfile and binary data

directory is modelname

imod.mf6.Drainage
=================
The Drain package is used to simulate head-dependent flux boundaries.
https://water.usgs.gov/ogw/modflow/mf6io.pdf#page=67

Parameters
----------
elevation: array of floats (xr.DataArray)
    elevation of the drain. (elev)
conductance: array of floats (xr.DataArray)
    is the conductance of the drain. (cond)
concentration: array of floats (xr.DataArray, optional)
    if this flow package is used in simulations also involving transport, then this array is used
    as the  concentration for inflow over this boundary.
concentration_boundary_type: ({"AUX", "AUXMIXED"}, optional)
    if this flow package is used in simulations also involving transport, then this keyword specifies
    how outflow over this boundary is computed.
print_input: ({True, False}, optional)
    keyword to indicate that the list of drain information will be written
    to the listing file immediately after it is read. Default is False.
print_flows: ({True, False}, optional)
    Indicates that the list of drain flow rates will be printed to the
    listing file for every stress period time step in which "BUDGET PRINT"
    is specified in Output Control. If there is no Output Control option and
    PRINT FLOWS is specified, then flow rates are printed for the last time
    step of each stress period.
    Default is False.
save_flows: ({True, False}, optional)
    Indicates that drain flow terms will be written to the file specified
    with "BUDGET FILEOUT" in Output Control. Default is False.
observations: [Not yet supported.]
    Default is None.
validate: {True, False}
    Flag to indicate whether the package should be validated upon
    initialization. This raises a ValidationError if package input is
    provided in the wrong manner. Defaults to True.
repeat_stress: Optional[xr.DataArray] of datetimes
    Used to repeat data for e.g. repeating stress periods such as
    seasonality without duplicating the values. The DataArray should have
    dimensions ``("repeat", "repeat_items")``. The ``repeat_items``
    dimension should have size 2: the first value is the "key", the second
    value is the "value". For the "key" datetime, the data of the "value"
    datetime will be used. Can also be set with a dictionary using the
    ``set_repeat_stress`` method.

imod.mf6.Drainage Class Members
===============================
   * imod.mf6.Drainage.clip_box
   * imod.mf6.Drainage.from_file
   * imod.mf6.Drainage.get_non_grid_data
   * imod.mf6.Drainage.is_empty
   * imod.mf6.Drainage.mask
   * imod.mf6.Drainage.regrid_like
   * imod.mf6.Drainage.render
   * imod.mf6.Drainage.set_repeat_stress
   * imod.mf6.Drainage.to_netcdf
   * imod.mf6.Drainage.write

imod.mf6.Drainage.from_file
===========================
Loads an imod mf6 package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.mf6.Package.dataset.to_netcdf(),
as the checks upon package initialization are not done again!

Parameters
----------
path : str, pathlib.Path
    Path to the file.
**kwargs : keyword arguments
    Arbitrary keyword arguments forwarded to ``xarray.open_dataset()``, or
    ``xarray.open_zarr()``.
Refer to the examples.

Returns
-------
package : imod.mf6.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.mf6.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.mf6.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.mf6.Evapotranspiration
===========================
Evapotranspiration (EVT) Package.
Any number of EVT Packages can be specified for a single groundwater flow
model. All single-valued variables are free format.
https://water.usgs.gov/water-resources/software/MODFLOW-6/mf6io_6.0.4.pdf#page=86

Parameters
----------
surface: array of floats (xr.DataArray)
    is the elevation of the ET surface (L). A time-series name may be
    specified.
rate: array of floats (xr.DataArray)
    is the maximum ET flux rate (LT −1). A time-series name may be
    specified.
depth: array of floats (xr.DataArray)
    is the ET extinction depth (L). A time-series name may be specified.
proportion_rate: array of floats (xr.DataArray)
    is the proportion of the maximum ET flux rate at the bottom of a segment
    (dimensionless). A time-series name may be specified. (petm)
proportion_depth: array of floats (xr.DataArray)
    is the proportion of the ET extinction depth at the bottom of a segment
    (dimensionless). A timeseries name may be specified. (pxdp)
concentration: array of floats (xr.DataArray, optional)
    if this flow package is used in simulations also involving transport, then this array is used
    as the  concentration for inflow over this boundary.
concentration_boundary_type: ({"AUX", "AUXMIXED"}, optional)
    if this flow package is used in simulations also involving transport, then this keyword specifies
    how outflow over this boundary is computed.
fixed_cell: array of floats (xr.DataArray)
    indicates that evapotranspiration will not be reassigned to a cell
    underlying the cell specified in the list if the specified cell is
    inactive.
print_input: ({True, False}, optional)
    keyword to indicate that the list of evapotranspiration information will
    be written to the listing file immediately after it is read.
    Default is False.
print_flows: ({True, False}, optional)
    Indicates that the list of evapotranspiration flow rates will be printed
    to the listing file for every stress period time step in which "BUDGET
    PRINT" is specified in Output Control. If there is no Output Control
    option and PRINT FLOWS is specified, then flow rates are printed for the
    last time step of each stress period.
    Default is False.
save_flows: ({True, False}, optional)
    Indicates that evapotranspiration flow terms will be written to the file
    specified with "BUDGET FILEOUT" in Output Control.
    Default is False.
observations: [Not yet supported.]
    Default is None.
validate: {True, False}
    Flag to indicate whether the package should be validated upon
    initialization. This raises a ValidationError if package input is
    provided in the wrong manner. Defaults to True.
repeat_stress: Optional[xr.DataArray] of datetimes
    Used to repeat data for e.g. repeating stress periods such as
    seasonality without duplicating the values. The DataArray should have
    dimensions ``("repeat", "repeat_items")``. The ``repeat_items``
    dimension should have size 2: the first value is the "key", the second
    value is the "value". For the "key" datetime, the data of the "value"
    datetime will be used. Can also be set with a dictionary using the
    ``set_repeat_stress`` method.

imod.mf6.Evapotranspiration Class Members
=========================================
   * imod.mf6.Evapotranspiration.clip_box
   * imod.mf6.Evapotranspiration.from_file
   * imod.mf6.Evapotranspiration.get_non_grid_data
   * imod.mf6.Evapotranspiration.is_empty
   * imod.mf6.Evapotranspiration.mask
   * imod.mf6.Evapotranspiration.regrid_like
   * imod.mf6.Evapotranspiration.render
   * imod.mf6.Evapotranspiration.set_repeat_stress
   * imod.mf6.Evapotranspiration.to_netcdf
   * imod.mf6.Evapotranspiration.write

imod.mf6.Evapotranspiration.from_file
=====================================
Loads an imod mf6 package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.mf6.Package.dataset.to_netcdf(),
as the checks upon package initialization are not done again!

Parameters
----------
path : str, pathlib.Path
    Path to the file.
**kwargs : keyword arguments
    Arbitrary keyword arguments forwarded to ``xarray.open_dataset()``, or
    ``xarray.open_zarr()``.
Refer to the examples.

Returns
-------
package : imod.mf6.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.mf6.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.mf6.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.mf6.GeneralHeadBoundary
============================
The General-Head Boundary package is used to simulate head-dependent flux
boundaries.
https://water.usgs.gov/water-resources/software/MODFLOW-6/mf6io_6.0.4.pdf#page=75

Parameters
----------
head: array of floats (xr.DataArray)
    is the boundary head. (bhead)
conductance: array of floats (xr.DataArray)
    is the hydraulic conductance of the interface between the aquifer cell and
    the boundary.(cond)
concentration: array of floats (xr.DataArray, optional)
    if this flow package is used in simulations also involving transport, then this array is used
    as the  concentration for inflow over this boundary.
concentration_boundary_type: ({"AUX", "AUXMIXED"}, optional)
    if this flow package is used in simulations also involving transport, then this keyword specifies
    how outflow over this boundary is computed.
print_input: ({True, False}, optional)
    keyword to indicate that the list of general head boundary information
    will be written to the listing file immediately after it is read.
    Default is False.
print_flows: ({True, False}, optional)
    Indicates that the list of general head boundary flow rates will be
    printed to the listing file for every stress period time step in which
    "BUDGET PRINT" is specified in Output Control. If there is no Output
    Control option and PRINT FLOWS is specified, then flow rates are printed
    for the last time step of each stress period.
    Default is False.
save_flows: ({True, False}, optional)
    Indicates that general head boundary flow terms will be written to the
    file specified with "BUDGET FILEOUT" in Output Control.
    Default is False.
observations: [Not yet supported.]
    Default is None.
validate: {True, False}
    Flag to indicate whether the package should be validated upon
    initialization. This raises a ValidationError if package input is
    provided in the wrong manner. Defaults to True.
repeat_stress: Optional[xr.DataArray] of datetimes
    Used to repeat data for e.g. repeating stress periods such as
    seasonality without duplicating the values. The DataArray should have
    dimensions ``("repeat", "repeat_items")``. The ``repeat_items``
    dimension should have size 2: the first value is the "key", the second
    value is the "value". For the "key" datetime, the data of the "value"
    datetime will be used. Can also be set with a dictionary using the
    ``set_repeat_stress`` method.

imod.mf6.GeneralHeadBoundary Class Members
==========================================
   * imod.mf6.GeneralHeadBoundary.clip_box
   * imod.mf6.GeneralHeadBoundary.from_file
   * imod.mf6.GeneralHeadBoundary.get_non_grid_data
   * imod.mf6.GeneralHeadBoundary.is_empty
   * imod.mf6.GeneralHeadBoundary.mask
   * imod.mf6.GeneralHeadBoundary.regrid_like
   * imod.mf6.GeneralHeadBoundary.render
   * imod.mf6.GeneralHeadBoundary.set_repeat_stress
   * imod.mf6.GeneralHeadBoundary.to_netcdf
   * imod.mf6.GeneralHeadBoundary.write

imod.mf6.GeneralHeadBoundary.from_file
======================================
Loads an imod mf6 package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.mf6.Package.dataset.to_netcdf(),
as the checks upon package initialization are not done again!

Parameters
----------
path : str, pathlib.Path
    Path to the file.
**kwargs : keyword arguments
    Arbitrary keyword arguments forwarded to ``xarray.open_dataset()``, or
    ``xarray.open_zarr()``.
Refer to the examples.

Returns
-------
package : imod.mf6.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.mf6.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.mf6.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.mf6.HorizontalFlowBarrierHydraulicCharacteristic
=====================================================
 Horizontal Flow Barrier (HFB) package

Input to the Horizontal Flow Barrier (HFB) Package is read from the file
that has type "HFB6" in the Name File. Only one HFB Package can be
specified for a GWF model.
https://water.usgs.gov/water-resources/software/MODFLOW-6/mf6io_6.2.2.pdf

Parameters
----------
geometry: gpd.GeoDataFrame
    Dataframe that describes:
     - geometry: the geometries of the barriers,
     - hydraulic_characteristic: the hydraulic characteristic of the barriers
     - ztop: the top z-value of the barriers
     - zbottom: the bottom z-value of the barriers
print_input: bool

Examples
--------

>>> barrier_x = [-1000.0, 0.0, 1000.0]
>>> barrier_y = [500.0, 250.0, 500.0]
>>> barrier_gdf = gpd.GeoDataFrame(
>>>     geometry=[shapely.linestrings(barrier_x, barrier_y),],
>>>     data={
>>>         "hydraulic_characteristic": [1e-3,],
>>>         "ztop": [10.0,],
>>>         "zbottom": [0.0,],
>>>     },
>>> )
>>> hfb = imod.mf6.HorizontalFlowBarrierHydraulicCharacteristic(barrier_gdf)

imod.mf6.HorizontalFlowBarrierHydraulicCharacteristic Class Members
===================================================================
   * imod.mf6.HorizontalFlowBarrierHydraulicCharacteristic.clip_box
   * imod.mf6.HorizontalFlowBarrierHydraulicCharacteristic.from_file
   * imod.mf6.HorizontalFlowBarrierHydraulicCharacteristic.get_non_grid_data
   * imod.mf6.HorizontalFlowBarrierHydraulicCharacteristic.is_empty
   * imod.mf6.HorizontalFlowBarrierHydraulicCharacteristic.mask
   * imod.mf6.HorizontalFlowBarrierHydraulicCharacteristic.regrid_like
   * imod.mf6.HorizontalFlowBarrierHydraulicCharacteristic.render
   * imod.mf6.HorizontalFlowBarrierHydraulicCharacteristic.set_repeat_stress
   * imod.mf6.HorizontalFlowBarrierHydraulicCharacteristic.to_mf6_pkg
   * imod.mf6.HorizontalFlowBarrierHydraulicCharacteristic.to_netcdf
   * imod.mf6.HorizontalFlowBarrierHydraulicCharacteristic.write

imod.mf6.HorizontalFlowBarrierHydraulicCharacteristic.clip_box
==============================================================
Clip a package by a bounding box (time, layer, y, x).

Slicing intervals may be half-bounded, by providing None:

* To select 500.0 <= x <= 1000.0:
  ``clip_box(x_min=500.0, x_max=1000.0)``.
* To select x <= 1000.0: ``clip_box(x_min=None, x_max=1000.0)``
  or ``clip_box(x_max=1000.0)``.
* To select x >= 500.0: ``clip_box(x_min = 500.0, x_max=None.0)``
  or ``clip_box(x_min=1000.0)``.

Parameters
----------
time_min: optional
time_max: optional
layer_min: optional, int
layer_max: optional, int
x_min: optional, float
x_max: optional, float
y_min: optional, float
y_max: optional, float
top: optional, GridDataArray
bottom: optional, GridDataArray

Returns
-------
sliced : Package

imod.mf6.HorizontalFlowBarrierHydraulicCharacteristic.from_file
===============================================================
Loads an imod mf6 package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.mf6.Package.dataset.to_netcdf(),
as the checks upon package initialization are not done again!

Parameters
----------
path : str, pathlib.Path
    Path to the file.
**kwargs : keyword arguments
    Arbitrary keyword arguments forwarded to ``xarray.open_dataset()``, or
    ``xarray.open_zarr()``.
Refer to the examples.

Returns
-------
package : imod.mf6.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.mf6.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.mf6.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.mf6.HorizontalFlowBarrierHydraulicCharacteristic.is_empty
==============================================================
Returns True if the package is empty- for example if it contains only no-data values.

imod.mf6.HorizontalFlowBarrierHydraulicCharacteristic.mask
==========================================================
The mask method is irrelevant for this package as it is
grid-agnostic, instead this method retuns a copy of itself.

imod.mf6.HorizontalFlowBarrierHydraulicCharacteristic.render
============================================================
Render fills in the template only, doesn't write binary data

imod.mf6.HorizontalFlowBarrierHydraulicCharacteristic.to_mf6_pkg
================================================================
Write package to Modflow 6 package.

Based on the model grid, top and bottoms, the layers in which the barrier belong are computed. If the
barrier only partially occupies a layer an effective resistance or hydraulic conductivity for that layer is
calculated. This calculation is skipped for the Multiplier type.

Parameters
----------
idomain: GridDataArray
     Grid with active cells.
top: GridDataArray
    Grid with top of model layers.
bottom: GridDataArray
    Grid with bottom of model layers.
k: GridDataArray
    Grid with hydraulic conductivities.
validate: bool
    Run validation before converting

Returns
-------

imod.mf6.HorizontalFlowBarrierMultiplier
========================================
 Horizontal Flow Barrier (HFB) package

Input to the Horizontal Flow Barrier (HFB) Package is read from the file
that has type "HFB6" in the Name File. Only one HFB Package can be
specified for a GWF model.
https://water.usgs.gov/water-resources/software/MODFLOW-6/mf6io_6.2.2.pdf

If parts of the barrier overlap a layer the multiplier is applied to the entire layer.

Parameters
----------
geometry: gpd.GeoDataFrame
    Dataframe that describes:
     - geometry: the geometries of the barriers,
     - multiplier: the multiplier of the barriers
     - ztop: the top z-value of the barriers
     - zbottom: the bottom z-value of the barriers
print_input: bool

Examples
--------

>>> barrier_x = [-1000.0, 0.0, 1000.0]
>>> barrier_y = [500.0, 250.0, 500.0]
>>> barrier_gdf = gpd.GeoDataFrame(
>>>     geometry=[shapely.linestrings(barrier_x, barrier_y),],
>>>     data={
>>>         "multiplier": [1.5,],
>>>         "ztop": [10.0,],
>>>         "zbottom": [0.0,],
>>>     },
>>> )
>>> hfb = imod.mf6.HorizontalFlowBarrierMultiplier(barrier_gdf)

imod.mf6.HorizontalFlowBarrierMultiplier Class Members
======================================================
   * imod.mf6.HorizontalFlowBarrierMultiplier.clip_box
   * imod.mf6.HorizontalFlowBarrierMultiplier.from_file
   * imod.mf6.HorizontalFlowBarrierMultiplier.get_non_grid_data
   * imod.mf6.HorizontalFlowBarrierMultiplier.is_empty
   * imod.mf6.HorizontalFlowBarrierMultiplier.mask
   * imod.mf6.HorizontalFlowBarrierMultiplier.regrid_like
   * imod.mf6.HorizontalFlowBarrierMultiplier.render
   * imod.mf6.HorizontalFlowBarrierMultiplier.set_repeat_stress
   * imod.mf6.HorizontalFlowBarrierMultiplier.to_mf6_pkg
   * imod.mf6.HorizontalFlowBarrierMultiplier.to_netcdf
   * imod.mf6.HorizontalFlowBarrierMultiplier.write

imod.mf6.HorizontalFlowBarrierMultiplier.from_file
==================================================
Loads an imod mf6 package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.mf6.Package.dataset.to_netcdf(),
as the checks upon package initialization are not done again!

Parameters
----------
path : str, pathlib.Path
    Path to the file.
**kwargs : keyword arguments
    Arbitrary keyword arguments forwarded to ``xarray.open_dataset()``, or
    ``xarray.open_zarr()``.
Refer to the examples.

Returns
-------
package : imod.mf6.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.mf6.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.mf6.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.mf6.HorizontalFlowBarrierResistance
========================================
Horizontal Flow Barrier (HFB) package

Input to the Horizontal Flow Barrier (HFB) Package is read from the file
that has type "HFB6" in the Name File. Only one HFB Package can be
specified for a GWF model.
https://water.usgs.gov/water-resources/software/MODFLOW-6/mf6io_6.2.2.pdf

Parameters
----------
geometry: gpd.GeoDataFrame
    Dataframe that describes:
     - geometry: the geometries of the barriers,
     - resistance: the resistance of the barriers
     - ztop: the top z-value of the barriers
     - zbottom: the bottom z-value of the barriers
print_input: bool

Examples
--------

>>> barrier_x = [-1000.0, 0.0, 1000.0]
>>> barrier_y = [500.0, 250.0, 500.0]
>>> barrier_gdf = gpd.GeoDataFrame(
>>>     geometry=[shapely.linestrings(barrier_x, barrier_y),],
>>>     data={
>>>         "resistance": [1e3,],
>>>         "ztop": [10.0,],
>>>         "zbottom": [0.0,],
>>>     },
>>> )
>>> hfb = imod.mf6.HorizontalFlowBarrierResistance(barrier_gdf)

imod.mf6.HorizontalFlowBarrierResistance Class Members
======================================================
   * imod.mf6.HorizontalFlowBarrierResistance.clip_box
   * imod.mf6.HorizontalFlowBarrierResistance.from_file
   * imod.mf6.HorizontalFlowBarrierResistance.get_non_grid_data
   * imod.mf6.HorizontalFlowBarrierResistance.is_empty
   * imod.mf6.HorizontalFlowBarrierResistance.mask
   * imod.mf6.HorizontalFlowBarrierResistance.regrid_like
   * imod.mf6.HorizontalFlowBarrierResistance.render
   * imod.mf6.HorizontalFlowBarrierResistance.set_repeat_stress
   * imod.mf6.HorizontalFlowBarrierResistance.to_mf6_pkg
   * imod.mf6.HorizontalFlowBarrierResistance.to_netcdf
   * imod.mf6.HorizontalFlowBarrierResistance.write

imod.mf6.HorizontalFlowBarrierResistance.from_file
==================================================
Loads an imod mf6 package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.mf6.Package.dataset.to_netcdf(),
as the checks upon package initialization are not done again!

Parameters
----------
path : str, pathlib.Path
    Path to the file.
**kwargs : keyword arguments
    Arbitrary keyword arguments forwarded to ``xarray.open_dataset()``, or
    ``xarray.open_zarr()``.
Refer to the examples.

Returns
-------
package : imod.mf6.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.mf6.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.mf6.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.mf6.InitialConditions
==========================
Initial Conditions (IC) Package information is read from the file that is
specified by "IC6" as the file type. Only one IC Package can be specified
for a GWF model.
https://water.usgs.gov/water-resources/software/MODFLOW-6/mf6io_6.0.4.pdf#page=46

Parameters
----------
head: array of floats (xr.DataArray)
    for backwards compatibility this argument is maintained, but please use
    the start-argument instead.
start: array of floats (xr.DataArray)
    is the initial (starting) head or concentration—that is, the simulation's
    initial state.
    STRT must be specified for all simulations, including steady-state simulations.
    One value is read for every model cell. For
    simulations in which the first stress period is steady state, the values
    used for STRT generally do not affect the simulation (exceptions may
    occur if cells go dry and (or) rewet). The execution time, however, will
    be less if STRT includes hydraulic heads that are close to the
    steadystate solution. A head value lower than the cell bottom can be
    provided if a cell should start as dry. (strt)
validate: {True, False}
    Flag to indicate whether the package should be validated upon
    initialization. This raises a ValidationError if package input is
    provided in the wrong manner. Defaults to True.

imod.mf6.InitialConditions Class Members
========================================
   * imod.mf6.InitialConditions.clip_box
   * imod.mf6.InitialConditions.from_file
   * imod.mf6.InitialConditions.get_non_grid_data
   * imod.mf6.InitialConditions.is_empty
   * imod.mf6.InitialConditions.mask
   * imod.mf6.InitialConditions.regrid_like
   * imod.mf6.InitialConditions.to_netcdf

imod.mf6.InitialConditions.from_file
====================================
Loads an imod mf6 package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.mf6.Package.dataset.to_netcdf(),
as the checks upon package initialization are not done again!

Parameters
----------
path : str, pathlib.Path
    Path to the file.
**kwargs : keyword arguments
    Arbitrary keyword arguments forwarded to ``xarray.open_dataset()``, or
    ``xarray.open_zarr()``.
Refer to the examples.

Returns
-------
package : imod.mf6.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.mf6.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.mf6.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.mf6.NodePropertyFlow
=========================
Node Property Flow package.

In this package the hydraulic conductivity and rewetting in the model is
specified. A single NPF Package is required for each GWF model.
https://water.usgs.gov/water-resources/software/MODFLOW-6/mf6io_6.0.4.pdf#page=51

A note about regridding: the fields k, k22, k33 define the principal
components of an anisotropic conductivity tensor. By default, k and k22 are
in the horizontal plane and k33 is vertical. Angle1, angle2 and angle3
define the rotation of this tensor. The regridding methods associated by
default are chosen based on the assumption that k and k22 are horizontal and
k33 is vertical. If this is not the case, it is up to the user to regrid the
npf package using other regridding methods. This may be recommended if for
example the rotation is such that k has become vertical and k33 horizontal.

Parameters
----------
icelltype: array of int (xr.DataArray)
    flag for each cell that specifies how saturated thickness is treated. 0
    means saturated thickness is held constant; >0 means saturated thickness
    varies with computed head when head is below the cell top; <0 means
    saturated thickness varies with computed head unless the
    starting_head_as_confined_thickness option is in effect. When
    starting_head_as_confined_thickness is in effect, a negative value of
    icelltype indicates that saturated thickness will be computed as
    strt-bot and held constant.
k: array of floats (xr.DataArray)
    is the hydraulic conductivity. For the common case in which the user
    would like to specify the horizontal hydraulic conductivity and the
    vertical hydraulic conductivity, then K should be assigned as the
    horizontal hydraulic conductivity, K33 should be assigned as the
    vertical hydraulic conductivity, and K22 and the three rotation
    angles should not be specified. When more sophisticated anisotropy is
    required, then K corresponds to the K11 hydraulic conductivity axis. All
    included cells (idomain > 0) must have a K value greater than zero
rewet: ({True, False}, optional)
    activates model rewetting.
    Default is False.
rewet_layer: float
    is a combination of the wetting threshold and a flag to indicate which
    neighboring cells can cause a cell to become wet. If rewet_layer < 0,
    only a cell below a dry cell can cause the cell to become wet. If
    rewet_layer > 0, the cell below a dry cell and horizontally adjacent
    cells can cause a cell to become wet. If rewet_layer is 0, the cell
    cannot be wetted. The absolute value of rewet_layer is the wetting
    threshold. When the sum of BOT and the absolute value of rewet_layer at
    a dry cell is equaled or exceeded by the head at an adjacent cell, the
    cell is wetted. rewet_layer must be specified if "rewet" is specified in
    the OPTIONS block. If "rewet" is not specified in the options block,
    then rewet_layer can be entered, and memory will be allocated for it,
    even though it is not used. (WETDRY)
    Default is None.
rewet_factor:
    is a keyword and factor that is included in the calculation of the head
    that is initially established at a cell when that cell is converted from
    dry to wet. (WETFCT)
    Default is None.
rewet_iterations:
    is a keyword and iteration interval for attempting to wet cells. Wetting
    is attempted every rewet_iterations iteration. This applies to outer
    iterations and not inner iterations. If rewet_iterations is specified as
    zero or less, then the value is changed to 1. (IWETIT)
    Default is None.
rewet_method:
    is a keyword and integer flag that determines which equation is used to
    define the initial head at cells that become wet. If rewet_method is 0,
    h = BOT + rewet_factor (hm - BOT). If rewet_method is not 0, h = BOT +
    rewet_factor (THRESH). (IHDWET)
    Default is None.
k22: array of floats (xr.DataArray)
    is the hydraulic conductivity of the second ellipsoid axis; for an
    unrotated case this is the hydraulic conductivity in the y direction. If
    K22 is not included, then K22 is set equal to K.
    For a regular MODFLOW grid (DIS Package is used) in which no rotation
    angles are specified, K22 is the hydraulic conductivity along columns in
    the y direction. For an unstructured DISU grid, the user must assign
    principal x and y axes and provide the angle for each cell face relative
    to the assigned x direction. All included cells (idomain > 0) must have
    a K22 value greater than zero.
    Default is None.
k33: array of floats (xr.DataArray)
    is the hydraulic conductivity of the third ellipsoid axis; for an
    unrotated case, this is the vertical hydraulic conductivity. When
    anisotropy is applied, K33 corresponds to the K33 tensor component. All
    included cells (idomain > 0) must have a K33 value greater than zero.
    Default is None.
angle1: float
    is a rotation angle of the hydraulic conductivity tensor in degrees. The
    angle represents the first of three sequential rotations of the
    hydraulic conductivity ellipsoid. With the K11, K22, and K33 axes of the
    ellipsoid initially aligned with the x, y, and z coordinate axes,
    respectively, angle1 rotates the ellipsoid about its K33 axis (within
    the x - y plane). A positive value represents counter-clockwise rotation
    when viewed from any point on the positive K33 axis, looking toward the
    center of the ellipsoid. A value of zero indicates that the K11 axis
    lies within the x - z plane. If angle1 is not specified, default values
    of zero are assigned to angle1, angle2, and angle3, in which case the
    K11, K22, and K33 axes are aligned with the x, y, and z axes,
    respectively.
    Default is None.
angle2: float
    is a rotation angle of the hydraulic conductivity tensor in degrees. The
    angle represents the second of three sequential rotations of the
    hydraulic conductivity ellipsoid. Following the rotation by angle1
    described above, angle2 rotates the ellipsoid about its K22 axis (out of
    the x - y plane). An array can be specified for angle2 only if angle1 is
    also specified. A positive value of angle2 represents clockwise rotation
    when viewed from any point on the positive K22 axis, looking toward the
    center of the ellipsoid. A value of zero indicates that the K11 axis
    lies within the x - y plane. If angle2 is not specified, default values
    of zero are assigned to angle2 and angle3; connections that are not
    user-designated as vertical are assumed to be strictly horizontal (that
    is, to have no z component to their orientation); and connection lengths
    are based on horizontal distances.
    Default is None.
angle3: float
    is a rotation angle of the hydraulic conductivity tensor in degrees. The
    angle represents the third of three sequential rotations of the
    hydraulic conductivity ellipsoid. Following the rotations by angle1 and
    angle2 described above, angle3 rotates the ellipsoid about its K11 axis.
    An array can be specified for angle3 only if angle1 and angle2 are also
    specified. An array must be specified for angle3 if angle2 is specified.
    A positive value of angle3 represents clockwise rotation when viewed
    from any point on the positive K11 axis, looking toward the center of
    the ellipsoid. A value of zero indicates that the K22 axis lies within
    the x - y plane.
    Default is None.
alternative_cell_averaging : str
    Method calculating horizontal cell connection conductance.
    Options: {"LOGARITHMIC", "AMT-LMK", or "AMT-HMK"}
    Default: uses harmonic mean for averaging
save_flows: ({True, False}, optional)
    keyword to indicate that cell-by-cell flow terms will be written to the
    file specified with "budget save file" in Output Control.
    Default is False.
starting_head_as_confined_thickness: ({True, False}, optional)
    indicates that cells having a negative icelltype are confined, and their
    cell thickness for conductance calculations will be computed as strt-bot
    rather than top-bot.
    (THICKSTRT)
    Default is False.
variable_vertical_conductance: ({True, False}, optional)
    keyword to indicate that the vertical conductance will be calculated
    using the saturated thickness and properties of the overlying cell and
    the thickness and properties of the underlying cell. if the dewatered
    keyword is also specified, then the vertical conductance is calculated
    using only the saturated thickness and properties of the overlying cell
    if the head in the underlying cell is below its top. if these keywords
    are not specified, then the default condition is to calculate the
    vertical conductance at the start of the simulation using the initial
    head and the cell properties. the vertical conductance remains constant
    for the entire simulation.
    (VARIABLECV)
    Default is False.
dewatered: ({True, False}, optional)
    If the dewatered keyword is specified, then the vertical conductance is
    calculated using only the saturated thickness and properties of the
    overlying cell if the head in the underlying cell is below its top.
    Default is False.
perched: ({True, False}, optional)
    keyword to indicate that when a cell is overlying a dewatered
    convertible cell, the head difference used in Darcy’s Law is equal to
    the head in the overlying cell minus the bottom elevation of the
    overlying cell. If not specified, then the default is to use the head
    difference between the two cells.
    Default is False.
save_specific_discharge: ({True, False}, optional)
    keyword to indicate that x, y, and z components of specific discharge
    will be calculated at cell centers and written to the cell-by-cell flow
    file, which is specified with"budget save file" in Output Control. If
    this option is activated, then additional information may be required in
    the discretization packages and the GWF Exchange package (if GWF models
    are coupled). Specifically, angldegx must be specified in the
    connectiondata block of the disu package; angldegx must also be
    specified for the GWF Exchange as an auxiliary variable. disu package
    has not been implemented yet.
    Default is False.
save_saturation: ({True, False}, optional)
    keyword to indicate that cell saturation will be written to the budget
    file, which is specified with "BUDGET SAVE FILE" in Output Control.
    Saturation will be saved to the budget file as an auxiliary variable
    saved with the DATA-SAT text label. Saturation is a cell variable that
    ranges from zero to one and can be used by post processing programs to
    determine how much of a cell volume is saturated. If ICELLTYPE is 0,
    then saturation is always one.
xt3d_option:  ({True, False}, optional)
    If True, the XT3D formulation will be used. By default False.
rhs_option: ({True, False}, optional)
    If True, then the XT3D additional terms will be added to the right-hand
    side. If False, then the XT3D terms will be put into the coefficient
    matrix. By default False.
validate: {True, False}
    Flag to indicate whether the package should be validated upon
    initialization. This raises a ValidationError if package input is
    provided in the wrong manner. Defaults to True.

imod.mf6.NodePropertyFlow Class Members
=======================================
   * imod.mf6.NodePropertyFlow.clip_box
   * imod.mf6.NodePropertyFlow.from_file
   * imod.mf6.NodePropertyFlow.get_non_grid_data
   * imod.mf6.NodePropertyFlow.get_xt3d_option
   * imod.mf6.NodePropertyFlow.is_dewatered
   * imod.mf6.NodePropertyFlow.is_empty
   * imod.mf6.NodePropertyFlow.is_variable_vertical_conductance
   * imod.mf6.NodePropertyFlow.mask
   * imod.mf6.NodePropertyFlow.regrid_like
   * imod.mf6.NodePropertyFlow.set_xt3d_option
   * imod.mf6.NodePropertyFlow.to_netcdf

imod.mf6.NodePropertyFlow.from_file
===================================
Loads an imod mf6 package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.mf6.Package.dataset.to_netcdf(),
as the checks upon package initialization are not done again!

Parameters
----------
path : str, pathlib.Path
    Path to the file.
**kwargs : keyword arguments
    Arbitrary keyword arguments forwarded to ``xarray.open_dataset()``, or
    ``xarray.open_zarr()``.
Refer to the examples.

Returns
-------
package : imod.mf6.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.mf6.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.mf6.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.mf6.NodePropertyFlow.get_xt3d_option
=========================================
Returns the xt3d option value for this object.

imod.mf6.NodePropertyFlow.is_dewatered
======================================
Returns the "dewatered" option value for this object. Used only when variable_vertical_conductance is true

imod.mf6.NodePropertyFlow.is_variable_vertical_conductance
==========================================================
Returns the VariableCV option value for this object.

imod.mf6.NodePropertyFlow.set_xt3d_option
=========================================
Returns the xt3d option value for this object.

imod.mf6.Recharge
=================
Recharge Package.
Any number of RCH Packages can be specified for a single groundwater flow
model.
https://water.usgs.gov/water-resources/software/MODFLOW-6/mf6io_6.0.4.pdf#page=79

Parameters
----------
rate: array of floats (xr.DataArray)
    is the recharge flux rate (LT −1). This rate is multiplied inside the
    program by the surface area of the cell to calculate the volumetric
    recharge rate. A time-series name may be specified.
concentration: array of floats (xr.DataArray, optional)
    if this flow package is used in simulations also involving transport, then this array is used
    as the  concentration for inflow over this boundary.
concentration_boundary_type: ({"AUX", "AUXMIXED"}, optional)
    if this flow package is used in simulations also involving transport, then this keyword specifies
    how outflow over this boundary is computed.
print_input: ({True, False}, optional)
    keyword to indicate that the list of recharge information will be
    written to the listing file immediately after it is read.
    Default is False.
print_flows: ({True, False}, optional)
    Indicates that the list of recharge flow rates will be printed to the
    listing file for every stress period time step in which "BUDGET PRINT"is
    specified in Output Control. If there is no Output Control option and
    PRINT FLOWS is specified, then flow rates are printed for the last time
    step of each stress period.
    Default is False.
save_flows: ({True, False}, optional)
    Indicates that recharge flow terms will be written to the file specified
    with "BUDGET FILEOUT" in Output Control.
    Default is False.
observations: [Not yet supported.]
    Default is None.
validate: {True, False}
    Flag to indicate whether the package should be validated upon
    initialization. This raises a ValidationError if package input is
    provided in the wrong manner. Defaults to True.
repeat_stress: Optional[xr.DataArray] of datetimes
    Used to repeat data for e.g. repeating stress periods such as
    seasonality without duplicating the values. The DataArray should have
    dimensions ``("repeat", "repeat_items")``. The ``repeat_items``
    dimension should have size 2: the first value is the "key", the second
    value is the "value". For the "key" datetime, the data of the "value"
    datetime will be used. Can also be set with a dictionary using the
    ``set_repeat_stress`` method.
fixed_cell: ({True, False}, optional)
    indicates that recharge will not be reassigned to a cell underlying the
    cell specified in the list if the specified cell is inactive.

imod.mf6.Recharge Class Members
===============================
   * imod.mf6.Recharge.clip_box
   * imod.mf6.Recharge.from_file
   * imod.mf6.Recharge.get_non_grid_data
   * imod.mf6.Recharge.is_empty
   * imod.mf6.Recharge.mask
   * imod.mf6.Recharge.regrid_like
   * imod.mf6.Recharge.render
   * imod.mf6.Recharge.set_repeat_stress
   * imod.mf6.Recharge.to_netcdf
   * imod.mf6.Recharge.write

imod.mf6.Recharge.from_file
===========================
Loads an imod mf6 package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.mf6.Package.dataset.to_netcdf(),
as the checks upon package initialization are not done again!

Parameters
----------
path : str, pathlib.Path
    Path to the file.
**kwargs : keyword arguments
    Arbitrary keyword arguments forwarded to ``xarray.open_dataset()``, or
    ``xarray.open_zarr()``.
Refer to the examples.

Returns
-------
package : imod.mf6.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.mf6.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.mf6.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.mf6.River
==============
River package.
Any number of RIV Packages can be specified for a single groundwater flow
model.
https://water.usgs.gov/water-resources/software/MODFLOW-6/mf6io_6.0.4.pdf#page=71

Parameters
----------
stage: array of floats (xr.DataArray)
    is the head in the river.
conductance: array of floats (xr.DataArray)
    is the riverbed hydraulic conductance.
bottom_elevation: array of floats (xr.DataArray)
    is the elevation of the bottom of the riverbed.
concentration: array of floats (xr.DataArray, optional)
    if this flow package is used in simulations also involving transport, then this array is used
    as the  concentration for inflow over this boundary.
concentration_boundary_type: ({"AUX", "AUXMIXED"}, optional)
    if this flow package is used in simulations also involving transport, then this keyword specifies
    how outflow over this boundary is computed.
print_input: ({True, False}, optional)
    keyword to indicate that the list of river information will be written
    to the listing file immediately after it is read. Default is False.
print_flows: ({True, False}, optional)
    Indicates that the list of river flow rates will be printed to the
    listing file for every stress period time step in which "BUDGET PRINT"
    is specified in Output Control. If there is no Output Control option and
    PRINT FLOWS is specified, then flow rates are printed for the last time
    step of each stress period. Default is False.
save_flows: ({True, False}, optional)
    Indicates that river flow terms will be written to the file specified
    with "BUDGET FILEOUT" in Output Control. Default is False.
observations: [Not yet supported.]
    Default is None.
validate: {True, False}
    Flag to indicate whether the package should be validated upon
    initialization. This raises a ValidationError if package input is
    provided in the wrong manner. Defaults to True.
repeat_stress: Optional[xr.DataArray] of datetimes
    Used to repeat data for e.g. repeating stress periods such as
    seasonality without duplicating the values. The DataArray should have
    dimensions ``("repeat", "repeat_items")``. The ``repeat_items``
    dimension should have size 2: the first value is the "key", the second
    value is the "value". For the "key" datetime, the data of the "value"
    datetime will be used. Can also be set with a dictionary using the
    ``set_repeat_stress`` method.

imod.mf6.River Class Members
============================
   * imod.mf6.River.clip_box
   * imod.mf6.River.from_file
   * imod.mf6.River.get_non_grid_data
   * imod.mf6.River.is_empty
   * imod.mf6.River.mask
   * imod.mf6.River.regrid_like
   * imod.mf6.River.render
   * imod.mf6.River.set_repeat_stress
   * imod.mf6.River.to_netcdf
   * imod.mf6.River.write

imod.mf6.River.from_file
========================
Loads an imod mf6 package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.mf6.Package.dataset.to_netcdf(),
as the checks upon package initialization are not done again!

Parameters
----------
path : str, pathlib.Path
    Path to the file.
**kwargs : keyword arguments
    Arbitrary keyword arguments forwarded to ``xarray.open_dataset()``, or
    ``xarray.open_zarr()``.
Refer to the examples.

Returns
-------
package : imod.mf6.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.mf6.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.mf6.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.mf6.SpecificStorage
========================
Storage Package with specific storage.

From wikipedia (https://en.wikipedia.org/wiki/Specific_storage):

"The specific storage is the amount of water that a portion of an aquifer
releases from storage, per unit mass or volume of aquifer, per unit change
in hydraulic head, while remaining fully saturated."

If the STO Package is not included for a model, then storage changes will
not be calculated, and thus, the model will be steady state. Only one STO
Package can be specified for a GWF model.

Parameters
----------
specific_storage: array of floats (xr.DataArray)
    Is specific storage. Specific storage values must be greater than
    or equal to 0. (ss)
specific_yield: array of floats (xr.DataArray)
    Is specific yield. Specific yield values must be greater than or
    equal to 0. Specific yield does not have to be specified if there are no
    convertible cells (convertible=0 in every cell). (sy)
transient: ({True, False}), or a DataArray with a time coordinate and dtype Bool
    Boolean to indicate if the model is transient or steady-state.
convertible: array of int (xr.DataArray)
    Is a flag for each cell that specifies whether or not a cell is
    convertible for the storage calculation. 0 indicates confined storage is
    used. >0 indicates confined storage is used when head is above cell top
    and a mixed formulation of unconfined and confined storage is used when
    head is below cell top. (iconvert)
save_flows: ({True, False}, optional)
    Indicates that storage flow terms will be written to the file specified
    with "BUDGET FILEOUT" in Output Control. Default is False.
validate: {True, False}
    Flag to indicate whether the package should be validated upon
    initialization. This raises a ValidationError if package input is
    provided in the wrong manner. Defaults to True.

imod.mf6.SpecificStorage Class Members
======================================
   * imod.mf6.SpecificStorage.clip_box
   * imod.mf6.SpecificStorage.from_file
   * imod.mf6.SpecificStorage.get_non_grid_data
   * imod.mf6.SpecificStorage.is_empty
   * imod.mf6.SpecificStorage.mask
   * imod.mf6.SpecificStorage.regrid_like
   * imod.mf6.SpecificStorage.to_netcdf

imod.mf6.SpecificStorage.from_file
==================================
Loads an imod mf6 package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.mf6.Package.dataset.to_netcdf(),
as the checks upon package initialization are not done again!

Parameters
----------
path : str, pathlib.Path
    Path to the file.
**kwargs : keyword arguments
    Arbitrary keyword arguments forwarded to ``xarray.open_dataset()``, or
    ``xarray.open_zarr()``.
Refer to the examples.

Returns
-------
package : imod.mf6.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.mf6.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.mf6.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.mf6.StorageCoefficient
===========================
Storage Package with a storage coefficient.  Be careful,
this is not the same as the specific storage.

From wikipedia (https://en.wikipedia.org/wiki/Specific_storage):

"Storativity or the storage coefficient is the volume of water released
from storage per unit decline in hydraulic head in the aquifer, per
unit area of the aquifer.  Storativity is a dimensionless quantity, and
is always greater than 0.

Under confined conditions:

S = Ss * b, where S is the storage coefficient,
Ss the specific storage, and b the aquifer thickness.

Under unconfined conditions:

S = Sy, where Sy is the specific yield"

If the STO Package is not included for a model, then storage changes will
not be calculated, and thus, the model will be steady state. Only one STO
Package can be specified for a GWF model.

Parameters
----------
storage_coefficient: array of floats (xr.DataArray)
    Is storage coefficient. Storage coefficient values must be greater than
    or equal to 0. (ss)
specific_yield: array of floats (xr.DataArray)
    Is specific yield. Specific yield values must be greater than or
    equal to 0. Specific yield does not have to be specified if there are no
    convertible cells (convertible=0 in every cell). (sy)
transient: ({True, False})
    Boolean to indicate if the model is transient or steady-state.
convertible: array of int (xr.DataArray)
    Is a flag for each cell that specifies whether or not a cell is
    convertible for the storage calculation. 0 indicates confined storage is
    used. >0 indicates confined storage is used when head is above cell top
    and a mixed formulation of unconfined and confined storage is used when
    head is below cell top. (iconvert)
save_flows: ({True, False}, optional)
    Indicates that storage flow terms will be written to the file specified
    with "BUDGET FILEOUT" in Output Control. Default is False.
validate: {True, False}
    Flag to indicate whether the package should be validated upon
    initialization. This raises a ValidationError if package input is
    provided in the wrong manner. Defaults to True.

imod.mf6.StorageCoefficient Class Members
=========================================
   * imod.mf6.StorageCoefficient.clip_box
   * imod.mf6.StorageCoefficient.from_file
   * imod.mf6.StorageCoefficient.get_non_grid_data
   * imod.mf6.StorageCoefficient.is_empty
   * imod.mf6.StorageCoefficient.mask
   * imod.mf6.StorageCoefficient.regrid_like
   * imod.mf6.StorageCoefficient.to_netcdf

imod.mf6.StorageCoefficient.from_file
=====================================
Loads an imod mf6 package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.mf6.Package.dataset.to_netcdf(),
as the checks upon package initialization are not done again!

Parameters
----------
path : str, pathlib.Path
    Path to the file.
**kwargs : keyword arguments
    Arbitrary keyword arguments forwarded to ``xarray.open_dataset()``, or
    ``xarray.open_zarr()``.
Refer to the examples.

Returns
-------
package : imod.mf6.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.mf6.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.mf6.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.mf6.UnsaturatedZoneFlow
============================
Unsaturated Zone Flow (UZF) package.

TODO: Support timeseries file? Observations? Water Mover?

Parameters
----------
surface_depression_depth: array of floats (xr.DataArray)
    is the surface depression depth of the UZF cell.
kv_sat: array of floats (xr.DataArray)
    is the vertical saturated hydraulic conductivity of the UZF cell.
    NOTE: the UZF package determines the location of inactive cells where kv_sat is np.nan
theta_res: array of floats (xr.DataArray)
    is the residual (irreducible) water content of the UZF cell.
theta_sat: array of floats (xr.DataArray)
    is the saturated water content of the UZF cell.
theta_init: array of floats (xr.DataArray)
    is the initial water content of the UZF cell.
epsilon: array of floats (xr.DataArray)
    is the epsilon exponent of the UZF cell.
infiltration_rate: array of floats (xr.DataArray)
    defines the applied infiltration rate of the UZF cell (LT -1).
et_pot: array of floats (xr.DataArray, optional)
    defines the potential evapotranspiration rate of the UZF cell and specified
    GWF cell. Evapotranspiration is first removed from the unsaturated zone and any remaining
    potential evapotranspiration is applied to the saturated zone. If IVERTCON is greater than zero
    then residual potential evapotranspiration not satisfied in the UZF cell is applied to the underlying
    UZF and GWF cells.
extinction_depth: array of floats (xr.DataArray, optional)
    defines the evapotranspiration extinction depth of the UZF cell. If
    IVERTCON is greater than zero and EXTDP extends below the GWF cell bottom then remaining
    potential evapotranspiration is applied to the underlying UZF and GWF cells. EXTDP is always
    specified, but is only used if SIMULATE ET is specified in the OPTIONS block.
extinction_theta: array of floats (xr.DataArray, optional)
    defines the evapotranspiration extinction water content of the UZF
    cell. If specified, ET in the unsaturated zone will be simulated either as a function of the
    specified PET rate while the water content (THETA) is greater than the ET extinction water content
air_entry_potential: array of floats (xr.DataArray, optional)
    defines the air entry potential (head) of the UZF cell. If specified, ET will be
    simulated using a capillary pressure based formulation.
    Capillary pressure is calculated using the Brooks-Corey retention function ("air_entry")
root_potential: array of floats (xr.DataArray, optional)
    defines the root potential (head) of the UZF cell. If specified, ET will be
    simulated using a capillary pressure based formulation.
    Capillary pressure is calculated using the Brooks-Corey retention function ("air_entry"
root_activity: array of floats (xr.DataArray, optional)
    defines the root activity function of the UZF cell. ROOTACT is
    the length of roots in a given volume of soil divided by that volume. Values range from 0 to about 3
    cm-2, depending on the plant community and its stage of development. If specified, ET will be
    simulated using a capillary pressure based formulation.
    Capillary pressure is calculated using the Brooks-Corey retention function ("air_entry"
groundwater_ET_function: ({"linear", "square"}, optional)
    keyword specifying that groundwater evapotranspiration will be simulated using either
    the original ET formulation of MODFLOW-2005 ("linear"). Or by assuming a constant ET
    rate for groundwater levels between land surface (TOP) and land surface minus the ET extinction
    depth (TOP-EXTDP) ("square"). In the latter case, groundwater ET is smoothly reduced
    from the PET rate to zero over a nominal interval at TOP-EXTDP.
simulate_seepage: ({True, False}, optional)
    keyword specifying that groundwater discharge (GWSEEP) to land surface will be
    simulated. Groundwater discharge is nonzero when groundwater head is greater than land surface.
print_input: ({True, False}, optional)
    keyword to indicate that the list of UZF information will be written to the listing file
    immediately after it is read.
    Default is False.
print_flows: ({True, False}, optional)
    keyword to indicate that the list of UZF flow rates will be printed to the listing file for
    every stress period time step in which "BUDGET PRINT" is specified in Output Control. If there is
    no Output Control option and "PRINT FLOWS" is specified, then flow rates are printed for the last
    time step of each stress period.
    Default is False.
save_flows: ({True, False}, optional)
    keyword to indicate that UZF flow terms will be written to the file specified with "BUDGET
    FILEOUT" in Output Control.
    Default is False.
budget_fileout: ({"str"}, optional)
    path to output cbc-file for UZF budgets
budgetcsv_fileout: ({"str"}, optional)
    path to output csv-file for summed budgets
observations: [Not yet supported.]
    Default is None.
water_mover: [Not yet supported.]
    Default is None.
timeseries: [Not yet supported.]
    Default is None.
    TODO: We could allow the user to either use xarray DataArrays to specify BCS or
    use a pd.DataFrame and use the MF6 timeseries files to read input. The latter could
    save memory for laterally large-scale models, through efficient use of the UZF cell identifiers.
validate: {True, False}
    Flag to indicate whether the package should be validated upon
    initialization. This raises a ValidationError if package input is
    provided in the wrong manner. Defaults to True.

imod.mf6.UnsaturatedZoneFlow Class Members
==========================================
   * imod.mf6.UnsaturatedZoneFlow.clip_box
   * imod.mf6.UnsaturatedZoneFlow.fill_stress_perioddata
   * imod.mf6.UnsaturatedZoneFlow.from_file
   * imod.mf6.UnsaturatedZoneFlow.get_non_grid_data
   * imod.mf6.UnsaturatedZoneFlow.is_empty
   * imod.mf6.UnsaturatedZoneFlow.mask
   * imod.mf6.UnsaturatedZoneFlow.regrid_like
   * imod.mf6.UnsaturatedZoneFlow.render
   * imod.mf6.UnsaturatedZoneFlow.set_repeat_stress
   * imod.mf6.UnsaturatedZoneFlow.to_netcdf
   * imod.mf6.UnsaturatedZoneFlow.write

imod.mf6.UnsaturatedZoneFlow.fill_stress_perioddata
===================================================
Modflow6 requires something to be filled in the stress perioddata,
even though the data is not used in the current configuration.
Only an infiltration rate is required,
the rest can be filled with dummy values if not provided.

imod.mf6.UnsaturatedZoneFlow.from_file
======================================
Loads an imod mf6 package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.mf6.Package.dataset.to_netcdf(),
as the checks upon package initialization are not done again!

Parameters
----------
path : str, pathlib.Path
    Path to the file.
**kwargs : keyword arguments
    Arbitrary keyword arguments forwarded to ``xarray.open_dataset()``, or
    ``xarray.open_zarr()``.
Refer to the examples.

Returns
-------
package : imod.mf6.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.mf6.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.mf6.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.mf6.UnsaturatedZoneFlow.render
===================================
Render fills in the template only, doesn't write binary data

imod.mf6.UnsaturatedZoneFlow.write
==================================
writes the blockfile and binary data

directory is modelname

imod.mf6.Well
=============
Agnostic WEL package, which accepts x, y and a top and bottom of the well screens.

This package can be written to any provided model grid.
Any number of WEL Packages can be specified for a single groundwater flow model.
https://water.usgs.gov/water-resources/software/MODFLOW-6/mf6io_6.0.4.pdf#page=63

Parameters
----------

y: float or list of floats
    is the y location of the well.
x: float or list of floats
    is the x location of the well.
screen_top: float or list of floats
    is the top of the well screen.
screen_bottom: float or list of floats
    is the bottom of the well screen.
rate: float, list of floats or xr.DataArray
    is the volumetric well rate. A positive value indicates well
    (injection) and a negative value indicates discharge (extraction) (q).
    If provided as DataArray, an ``"index"`` dimension is required and an
    optional ``"time"`` dimension and coordinate specify transient input.
    In the latter case, it is important that dimensions are in the order:
    ``("time", "index")``
concentration: array of floats (xr.DataArray, optional)
    if this flow package is used in simulations also involving transport, then this array is used
    as the  concentration for inflow over this boundary.
concentration_boundary_type: ({"AUX", "AUXMIXED"}, optional)
    if this flow package is used in simulations also involving transport, then this keyword specifies
    how outflow over this boundary is computed.
id: list of Any, optional
    assign an identifier code to each well. if not provided, one will be generated
    Must be convertible to string, and unique entries.
minimum_k: float, optional
    on creating point wells, no point wells will be placed in cells with a lower horizontal conductivity than this
minimum_thickness: float, optional
    on creating point wells, no point wells will be placed in cells with a lower thickness than this
print_input: ({True, False}, optional)
    keyword to indicate that the list of well information will be written to
    the listing file immediately after it is read.
    Default is False.
print_flows: ({True, False}, optional)
    Indicates that the list of well flow rates will be printed to the
    listing file for every stress period time step in which "BUDGET PRINT"
    is specified in Output Control. If there is no Output Control option
    and PRINT FLOWS is specified, then flow rates are printed for the last
    time step of each stress period.
    Default is False.
save_flows: ({True, False}, optional)
    Indicates that well flow terms will be written to the file specified
    with "BUDGET FILEOUT" in Output Control.
    Default is False.
observations: [Not yet supported.]
    Default is None.
validate: {True, False}
    Flag to indicate whether the package should be validated upon
    initialization. This raises a ValidationError if package input is
    provided in the wrong manner. Defaults to True.
repeat_stress: Optional[xr.DataArray] of datetimes
    Used to repeat data for e.g. repeating stress periods such as
    seasonality without duplicating the values. The DataArray should have
    dimensions ``("repeat", "repeat_items")``. The ``repeat_items``
    dimension should have size 2: the first value is the "key", the second
    value is the "value". For the "key" datetime, the data of the "value"
    datetime will be used. Can also be set with a dictionary using the
    ``set_repeat_stress`` method.

Examples
---------

>>> screen_top = [0.0, 0.0]
>>> screen_bottom = [-2.0, -2.0]
>>> y = [83.0, 77.0]
>>> x = [81.0, 82.0]
>>> rate = [1.0, 1.0]

>>> imod.mf6.Well(x, y, screen_top, screen_bottom, rate)

For a transient well:

>>> weltimes = pd.date_range("2000-01-01", "2000-01-03")

>>> rate_factor_time = xr.DataArray([0.5, 1.0], coords={"time": weltimes}, dims=("time",))
>>> rate_transient = rate_factor_time * xr.DataArray(rate, dims=("index",))

>>> imod.mf6.Well(x, y, screen_top, screen_bottom, rate_transient)

imod.mf6.Well Class Members
===========================
   * imod.mf6.Well.clip_box
   * imod.mf6.Well.from_file
   * imod.mf6.Well.get_non_grid_data
   * imod.mf6.Well.is_empty
   * imod.mf6.Well.mask
   * imod.mf6.Well.regrid_like
   * imod.mf6.Well.render
   * imod.mf6.Well.set_repeat_stress
   * imod.mf6.Well.to_mf6_pkg
   * imod.mf6.Well.to_netcdf
   * imod.mf6.Well.write

imod.mf6.Well.clip_box
======================
Clip a package by a bounding box (time, layer, y, x).

The well package doesn't use the layer attribute to describe its depth and length.
Instead, it uses the screen_top and screen_bottom parameters which corresponds with
the z-coordinates of the top and bottom of the well. To go from a layer_min and
layer_max to z-values used for clipping the well a top and bottom array have to be
provided as well.

Slicing intervals may be half-bounded, by providing None:

* To select 500.0 <= x <= 1000.0:
  ``clip_box(x_min=500.0, x_max=1000.0)``.
* To select x <= 1000.0: ``clip_box(x_min=None, x_max=1000.0)``
  or ``clip_box(x_max=1000.0)``.
* To select x >= 500.0: ``clip_box(x_min = 500.0, x_max=None.0)``
  or ``clip_box(x_min=1000.0)``.

Parameters
----------
time_min: optional
time_max: optional
layer_min: optional, int
layer_max: optional, int
x_min: optional, float
x_max: optional, float
y_min: optional, float
y_max: optional, float
top: optional, GridDataArray
bottom: optional, GridDataArray
state_for_boundary: optional, GridDataArray

Returns
-------
sliced : Package

imod.mf6.Well.from_file
=======================
Loads an imod mf6 package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.mf6.Package.dataset.to_netcdf(),
as the checks upon package initialization are not done again!

Parameters
----------
path : str, pathlib.Path
    Path to the file.
**kwargs : keyword arguments
    Arbitrary keyword arguments forwarded to ``xarray.open_dataset()``, or
    ``xarray.open_zarr()``.
Refer to the examples.

Returns
-------
package : imod.mf6.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.mf6.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.mf6.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.mf6.Well.mask
==================
Mask wells based on two-dimensional domain. For three-dimensional
masking: Wells falling in inactive cells are automatically removed in
the call to write to Modflow 6 package. You can verify this by calling
the ``to_mf6_pkg`` method.

imod.mf6.Well.render
====================
Render fills in the template only, doesn't write binary data

imod.mf6.Well.to_mf6_pkg
========================
Write package to Modflow 6 package.

Based on the model grid and top and bottoms, cellids are determined.
When well screens hit multiple layers, groundwater extractions are
distributed based on layer transmissivities. Wells located in inactive
cells are removed.

Note
----
The well distribution based on transmissivities assumes confined
aquifers. If wells fall dry (and the rate distribution has to be
recomputed at runtime), it is better to use the Multi-Aquifer Well
package.

Parameters
----------
is_partitioned: bool
validate: bool
    Run validation before converting
active: {xarry.DataArray, xugrid.UgridDataArray}
    Grid with active cells.
top: {xarry.DataArray, xugrid.UgridDataArray}
    Grid with top of model layers.
bottom: {xarry.DataArray, xugrid.UgridDataArray}
    Grid with bottom of model layers.
k: {xarry.DataArray, xugrid.UgridDataArray}
    Grid with hydraulic conductivities.
Returns
-------
Mf6Wel
    Object with wells as list based input.

imod.mf6.Well.write
===================
writes the blockfile and binary data

directory is modelname

imod.mf6.WellDisStructured
==========================
WEL package for structured discretization (DIS) models .
Any number of WEL Packages can be specified for a single groundwater flow model.
https://water.usgs.gov/water-resources/software/MODFLOW-6/mf6io_6.0.4.pdf#page=63

.. warning::
    This class is deprecated and will be deleted in a future release.
    Consider changing your code to use the ``imod.mf6.Well`` package.

Parameters
----------
layer: list of int
    Model layer in which the well is located.
row: list of int
    Row in which the well is located.
column: list of int
    Column in which the well is located.
rate: float or list of floats
    is the volumetric well rate. A positive value indicates well
    (injection) and a negative value indicates discharge (extraction) (q).
concentration: array of floats (xr.DataArray, optional)
    if this flow package is used in simulations also involving transport, then this array is used
    as the  concentration for inflow over this boundary.
concentration_boundary_type: ({"AUX", "AUXMIXED"}, optional)
    if this flow package is used in simulations also involving transport, then this keyword specifies
    how outflow over this boundary is computed.
print_input: ({True, False}, optional)
    keyword to indicate that the list of well information will be written to
    the listing file immediately after it is read.
    Default is False.
print_flows: ({True, False}, optional)
    Indicates that the list of well flow rates will be printed to the
    listing file for every stress period time step in which "BUDGET PRINT"
    is specified in Output Control. If there is no Output Control option
    and PRINT FLOWS is specified, then flow rates are printed for the last
    time step of each stress period.
    Default is False.
save_flows: ({True, False}, optional)
    Indicates that well flow terms will be written to the file specified
    with "BUDGET FILEOUT" in Output Control.
    Default is False.
observations: [Not yet supported.]
    Default is None.
validate: {True, False}
    Flag to indicate whether the package should be validated upon
    initialization. This raises a ValidationError if package input is
    provided in the wrong manner. Defaults to True.
repeat_stress: Optional[xr.DataArray] of datetimes
    Used to repeat data for e.g. repeating stress periods such as
    seasonality without duplicating the values. The DataArray should have
    dimensions ``("repeat", "repeat_items")``. The ``repeat_items``
    dimension should have size 2: the first value is the "key", the second
    value is the "value". For the "key" datetime, the data of the "value"
    datetime will be used. Can also be set with a dictionary using the
    ``set_repeat_stress`` method.

imod.mf6.WellDisStructured Class Members
========================================
   * imod.mf6.WellDisStructured.clip_box
   * imod.mf6.WellDisStructured.from_file
   * imod.mf6.WellDisStructured.get_non_grid_data
   * imod.mf6.WellDisStructured.is_empty
   * imod.mf6.WellDisStructured.mask
   * imod.mf6.WellDisStructured.regrid_like
   * imod.mf6.WellDisStructured.render
   * imod.mf6.WellDisStructured.set_repeat_stress
   * imod.mf6.WellDisStructured.to_netcdf
   * imod.mf6.WellDisStructured.write

imod.mf6.WellDisStructured.clip_box
===================================
Clip a package by a bounding box (time, layer, y, x).

Slicing intervals may be half-bounded, by providing None:

* To select 500.0 <= x <= 1000.0:
  ``clip_box(x_min=500.0, x_max=1000.0)``.
* To select x <= 1000.0: ``clip_box(x_min=None, x_max=1000.0)``
  or ``clip_box(x_max=1000.0)``.
* To select x >= 500.0: ``clip_box(x_min = 500.0, x_max=None.0)``
  or ``clip_box(x_min=1000.0)``.

Parameters
----------
time_min: optional
time_max: optional
layer_min: optional, int
layer_max: optional, int
x_min: optional, float
x_max: optional, float
y_min: optional, float
y_max: optional, float
top: optional, GridDataArray
bottom: optional, GridDataArray
state_for_boundary: optional, GridDataArray

Returns
-------
sliced : Package

imod.mf6.WellDisStructured.from_file
====================================
Loads an imod mf6 package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.mf6.Package.dataset.to_netcdf(),
as the checks upon package initialization are not done again!

Parameters
----------
path : str, pathlib.Path
    Path to the file.
**kwargs : keyword arguments
    Arbitrary keyword arguments forwarded to ``xarray.open_dataset()``, or
    ``xarray.open_zarr()``.
Refer to the examples.

Returns
-------
package : imod.mf6.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.mf6.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.mf6.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.mf6.WellDisVertices
========================
WEL package for discretization by vertices (DISV) models. Any number of WEL
Packages can be specified for a single groundwater flow model.
https://water.usgs.gov/water-resources/software/MODFLOW-6/mf6io_6.0.4.pdf#page=63

.. warning::
    This class is deprecated and will be deleted in a future release.
    Consider changing your code to use the ``imod.mf6.Well`` package.

Parameters
----------
layer: list of int
    Modellayer in which the well is located.
cell2d: list of int
    Cell in which the well is located.
rate: float or list of floats
    is the volumetric well rate. A positive value indicates well (injection)
    and a negative value indicates discharge (extraction) (q).
concentration: array of floats (xr.DataArray, optional)
    if this flow package is used in simulations also involving transport,
    then this array is used as the  concentration for inflow over this
    boundary.
concentration_boundary_type: ({"AUX", "AUXMIXED"}, optional)
    if this flow package is used in simulations also involving transport,
    then this keyword specifies how outflow over this boundary is computed.
print_input: ({True, False}, optional)
    keyword to indicate that the list of well information will be written to
    the listing file immediately after it is read. Default is False.
print_flows: ({True, False}, optional)
    Indicates that the list of well flow rates will be printed to the
    listing file for every stress period time step in which "BUDGET PRINT"
    is specified in Output Control. If there is no Output Control option and
    PRINT FLOWS is specified, then flow rates are printed for the last time
    step of each stress period. Default is False.
save_flows: ({True, False}, optional)
    Indicates that well flow terms will be written to the file specified
    with "BUDGET FILEOUT" in Output Control. Default is False.
observations: [Not yet supported.]
    Default is None.
validate: {True, False}
    Flag to indicate whether the package should be validated upon
    initialization. This raises a ValidationError if package input is
    provided in the wrong manner. Defaults to True.

imod.mf6.WellDisVertices Class Members
======================================
   * imod.mf6.WellDisVertices.clip_box
   * imod.mf6.WellDisVertices.from_file
   * imod.mf6.WellDisVertices.get_non_grid_data
   * imod.mf6.WellDisVertices.is_empty
   * imod.mf6.WellDisVertices.mask
   * imod.mf6.WellDisVertices.regrid_like
   * imod.mf6.WellDisVertices.render
   * imod.mf6.WellDisVertices.set_repeat_stress
   * imod.mf6.WellDisVertices.to_netcdf
   * imod.mf6.WellDisVertices.write

imod.mf6.WellDisVertices.clip_box
=================================
Clip a package by a bounding box (time, layer, y, x).

Slicing intervals may be half-bounded, by providing None:

* To select 500.0 <= x <= 1000.0:
  ``clip_box(x_min=500.0, x_max=1000.0)``.
* To select x <= 1000.0: ``clip_box(x_min=None, x_max=1000.0)``
  or ``clip_box(x_max=1000.0)``.
* To select x >= 500.0: ``clip_box(x_min = 500.0, x_max=None.0)``
  or ``clip_box(x_min=1000.0)``.

Parameters
----------
time_min: optional
time_max: optional
layer_min: optional, int
layer_max: optional, int
x_min: optional, float
x_max: optional, float
y_min: optional, float
y_max: optional, float
top: optional, GridDataArray
bottom: optional, GridDataArray
state_for_boundary: optional, GridDataArray

Returns
-------
clipped: Package

imod.mf6.WellDisVertices.from_file
==================================
Loads an imod mf6 package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.mf6.Package.dataset.to_netcdf(),
as the checks upon package initialization are not done again!

Parameters
----------
path : str, pathlib.Path
    Path to the file.
**kwargs : keyword arguments
    Arbitrary keyword arguments forwarded to ``xarray.open_dataset()``, or
    ``xarray.open_zarr()``.
Refer to the examples.

Returns
-------
package : imod.mf6.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.mf6.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.mf6.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.mf6.AdvectionCentral
=========================
The central-in-space weighting scheme is based on a simple
distance-weighted linear interpolation between the center of cell n and the
center of cell m to calculate solute concentration at the shared face
between cell n and cell m. Although central-in-space is a misnomer for
grids without equal spacing between connected cells, it is retained here
for consistency with nomenclature used by other MODFLOW-based transport
programs, such as MT3D.
Note: all constructor arguments will be ignored

imod.mf6.AdvectionCentral Class Members
=======================================
   * imod.mf6.AdvectionCentral.clip_box
   * imod.mf6.AdvectionCentral.from_file
   * imod.mf6.AdvectionCentral.get_non_grid_data
   * imod.mf6.AdvectionCentral.is_empty
   * imod.mf6.AdvectionCentral.mask
   * imod.mf6.AdvectionCentral.regrid_like
   * imod.mf6.AdvectionCentral.to_netcdf

imod.mf6.AdvectionCentral.from_file
===================================
Loads an imod mf6 package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.mf6.Package.dataset.to_netcdf(),
as the checks upon package initialization are not done again!

Parameters
----------
path : str, pathlib.Path
    Path to the file.
**kwargs : keyword arguments
    Arbitrary keyword arguments forwarded to ``xarray.open_dataset()``, or
    ``xarray.open_zarr()``.
Refer to the examples.

Returns
-------
package : imod.mf6.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.mf6.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.mf6.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.mf6.AdvectionCentral.mask
==============================
The mask method is irrelevant for this package , instead this method
retuns a copy of itself.

imod.mf6.AdvectionTVD
=====================
An implicit second order TVD scheme. More expensive than upstream
weighting but more robust.
Note: all constructor arguments will be ignored

imod.mf6.AdvectionTVD Class Members
===================================
   * imod.mf6.AdvectionTVD.clip_box
   * imod.mf6.AdvectionTVD.from_file
   * imod.mf6.AdvectionTVD.get_non_grid_data
   * imod.mf6.AdvectionTVD.is_empty
   * imod.mf6.AdvectionTVD.mask
   * imod.mf6.AdvectionTVD.regrid_like
   * imod.mf6.AdvectionTVD.to_netcdf

imod.mf6.AdvectionTVD.from_file
===============================
Loads an imod mf6 package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.mf6.Package.dataset.to_netcdf(),
as the checks upon package initialization are not done again!

Parameters
----------
path : str, pathlib.Path
    Path to the file.
**kwargs : keyword arguments
    Arbitrary keyword arguments forwarded to ``xarray.open_dataset()``, or
    ``xarray.open_zarr()``.
Refer to the examples.

Returns
-------
package : imod.mf6.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.mf6.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.mf6.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.mf6.AdvectionUpstream
==========================
The upstream weighting (first order upwind) scheme sets the concentration
at the cellface between two adjacent cells equal to the concentration in
the cell where the flow comes from. It surpresses oscillations.
Note: all constructor arguments will be ignored

imod.mf6.AdvectionUpstream Class Members
========================================
   * imod.mf6.AdvectionUpstream.clip_box
   * imod.mf6.AdvectionUpstream.from_file
   * imod.mf6.AdvectionUpstream.get_non_grid_data
   * imod.mf6.AdvectionUpstream.is_empty
   * imod.mf6.AdvectionUpstream.mask
   * imod.mf6.AdvectionUpstream.regrid_like
   * imod.mf6.AdvectionUpstream.to_netcdf

imod.mf6.AdvectionUpstream.from_file
====================================
Loads an imod mf6 package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.mf6.Package.dataset.to_netcdf(),
as the checks upon package initialization are not done again!

Parameters
----------
path : str, pathlib.Path
    Path to the file.
**kwargs : keyword arguments
    Arbitrary keyword arguments forwarded to ``xarray.open_dataset()``, or
    ``xarray.open_zarr()``.
Refer to the examples.

Returns
-------
package : imod.mf6.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.mf6.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.mf6.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.mf6.ConstantConcentration
==============================
Constant Concentration package.

Parameters
----------
concentration: array of floats (xr.DataArray)
    Concentration of the boundary.
print_input: ({True, False}, optional)
    keyword to indicate that the list of constant head information will
    be written to the listing file immediately after it is read. Default is
    False.
print_flows: ({True, False}, optional)
    Indicates that the list of constant head flow rates will be printed to
    the listing file for every stress period time step in which "BUDGET
    PRINT" is specified in Output Control. If there is no Output Control
    option and PRINT FLOWS is specified, then flow rates are printed for the
    last time step of each stress period.
    Default is False.
save_flows: ({True, False}, optional)
    Indicates that constant head flow terms will be written to the file
    specified with "BUDGET FILEOUT" in Output Control. Default is False.
observations: [Not yet supported.]
    Default is None.
validate: {True, False}
    Flag to indicate whether the package should be validated upon
    initialization. This raises a ValidationError if package input is
    provided in the wrong manner. Defaults to True.

imod.mf6.ConstantConcentration Class Members
============================================
   * imod.mf6.ConstantConcentration.clip_box
   * imod.mf6.ConstantConcentration.from_file
   * imod.mf6.ConstantConcentration.get_non_grid_data
   * imod.mf6.ConstantConcentration.is_empty
   * imod.mf6.ConstantConcentration.mask
   * imod.mf6.ConstantConcentration.regrid_like
   * imod.mf6.ConstantConcentration.render
   * imod.mf6.ConstantConcentration.set_repeat_stress
   * imod.mf6.ConstantConcentration.to_netcdf
   * imod.mf6.ConstantConcentration.write

imod.mf6.ConstantConcentration.from_file
========================================
Loads an imod mf6 package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.mf6.Package.dataset.to_netcdf(),
as the checks upon package initialization are not done again!

Parameters
----------
path : str, pathlib.Path
    Path to the file.
**kwargs : keyword arguments
    Arbitrary keyword arguments forwarded to ``xarray.open_dataset()``, or
    ``xarray.open_zarr()``.
Refer to the examples.

Returns
-------
package : imod.mf6.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.mf6.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.mf6.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.mf6.Dispersion
===================
Molecular Diffusion and Dispersion.

Parameters
----------
diffusion_coefficient: xr.DataArray
    effective molecular diffusion coefficient. (DIFFC)
longitudinal_horizontal: xr.DataArray
    longitudinal dispersivity in horizontal direction. If flow is strictly
    horizontal, then this is the longitudinal dispersivity that will be
    used. If flow is not strictly horizontal or strictly vertical, then the
    longitudinal dispersivity is a function of both ALH and ALV. If
    mechanical dispersion is represented (by specifying any dispersivity
    values) then this array is required. (ALH)
transverse_horizontal1: xr.DataArray
    transverse dispersivity in horizontal direction. This is the transverse
    dispersivity value for the second ellipsoid axis. If flow is strictly
    horizontal and directed in the x direction (along a row for a regular
    grid), then this value controls spreading in the y direction.
    If mechanical dispersion is represented (by specifying any dispersivity
    values) then this array is required. (ATH1)
longitudinal_vertical: xr.DataArray, optional
    longitudinal dispersivity in vertical direction. If flow is strictly
    vertical, then this is the longitudinal dispsersivity value that will be
    used. If flow is not strictly horizontal or strictly vertical, then the
    longitudinal dispersivity is a function of both ALH and ALV. If this
    value is not specified and mechanical dispersion is represented, then
    this array is set equal to ALH. (ALV)
transverse_horizontal2: xr.DataArray, optional
    transverse dispersivity in horizontal direction. This is the transverse
    dispersivity value for the third ellipsoid axis. If flow is strictly
    horizontal and directed in the x direction (along a row for a regular
    grid), then this value controls spreading in the z direction. If this
    value is not specified and mechanical dispersion is represented, then
    this array is set equal to ATH1. (ATH2)
tranverse_vertical: xr.DataArray, optional
    transverse dispersivity when flow is in vertical direction. If flow is
    strictly vertical and directed in the z direction, then this value
    controls spreading in the x and y directions. If this value is not
    specified and mechanical dispersion is represented, then this array is
    set equal to ATH2. (ATV)
xt3d_off: bool, optional
    deactivate the xt3d method and use the faster and less accurate
    approximation. (XT3D_OFF)
xt3d_rhs: bool, optional
    add xt3d terms to right-hand side, when possible. This option uses less
    memory, but may require more iterations. (XT3D_RHS)
validate: {True, False}
    Flag to indicate whether the package should be validated upon
    initialization. This raises a ValidationError if package input is
    provided in the wrong manner. Defaults to True.

imod.mf6.Dispersion Class Members
=================================
   * imod.mf6.Dispersion.clip_box
   * imod.mf6.Dispersion.from_file
   * imod.mf6.Dispersion.get_non_grid_data
   * imod.mf6.Dispersion.is_empty
   * imod.mf6.Dispersion.mask
   * imod.mf6.Dispersion.regrid_like
   * imod.mf6.Dispersion.to_netcdf

imod.mf6.Dispersion.from_file
=============================
Loads an imod mf6 package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.mf6.Package.dataset.to_netcdf(),
as the checks upon package initialization are not done again!

Parameters
----------
path : str, pathlib.Path
    Path to the file.
**kwargs : keyword arguments
    Arbitrary keyword arguments forwarded to ``xarray.open_dataset()``, or
    ``xarray.open_zarr()``.
Refer to the examples.

Returns
-------
package : imod.mf6.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.mf6.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.mf6.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.mf6.ImmobileStorageTransfer
================================
The Immobile Storage and Transfer (IST) package represents an immobile
fraction of groundwater. Any number of IST Packages can be specified for a
single GWT model. This allows the user to specify triple porosity systems,
or systems with as many immobile domains as necessary.

Parameters
----------
initial_immmobile_concentration : array of floats (xr.DataArray)
    initial concentration of the immobile domain in mass per length cubed
    (cim).
immobile_porosity : array of floats (xr.DataArray)
    porosity of the immobile domain specified as the volume of immobile
    pore space per total volume (dimensionless) (thetaim).
mobile_immobile_mass_transfer_rate: array of floats (xr.DataArray)
    mass transfer rate coefficient between the mobile and immobile domains,
    in dimensions of per time (zetaim).
decay: array of floats (xr.DataArray).
    is the rate coefficient for first or zero-order decay for the aqueous
    phase of the immobile domain. A negative value indicates solute
    production. The dimensions of decay for first-order decay is one over
    time. The dimensions of decay for zero-order decay is mass per length
    cubed per time. Decay will have no effect on simulation results unless
    either first- or zero-order decay is specified in the options block.
decay_sorbed: array of floats (xr.DataArray)
    is the rate coefficient for first or zero-order decay for the sorbed
    phase of the immobile domain. A negative value indicates solute
    production. The dimensions of decay_sorbed for first-order decay is one
    over time. The dimensions of decay_sorbed for zero-order decay is mass
    of solute per mass of aquifer per time. If decay_sorbed is not
    specified and both decay and sorption are active, then the program will
    terminate with an error. decay_sorbed will have no effect on simulation
    results unless the SORPTION keyword and either first- or zero-order
    decay are specified in the options block.
bulk_density: array of floats (xr.DataArray)
    is the bulk density of the aquifer in mass per length cubed.
    bulk_density will have no effect on simulation results unless the
    SORPTION keyword is specified in the options block.
distribution_coefficient: array of floats (xr.DataArray)
    is the distribution coefficient for the equilibrium-controlled linear
    sorption isotherm in dimensions of length cubed per mass. distcoef will
    have no effect on simulation results unless the SORPTION keyword is
    specified in the options block.
save_flows: ({True, False}, optional)
    Indicates that drain flow terms will be written to the file specified
    with "BUDGET FILEOUT" in Output Control. Default is False.
budgetbinfile:
    name of the binary output file to write budget information.
budgetcsvfile:
    name of the comma-separated value (CSV) output file to write budget
    summary information. A budget summary record will be written to this
    file for each time step of the simulation.
sorption: ({True, False}, optional)
    is a text keyword to indicate that sorption will be activated. Use of
    this keyword requires that BULK_DENSITY and DISTCOEF are specified in
    the GRIDDATA block. The linear sorption isotherm is the only isotherm
    presently supported in the IST Package.
first_order_decay: ({True, False}, optional)
    is a text keyword to indicate that first-order decay will occur. Use of
    this keyword requires that DECAY and DECAY_SORBED (if sorption is
    active) are specified in the GRIDDATA block.
zero_order_decay: ({True, False}, optional)
    is a text keyword to indicate that zero-order decay will occur. Use of
    this keyword requires that DECAY and DECAY_SORBED (if sorption is
    active) are specified in the GRIDDATA block.
cimfile: (str)
    name of the output file to write immobile concentrations. This file is
    a binary file that has the same format and structure as a binary head
    and concentration file. The value for the text variable written to the
    file is CIM. Immobile domain concentrations will be written to this
    file at the same interval as mobile domain concentrations are saved, as
    specified in the GWT Model Output Control file.
columns: (int, optional), default is 10
    number of columns for writing data.
width: (int, optional), default is 10
    width for writing each number.
digits: (int, optional), default is 7
    number of digits to use for writing a number.
format: (str, optional) default exponential
    One of {"exponential", "fixed", "general", "scientific"}.
validate: {True, False}
    Flag to indicate whether the package should be validated upon
    initialization. This raises a ValidationError if package input is
    provided in the wrong manner. Defaults to True.

imod.mf6.ImmobileStorageTransfer Class Members
==============================================
   * imod.mf6.ImmobileStorageTransfer.clip_box
   * imod.mf6.ImmobileStorageTransfer.from_file
   * imod.mf6.ImmobileStorageTransfer.get_non_grid_data
   * imod.mf6.ImmobileStorageTransfer.is_empty
   * imod.mf6.ImmobileStorageTransfer.mask
   * imod.mf6.ImmobileStorageTransfer.regrid_like
   * imod.mf6.ImmobileStorageTransfer.to_netcdf

imod.mf6.ImmobileStorageTransfer.from_file
==========================================
Loads an imod mf6 package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.mf6.Package.dataset.to_netcdf(),
as the checks upon package initialization are not done again!

Parameters
----------
path : str, pathlib.Path
    Path to the file.
**kwargs : keyword arguments
    Arbitrary keyword arguments forwarded to ``xarray.open_dataset()``, or
    ``xarray.open_zarr()``.
Refer to the examples.

Returns
-------
package : imod.mf6.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.mf6.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.mf6.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.mf6.MobileStorageTransfer
==============================
Mobile Storage.

Parameters
----------
porosity: array of floats (xr.DataArray)
    volume of interconnected voids per volume of rock (percentage).
decay : array of floats (xr.DataArray, optional)
    is the rate coefficient for first or zero-order decay for the aqueous phase of the mobile domain.
    A negative value indicates solute production. The dimensions of decay for first-order decay is one
    over time. The dimensions of decay for zero-order decay is mass per length cubed per time. decay will
    have no effect on simulation results unless either first- or zero-order decay is specified in the
    options block.
decay_sorbed : array of floats (xr.DataArray, optional)
    is the rate coefficient for first or zero-order decay for the sorbed phase of the mobile domain.
    A negative value indicates solute production. The dimensions of decay_sorbed for first-order decay
    is one over time. The dimensions of decay_sorbed for zero-order decay is mass of solute per mass of
    aquifer per time. If decay_sorbed is not specified and both decay and sorption are active, then the
    program will terminate with an error. decay_sorbed will have no effect on simulation results unless
    the SORPTION keyword and either first- or zero-order decay are specified in the options block.
bulk_density : array of floats (xr.DataArray, optional)
    is the bulk density of the aquifer in mass per length cubed. bulk_density is not required unless
    the SORPTION keyword is specified.
distcoef  : array of floats (xr.DataArray, optional)
    is the distribution coefficient for the equilibrium-controlled linear sorption isotherm in dimensions
    of length cubed per mass. distcoef is not required unless the SORPTION keyword is specified.
sp2  : array of floats (xr.DataArray, optional)
    is the exponent for the Freundlich isotherm and the sorption capacity for the Langmuir isotherm.
save_flows: ({True, False}, optional)
    Indicates that recharge flow terms will be written to the file specified
    with "BUDGET FILEOUT" in Output Control.
    Default is False.
zero_order_decay: bool, optional
    Requires decay to be specified
first_order_decay: bool, optional
    Requires decay to be specified
sorption: ({linear, freundlich, langmuir}, optional)
    Type of sorption, if any.
validate: {True, False}
    Flag to indicate whether the package should be validated upon
    initialization. This raises a ValidationError if package input is
    provided in the wrong manner. Defaults to True.

imod.mf6.MobileStorageTransfer Class Members
============================================
   * imod.mf6.MobileStorageTransfer.clip_box
   * imod.mf6.MobileStorageTransfer.from_file
   * imod.mf6.MobileStorageTransfer.get_non_grid_data
   * imod.mf6.MobileStorageTransfer.is_empty
   * imod.mf6.MobileStorageTransfer.mask
   * imod.mf6.MobileStorageTransfer.regrid_like
   * imod.mf6.MobileStorageTransfer.to_netcdf

imod.mf6.MobileStorageTransfer.from_file
========================================
Loads an imod mf6 package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.mf6.Package.dataset.to_netcdf(),
as the checks upon package initialization are not done again!

Parameters
----------
path : str, pathlib.Path
    Path to the file.
**kwargs : keyword arguments
    Arbitrary keyword arguments forwarded to ``xarray.open_dataset()``, or
    ``xarray.open_zarr()``.
Refer to the examples.

Returns
-------
package : imod.mf6.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.mf6.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.mf6.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.mf6.MassSourceLoading
==========================
Mass Source Loading (SRC) package for structured discretization (DIS)
models. Any number of SRC Packages can be specified for a single
groundwater flow model.
https://water.usgs.gov/water-resources/software/MODFLOW-6/mf6io_6.3.0.pdf#page=202

Parameters
----------
rate: xr.DataArray of floats
    is the mass source loading rate. A positive value indicates addition of
    solute mass and a negative value indicates removal of solute mass
    (smassrate).
print_input: ({True, False}, optional), default is False.
    keyword to indicate that the list of mass source information will be
    written to the listing file immediately after it is read.
print_flows: ({True, False}, optional), default is False.
    keyword to indicate that the list of mass source flow rates will be
    printed to the listing file for every stress period time step in which
    "BUDGET PRINT" is specified in Output Control. If there is no Output
    Control option and "PRINT FLOWS" is specified, then flow rates are
    printed for the last time step of each stress period.
save_flows: ({True, False}, optional)
    Indicates that the mass source flow terms will be written to the file specified
    with "BUDGET FILEOUT" in Output Control.
    Default is False.
observations: [Not yet supported.]
    Default is None.
validate: {True, False}
    Flag to indicate whether the package should be validated upon
    initialization. This raises a ValidationError if package input is
    provided in the wrong manner. Defaults to True.

imod.mf6.MassSourceLoading Class Members
========================================
   * imod.mf6.MassSourceLoading.clip_box
   * imod.mf6.MassSourceLoading.from_file
   * imod.mf6.MassSourceLoading.get_non_grid_data
   * imod.mf6.MassSourceLoading.is_empty
   * imod.mf6.MassSourceLoading.mask
   * imod.mf6.MassSourceLoading.regrid_like
   * imod.mf6.MassSourceLoading.render
   * imod.mf6.MassSourceLoading.set_repeat_stress
   * imod.mf6.MassSourceLoading.to_netcdf
   * imod.mf6.MassSourceLoading.write

imod.mf6.MassSourceLoading.from_file
====================================
Loads an imod mf6 package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.mf6.Package.dataset.to_netcdf(),
as the checks upon package initialization are not done again!

Parameters
----------
path : str, pathlib.Path
    Path to the file.
**kwargs : keyword arguments
    Arbitrary keyword arguments forwarded to ``xarray.open_dataset()``, or
    ``xarray.open_zarr()``.
Refer to the examples.

Returns
-------
package : imod.mf6.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.mf6.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.mf6.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.mf6.SourceSinkMixing
=========================
Parameters
----------
package_names: array_like of str
concentration_boundary_type: array_like of str
auxiliary_variable_name: array_like of str
print_flows: bool
save_flows: bool
validate: bool

imod.mf6.SourceSinkMixing Class Members
=======================================
   * imod.mf6.SourceSinkMixing.clip_box
   * imod.mf6.SourceSinkMixing.from_file
   * imod.mf6.SourceSinkMixing.from_flow_model
   * imod.mf6.SourceSinkMixing.get_non_grid_data
   * imod.mf6.SourceSinkMixing.is_empty
   * imod.mf6.SourceSinkMixing.mask
   * imod.mf6.SourceSinkMixing.regrid_like
   * imod.mf6.SourceSinkMixing.render
   * imod.mf6.SourceSinkMixing.set_repeat_stress
   * imod.mf6.SourceSinkMixing.to_netcdf
   * imod.mf6.SourceSinkMixing.write

imod.mf6.SourceSinkMixing.from_file
===================================
Loads an imod mf6 package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.mf6.Package.dataset.to_netcdf(),
as the checks upon package initialization are not done again!

Parameters
----------
path : str, pathlib.Path
    Path to the file.
**kwargs : keyword arguments
    Arbitrary keyword arguments forwarded to ``xarray.open_dataset()``, or
    ``xarray.open_zarr()``.
Refer to the examples.

Returns
-------
package : imod.mf6.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.mf6.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.mf6.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.mf6.SourceSinkMixing.from_flow_model
=========================================
Derive a Source and Sink Mixing package from a Groundwater Flow model's
boundary conditions (e.g. GeneralHeadBoundary). Boundary condition
packages which have the ``concentration`` variable set are included.

Parameters
----------
model: GroundwaterFlowModel
    Groundwater flow model from which sources & sinks have to be
    inferred.
species: str
    Name of species to create a transport model for. This name will be
    looked for in the ``species`` dimensions of the ``concentration``
    argument.
print_flows: ({True, False}, optional)
    Indicates that the list of general head boundary flow rates will be
    printed to the listing file for every stress period time step in which
    "BUDGET PRINT" is specified in Output Control. If there is no Output
    Control option and PRINT FLOWS is specified, then flow rates are printed
    for the last time step of each stress period.
    Default is False.
save_flows: ({True, False}, optional)
    Indicates that general head boundary flow terms will be written to the
    file specified with "BUDGET FILEOUT" in Output Control.
    Default is False.
validate: ({True, False}, optional)
    Flag to indicate whether the package should be validated upon
    initialization. This raises a ValidationError if package input is
    provided in the wrong manner. Defaults to True.
is_split:  ({True, False}, optional)
    Flag to indicate if the simulation has been split into partitions.

imod.mf6.SourceSinkMixing.render
================================
Render fills in the template only, doesn't write binary data

imod.mf6.regrid.ConstantHeadRegridMethod
========================================
Object containing regridder methods for the :class:`imod.mf6.ConstantHead`
package. This can be provided to the ``regrid_like`` method to regrid with
custom settings.

Parameters
----------
head: tuple, default (RegridderType.OVERLAP, "mean")
concentration: tuple, default (RegridderType.OVERLAP, "mean")

Examples
--------
Regrid with custom settings:

>>> regrid_method = ConstantHeadRegridMethod(head=(RegridderType.BARYCENTRIC,))
>>> chd.regrid_like(target_grid, RegridderWeightsCache(), regrid_method)

The RegridderType.OVERLAP and RegridderType.RELATIVEOVERLAP require an extra
method as string.

>>> regrid_method = ConstantHeadRegridMethod(head=(RegridderType.OVERLAP, "max",))

imod.mf6.regrid.DiscretizationRegridMethod
==========================================
Object containing regridder methods for the
:class:`imod.mf6.StructuredDiscretization` and
:class:`imod.mf6.VerticesDiscretization` packages. This can be provided to
the ``regrid_like`` method to regrid with custom settings.

Parameters
----------
top: tuple, default (RegridderType.OVERLAP, "mean")
bottom: tuple, default (RegridderType.OVERLAP, "mean")
idomain: tuple, default (RegridderType.OVERLAP, "mode")

Examples
--------
Regrid with custom settings:

>>> regrid_method = DiscretizationRegridMethod(top=(RegridderType.BARYCENTRIC,))
>>> dis.regrid_like(target_grid, RegridderWeightsCache(), regrid_method)

The RegridderType.OVERLAP and RegridderType.RELATIVEOVERLAP require an extra
method as string.

>>> regrid_method = DiscretizationRegridMethod(top=(RegridderType.OVERLAP, "max",))

imod.mf6.regrid.DispersionRegridMethod
======================================
Object containing regridder methods for the :class:`imod.mf6.Dispersion`
package. This can be provided to the ``regrid_like`` method to regrid with
custom settings.

Parameters
----------
diffusion_coefficient: tuple, default (RegridderType.OVERLAP, "mean")
longitudinal_horizontal: tuple, default (RegridderType.OVERLAP, "mean")
transversal_horizontal: tuple, default (RegridderType.OVERLAP, "mean")
longitudinal_vertical: tuple, default (RegridderType.OVERLAP, "mean")
transversal_horizontal2: tuple, default (RegridderType.OVERLAP, "mean")
transversal_vertical: tuple, default (RegridderType.OVERLAP, "mean")

Examples
--------
Regrid with custom settings:

>>> regrid_method = DispersionRegridMethod(longitudinal_horizontal=(RegridderType.BARYCENTRIC,))
>>> dsp.regrid_like(target_grid, RegridderWeightsCache(), regrid_method)

The RegridderType.OVERLAP and RegridderType.RELATIVEOVERLAP require an extra
method as string.

>>> regrid_method = DispersionRegridMethod(longitudinal_horizontal=(RegridderType.OVERLAP, "max",))

imod.mf6.regrid.DrainageRegridMethod
====================================
Object containing regridder methods for the :class:`imod.mf6.Drainage`
package. This can be provided to the ``regrid_like`` method to regrid with
custom settings.

Parameters
----------
elevation: tuple, default (RegridderType.OVERLAP, "mean")
conductance: tuple, default (RegridderType.RELATIVEOVERLAP,"conductance")
concentration: tuple, default (RegridderType.OVERLAP, "mean")

Examples
--------
Regrid with custom settings:

>>> regrid_method = DrainageRegridMethod(elevation=(RegridderType.BARYCENTRIC,))
>>> drn.regrid_like(target_grid, RegridderWeightsCache(), regrid_method)

The RegridderType.OVERLAP and RegridderType.RELATIVEOVERLAP require an extra
method as string.

>>> regrid_method = DrainageRegridMethod(elevation=(RegridderType.OVERLAP, "max",))

imod.mf6.regrid.EmptyRegridMethod
=================================
Base class for protocol classes.

Protocol classes are defined as::

    class Proto(Protocol):
        def meth(self) -> int:
            ...

Such classes are primarily used with static type checkers that recognize
structural subtyping (static duck-typing).

For example::

    class C:
        def meth(self) -> int:
            return 0

    def func(x: Proto) -> int:
        return x.meth()

    func(C())  # Passes static type check

See PEP 544 for details. Protocol classes decorated with
@typing.runtime_checkable act as simple-minded runtime protocols that check
only the presence of given attributes, ignoring their type signatures.
Protocol classes can be generic, they are defined as::

    class GenProto[T](Protocol):
        def meth(self) -> T:
            ...

imod.mf6.regrid.EvapotranspirationRegridMethod
==============================================
Object containing regridder methods for the
:class:`imod.mf6.Evapotranspiration` package. This can be provided to the
``regrid_like`` method to regrid with custom settings.

Parameters
----------
surface: tuple, default (RegridderType.OVERLAP, "mean")
rate: tuple, default (RegridderType.OVERLAP, "mean")
depth: tuple, default (RegridderType.OVERLAP, "mean")
proportion_rate: tuple, default (RegridderType.OVERLAP, "mean")
proportion_depth: tuple, default (RegridderType.OVERLAP, "mean")

Examples
--------
Regrid with custom settings:

>>> regrid_method = EvapotranspirationRegridMethod(surface=(RegridderType.BARYCENTRIC,))
>>> evt.regrid_like(target_grid, RegridderWeightsCache(), regrid_method)

The RegridderType.OVERLAP and RegridderType.RELATIVEOVERLAP require an extra
method as string.

>>> regrid_method = EvapotranspirationRegridMethod(surface=(RegridderType.OVERLAP, "max",))

imod.mf6.regrid.GeneralHeadBoundaryRegridMethod
===============================================
Object containing regridder methods for the
:class:`imod.mf6.GeneralHeadBoundary` package. This can be provided to the
``regrid_like`` method to regrid with custom settings.

Parameters
----------
head: tuple, default (RegridderType.OVERLAP, "mean")
conductance: tuple, default (RegridderType.RELATIVEOVERLAP,"conductance")
concentration: tuple, default (RegridderType.OVERLAP, "mean")

Examples
--------
Regrid with custom settings:

>>> regrid_method = GeneralHeadBoundaryRegridMethod(head=(RegridderType.BARYCENTRIC,))
>>> ghb.regrid_like(target_grid, RegridderWeightsCache(), regrid_method)

The RegridderType.OVERLAP and RegridderType.RELATIVEOVERLAP require an extra
method as string.

>>> regrid_method = GeneralHeadBoundaryRegridMethod(head=(RegridderType.OVERLAP, "max",))

imod.mf6.regrid.InitialConditionsRegridMethod
=============================================
Object containing regridder methods for the
:class:`imod.mf6.InitialConditions` package. This can be provided to the
``regrid_like`` method to regrid with custom settings.

Parameters
----------
start: tuple, default (RegridderType.OVERLAP, "mean")

Examples
--------
Regrid with custom settings:

>>> regrid_method = InitialConditionsRegridMethod(start=(RegridderType.BARYCENTRIC,))
>>> ic.regrid_like(target_grid, RegridderWeightsCache(), regrid_method)

The RegridderType.OVERLAP and RegridderType.RELATIVEOVERLAP require an extra
method as string.

>>> regrid_method = InitialConditionsRegridMethod(start=(RegridderType.OVERLAP, "max",))

imod.mf6.regrid.MobileStorageTransferRegridMethod
=================================================
Object containing regridder methods for the
:class:`imod.mf6.MobileStorageTransfer` package. This can be provided to the
``regrid_like`` method to regrid with custom settings.

Parameters
----------
porosity: tuple, default (RegridderType.OVERLAP, "mean")
decay: tuple, default (RegridderType.OVERLAP, "mean")
decay_sorbed: tuple, default (RegridderType.OVERLAP, "mean")
bulk_density: tuple, default (RegridderType.OVERLAP, "mean")
distcoef: tuple, default (RegridderType.OVERLAP, "mean")
sp2: tuple, default (RegridderType.OVERLAP, "mean")

Examples
--------
Regrid with custom settings:

>>> regrid_method = MobileStorageTransferRegridMethod(porosity=(RegridderType.BARYCENTRIC,))
>>> mst.regrid_like(target_grid, RegridderWeightsCache(), regrid_method)

The RegridderType.OVERLAP and RegridderType.RELATIVEOVERLAP require an extra
method as string.

>>> regrid_method = MobileStorageTransferRegridMethod(porosity=(RegridderType.OVERLAP, "max",))

imod.mf6.regrid.NodePropertyFlowRegridMethod
============================================
Object containing regridder methods for the
:class:`imod.mf6.NodePropertyFlow` package. This can be provided to the
``regrid_like`` method to regrid with custom settings.

Parameters
----------
icelltype: tuple, defaults (RegridderType.OVERLAP, "mean")
k: tuple, defaults ( RegridderType.OVERLAP,"geometric_mean")
k22: tuple, defaults (RegridderType.OVERLAP,"geometric_mean")
k33: tuple, defaults (RegridderType.OVERLAP,"harmonic_mean")
angle1: tuple, defaults (RegridderType.OVERLAP, "mean")
angle2: tuple, defaults (RegridderType.OVERLAP, "mean")
angle3: tuple, defaults (RegridderType.OVERLAP, "mean")
rewet_layer: tuple, defaults (RegridderType.OVERLAP, "mean")

Examples
--------
Regrid with custom settings:

>>> regrid_method = NodePropertyFlowRegridMethod(k=(RegridderType.BARYCENTRIC,))
>>> npf.regrid_like(target_grid, RegridderWeightsCache(), regrid_method)

The RegridderType.OVERLAP and RegridderType.RELATIVEOVERLAP require an extra
method as string.

>>> regrid_method = NodePropertyFlowRegridMethod(k=(RegridderType.OVERLAP, "max",))

imod.mf6.regrid.RechargeRegridMethod
====================================
Object containing regridder methods for the :class:`imod.mf6.Recharge`
package. This can be provided to the ``regrid_like`` method to regrid with
custom settings.

Parameters
----------
rate: tuple, defaults (RegridderType.OVERLAP, "mean")
concentration: tuple, defaults (RegridderType.OVERLAP, "mean")

Examples
--------
Regrid with custom settings:

>>> regrid_method = RechargeRegridMethod(rate=(RegridderType.BARYCENTRIC,))
>>> rch.regrid_like(target_grid, RegridderWeightsCache(), regrid_method)

The RegridderType.OVERLAP and RegridderType.RELATIVEOVERLAP require an extra
method as string.

>>> regrid_method = RechargeRegridMethod(rate=(RegridderType.OVERLAP, "max",))

imod.mf6.regrid.RiverRegridMethod
=================================
Object containing regridder methods for the :class:`imod.mf6.River` package.
This can be provided to the ``regrid_like`` method to regrid with custom
settings.

Parameters
----------
stage: tuple, default (RegridderType.OVERLAP, "mean")
conductance: tuple, default (RegridderType.RELATIVEOVERLAP, "conductance")
bottom_elevation: tuple, default (RegridderType.OVERLAP, "mean")
concentration: tuple, default (RegridderType.OVERLAP, "mean")

Examples
--------
Regrid with custom settings:

>>> regrid_method = RiverRegridMethod(stage=(RegridderType.BARYCENTRIC,))
>>> riv.regrid_like(target_grid, RegridderWeightsCache(), regrid_method)

The RegridderType.OVERLAP and RegridderType.RELATIVEOVERLAP require an extra
method as string.

>>> regrid_method = RiverRegridMethod(stage=(RegridderType.OVERLAP, "max",))

imod.mf6.regrid.SpecificStorageRegridMethod
===========================================
Object containing regridder methods for the
:class:`imod.mf6.SpecificStorage` package. This can be provided to the
``regrid_like`` method to regrid with custom settings.

Parameters
----------
convertible: tuple, default (RegridderType.OVERLAP, "mode")
specific_storage: tuple, default (RegridderType.OVERLAP, "mean")
specific_yield: tuple, default (RegridderType.OVERLAP, "mean")

Examples
--------
Regrid with custom settings:

>>> regrid_method = SpecificStorageRegridMethod(specific_storage=(RegridderType.BARYCENTRIC,))
>>> sto.regrid_like(target_grid, RegridderWeightsCache(), regrid_method)

The RegridderType.OVERLAP and RegridderType.RELATIVEOVERLAP require an extra
method as string.

>>> regrid_method = SpecificStorageRegridMethod(specific_storage=(RegridderType.OVERLAP, "max",))

imod.mf6.regrid.StorageCoefficientRegridMethod
==============================================
Object containing regridder methods for the
:class:`imod.mf6.StorageCoefficient` package. This can be provided to the
``regrid_like`` method to regrid with custom settings.

Parameters
----------
convertible: tuple, default (RegridderType.OVERLAP, "mode")
storage_coefficient: tuple, default (RegridderType.OVERLAP, "mean")
specific_yield: tuple, default (RegridderType.OVERLAP, "mean")

Examples
--------
Regrid with custom settings:

>>> regrid_method = StorageCoefficientRegridMethod(storage_coefficient=(RegridderType.BARYCENTRIC,))
>>> sto.regrid_like(target_grid, RegridderWeightsCache(), regrid_method)

The RegridderType.OVERLAP and RegridderType.RELATIVEOVERLAP require an extra
method as string.

>>> regrid_method = StorageCoefficientRegridMethod(storage_coefficient=(RegridderType.OVERLAP, "max",))

