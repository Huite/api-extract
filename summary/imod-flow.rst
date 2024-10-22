imod.flow.ImodflowModel
=======================
Class representing iMODFLOW model input. Running it requires iMOD5.

`Download iMOD5 here <https://oss.deltares.nl/web/imod/download-imod5>`_

Attributes
----------
modelname : str check : str, optional
    When to perform model checks {None, "defer", "eager"}. Defaults to
    "defer".

Examples
--------

>>> m = Imodflow("example")
>>> m["riv"] = River(...)
>>> # ...etc.
>>> m.create_time_discretization(endtime)
>>> m.write()

imod.flow.ImodflowModel Class Members
=====================================
   * imod.flow.ImodflowModel.clear
   * imod.flow.ImodflowModel.create_time_discretization
   * imod.flow.ImodflowModel.get
   * imod.flow.ImodflowModel.items
   * imod.flow.ImodflowModel.keys
   * imod.flow.ImodflowModel.pop
   * imod.flow.ImodflowModel.popitem
   * imod.flow.ImodflowModel.render
   * imod.flow.ImodflowModel.setdefault
   * imod.flow.ImodflowModel.update
   * imod.flow.ImodflowModel.values
   * imod.flow.ImodflowModel.write

imod.flow.ImodflowModel.clear
=============================
D.clear() -> None.  Remove all items from D.

imod.flow.ImodflowModel.create_time_discretization
==================================================
Collect all unique times from model packages and additional given `times`. These
unique times are used as stress periods in the model. All stress packages must
have the same starting time.

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


with the stress periods defined between these dates. I.e. the model times are the set of all times you include in the model.

Parameters
----------
times : str, datetime; or iterable of str, datetimes.
    Times to add to the time discretization. At least one single time
    should be given, which will be used as the ending time of the
    simulation.

Examples
--------
Add a single time:

>>> m.create_time_discretization("2001-01-01")

Add a daterange:

>>> m.create_time_discretization(pd.daterange("2000-01-01", "2001-01-01"))

Add a list of times:

>>> m.create_time_discretization(["2000-01-01", "2001-01-01"])

imod.flow.ImodflowModel.get
===========================
D.get(k[,d]) -> D[k] if k in D, else d.  d defaults to None.

imod.flow.ImodflowModel.items
=============================
D.items() -> a set-like object providing a view on D's items

imod.flow.ImodflowModel.keys
============================
D.keys() -> a set-like object providing a view on D's keys

imod.flow.ImodflowModel.pop
===========================
D.pop(k[,d]) -> v, remove specified key and return the corresponding value.
If key is not found, d is returned if given, otherwise KeyError is raised.

imod.flow.ImodflowModel.popitem
===============================
D.popitem() -> (k, v), remove and return some (key, value) pair
as a 2-tuple; but raise KeyError if D is empty.

imod.flow.ImodflowModel.render
==============================
Render the runfile as a string, package by package.

imod.flow.ImodflowModel.setdefault
==================================
D.setdefault(k[,d]) -> D.get(k,d), also set D[k]=d if k not in D

imod.flow.ImodflowModel.update
==============================
D.update([E, ]**F) -> None.  Update D from mapping/iterable E and F.
If E present and has a .keys() method, does:     for k in E: D[k] = E[k]
If E present and lacks .keys() method, does:     for (k, v) in E: D[k] = v
In either case, this is followed by: for k, v in F.items(): D[k] = v

imod.flow.ImodflowModel.values
==============================
D.values() -> an object providing a view on D's values

imod.flow.ImodflowModel.write
=============================
Writes model input files.

Parameters
----------
directory : str, pathlib.Path
    Directory into which the model input will be written. The model
    input will be written into a directory called modelname.
result_dir : str, pathlib.Path
    Path to directory in which output will be written when running the
    model. Is written as the value of the ``result_dir`` key in the
    runfile. See the examples.
resultdir_is_workdir: boolean, optional
    Wether the set all input paths in the runfile relative to the output
    directory. Because iMOD-wq generates a number of files in its
    working directory, it may be advantageous to set the working
    directory to a different path than the runfile location.
convert_to: str
    The type of object to convert the projectfile to in the
    configuration ini file. Should be one of ``["mf2005_namfile",
    "mf6_namfile", "runfile"]``.

Returns
-------
None

Examples
--------
Say we wish to write the model input to a file called input, and we
desire that when running the model, the results end up in a directory
called output. We may run:

>>> model.write(directory="input", result_dir="output")

And in the ``config_run.ini``, a value of ``../../output`` will be
written for ``result_dir``. This ``config_run.ini`` has to be called
with iMOD 5 to convert the model projectfile to a Modflow 2005 namfile.
To specify a conversion to a runfile, run:

>>> model.write(directory="input", convert_to="runfile")

You can then run the following command to convert the projectfile to a runfile:

>>> path/to/iMOD5.exe ./input/config_run.ini

`Download iMOD5 here <https://oss.deltares.nl/web/imod/download-imod5>`_

imod.flow.Boundary
==================
Specify the locations of active, inactive, and specified head in cells

Parameters
----------
ibound: xr.DataArray of ints
    is the boundary variable with dimensions ``("layer", "y", "x")``.

    * If IBOUND(J,I,K) < 0, cell J,I,K has a constant head.
    * If IBOUND(J,I,K) = 0, cell J,I,K is inactive.
    * If IBOUND(J,I,K) > 0, cell J,I,K is active.

imod.flow.Boundary Class Members
================================
   * imod.flow.Boundary.compose
   * imod.flow.Boundary.from_file

imod.flow.Boundary.compose
==========================
Composes package, not useful for boundary conditions

Parameters
----------
directory : str
    Path to working directory, where files will be written.
    Necessary to generate the paths for the projectfile.
globaltimes : list #TODO make this an *arg, change order.
    Not used, only included to comply with BoundaryCondition.compose
nlayer : int
    Number of layers
**ignored
    Contains keyword arguments unused for packages

imod.flow.Boundary.from_file
============================
Loads an imod-flow package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.flow.Package.save(),
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
package : imod.flow.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.flow.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.flow.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.flow.Top
=============
The top of the aquifers

Parameters
----------
top: xr.DataArray of floats
    is the top elevation with dimensions ``("layer", "y", "x")``. For the
    common situation in which the top layer represents a water-table
    aquifer, it may be reasonable to set`top` equal to land-surface
    elevation.  The DataArray should at least include the `layer`
    dimension.

imod.flow.Top Class Members
===========================
   * imod.flow.Top.compose
   * imod.flow.Top.from_file

imod.flow.Top.from_file
=======================
Loads an imod-flow package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.flow.Package.save(),
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
package : imod.flow.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.flow.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.flow.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.flow.Bottom
================
The bottom of the aquifers

Parameters
----------
bottom: xr.DataArray of floats
    is the bottom elevation of model layers or Quasi-3d confining beds,
    with dimensions ``("layer", "y", "x")``. The DataArray should at least
    include the `layer` dimension.

imod.flow.Bottom Class Members
==============================
   * imod.flow.Bottom.compose
   * imod.flow.Bottom.from_file

imod.flow.Bottom.from_file
==========================
Loads an imod-flow package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.flow.Package.save(),
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
package : imod.flow.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.flow.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.flow.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.flow.PreconditionedConjugateGradientSolver
===============================================
The Preconditioned Conjugate Gradient Solver is used to solve the finite
difference equations in each step of a MODFLOW stress period.

Parameters
----------
max_iter: int
    is the maximum number of outer iterations - that is, calss to the
    solutions routine (MXITER). For a linear problem max_iter should be 1, unless
    more than 50 inner iterations are required, when max_iter could be as
    large as 10. A larger number (generally less than 100) is required for a
    nonlinear problem.
inner_iter: int
    is the number of inner iterations (iter1). For nonlinear problems,
    inner_iter usually ranges from 10 to 30; a value of 30 will be
    sufficient for most linear problems.
rclose: float
    is the residual criterion for convergence, in units of cubic length per
    time. The units for length and time are the same as established for all
    model data. When the maximum absolute value of the residual at all nodes
    during an iteration is less than or equal to RCLOSE, and the criterion
    for HCLOSE is also satisfied (see below), iteration stops.

    Default value: 100.0. **Nota bene**: this is aimed at regional modelling.
    For detailed studies (e.g. lab experiments) much smaller values can be
    required.
    Very general rule of thumb: should be less than 10% of smallest cell volume.
hclose: float
    is the head change criterion for convergence, in units of length. When
    the maximum absolute value of head change from all nodes during an
    iteration is less than or equal to HCLOSE, and the criterion for RCLOSE
    is also satisfied (see above), iteration stops.
    Default value: 1.0e-4. **Nota bene**: This is aimed at regional modelling, `
    for detailed studies (e.g. lab experiments) much smaller values can be
    required.
relax: float, optional
    is the relaxation parameter used. Usually, RELAX = 1.0, but for some
    problems a value of 0.99, 0.98, or 0.97 will reduce the number of
    iterations required for convergence.
    Default value: 0.98.
matrix_conditioning_method: int, optional
    the flag used to select the matrix conditioning method
        1 is for Modified Incomplete Cholesky (for use on scalar computers)
        2 is for Polynomial (for use on vector computers or to conserve computer memory)
damp: float, optional
    the damping factor.
    It is typically set equal to one, which indicates
    no damping. A value less than 1 and greater than 0 causes damping. DAMP
    does not affect inner iterations; instead, it affects outer iterations.
    Default value: 1.0.
damp_transient: float, optional
    the damping factor for transient stress periods.
    it is read only when `damp` is specified as a negative value.
    If damp_transient is not read, then the single damping factor,
    `damp`, is used for both transient and steady-state stress periods.
printout_interval: int, optional
    is the printout interval for PCG.
    If equal to zero, it is changed to 999.
    The maximum head change (positive or negative) and
    residual change are printed for each iteration of a time step
    whenever the time step is an even multiple of printout_interval.
    This printout also occurs at the end of each stress period
    regardless of the value of printout_interval.
print_convergence_info: int, optional
    a flag that controls printing of convergence information from the solver:
        0 is for printing tables of maximum head change and residual each iteration
        1 is for printing only the total number of iterations
        2 is for no printing
        3 is for printing only if convergence fails

imod.flow.PreconditionedConjugateGradientSolver Class Members
=============================================================
   * imod.flow.PreconditionedConjugateGradientSolver.compose
   * imod.flow.PreconditionedConjugateGradientSolver.from_file

imod.flow.PreconditionedConjugateGradientSolver.from_file
=========================================================
Loads an imod-flow package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.flow.Package.save(),
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
package : imod.flow.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.flow.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.flow.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.flow.ConstantHead
======================
The Constant Head package. The Time-Variant Specified-Head package is used
to simulate specified head boundaries that can change between
stress periods.

Parameters
----------
head: xr.DataArray of floats
    is the head at the boundary

imod.flow.ConstantHead Class Members
====================================
   * imod.flow.ConstantHead.compose
   * imod.flow.ConstantHead.from_file
   * imod.flow.ConstantHead.periodic_stress
   * imod.flow.ConstantHead.repeat_stress

imod.flow.ConstantHead.compose
==============================
Composes all variables for one system.

Parameters
----------
globaltimes : list, np.array
    Holds the global times, i.e. the combined unique times of every
    boundary condition that are used to define the stress periods.  The
    times of the BoundaryCondition do not have to match all the global
    times. When a globaltime is not present in the BoundaryCondition,
    the value of the first previous available time is filled in. The
    effective result is a forward fill in time.
directory : str
    Path to working directory, where files will be written.  Necessary
    to generate the paths for the projectfile.
nlayer : int
    Number of layers
system_index : int
    System number. Defaults to 1, but for package groups it is used
pkggroup_times : optional, list, np.array
    Holds the package_group times.  Packages in one group need to be
    synchronized for iMODFLOW.

Returns
-------
composition : defaultdict
    A nested dictionary containing following the projectfile hierarchy:
    ``{_pkg_id : {stress_period : {varname : {system_index : {lay_nr : value}}}}}``
    where 'value' can be a scalar or a path to a file.
    The stress period number may be the wildcard '?'.

imod.flow.ConstantHead.from_file
================================
Loads an imod-flow package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.flow.Package.save(),
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
package : imod.flow.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.flow.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.flow.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.flow.ConstantHead.periodic_stress
======================================
Periodic stress periods.

Adds periodic stresses to each variable in the package.  iMODFLOW will
then implicitly repeat these.

The iMOD manual currently states: 'A PERIOD repeats until another time
definition is more close to the current time step'.

Parameters
----------
periods: dict of datetime - string pairs
    contains a datetime as key which maps to a period label.  This will
    be used for the entire package.
use_cftime: bool
    Whether to force datetimes to cftime or not.

Examples
--------
The following example assigns a higher head to the summer period than
winter period.  iMODFLOW will switch to period "summer" once
'xxxx-04-01' has passed, and "winter" once 'xxxx-10-01' has passed.

>>> times = [np.datetime64('2000-04-01'), np.datetime64('2000-10-01')]

>>> head_periodic = xr.DataArray([2., 1.], coords={"time": times}, dims=["time"])

>>> timemap = {times[0]: "summer", times[1]: "winter"}

>>> ghb = GeneralHeadBoundary(head = head_periodic, conductance = 10.)
>>> ghb.periodic_stress(timemap)

imod.flow.ConstantHead.repeat_stress
====================================
Repeat stress periods.

Parameters
----------
use_cftime: bool
    Whether to force datetimes to cftime or not.
**repeats: dict of datetime - datetime pairs
    keyword argument with variable name as keyword and
    a dict as value. This dict contains a datetime as key
    which maps to another already existing datetime in the
    BoundaryCondition, for which data should be repeated.

imod.flow.Drain
===============
The Drain package is used to simulate head-dependent flux boundaries. In
the Drain package if the head in the cell falls below a certain threshold,
the flux from the drain to the model cell drops to zero.

Parameters
----------
elevation: float or xr.DataArray of floats
    elevation of the drain, dims ``("layer", "y", "x")``.
conductance: float or xr.DataArray of floats
    is the conductance of the drain, dims ``("layer", "y", "x")``.

imod.flow.Drain Class Members
=============================
   * imod.flow.Drain.compose
   * imod.flow.Drain.from_file
   * imod.flow.Drain.periodic_stress
   * imod.flow.Drain.repeat_stress

imod.flow.Drain.from_file
=========================
Loads an imod-flow package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.flow.Package.save(),
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
package : imod.flow.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.flow.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.flow.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.flow.EvapoTranspiration
============================
Recharge provides a fixed flux boundary condition to the top layer of the
groundwater system. Note that unlike in iMOD-WQ, there is only the option
in iMODFLOW to apply the recharge package to the top layer.

Parameters
----------
rate: float or xr.DataArray of floats
    evaporation rate in mm/day (NOTA BENE!), dims ``("time", "y", "x")``.
top_elevation: floats or xr.DataArray of floats
    Top elevation in m+MSL for maximal evapotranspiration strength.
extinction_depth: float or xr.Datarray of floats
    Depth [m] in which evapotranspiration strength reduced to zero.

imod.flow.EvapoTranspiration Class Members
==========================================
   * imod.flow.EvapoTranspiration.compose
   * imod.flow.EvapoTranspiration.from_file
   * imod.flow.EvapoTranspiration.periodic_stress
   * imod.flow.EvapoTranspiration.repeat_stress

imod.flow.EvapoTranspiration.compose
====================================
Composes all variables for one system.

Parameters
----------
globaltimes : list, np.array
    Holds the global times, i.e. the combined unique times of every
    boundary condition that are used to define the stress periods.  The
    times of the BoundaryCondition do not have to match all the global
    times. When a globaltime is not present in the BoundaryCondition,
    the value of the first previous available time is filled in. The
    effective result is a forward fill in time.
directory : str
    Path to working directory, where files will be written.  Necessary
    to generate the paths for the projectfile.
nlayer : int
    Number of layers
system_index : int
    System number. Defaults to 1, but for package groups it is used
pkggroup_times : optional, list, np.array
    Holds the package_group times.  Packages in one group need to be
    synchronized for iMODFLOW.

Returns
-------
composition : defaultdict
    A nested dictionary containing following the projectfile hierarchy:
    ``{_pkg_id : {stress_period : {varname : {system_index : {lay_nr : value}}}}}``
    where 'value' can be a scalar or a path to a file.
    The stress period number may be the wildcard '?'.

imod.flow.EvapoTranspiration.from_file
======================================
Loads an imod-flow package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.flow.Package.save(),
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
package : imod.flow.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.flow.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.flow.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.flow.GeneralHeadBoundary
=============================
The General-Head Boundary package is used to simulate head-dependent flux
boundaries. In the General-Head Boundary package the flux is always
proportional to the difference in head.

Parameters
----------
head: float or xr.DataArray of floats
    head value for the GHB (BHEAD), dims ``("layer", "y", "x")``.
conductance: float or xr.DataArray of floats
    the conductance of the GHB (COND), dims ``("layer", "y", "x")``.

imod.flow.GeneralHeadBoundary Class Members
===========================================
   * imod.flow.GeneralHeadBoundary.compose
   * imod.flow.GeneralHeadBoundary.from_file
   * imod.flow.GeneralHeadBoundary.periodic_stress
   * imod.flow.GeneralHeadBoundary.repeat_stress

imod.flow.GeneralHeadBoundary.from_file
=======================================
Loads an imod-flow package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.flow.Package.save(),
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
package : imod.flow.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.flow.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.flow.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.flow.HorizontalHydraulicConductivity
=========================================
Horizontal hydraulic conductivity [L/T] of the aquifers, between TOP and
BOT.

This variable behaves somewhat similar to the horizontal hydraulic
conductivity in MODFLOW 2005's "Layer Property Flow" schematization.

Note however that this does not hold for the vertical hydraulic
conductivity: iMODFLOW uses the vertical hydraulic conductivity to specify
the hydraulic conductivity of aquitards (between BOT and TOP)

Parameters
----------
k_horizontal : xr.DataArray
    Horizontal hydraulic conductivity, dims ``("layer", "y", "x")``.

imod.flow.HorizontalHydraulicConductivity Class Members
=======================================================
   * imod.flow.HorizontalHydraulicConductivity.compose
   * imod.flow.HorizontalHydraulicConductivity.from_file

imod.flow.HorizontalHydraulicConductivity.from_file
===================================================
Loads an imod-flow package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.flow.Package.save(),
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
package : imod.flow.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.flow.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.flow.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.flow.HorizontalAnisotropy
==============================
Horizontal anisotropy is a phenomenon in which the horizontal hydraulic
conductivity is not equal along the x and y Cartesian axes. iMODFLOW can
calculate this anisotropy based on a anisotropy factor and an anisotropy
angle. iMODFLOW also accounts for the cross-terms in the horizontal
hydraulic conductivity tensor.

See also section 12.14 "ANI Horizontal anisotropy module" in the iMOD v5.2
manual for further explanation.

Parameters
----------
anisotropy_factor : xr.DataArray
    The anisotropy factor is defined perpendicular to the main principal
    axis. The factor is between 0.0 (full anisotropic) and 1.0 (full isotropic)
anisotropy_angle : xr.DataArray
    The anisotropy angle denotes the angle along the main principal axis
    (highest permeability k) measured in degrees from
    north (0°), east (90°), south (180°) and west (270°).

imod.flow.HorizontalAnisotropy Class Members
============================================
   * imod.flow.HorizontalAnisotropy.compose
   * imod.flow.HorizontalAnisotropy.from_file

imod.flow.HorizontalAnisotropy.from_file
========================================
Loads an imod-flow package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.flow.Package.save(),
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
package : imod.flow.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.flow.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.flow.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.flow.HorizontalFlowBarrier
===============================
Horizontal barriers obstructing flow such as semi- or impermeable fault
zone or a sheet pile wall are defined for each model layer by a `*.GEN` line
file.

Parameters
----------
id_name: str or list of str
    name of the barrier
geometry: object array of shapely LineStrings
    geometry of barriers, should be lines
resistance: float or list of floats
    resistance of the barrier (d).
layer: Optional, int
    layer where barrier is located. Defaults to None.

imod.flow.HorizontalFlowBarrier Class Members
=============================================
   * imod.flow.HorizontalFlowBarrier.compose
   * imod.flow.HorizontalFlowBarrier.from_file

imod.flow.HorizontalFlowBarrier.from_file
=========================================
Loads an imod-flow package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.flow.Package.save(),
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
package : imod.flow.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.flow.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.flow.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.flow.MetaSwap
==================
The MetaSWAP package (CAP), provides the input to be converted to a
MetaSWAP model, which is an external model code used to simulate the
unsaturated zone.

Note that only two-dimensional DataArrays with dimensions ``("y", "x")``
should be supplied to this package.  In the current implementation
time-related files are provided as external files ("extra files"). Similar
to the iMODFLOW implementation of the projectfile. For now these need to be
provided as a path.

MetaSWAP is developed by Alterra, Wageningen as part of the SIMGRO model
code. The SIMGRO framework is intended for regions with an undulating
topography and unconsolidated sediments in the (shallow) subsoil. Both
shallow and deep groundwater levels can be modelled by MetaSWAP. This
model is based on a simplification of ‘straight Richards’, meaning that no
special processes like hysteresis, preferential flow and bypass flow are
modelled. Snow is not modelled, and neither the influence of frost on the
soil water conductivity. A perched watertable can be present in the SVAT
column model, but interflow is not modelled. Processes that are typical for
steep slopes are not included. The code contains several parameterized
water management schemes, including irrigation and water level management.

References:

* Van Walsum, P. E. V., 2017a. SIMGRO V7.3.3.2, Input and Output reference
  manual. Tech. Rep.  Alterra-Report 913.3, Alterra, Wageningen. 98 pp.

* Van Walsum, P. E. V., 2017b. SIMGRO V7.3.3.2, Users Guide. Tech. Rep.
  Alterra-Report 913.2, Alterra, Wageningen. 111 pp.

* Van Walsum, P. E. V. and P. Groenendijk, 2008. "Quasi Steady-State
  Simulation on of the Unsaturated Zone in Groundwater Modeling of Lowland
  Regions." Vadose Zone Journal 7: 769-778.

* Van Walsum, P. E. V., A. A. Veldhuizen and P. Groenendijk, 2016. SIMGRO
  V7.2.27, Theory and model implementation. Tech. Rep. Alterra-Report 913.1,
  Alterra, Wageningen. 93 pp 491.

Parameters
----------
boundary : int or xr.DataArray of ints
    2D boundary, used to specify active MetaSWAP elements, similar to
    ibound in the Boundary package

landuse : int or xr.DataArray of ints
    Landuse codes, referred to in the lookup table file luse_mswp.inp

rootzone_thickness : float or xr.DataArray of floats
    Rootzone thickness in cm (min. value is 10 centimeter).

soil_physical_unit : int or xr.DataArray of ints
    Soil Physical Unit, referred to in the lookup table file fact_mswp.inp.

meteo_station_number : float or xr.DataArray of ints
    Meteo station number, referred to by mete_mswp.inp.

surface_elevation : float or xr.DataArray of floats
    Surface Elevation (m+MSL)

sprinkling_type : int or xr.DataArray of ints
    Sprinkling type ("Artificial Recharge Type" in iMOD manual):

    * 0 = no occurrence
    * 1 = from groundwater
    * 2 = from surface water

sprinkling_layer : int or xr.DataArray of ints
    Number of modellayer from which water is extracted ("Artificial
            Recharge Location" in iMOD manual)

sprinkling_capacity : float or xr.DataArray of floats
    Sprinkling capacity (mm/d) sets the maximum amount extracted for
    sprinkling ("Artificial Recharge Capacity" in iMOD manual)

wetted_area : float or xr.DataArray of floats
    Total area (m2) occupied by surface water elements.  Values will be
    truncated by maximum cellsize.

urban_area : float or xr.DataArray of floats
    Total area (m2) occupied by urban area.  Values will be truncated by
    maximum cellsize.

ponding_depth_urban : float or xr.DataArray of floats
    Ponding Depth Urban Area (m), specifying the acceptable depth of the
    ponding of water on the surface in the urban area before surface runoff
    occurs.

ponding_depth_rural : float or xr.DataArray of floats
    Ponding Depth Rural Area (m), specifying the acceptable depth of the
    ponding of water on the surface in the rural area before surface runoff
    occurs.

runoff_resistance_urban : float or xr.DataArray of floats
    Runoff Resistance Urban Area (day), specifying the resistance surface
    flow encounters in the urban area. The minimum value is equal to the
    model time period.

runoff_resistance_rural : float or xr.DataArray of floats
    Runoff Resistance Rural Area (day), specifying the resistance surface
    flow encounters in the rural area. The minimum value is equal to the
    model time period.

runon_resistance_urban : float or xr.DataArray of floats
    Runon Resistance Urban Area (day), specifying the resistance surface
    flow encounters to a model cell from an adjacent cell in the urban
    area. The minimum value is equal to the model time period.

runon_resistance_rural : float or xr.DataArray of floats
    Runon Resistance Rural Area (day), specifying the resistance surface
    flow encounters to a model cell from an adjacent cell in the rural
    area. The minimum value is equal to the model time period.

infiltration_capacity_urban : float or xr.DataArray of floats
    the infiltration capacity (m/d) of the soil surface in the urban area.
    The range is 0-1000 m/d. The NoDataValue -9999 indicates unlimited
    infiltration is possible.

infiltration_capacity_rural : float or xr.DataArray of floats
    the infiltration capacity (m/d) of the soil surface in the urban area.
    The range is 0-1000 m/d. The NoDataValue -9999 indicates unlimited
    infiltration is possible.

perched_water_table : float or xr.DataArray of floats
    Depth of the perched water table level (m)

soil_moisture_factor : float
    Soil Moisture Factor to adjust the soil moisture coefficient. This
    factor may be used during calibration. Default value is 1.0.

conductivity_factor : float
    Conductivity Factor to adjust the vertical conductivity. This factor
    may be used during calibration. Default value is 1.0.

lookup_and_forcing_files : list of pathlib.Path or str
    List of paths to files required by MetaSWAP. This a list of
    lookup tables and meteorological information that is required by
    MetaSwap. Note that MetaSwap looks for files with a specific name, so
    calling "luse_svat.inp" something else will result in errors. To view
    the files required, you can call: ``print(MetaSwap()._required_extra)``

imod.flow.MetaSwap Class Members
================================
   * imod.flow.MetaSwap.check_lookup_and_forcing_files
   * imod.flow.MetaSwap.compose
   * imod.flow.MetaSwap.from_file

imod.flow.MetaSwap.check_lookup_and_forcing_files
=================================================
Check for presence of required MetaSWAP input files.

imod.flow.MetaSwap.from_file
============================
Loads an imod-flow package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.flow.Package.save(),
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
package : imod.flow.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.flow.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.flow.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.flow.Recharge
==================
Recharge provides a fixed flux boundary condition to the top layer of the
groundwater system.  Note that unlike in iMOD-WQ, there is only the option
in iMODFLOW to apply the recharge package to the top layer.

Parameters
----------
rate: float or xr.DataArray of floats
    recharge rate in mm/day (NOTA BENE!), dims ``("time", "y", "x")``.

imod.flow.Recharge Class Members
================================
   * imod.flow.Recharge.compose
   * imod.flow.Recharge.from_file
   * imod.flow.Recharge.periodic_stress
   * imod.flow.Recharge.repeat_stress

imod.flow.Recharge.from_file
============================
Loads an imod-flow package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.flow.Package.save(),
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
package : imod.flow.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.flow.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.flow.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.flow.River
===============
The River package is used to simulate head-dependent flux boundaries. In
the River package if the head in the cell falls below a certain threshold,
the flux from the river to the model cell is set to a specified lower
bound.

Parameters
----------
stage: float or xr.DataArray of floats
    is the head in the river (STAGE), dims ``("layer", "y", "x")``.
bottom_elevation: float or xr.DataArray of floats
    is the bottom of the riverbed (RBOT), dims ``("layer", "y", "x")``.
conductance: float or xr.DataArray of floats
    is the conductance of the river, dims ``("layer", "y", "x")``.
infiltration_factor: float or xr.DataArray of floats
    is the infiltration factor, dims ``("layer", "y", "x")``. This factor
    defines the reduces the conductance for infiltrating water compared to
    exfiltrating water:

    ``cond = A / (c * iff)``

    where ``A`` [L2] is the area where surface water and groundwater
    interact, ``c``  [L] is the resistance, and ``iff`` is the infiltration
    factor.

    The infiltration factor is thus equal or larger than 1.

imod.flow.River Class Members
=============================
   * imod.flow.River.compose
   * imod.flow.River.from_file
   * imod.flow.River.periodic_stress
   * imod.flow.River.repeat_stress

imod.flow.River.from_file
=========================
Loads an imod-flow package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.flow.Package.save(),
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
package : imod.flow.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.flow.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.flow.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.flow.StartingHead
======================
The initial head in all cells

Parameters
----------
starting_head: float or xr.DataArray of floats
    is initial (starting) head—that is, head at the beginning of the
    simulation (SHD). starting_head must be specified for all simulations,
    including steady-state simulations. One value is read for every model
    cell. Usually, these values are read a layer at a time.

imod.flow.StartingHead Class Members
====================================
   * imod.flow.StartingHead.compose
   * imod.flow.StartingHead.from_file

imod.flow.StartingHead.from_file
================================
Loads an imod-flow package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.flow.Package.save(),
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
package : imod.flow.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.flow.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.flow.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.flow.VerticalHydraulicConductivity
=======================================
Vertical hydraulic conductivity [L/T] for aquitards (between BOT and TOP).

Note that this is different from MODFLOW 2005's "Layer Property Flow"
schematization.  To specify the vertical hydraulic conductivity for
aquifers, use VerticalAnisotropy in combination with
HorizontalHydraulicConductivity.

Parameters
----------
k_vertical : xr.DataArray
    Vertical hydraulic conductivity, dims ``("layer", "y", "x")``.

imod.flow.VerticalHydraulicConductivity Class Members
=====================================================
   * imod.flow.VerticalHydraulicConductivity.compose
   * imod.flow.VerticalHydraulicConductivity.from_file

imod.flow.VerticalHydraulicConductivity.from_file
=================================================
Loads an imod-flow package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.flow.Package.save(),
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
package : imod.flow.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.flow.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.flow.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.flow.VerticalAnisotropy
============================
Vertical anisotropy for aquifers [-], defined as the horizontal hydraulic
conductivity over the vertical hydraulic conductivity.

Use this package in combination with HorizontalHydraulicConductivity to
specify the vertical hydraulic conductivity.

Parameters
----------
vertical_anisotropy : xr.DataArray
    Vertical anisotropy factor (Kv/Kh), dims ``("layer", "y", "x")``.

imod.flow.VerticalAnisotropy Class Members
==========================================
   * imod.flow.VerticalAnisotropy.compose
   * imod.flow.VerticalAnisotropy.from_file

imod.flow.VerticalAnisotropy.from_file
======================================
Loads an imod-flow package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.flow.Package.save(),
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
package : imod.flow.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.flow.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.flow.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.flow.StorageCoefficient
============================
Storage coefficient [-].  Be careful, this is not the same as the specific
storage.

From wikipedia (https://en.wikipedia.org/wiki/Specific_storage):

Storativity or the storage coefficient is the volume of water released
from storage per unit decline in hydraulic head in the aquifer, per
unit area of the aquifer.  Storativity is a dimensionless quantity, and
is always greater than 0.

Under confined conditions:

S = Ss * b, where S is the storage coefficient,
Ss the specific storage, and b the aquifer thickness.

Under unconfined conditions:

S = Sy, where Sy is the specific yield

Parameters
----------
storage_coefficient : float or xr.DataArray
    Storage coefficient, dims = ("layer", "y", "x").

imod.flow.StorageCoefficient Class Members
==========================================
   * imod.flow.StorageCoefficient.compose
   * imod.flow.StorageCoefficient.from_file

imod.flow.StorageCoefficient.from_file
======================================
Loads an imod-flow package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.flow.Package.save(),
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
package : imod.flow.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.flow.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.flow.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.flow.SpecificStorage
=========================
Specific storage [L-1].  Be careful, this is not the same as the storage
coefficient.

From wikipedia (https://en.wikipedia.org/wiki/Specific_storage):

The specific storage is the amount of water that a portion of an aquifer
releases from storage, per unit mass or volume of aquifer, per unit change
in hydraulic head, while remaining fully saturated.

Parameters
----------
specific_storage : float or xr.DataArray
    Specific storage, dims ``("layer", "y", "x")``.

imod.flow.SpecificStorage Class Members
=======================================
   * imod.flow.SpecificStorage.compose
   * imod.flow.SpecificStorage.from_file

imod.flow.SpecificStorage.from_file
===================================
Loads an imod-flow package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.flow.Package.save(),
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
package : imod.flow.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.flow.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.flow.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.flow.Well
==============
The Well package is used to simulate a specified flux to individual cells
and specified in units of length3/time.

Parameters
----------
id_name: str or list of str
    name of the well(s).
x: float or list of floats
    x coordinate of the well(s).
y: float or list of floats
    y coordinate of the well(s).
rate: float or list of floats.
    pumping rate in the well(s).
layer: "None" or int, optional
    layer from which the pumping takes place.
time: "None" or listlike of np.datetime64, datetime.datetime, pd.Timestamp,
    cftime.datetime
    time during which the pumping takes place. Only need to specify if
    model is transient.

imod.flow.Well Class Members
============================
   * imod.flow.Well.compose
   * imod.flow.Well.from_file
   * imod.flow.Well.periodic_stress
   * imod.flow.Well.repeat_stress

imod.flow.Well.from_file
========================
Loads an imod-flow package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.flow.Package.save(),
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
package : imod.flow.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.flow.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.flow.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

