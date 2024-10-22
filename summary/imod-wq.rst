imod.wq.SeawatModel
===================
iMOD-WQ SEAWAT model.

Attributes
----------
modelname : str
check : str, optional
    When to perform model checks {None, "defer", "eager"}.
    Defaults to "defer".

Examples
--------

>>> m = SeawatModel("example")
>>> m["riv"] = River(...)
>>> # ...etc.
>>> m.create_time_discretization(endtime)
>>> m.write()

imod.wq.SeawatModel Class Members
=================================
   * imod.wq.SeawatModel.clear
   * imod.wq.SeawatModel.clip
   * imod.wq.SeawatModel.create_time_discretization
   * imod.wq.SeawatModel.get
   * imod.wq.SeawatModel.items
   * imod.wq.SeawatModel.keys
   * imod.wq.SeawatModel.pop
   * imod.wq.SeawatModel.popitem
   * imod.wq.SeawatModel.render
   * imod.wq.SeawatModel.setdefault
   * imod.wq.SeawatModel.to_netcdf
   * imod.wq.SeawatModel.update
   * imod.wq.SeawatModel.values
   * imod.wq.SeawatModel.write

imod.wq.SeawatModel.clear
=========================
D.clear() -> None.  Remove all items from D.

imod.wq.SeawatModel.clip
========================
Method to clip the model to a certain `extent`. The spatial resolution of the clipped model is unchanged.
Boundary conditions of clipped model can be derived from parent model calculation results and are applied
along the edge of `extent` (CHD and TVC). Packages from parent that have no data within extent are optionally removed.

Parameters
----------
extent : tuple, geopandas.GeoDataFrame, xarray.DataArray
    Extent of the clipped model. Tuple must be in the form of (`xmin`,`xmax`,`ymin`,`ymax`). If a GeoDataFrame, all
    polygons are included in the model extent. If a DataArray, non-null/non-zero values are taken as the new extent.

heads_boundary : xarray.DataArray, optional.
    Heads to be applied as a Constant Head boundary condition along the edge of the model extent. These heads are assumed
    to be derived from calculations with the parent model. Timestamp of boundary condition is shifted to correct for difference
    between 'end of period' timestamp of results and 'start of period' timestamp of boundary condition.
    If None (default), no constant heads boundary condition is applied.

concentration_boundary : xarray.DataArray, optional.
    Concentration to be applied as a Time Varying Concentration boundary condition along the edge of the model extent.
    These concentrations can be derived from calculations with the parent model. Timestamp of boundary condition is shifted
    to correct for difference between 'end of period' timestamp of results and 'start of period' timestamp of boundary condition.
    If None (default), no time varying concentration boundary condition is applied.

    *Note that the Time Varying Concentration boundary sets a constant concentration for the entire stress period,
    unlike the linearly varying Constant Head. This will inevitably cause a time shift in concentrations along the boundary.
    This shift becomes more significant when stress periods are longer. If necessary, consider interpolating concentrations
    along the time axis, to reduce the length of stress periods (see examples).*

delete_empty_pkg : bool, optional.
    Set to True to delete packages that contain no data in the clipped model. Defaults to False.

Examples
--------
Given a full model, clip a 1 x 1 km rectangular submodel without boundary conditions along its edge:

>>> extent = (1000., 2000., 5000., 6000.)  # xmin, xmax, ymin, ymax
>>> clipped = ml.clip(extent)

Load heads and concentrations from full model results:

>>> heads = imod.idf.open("head/head_*.idf")
>>> conc = imod.idf.open("conc/conc_*.idf")
>>> clipped = ml.clip(extent, heads, conc)

Use a shape of a model area:

>>> extent = geopandas.read_file("clipped_model_area.shp")
>>> clipped = ml.clip(extent, heads, conc)

Interpolate concentration results to annual results using xarray.interp(), to improve time resolution of concentration boundary:

>>> conc = imod.idf.open("conc/conc_*.idf")
>>> dates = pd.date_range(conc.time.values[0], conc.time.values[-1], freq="AS")
>>> conc_interpolated = conc.load().interp(time=dates, method="linear")
>>> clipped = ml.clip(extent, heads, conc_interpolated)

imod.wq.SeawatModel.create_time_discretization
==============================================
Collect all unique times from model packages and additional given `times`. These
unique times are used as stress periods in the model. All stress packages must
have the same starting time.

The time discretization in imod-python works as follows:

- The datetimes of all packages you send in are always respected
- Subsequently, the input data you use is always included fully as well
- All times are treated as starting times for the stress: a stress is always applied until the next specified date
- For this reason, a final time is required to determine the length of the last stress period
- Additional times can be provided to force shorter stress periods & more detailed output
- Every stress has to be defined on the first stress period (this is a modflow requirement)

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

Note
----
To set the other parameters of the TimeDiscretization object, you have
to set these to the object after calling this function.

Examples
--------
Add a single time:

>>> m.create_time_discretization("2001-01-01")

Add a daterange:

>>> m.create_time_discretization(pd.daterange("2000-01-01", "2001-01-01"))

Add a list of times:

>>> m.create_time_discretization(["2000-01-01", "2001-01-01"])

imod.wq.SeawatModel.get
=======================
D.get(k[,d]) -> D[k] if k in D, else d.  d defaults to None.

imod.wq.SeawatModel.items
=========================
D.items() -> a set-like object providing a view on D's items

imod.wq.SeawatModel.keys
========================
D.keys() -> a set-like object providing a view on D's keys

imod.wq.SeawatModel.pop
=======================
D.pop(k[,d]) -> v, remove specified key and return the corresponding value.
If key is not found, d is returned if given, otherwise KeyError is raised.

imod.wq.SeawatModel.popitem
===========================
D.popitem() -> (k, v), remove and return some (key, value) pair
as a 2-tuple; but raise KeyError if D is empty.

imod.wq.SeawatModel.render
==========================
Render the runfile as a string, package by package.

imod.wq.SeawatModel.setdefault
==============================
D.setdefault(k[,d]) -> D.get(k,d), also set D[k]=d if k not in D

imod.wq.SeawatModel.to_netcdf
=============================
Convenience function to write all model packages
to netcdf files.

Parameters
----------
directory : str, pathlib.Path
    Directory into which the different model packages will be written.

pattern : str, optional.
    Pattern for filename of each package, in which `pkgname`
    signifies the package name. Default is `"{pkgname}.nc"`,
    so `model["river"]` would get written to `path / river.nc`.

kwargs :
    Additional kwargs to be forwarded to `xarray.Dataset.to_netcdf`.

imod.wq.SeawatModel.update
==========================
D.update([E, ]**F) -> None.  Update D from mapping/iterable E and F.
If E present and has a .keys() method, does:     for k in E: D[k] = E[k]
If E present and lacks .keys() method, does:     for (k, v) in E: D[k] = v
In either case, this is followed by: for k, v in F.items(): D[k] = v

imod.wq.SeawatModel.values
==========================
D.values() -> an object providing a view on D's values

imod.wq.SeawatModel.write
=========================
Writes model input files.

Parameters
----------
directory : str, pathlib.Path
    Directory into which the model input will be written. The model
    input will be written into a directory called modelname.
result_dir : str, pathlib.Path
    Path to directory in which output will be written when running the
    model. Is written as the value of the ``result_dir`` key in the
    runfile.

    See the examples.
resultdir_is_workdir: boolean, optional
    Wether the set all input paths in the runfile relative to the output
    directory. Because iMOD-wq generates a number of files in its working
    directory, it may be advantageous to set the working directory to
    a different path than the runfile location.

Returns
-------
None

Examples
--------
Say we wish to write the model input to a file called input, and we
desire that when running the model, the results end up in a directory
called output. We may run:

>>> model.write(directory="input", result_dir="output")

And in the runfile, a value of ``../../output`` will be written for
result_dir.

imod.wq.TimeDiscretization
==========================
Time discretisation package class.

Parameters
----------
timestep_duration: float
    is the length of the current stress period (PERLEN). If the flow
    solution is transient, timestep_duration specified here must be equal to
    that specified for the flow model. If the flow solution is steady-state,
    timestep_duration can be set to any desired length.
n_timesteps: int, optional
    is the number of time steps for the transient flow solution in the
    current stress period (NSTP). If the flow solution is steady-state,
    n_timestep=1. Default value is 1.
transient: bool, optional
    Flag indicating wether the flow simulation is transient (True) or False
    (Steady State).
    Default is True.
timestep_multiplier: float, optional
    is the multiplier for the length of successive time steps used in the
    transient flow solution (TSMULT); it is used only if n_timesteps>1.
    timestep_multiplier>0, the length of each flow time step within the
    current stress period is calculated using the geometric progression as
    in MODFLOW. Note that both n_timesteps and timestep_multiplier specified
    here must be identical to those specified in the flow model if the flow
    model is transient.
    timestep_multiplier ≤ 0, the length of each flow time step within the
    current stress period is read from the record TSLNGH. This option is
    needed in case the length of time steps for the flow solution is not
    based on a geometric progression in a flow model, unlike MODFLOW.
    Default is 1.0.
max_n_transport_timestep: int, optional
    is the maximum number of transport steps allowed within one time step of
    the flow solution (mxstrn). If the number of transport steps within a
    flow time step exceeds max_n_transport_timestep, the simulation is
    terminated.
    Default is 50_000.
transport_timestep_multiplier: float or {"None"}, optional
    is the multiplier for successive transport steps within a flow time step
    (TTSMULT).
    If the Generalized Conjugate Gradient (GCG) solver is used and the
    solution option for the advection term is the standard finite difference
    method. A value between 1.0 and 2.0 is generally adequate. If the GCG
    package is not used, the transport solution is solved explicitly as in
    the original MT3D code, and transport_timestep_multiplier is always set
    to 1.0 regardless of the user-specified input. Note that for the
    particle tracking based solution options and the 3rd-order TVD scheme,
    transport_timestep_multiplier does not apply.
    Default is {"None"}.
transport_initial_timestep: int, optional
    is the user-specified transport stepsize within each time step of the
    flow solution (DT0).
    transport_initial_timestep is interpreted differently depending on
    whether the solution option chosen is explicit or implicit: For explicit
    solutions (i.e., the GCG solver is not used), the program will always
    calculate a maximum transport stepsize which meets the various stability
    criteria. Setting transport_initial_timestep to zero causes the model
    calculated transport stepsize to be used in the simulation. However, the
    model-calculated transport_initial_timestep may not always be optimal.
    In this situation, transport_initial_timestep should be adjusted to find
    a value that leads to the best results. If transport_initial_timestep is
    given a value greater than the model-calculated stepsize, the
    model-calculated stepsize, instead of the user-specified value, will be
    used in the simulation.
    For implicit solutions (i.e., the GCG solver is used),
    transport_initial_timestep is the initial transport stepsize. If it is
    specified as zero, the model-calculated value of
    transport_initial_timestep, based on the user-specified Courant number
    in the Advection Package, will be used. The subsequent transport
    stepsize may increase or remain constant depending on the userspecified
    transport stepsize multiplier transport_timestep_multiplier and the
    solution scheme for the advection term.
    Default is 0.

imod.wq.TimeDiscretization Class Members
========================================
   * imod.wq.TimeDiscretization.from_file

imod.wq.TimeDiscretization.from_file
====================================
Loads an imod-wq package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.wq.Package.save(),
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
package : imod.wq.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.wq.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.wq.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.wq.OutputControl
=====================
The Output Control Option is used to specify if head, drawdown, or budget
data should be saved and in which format.

Parameters
----------
save_head_idf: bool, optional
    Save calculated head values in IDF format.
    Default value is False.
save_concentration_idf: bool, optional
    Save calculated concentration values in IDF format.
    Default value is False.
save_budget_idf: bool, optional
    Save calculated budget in IDF format.
    Default value is False.
save_head_tec: bool, optional
    Save calculated head values in a format compatible with Tecplot.
    Default value is False.
save_concentration_tec: bool, optional
    Save calculated concentration values in a format compatible with
    Tecplot.
    Default value is False.
save_budget_tec: bool, optional
    Save calculated budget in a format compatible with Tecplot.
    Default value is False.
save_head_vtk: bool, optional
    Save calculated head values in a format compatible with ParaView (VTK).
    Default value is False.
save_concentration_vtk: bool, optional
    Save calculated concentration values in a format compatible with
    ParaView (VTK).
    Default value is False.
save_budget_vtk: bool, optional
    Save calculated budget in a format compatible with ParaView (VTK).
    Default value is False.

imod.wq.OutputControl Class Members
===================================
   * imod.wq.OutputControl.from_file

imod.wq.OutputControl.from_file
===============================
Loads an imod-wq package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.wq.Package.save(),
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
package : imod.wq.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.wq.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.wq.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.wq.PreconditionedConjugateGradientSolver
=============================================
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
damp: float, optional
    is the damping factor. It is typically set equal to one, which indicates
    no damping. A value less than 1 and greater than 0 causes damping. DAMP
    does not affect inner iterations; instead, it affects outer iterations.
    Default value: 1.0.

imod.wq.PreconditionedConjugateGradientSolver Class Members
===========================================================
   * imod.wq.PreconditionedConjugateGradientSolver.from_file

imod.wq.PreconditionedConjugateGradientSolver.from_file
=======================================================
Loads an imod-wq package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.wq.Package.save(),
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
package : imod.wq.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.wq.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.wq.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.wq.GeneralizedConjugateGradientSolver
==========================================
The Generalized Conjugate Gradient Solver solves the matrix equations
resulting from the implicit solution of the transport equation.

Parameters
----------
max_iter: int
    is the maximum number of outer iterations (MXITER); it should be set to an
    integer greater than one (1) only when a nonlinear sorption isotherm is
    included in simulation.
iter1: int
    is the maximum number of inner iterations (iter1); a value of 30-50
    should be adequate for most problems.
isolve: int
    is the type of preconditioners to be used with the Lanczos/ORTHOMIN
    acceleration scheme:
    isolve = 1: Jacobi
    isolve = 2: SSOR
    isolve = 3: Modified Incomplete Cholesky (MIC)
    (MIC usually converges faster, but it needs significantly more memory)
lump_dispersion: bool
    is an integer flag for treatment of dispersion tensor cross terms:
    ncrs = 0: lump all dispersion cross terms to the right-hand-side
    (approximate but highly efficient).
    ncrs = 1: include full dispersion tensor (memory intensive).
cclose: float
    is the convergence criterion in terms of relative concentration; a real
    value between 10-4 and 10-6 is generally adequate.

imod.wq.GeneralizedConjugateGradientSolver Class Members
========================================================
   * imod.wq.GeneralizedConjugateGradientSolver.from_file

imod.wq.GeneralizedConjugateGradientSolver.from_file
====================================================
Loads an imod-wq package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.wq.Package.save(),
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
package : imod.wq.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.wq.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.wq.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.wq.ParallelKrylovFlowSolver
================================
The Parallel Krylov Flow Solver is used for parallel solving of the flow
model.

Parameters
----------
max_iter: int
    is the maximum number of outer iterations (MXITER); it should be set to
    an integer greater than one (1) only when a nonlinear sorption isotherm
    is included in simulation.
inner_iter: int
    is the maximum number of inner iterations (INNERIT); a value of 30-50
    should be adequate for most problems.
hclose: float
    is the head change criterion for convergence (HCLOSEPKS), in units of
    length. When the maximum absolute value of head change from all nodes
    during an iteration is less than or equal to HCLOSE, and the criterion
    for RCLOSE is also satisfied (see below), iteration stops.
rclose: float
    is the residual criterion for convergence (RCLOSEPKS), in units of cubic
    length per time. The units for length and time are the same as
    established for all model data. When the maximum absolute value of the
    residual at all nodes during an iteration is less than or equal to
    RCLOSE, and the criterion for HCLOSE is also satisfied (see above),
    iteration stops.
relax: float
    is the relaxation parameter used. Usually, RELAX = 1.0, but for some
    problems a value of 0.99, 0.98, or 0.97 will reduce the number of
    iterations required for convergence.
h_fstrict: float, optional
    is a factor to apply to HCLOSE to set a stricter hclose for the linear
    inner iterations (H_FSTRICTPKS). HCLOSE_inner is calculated as follows:
    HCLOSEPKS * H_FSTRICTPKS.
r_fstrict: float, optional
    is a factor to apply to RCLOSE to set a stricter rclose for the linear
    inner iterations (R_FSTRICTPKS). RCLOSE_inner is calculated as follows:
    RCLOSEPKS * R_FSTRICTPKS.
partition: {"uniform", "rcb"}, optional
    Partitioning option (PARTOPT). "uniform" partitions the model domain
    into equally sized subdomains. "rcb" (Recursive Coordinate Bisection)
    uses a 2D pointer grid with weights to partition the model domain.
    Default value: "uniform"
solver: {"pcg"}, optional
    Flag indicating the linear solver to be used (ISOLVER).
    Default value: "pcg"
preconditioner: {"ilu"}, optional
    Flag inicating the preconditioner to be used (NPC).
    Devault value: "ilu"
deflate: {True, False}, optional
    Flag for deflation preconditioner.
    Default value: False
debug: {True, False}, optional
    Debug option.
    Default value: False
load_balance_weight: xarray.DataArray, optional
    2D grid with load balance weights, used when partition = "rcb"
    (Recursive Coordinate Bisection). If None (default), then the module
    will create a load balance grid by summing active cells over layers:
    `(ibound != 0).sum("layer")`

    Note that even though the iMOD-SEAWAT helpfile states .idf is
    accepted, it is not. This load balance grid should be a .asc file
    (without a header). Formatting is done as follows:
    `pd.DataFrame(load_balance_weight.values).to_csv(path, sep='\t',
    header=False, index=False, float_format = "%8.2f")`

imod.wq.ParallelKrylovFlowSolver Class Members
==============================================
   * imod.wq.ParallelKrylovFlowSolver.from_file
   * imod.wq.ParallelKrylovFlowSolver.save

imod.wq.ParallelKrylovFlowSolver.from_file
==========================================
Loads an imod-wq package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.wq.Package.save(),
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
package : imod.wq.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.wq.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.wq.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.wq.ParallelKrylovFlowSolver.save
=====================================
Overloaded method to write .asc instead of .idf.
(This is an idiosyncracy of the parallel iMODwq code.)

imod.wq.ParallelKrylovTransportSolver
=====================================
The Parallel Krylov Transport Solver is used for parallel solving of the
transport model.

Parameters
----------
max_iter: int
    is the maximum number of outer iterations (MXITER); it should be set to
    an integer greater than one (1) only when a nonlinear sorption isotherm
    is included in simulation.
inner_iter: int
    is the maximum number of inner iterations (INNERIT); a value of 30-50
    should be adequate for most problems.
cclose: float, optional
    is the convergence criterion in terms of relative concentration; a real
    value between 10-4 and 10-6 is generally adequate.
    Default value: 1.0e-6.
relax: float, optional
    is the relaxation parameter used. Usually, RELAX = 1.0, but for some
    problems a value of 0.99, 0.98, or 0.97 will reduce the number of
    iterations required for convergence.
    Default value: 0.98.
partition: {"uniform", "rcb"}, optional
    Partitioning option (PARTOPT). "uniform" partitions the model domain
    into equally sized subdomains. "rcb" (Recursive Coordinate Bisection)
    uses a 2D pointer grid with weights to partition the model domain.
    Default value: "uniform".
solver: {"bicgstab", "gmres", "gcr"}, optional
    Flag indicating the linear solver to be used (ISOLVER).
    Default value: "bicgstab"
preconditioner: {"ilu"}, optional
    Flag inicating the preconditioner to be used (NPC).
    Devault value: "ilu".
debug: {True, False}, optional
    Debug option.
    Default value: False
load_balance_weight: xarray.DataArray, optional
    2D grid with load balance weights, used when partition = "rcb"
    (Recursive Coordinate Bisection). If None (default), then the module
    will create a load balance grid by summing active cells over layers:
    `(ibound != 0).sum("layer")`

    Note that even though the iMOD-SEAWAT helpfile states .idf is
    accepted, it is not. This load balance grid should be a .asc file
    (without a header). Formatting is done as follows:
    `pd.DataFrame(load_balance_weight.values).to_csv(path, sep='\t',
    header=False, index=False, float_format = "%8.2f")`

imod.wq.ParallelKrylovTransportSolver Class Members
===================================================
   * imod.wq.ParallelKrylovTransportSolver.from_file
   * imod.wq.ParallelKrylovTransportSolver.save

imod.wq.ParallelKrylovTransportSolver.from_file
===============================================
Loads an imod-wq package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.wq.Package.save(),
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
package : imod.wq.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.wq.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.wq.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.wq.BasicFlow
=================
The Basic package is used to specify certain data used in all models.
These include:

1. the locations of acitve, inactive, and specified head in cells,
2. the head stored in inactive cells,
3. the initial head in all cells, and
4. the top and bottom of the aquifer

The number of layers (NLAY) is automatically calculated using the IBOUND.
Thickness is calculated using the specified tops en bottoms.
The Basic package input file is required in all models.

Parameters
----------
ibound: xr.DataArray of integers
    is the boundary variable.
    If IBOUND(J,I,K) < 0, cell J,I,K has a constant head.
    If IBOUND(J,I,K) = 0, cell J,I,K is inactive.
    If IBOUND(J,I,K) > 0, cell J,I,K is active.
top: float or xr.DataArray of floats
    is the top elevation of layer 1. For the common situation in which the
    top layer represents a water-table aquifer, it may be reasonable to set
    `top` equal to land-surface elevation.
bottom: xr.DataArray of floats
    is the bottom elevation of model layers or Quasi-3d confining beds. The
    DataArray should at least include the `layer` dimension.
starting_head: float or xr.DataArray of floats
    is initial (starting) head—that is, head at the beginning of the
    simulation (STRT). starting_head must be specified for all simulations,
    including steady-state simulations. One value is read for every model
    cell. Usually, these values are read a layer at a time.
inactive_head: float, optional
    is the value of head to be assigned to all inactive (no flow) cells
    (IBOUND = 0) throughout the simulation (HNOFLO). Because head at
    inactive cells is unused in model calculations, this does not affect
    model results but serves to identify inactive cells when head is
    printed. This value is also used as drawdown at inactive cells if the
    drawdown option is used. Even if the user does not anticipate having
    inactive cells, a value for inactive_head must be entered. Default
    value is 1.0e30.

imod.wq.BasicFlow Class Members
===============================
   * imod.wq.BasicFlow.from_file
   * imod.wq.BasicFlow.thickness

imod.wq.BasicFlow.from_file
===========================
Loads an imod-wq package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.wq.Package.save(),
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
package : imod.wq.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.wq.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.wq.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.wq.BasicFlow.thickness
===========================
Computes layer thickness from top and bottom data.

Returns
-------
thickness : xr.DataArray

imod.wq.ConstantHead
====================
The Constant Head package. The Time-Variant Specified-Head package is used
to simulate specified head boundaries that can change within or between
stress periods.

Parameters
----------
head_start: xr.DataArray of floats
    is the head at the boundary at the start of the stress period.
head_end: xr.DataArray of floats
    is the head at the boundary at the end of the stress period.
concentration: xr.DataArray of floats
    concentrations for the constant heads. It gets automatically written to
    the SSM package.
save_budget: bool, optional
    is a flag indicating if the budget should be saved (ICHDCB).
    Default is False.

imod.wq.ConstantHead Class Members
==================================
   * imod.wq.ConstantHead.from_file

imod.wq.ConstantHead.from_file
==============================
Loads an imod-wq package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.wq.Package.save(),
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
package : imod.wq.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.wq.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.wq.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.wq.Drainage
================
The Drain package is used to simulate head-dependent flux boundaries. In the
Drain package if the head in the cell falls below a certain threshold, the
flux from the drain to the model cell drops to zero.

Parameters
----------
elevation: float or xr.DataArray of floats
    elevation of the drain.
conductance: float or xr.DataArray of floats
    is the conductance of the drain.
save_budget: bool, optional
    A flag that is used to determine if cell-by-cell budget data should be
    saved. If save_budget is True cell-by-cell budget data will be saved.
    Default is False.

imod.wq.Drainage Class Members
==============================
   * imod.wq.Drainage.from_file

imod.wq.Drainage.from_file
==========================
Loads an imod-wq package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.wq.Package.save(),
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
package : imod.wq.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.wq.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.wq.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.wq.EvapotranspirationTopLayer
==================================
Base package for (transient) boundary conditions:
* recharge
* general head boundary
* constant head
* river
* drainage

imod.wq.EvapotranspirationTopLayer Class Members
================================================
   * imod.wq.EvapotranspirationTopLayer.from_file

imod.wq.EvapotranspirationTopLayer.from_file
============================================
Loads an imod-wq package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.wq.Package.save(),
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
package : imod.wq.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.wq.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.wq.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.wq.EvapotranspirationLayers
================================
Base package for (transient) boundary conditions:
* recharge
* general head boundary
* constant head
* river
* drainage

imod.wq.EvapotranspirationLayers Class Members
==============================================
   * imod.wq.EvapotranspirationLayers.from_file

imod.wq.EvapotranspirationLayers.from_file
==========================================
Loads an imod-wq package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.wq.Package.save(),
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
package : imod.wq.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.wq.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.wq.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.wq.EvapotranspirationHighestActive
=======================================
Base package for (transient) boundary conditions:
* recharge
* general head boundary
* constant head
* river
* drainage

imod.wq.EvapotranspirationHighestActive Class Members
=====================================================
   * imod.wq.EvapotranspirationHighestActive.from_file

imod.wq.EvapotranspirationHighestActive.from_file
=================================================
Loads an imod-wq package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.wq.Package.save(),
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
package : imod.wq.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.wq.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.wq.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.wq.GeneralHeadBoundary
===========================
The General-Head Boundary package is used to simulate head-dependent flux
boundaries. In the General-Head Boundary package the flux is always
proportional to the difference in head.

Parameters
----------
head: float or xr.DataArray of floats
    head value for the GHB (BHEAD).
conductance: float or xr.DataArray of floats
    the conductance of the GHB (COND).
density: float or xr.DataArray of floats
    is the density used to convert the point head to the freshwater head
    (GHBSSMDENS).
concentration: "None" or xr.DataArray of floats, optional
    concentration of the GHB (CGHB), get automatically inserted into the SSM
    package.
    Default is "None".
save_budget: bool, optional
    is a flag indicating if the budget should be saved (IGHBCB).
    Default is False.

imod.wq.GeneralHeadBoundary Class Members
=========================================
   * imod.wq.GeneralHeadBoundary.from_file

imod.wq.GeneralHeadBoundary.from_file
=====================================
Loads an imod-wq package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.wq.Package.save(),
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
package : imod.wq.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.wq.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.wq.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.wq.LayerPropertyFlow
=========================
The Layer-Property Flow (LPF) package is used to specify properties
controlling flow between cells.

Parameters
----------
k_horizontal: float or xr.DataArray of floats
    is the hydraulic conductivity along rows (HK). HK is multiplied by
    horizontal anisotropy (see horizontal_anisotropy) to obtain hydraulic
    conductivity along columns.
k_vertical: float or xr.DataArray of floats
    is the vertical hydraulic conductivity (VKA).
horizontal_anisotropy: float or xr.DataArray of floats
    contains a value for each layer that is the horizontal anisotropy
    (CHANI). Use as many records as needed to enter a value of CHANI for
    each layer. The horizontal anisotropy is the ratio of the hydraulic
    conductivity along columns (the Y direction) to the hydraulic
    conductivity along rows (the X direction).
interblock: int
    contains a flag for each layer that defines the method of calculating
    interblock transmissivity (LAYAVG). Use as many records needed to enter
    a value for each layer.
    0 = harmonic mean (This is most appropriate for confined and unconfined
    aquifers with abrupt boundaries in transmissivity at the cell boundaries
    or for confined aquifers with uniform hydraulic conductivity).
    1 = logarithmic mean (This is most appropriate for confined aquifers
    with gradually varying transmissivities).
    2 = arithmetic mean of saturated thickness and logarithmic-mean
    hydraulic conductivity. (This is most appropriate for unconfined
    aquifers with gradually varying transmissivities).
layer_type: int
    contains a flag for each layer that specifies the layer type (LAYTYP).
    Use as many records needed to enter a value for each layer.
    0 = confined
    not 0 = convertible
specific_storage: float or xr.DataArray of floats
    is specific storage (SS). Read only for a transient simulation (at least
    one transient stress period). Include only if at least one stress period
    is transient.
    Specific storage is the amount of water released when the head in an aquifer
    drops by 1 m, in one meter of the aquifer (or model layer).
    The unit is: ((m3 / m2) / m head change) / m aquifer = m-1
specific_yield: float or xr.DataArray of floats
    is specific yield (SY). Read only for a transient simulation (at least
    one transient stress period) and if the layer is convertible (layer_type
    is not 0). Include only if at least one stress period is transient.
    The specific yield is the volume of water released from (or added to) the
    pore matrix for one meter of head change.
    The unit is: (m3 / m2) / m head change = dimensionless
save_budget: int
    is a flag and a unit number (ILPFCB).
    If save_budget > 0, it is the unit number to which cell-by-cell flow
    terms will be written when "SAVE BUDGET" or a non-zero value for
    save_budget is specified in Output Control. The terms that are saved are
    storage, constant-head flow, and flow between adjacent cells.
    If save_budget = 0, cell-by-cell flow terms will not be written.
    If save_budget < 0, cell-by-cell flow for constant-head cells will be
    written in the listing file when "SAVE BUDGET" or a non-zero value for
    ICBCFL is specified in Output Control. Cell-by-cell flow to storage and
    between adjacent cells will not be written to any file. The flow terms
    that will be saved are the flows through the right, front, and lower
    cell face. Positive values represent flows toward higher column, row, or
    layer numbers.
layer_wet: int
    contains a flag for each layer that indicates if wetting is active. Use
    as many records as needed to enter a value for each layer.
    0 = wetting is inactive
    not 0 = wetting is active
interval_wet: int
    is the iteration interval for attempting to wet cells. Wetting is
    attempted every interval_wet iteration (IWETIT). If using the PCG solver
    (Hill, 1990), this applies to outer iterations, not inner iterations. If
    interval_wet less than or equal to 0, it is changed to 1.
method_wet: int
    is a flag that determines which equation is used to define the initial
    head at cells that become wet (IHDWET).
    If method_wet = 0, this equation is used:
    h = BOT + WETFCT (hn - BOT).
    (hn is the head in the neighboring cell that is causing the dry cell to
    convert to an active cell.)
    If method_wet is not 0, this equation is used:
    h = BOT + WETFCT(THRESH).
    WETFCT is a factor that is included in the calculation of the head that
    is initially established at a cell when it is converted from dry to wet.
head_dry: float, optional
    is the head that is assigned to cells that are converted to dry during a
    simulation (HDRY). Although this value plays no role in the model calculations,
    it is useful as an indicator when looking at the resulting heads that
    are output from the model. HDRY is thus similar to HNOFLO in the Basic
    Package, which is the value assigned to cells that are no-flow cells at
    the start of a model simulation.
    Default value: 1.0e20.

imod.wq.LayerPropertyFlow Class Members
=======================================
   * imod.wq.LayerPropertyFlow.from_file

imod.wq.LayerPropertyFlow.from_file
===================================
Loads an imod-wq package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.wq.Package.save(),
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
package : imod.wq.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.wq.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.wq.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.wq.RechargeTopLayer
========================
The Recharge package is used to simulate a specified flux distributed over
the top of the model and specified in units of length/time (usually m/d). Within MODFLOW,
these rates are multiplied by the horizontal area of the cells to which they
are applied to calculate the volumetric flux rates. In this class the
Recharge gets applied to the top grid layer (NRCHOP=1).

Parameters
----------
rate: float or xr.DataArray of floats
    is the amount of recharge.
concentration: float or xr.DataArray of floats
    is the concentration of the recharge
save_budget: bool, optional
    flag indicating if the budget needs to be saved.
    Default is False.

imod.wq.RechargeTopLayer Class Members
======================================
   * imod.wq.RechargeTopLayer.from_file

imod.wq.RechargeTopLayer.from_file
==================================
Loads an imod-wq package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.wq.Package.save(),
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
package : imod.wq.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.wq.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.wq.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.wq.RechargeLayers
======================
The Recharge package is used to simulate a specified flux distributed over
the top of the model and specified in units of length/time (usually m/d). Within MODFLOW,
these rates are multiplied by the horizontal area of the cells to which they
are applied to calculate the volumetric flux rates. In this class the
Recharge gets applied to a specific, specified, layer (NRCHOP=2).

Parameters
----------
rate: float or xr.DataArray of floats
    is the amount of recharge.
recharge_layer: int or xr.DataArray of integers
    layer number variable that defines the layer in each vertical column
    where recharge is applied
concentration: float or xr.DataArray of floats
    is the concentration of the recharge
save_budget: bool, optional
    flag indicating if the budget needs to be saved.
    Default is False.

imod.wq.RechargeLayers Class Members
====================================
   * imod.wq.RechargeLayers.from_file

imod.wq.RechargeLayers.from_file
================================
Loads an imod-wq package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.wq.Package.save(),
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
package : imod.wq.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.wq.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.wq.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.wq.RechargeHighestActive
=============================
The Recharge package is used to simulate a specified flux distributed over
the top of the model and specified in units of length/time (usually m/d). Within MODFLOW,
these rates are multiplied by the horizontal area of the cells to which they
are applied to calculate the volumetric flux rates. In this class the
Recharge gets applied to the highest active cell in each vertical column
(NRCHOP=3).

Parameters
----------
rate: float or xr.DataArray of floats
    is the amount of recharge.
concentration: float or xr.DataArray of floats
    is the concentration of the recharge
save_budget: bool, optional
    flag indicating if the budget needs to be saved.
    Default is False.

imod.wq.RechargeHighestActive Class Members
===========================================
   * imod.wq.RechargeHighestActive.from_file

imod.wq.RechargeHighestActive.from_file
=======================================
Loads an imod-wq package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.wq.Package.save(),
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
package : imod.wq.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.wq.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.wq.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.wq.River
=============
The River package is used to simulate head-dependent flux boundaries. In the
River package if the head in the cell falls below a certain threshold, the
flux from the river to the model cell is set to a specified lower bound.

Parameters
----------
stage: float or xr.DataArray of floats
    is the head in the river (STAGE).
bottom_elevation: float or xr.DataArray of floats
    is the bottom of the riverbed (RBOT).
conductance: float or xr.DataArray of floats
    is the conductance of the river.
density: float or xr.DataArray of floats
    is the density used to convert the point head to the freshwater head
    (RIVSSMDENS).
concentration: "None", float or xr.DataArray of floats, optional
    is the concentration in the river.
    Default is None.
save_budget: bool, optional
    is a flag indicating if the budget should be saved (IRIVCB).
    Default is False.

imod.wq.River Class Members
===========================
   * imod.wq.River.from_file

imod.wq.River.from_file
=======================
Loads an imod-wq package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.wq.Package.save(),
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
package : imod.wq.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.wq.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.wq.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.wq.Well
============
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
    time during which the pumping takes place. Only need to specify if model
    is transient.
save_budget: bool, optional
    is a flag indicating if the budget should be saved (IRIVCB).
    Default is False.

imod.wq.Well Class Members
==========================
   * imod.wq.Well.from_file

imod.wq.Well.from_file
======================
Loads an imod-wq package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.wq.Package.save(),
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
package : imod.wq.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.wq.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.wq.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.wq.VariableDensityFlow
===========================
Variable Density Flow package.

Parameters
----------
density_species: int
    is the MT3DMS species number that will be used in the equation of state
    to compute fluid density (mtdnconc).
    If density_species = 0, fluid density is specified using items 6 and 7,
    and flow will be uncoupled with transport if the IMT Process is active.
    If density_species > 0, fluid density is calculated using the MT3DMS
    species number that corresponds with density_species. A value for
    density_species greater than zero indicates that flow will be coupled
    with transport.
    If density_species = -1, fluid density is calculated using one or more
    MT3DMS species. Items 4a, 4b, and 4c will be read instead of item 4.
density_min: float
    is the minimum fluid density (DENSEMIN). If the resulting density value
    calculated with the equation of state is less than density_min, the
    density value is set to density_min.
    If density_min = 0, the computed fluid density is not limited by
    density_min (this is the option to use for most simulations).
    If density_min > 0, a computed fluid density less than density_min is
    automatically reset to density_min.
density_max: float
    is the maximum fluid density (DENSEMAX). If the resulting density value
    calculated with the equation of state is greater than density_max, the
    density value is set to density_max.
    If density_max = 0, the computed fluid density is not limited by
    density_max (this is the option to use for most simulations).
    If density_max > 0, a computed fluid density larger than density_max is
    automatically reset to density_max.
density_ref: float
    is the fluid density at the reference concentration, temperature, and
    pressure (DENSEREF). For most simulations, density_ref is specified as
    the density of freshwater at 25 °C and at a reference pressure of zero.
    Value of 1000 is often used.
density_concentration_slope: float
    is the slope d(rho)/d(C) of the linear equation of state that relates
    fluid density to solute concentration (denseslp). Value of 0.7143 is
    often used.
density_criterion: float
    is the convergence parameter for the coupling between flow and transport
    and has units of fluid density (DNSCRIT). If the maximum density
    difference between two consecutive coupling iterations is not less than
    density_criterion, the program will continue to iterate on the flow and
    transport equations or will terminate if 'coupling' is exceeded.
coupling: int
    is a flag used to determine the flow and transport coupling procedure
    (nswtcpl).
    If coupling = 0 or 1, flow and transport will be explicitly coupled
    using a one-timestep lag. The explicit coupling option is normally much
    faster than the iterative option and is recommended for most
    applications.
    If coupling > 1, coupling is the maximum number of non-linear coupling
    iterations for the flow and transport solutions. SEAWAT-2000 will stop
    execution after coupling iterations if convergence between flow and
    transport has not occurred.
    If coupling = -1, the flow solution will be recalculated only for: The
    first transport step of the simulation, or The last transport step of
    the MODFLOW timestep, or The maximum density change at a cell is greater
    than density_criterion.
correct_water_table: bool
    is a flag used to activate the variable-density water-table corrections
    (IWTABLE).
    If correct_water_table = False, the water-table correction will not be
    applied.
    If correct_water_table = True, the water-table correction will be
    applied.
internodal: str, {"upstream", "central"}
    is a flag that determines the method for calculating the internodal
    density values used to conserve fluid mass (MFNADVFD).
    If internodal = "central", internodal conductance values used to
    conserve fluid mass are calculated using a central-in-space algorithm.
    If internodal = "upstream", internodal conductance values used to
    conserve fluid mass are calculated using an upstream-weighted algorithm.

imod.wq.VariableDensityFlow Class Members
=========================================
   * imod.wq.VariableDensityFlow.from_file

imod.wq.VariableDensityFlow.from_file
=====================================
Loads an imod-wq package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.wq.Package.save(),
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
package : imod.wq.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.wq.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.wq.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.wq.AdvectionTVD
====================
Total Variation Diminishing (TVD) formulation (ULTIMATE, MIXELM = -1).

Attributes
----------
courant : float
    Courant number (PERCEL) is the number of cells (or a fraction of a cell)
    advection will be allowed in any direction in one transport step. For
    implicit finite-difference or particle tracking based schemes, there is
    no limit on PERCEL, but for accuracy reasons, it is generally not set
    much greater than one. Note, however, that the PERCEL limit is checked
    over the entire model grid. Thus, even if PERCEL > 1, advection may not
    be more than one cell’s length at most model locations. For the explicit
    finite-difference, PERCEL is also a stability constraint, which must not
    exceed one and will be automatically reset to one if a value greater
    than one is specified.

imod.wq.AdvectionTVD Class Members
==================================
   * imod.wq.AdvectionTVD.from_file

imod.wq.AdvectionTVD.from_file
==============================
Loads an imod-wq package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.wq.Package.save(),
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
package : imod.wq.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.wq.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.wq.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.wq.AdvectionMOC
====================
Solve the advection term using the Method of Characteristics (MIXELM = 1)

Nota bene: number of particles settings have not been tested. The defaults
here are chosen conservatively, with many particles. This increases both
memory usage and computational effort.

Attributes
-----------
courant: float
    Courant number (PERCEL) is the number of cells (or a fraction of a cell)
    advection will be allowed in any direction in one transport step. For
    implicit finite-difference or particle tracking based schemes, there is
    no limit on PERCEL, but for accuracy reasons, it is generally not set
    much greater than one. Note, however, that the PERCEL limit is checked
    over the entire model grid. Thus, even if PERCEL > 1, advection may not
    be more than one cell’s length at most model locations. For the explicit
    finite-difference, PERCEL is also a stability constraint, which must not
    exceed one and will be automatically reset to one if a value greater
    than one is specified.
max_nparticles: int
    is the maximum total number of moving particles allowed (MXPART).
tracking: string {"euler", "runge-kutta", "hybrid"}, optional
    indicates which particle tracking algorithm is selected for the
    Eulerian-Lagrangian methods. ITRACK = 1, the first-order Euler algorithm is
    used; ITRACK = 2, the fourth-order Runge-Kutta algorithm is used; this
    option is computationally demanding and may be needed only when PERCEL is
    set greater than one. ITRACK = 3, the hybrid 1st and 4th order algorithm is
    used; the Runge- Kutta algorithm is used in sink/source cells and the cells
    next to sinks/sources while the Euler algorithm is used elsewhere.
    Default value is "hybrid".
weighting_factor: float, optional
    is a concentration weighting factor (WD) between 0.5 and 1. It is used for
    operator splitting in the particle tracking based methods. The value of
    0.5 is generally adequate. The value may be adjusted to achieve better
    mass balance. Generally, it can be increased toward 1.0 as advection
    becomes more dominant.
    Default value: 0.5.
dconcentration_epsilon: float, optional
    is a small Relative Cell Concentration Gradient below which advective
    transport is considered negligible. A value around 10-5 is generally
    adequate.
    Default value: 1.0e-5.
nplane: int, optional
    is a flag indicating whether the random or fixed pattern is selected for
    initial placement of moving particles. NPLANE = 0, the random pattern is
    selected for initial placement. Particles are distributed randomly in
    both the horizontal and vertical directions by calling a random number
    generator. This option is usually preferred and leads to smaller mass
    balance discrepancy in nonuniform or diverging/converging flow fields.
    NPLANE > 0, the fixed pattern is selected for initial placement. The
    value of NPLANE serves as the number of vertical "planes" on which
    initial particles are placed within each cell block. The fixed pattern
    may work better than the random pattern only in relatively uniform flow
    fields. For two-dimensional simulations in plan view, set NPLANE = 1.
    For cross sectional or three-dimensional simulations, NPLANE = 2 is
    normally adequate. Increase NPLANE if more resolution in the vertical
    direction is desired.
    Default value: 2.
nparticles_no_advection: int, optional
    is number of initial particles per cell to be placed at cells where the
    Relative Cell Concentration Gradient is less than or equal to DCEPS.
    Generally, NPL can be set to zero since advection is considered
    insignificant when the Relative Cell Concentration Gradient is less than
    or equal to DCEPS. Setting NPL equal to NPH causes a uniform number of
    particles to be placed in every cell over the entire grid (i.e., the
    uniform approach).
    Default value: 10.
nparticles_advection: int, optional
    is number of initial particles per cell to be placed at cells where the
    Relative Cell Concentration Gradient is greater than DCEPS. The
    selection of NPH depends on the nature of the flow field and also the
    computer memory limitation. Generally, use a smaller number in
    relatively uniform flow fields and a larger number in relatively
    nonuniform flow fields. However, values exceeding 16 in twodimensional
    simulation or 32 in three-dimensional simulation are rarely necessary.
    If the random pattern is chosen, NPH particles are randomly distributed
    within the cell block. If the fixed pattern is chosen, NPH is divided by
    NPLANE to yield the number of particles to be placed per vertical plane.
    Default value: 40.
cell_min_nparticles: int, optional
    is the minimum number of particles allowed per cell. If the number of
    particles in a cell at the end of a transport step is fewer than NPMIN,
    new particles are inserted into that cell to maintain a sufficient
    number of particles. NPMIN can be set to zero in relatively uniform flow
    fields, and a number greater than zero in diverging/converging flow
    fields. Generally, a value between zero and four is adequate.
    Default value is 5.
cell_max_nparticles: int, optional
    is the maximum number of particles allowed per cell. If the number of
    particles in a cell exceeds NPMAX, all particles are removed from that
    cell and replaced by a new set of particles equal to NPH to maintain
    mass balance. Generally, NPMAX can be set to approximately twice of NPH.
    Default value: 80.

imod.wq.AdvectionMOC Class Members
==================================
   * imod.wq.AdvectionMOC.from_file

imod.wq.AdvectionMOC.from_file
==============================
Loads an imod-wq package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.wq.Package.save(),
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
package : imod.wq.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.wq.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.wq.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.wq.AdvectionModifiedMOC
============================
Solve the advention term using the Modified Method of Characteristics (MIXELM = 2)
Courant number (PERCEL) is the number of cells (or a fraction of a
cell) advection will be allowed in any direction in one transport step.

Attributes
----------
courant: float
    Courant number (PERCEL) is the number of cells (or a fraction of a cell)
    advection will be allowed in any direction in one transport step. For
    implicit finite-difference or particle tracking based schemes, there is
    no limit on PERCEL, but for accuracy reasons, it is generally not set
    much greater than one. Note, however, that the PERCEL limit is checked
    over the entire model grid. Thus, even if PERCEL > 1, advection may not
    be more than one cell’s length at most model locations. For the explicit
    finite-difference, PERCEL is also a stability constraint, which must not
    exceed one and will be automatically reset to one if a value greater
    than one is specified.
tracking: string, {"euler", "runge-kutta", "hybrid"}
    indicates which particle tracking algorithm is selected for the
    Eulerian-Lagrangian methods. ITRACK = 1, the first-order Euler algorithm is
    used; ITRACK = 2, the fourth-order Runge-Kutta algorithm is used; this
    option is computationally demanding and may be needed only when PERCEL is
    set greater than one. ITRACK = 3, the hybrid 1st and 4th order algorithm is
    used; the Runge- Kutta algorithm is used in sink/source cells and the cells
    next to sinks/sources while the Euler algorithm is used elsewhere.
weighting_factor: float
    is a concentration weighting factor (WD) between 0.5 and 1. It is used for
    operator splitting in the particle tracking based methods. The value of
    0.5 is generally adequate. The value may be adjusted to achieve better
    mass balance. Generally, it can be increased toward 1.0 as advection
    becomes more dominant.
dconcentration_epsilon: float, optional
    is a small Relative Cell Concentration Gradient (DCEPS) below which advective
    transport is considered negligible. A value around 1.0e-5 is generally
    adequate.
    Default value: 1.0e-5.
sink_particle_placement: int
    indicates whether the random or fixed pattern is selected for initial
    placement of particles to approximate sink cells in the MMOC scheme.
    (NLSINK)
sink_nparticles: int
    is the number of particles used to approximate sink cells in the MMOC
    scheme. (NPSINK)

imod.wq.AdvectionModifiedMOC Class Members
==========================================
   * imod.wq.AdvectionModifiedMOC.from_file

imod.wq.AdvectionModifiedMOC.from_file
======================================
Loads an imod-wq package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.wq.Package.save(),
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
package : imod.wq.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.wq.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.wq.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.wq.AdvectionHybridMOC
==========================
Hybrid Method of Characteristics and Modified Method of Characteristics with
MOC or MMOC automatically and dynamically selected (MIXELM = 3)

Attributes
----------
courant: float
    Courant number (PERCEL) is the number of cells (or a fraction of a cell)
    advection will be allowed in any direction in one transport step. For
    implicit finite-difference or particle tracking based schemes, there is
    no limit on PERCEL, but for accuracy reasons, it is generally not set
    much greater than one. Note, however, that the PERCEL limit is checked
    over the entire model grid. Thus, even if PERCEL > 1, advection may not
    be more than one cell’s length at most model locations. For the explicit
    finite-difference, PERCEL is also a stability constraint, which must not
    exceed one and will be automatically reset to one if a value greater
    than one is specified.
max_particles: int
    is the maximum total number of moving particles allowed (MXPART).
tracking: int
    indicates which particle tracking algorithm is selected for the
    Eulerian-Lagrangian methods. ITRACK = 1, the first-order Euler algorithm is
    used; ITRACK = 2, the fourth-order Runge-Kutta algorithm is used; this
    option is computationally demanding and may be needed only when PERCEL is
    set greater than one. ITRACK = 3, the hybrid 1st and 4th order algorithm is
    used; the Runge- Kutta algorithm is used in sink/source cells and the cells
    next to sinks/sources while the Euler algorithm is used elsewhere.
weighting_factor: float
    is a concentration weighting factor (WD) between 0.5 and 1. It is used for
    operator splitting in the particle tracking based methods. The value of
    0.5 is generally adequate. The value may be adjusted to achieve better
    mass balance. Generally, it can be increased toward 1.0 as advection
    becomes more dominant.
dceps: float
    is a small Relative Cell Concentration Gradient below which advective
    transport is considered negligible. A value around 10-5 is generally
    adequate.
nplane: int
    is a flag indicating whether the random or fixed pattern is selected for
    initial placement of moving particles. NPLANE = 0, the random pattern is
    selected for initial placement. Particles are distributed randomly in
    both the horizontal and vertical directions by calling a random number
    generator. This option is usually preferred and leads to smaller mass
    balance discrepancy in nonuniform or diverging/converging flow fields.
    NPLANE > 0, the fixed pattern is selected for initial placement. The
    value of NPLANE serves as the number of vertical "planes" on which
    initial particles are placed within each cell block. The fixed pattern
    may work better than the random pattern only in relatively uniform flow
    fields. For two-dimensional simulations in plan view, set NPLANE = 1.
    For cross sectional or three-dimensional simulations, NPLANE = 2 is
    normally adequate. Increase NPLANE if more resolution in the vertical
    direction is desired.
npl: int
    is number of initial particles per cell to be placed at cells where the
    Relative Cell Concentration Gradient is less than or equal to DCEPS.
    Generally, NPL can be set to zero since advection is considered
    insignificant when the Relative Cell Concentration Gradient is less than
    or equal to DCEPS. Setting NPL equal to NPH causes a uniform number of
    particles to be placed in every cell over the entire grid (i.e., the
    uniform approach).
nph: int
    is number of initial particles per cell to be placed at cells where the
    Relative Cell Concentration Gradient is greater than DCEPS. The
    selection of NPH depends on the nature of the flow field and also the
    computer memory limitation. Generally, use a smaller number in
    relatively uniform flow fields and a larger number in relatively
    nonuniform flow fields. However, values exceeding 16 in twodimensional
    simulation or 32 in three-dimensional simulation are rarely necessary.
    If the random pattern is chosen, NPH particles are randomly distributed
    within the cell block. If the fixed pattern is chosen, NPH is divided by
    NPLANE to yield the number of particles to be placed per vertical plane.
npmin: int
    is the minimum number of particles allowed per cell. If the number of
    particles in a cell at the end of a transport step is fewer than NPMIN,
    new particles are inserted into that cell to maintain a sufficient
    number of particles. NPMIN can be set to zero in relatively uniform flow
    fields, and a number greater than zero in diverging/converging flow
    fields. Generally, a value between zero and four is adequate.
npmax: int
    is the maximum number of particles allowed per cell. If the number of
    particles in a cell exceeds NPMAX, all particles are removed from that
    cell and replaced by a new set of particles equal to NPH to maintain
    mass balance. Generally, NPMAX can be set to approximately twice of NPH.
dchmoc: real
    is the critical Relative Concentration Gradient for controlling the
    selective use of either MOC or MMOC in the HMOC solution scheme. The MOC
    solution is selected at cells where the Relative Concentration Gradient
    is greater than DCHMOC; The MMOC solution is selected at cells where the
    Relative Concentration Gradient is less than or equal to DCHMOC

imod.wq.AdvectionHybridMOC Class Members
========================================
   * imod.wq.AdvectionHybridMOC.from_file

imod.wq.AdvectionHybridMOC.from_file
====================================
Loads an imod-wq package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.wq.Package.save(),
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
package : imod.wq.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.wq.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.wq.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.wq.AdvectionFiniteDifference
=================================
Solve the advection term using the explicit Finite Difference method
(MIXELM = 0) with upstream weighting

Attributes
----------
courant: float
    Courant number (PERCEL) is the number of cells (or a fraction of a cell)
    advection will be allowed in any direction in one transport step. For
    implicit finite-difference or particle tracking based schemes, there is
    no limit on PERCEL, but for accuracy reasons, it is generally not set
    much greater than one. Note, however, that the PERCEL limit is checked
    over the entire model grid. Thus, even if PERCEL > 1, advection may not
    be more than one cell’s length at most model locations. For the explicit
    finite-difference, PERCEL is also a stability constraint, which must not
    exceed one and will be automatically reset to one if a value greater
    than one is specified.
weighting : string {"upstream", "central"}, optional
    Indication of which weighting scheme should be used, set to default
    value "upstream" (NADVFD = 0 or 1)
    Default value: "upstream"

imod.wq.AdvectionFiniteDifference Class Members
===============================================
   * imod.wq.AdvectionFiniteDifference.from_file

imod.wq.AdvectionFiniteDifference.from_file
===========================================
Loads an imod-wq package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.wq.Package.save(),
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
package : imod.wq.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.wq.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.wq.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.wq.BasicTransport
======================
Handles basic tasks that are required by the entire transport model. Among
these tasks are definition of the problem, specification of the boundary and
initial conditions, determination of the stepsize, preparation of mass
balance information, and printout of the simulation results.

Parameters
----------
icbund: xr.DataArray of int
    is an integer array specifying the boundary condition type (inactive,
    constant-concentration, or active) for every model cell. For
    multi-species simulation, ICBUND defines the boundary condition type
    shared by all species. Note that different species are allowed to have
    different constant-concentration conditions through an option in the
    Source and Sink Mixing Package.
    ICBUND=0, the cell is an inactive concentration cell for all species.
    Note that no-flow or "dry" cells are automatically converted into
    inactive concentration cells. Furthermore, active cells in terms of flow
    can be treated as inactive concentration cells to minimize the area
    needed for transport simulation, as long as the solute transport is
    insignificant near those cells.
    ICBUND<0, the cell is a constant-concentration cell for all species. The
    starting concentration of each species remains the same at the cell
    throughout the simulation. (To define different constantconcentration
    conditions for different species at the same cell location, refer to the
    Sink/Source Mixing Package.) Also note that unless explicitly defined as
    a constant-concentration cell, a constant-head cell in the flow model is
    not treated as a constantconcentration cell.
    If ICBUND>0, the cell is an active (variable) concentration cell where
    the concentration value will be calculated.
starting_concentration: float or xr.DataArray of floats
    is the starting concentration (initial condition) at the beginning of
    the simulation (unit: ML-3) (SCONC). For multispecies simulation, the
    starting concentration must be specified for all species, one species at
    a time.
porosity: float, optional
    is the "effective" porosity of the porous medium in a single porosity
    system (PRSITY).
    Default value is 0.35.
n_species: int, optional
    is the total number of chemical species included in the current
    simulation (NCOMP). For single-species simulation, set n_species = 1.
    Default value is 1.
inactive_concentration: float, optional
    is the value for indicating an inactive concentration cell (ICBUND=0)
    (CINACT). Even if it is not anticipated to have inactive cells in the
    model, a value for inactive_concentration still must be submitted.
    Default value is 1.0e30
minimum_active_thickness: float, optional
    is the minimum saturated thickness in a cell (THKMIN), expressed as the
    decimal fraction of the model layer thickness, below which the cell is
    considered inactive.
    Default value is 0.01 (i.e., 1% of the model layer thickness).

imod.wq.BasicTransport Class Members
====================================
   * imod.wq.BasicTransport.from_file

imod.wq.BasicTransport.from_file
================================
Loads an imod-wq package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.wq.Package.save(),
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
package : imod.wq.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.wq.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.wq.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.wq.Dispersion
==================
Solves the concentration change due to dispersion explicitly or formulates
the coefficient matrix of the dispersion term for the matrix solver.

Parameters
----------
longitudinal: float
    is the longitudinal dispersivity (AL), for every cell of the model grid
    (unit: L).
    Default value is 1.0 m. Nota bene: this is for regional applications.
ratio_horizontal: float
    is a 1D real array defining the ratio of the horizontal transverse
    dispersivity (TRPT), to the longitudinal dispersivity. Each value in the
    array corresponds to one model layer. Some recent field studies suggest
    that ratio_horizontal is generally not greater than 0.1.
ratio_vertical: float
    (TRPV) is the ratio of the vertical transverse dispersivity to the
    longitudinal dispersivity. Each value in the array corresponds to one
    model layer.
    Some recent field studies suggest that ratio_vertical is generally not
    greater than 0.01.
    Set ratio_vertical equal to ratio_horizontal to use the standard
    isotropic dispersion model. Otherwise, the modified isotropic dispersion
    model is used.
diffusion_coefficient: float
    is the effective molecular diffusion coefficient (unit: L2T-1). Set
    diffusion_coefficient = 0 if the effect of molecular diffusion is
    considered unimportant. Each value in the array corresponds to one model
    layer.

    iMOD-wq always uses meters and days.

imod.wq.Dispersion Class Members
================================
   * imod.wq.Dispersion.from_file

imod.wq.Dispersion.from_file
============================
Loads an imod-wq package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.wq.Package.save(),
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
package : imod.wq.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.wq.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.wq.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.wq.MassLoading
===================
Mass loading package. Has no direct effect on groundwater flow, is only
included via MT3DMS source and sinks. (SSM ITYPE 15)

Parameters
----------
concentration: xr.DataArray of floats

imod.wq.MassLoading Class Members
=================================
   * imod.wq.MassLoading.from_file

imod.wq.MassLoading.from_file
=============================
Loads an imod-wq package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.wq.Package.save(),
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
package : imod.wq.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.wq.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.wq.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

imod.wq.TimeVaryingConstantConcentration
========================================
Time varying constant concentration package. Has no direct effect on
groundwater flow, is only included via MT3DMS source and sinks. (SSM ITYPE
-1)

Parameters
----------
concentration: xr.DataArray of floats

imod.wq.TimeVaryingConstantConcentration Class Members
======================================================
   * imod.wq.TimeVaryingConstantConcentration.from_file

imod.wq.TimeVaryingConstantConcentration.from_file
==================================================
Loads an imod-wq package from a file (currently only netcdf and zarr are supported).
Note that it is expected that this file was saved with imod.wq.Package.save(),
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
package : imod.wq.Package
    Returns a package with data loaded from file.

Examples
--------

To load a package from a file, e.g. a River package:

>>> river = imod.wq.River.from_file("river.nc")

For large datasets, you likely want to process it in chunks. You can
forward keyword arguments to ``xarray.open_dataset()`` or
``xarray.open_zarr()``:

>>> river = imod.wq.River.from_file("river.nc", chunks={"time": 1})

Refer to the xarray documentation for the possible keyword arguments.

