imod.evaluate.calculate_gxg
===========================
Calculate GxG groundwater characteristics from head time series.

GLG and GHG (average lowest and average highest groundwater level respectively) are
calculated as the average of the three lowest (GLG) or highest (GHG) head values per
Dutch hydrological year (april - april), for head values measured at a semi-monthly frequency
(14th and 28th of every month). GVG (average spring groundwater level) is calculated as
the average of groundwater level on 14th and 28th of March, and 14th of April. Supplied head
values are resampled (nearest) to the 14/28 frequency.

Hydrological years without all 24 14/28 dates present are discarded for glg and ghg.
Years without the 3 dates for gvg are discarded.

Parameters
----------
head : xr.DataArray of floats
    Head relative to sea level, in m, or m below surface level if `below_surfacelevel` is
    set to True. Must be of dimensions ``("time", "y", "x")``.
below_surfacelevel : boolean, optional, default: False.
    False (default) if heads are relative to a datum (e.g. sea level). If
    True, heads are taken as m below surface level.
tolerance: pd.Timedelta, default: 7 days.
    Maximum time window allowed when searching for dates around the 14th
    and 28th of every month.

Returns
-------
gxg : xr.Dataset
    Dataset containing ``glg``: average lowest head, ``ghg``: average
    highest head, ``gvg``: average spring head, ``n_years_gvg``: numbers of
    years used for gvg, ``n_years_gxg``: numbers of years used for glg and
    ghg.

Examples
--------
Load the heads, and calculate groundwater characteristics after the year 2000:

>>> import imod
>>> heads = imod.idf.open("head*.idf")
>>> heads = heads.sel(time=heads.time.dt.year >= 2000, layer=1)
>>> gxg = imod.evaluate.calculate_gxg(heads)

Transform to meters below surface level by substracting from surface level:

>>> surflevel = imod.idf.open("surfacelevel.idf")
>>> gxg = surflevel - gxg

Or calculate from groundwater level relative to surface level directly:

>>> gwl = surflevel - heads
>>> gxg = imod.evaluate.calculate_gxg(gwl, below_surfacelevel=True)

imod.evaluate.convert_pointwaterhead_freshwaterhead
===================================================
Function to convert point water head (as outputted by seawat)
into freshwater head, using Eq.3 from Guo, W., & Langevin, C. D. (2002):

.. math:: h_{f}=\frac{\rho}{\rho_{f}}h-\frac{\rho-\rho_{f}}{\rho_{f}}Z

An edge case arises when the head is below the cell centre, or entirely below
the cell. Strictly applying Eq.3 would result in freshwater heads that are
lower than the original point water head, which is physically impossible. This
function then outputs the freshwaterhead for the uppermost underlying cell where
the original point water head exceeds the cell centre.

Parameters
----------
pointwaterhead : float or xr.DataArray of floats
    the point water head as outputted by SEAWAT, in m.
density : float or xr.DataArray of floats
    the water density at the same locations as `pointwaterhead`.
elevation : float or xr.DataArray of floats
    elevation at the same locations as `pointwaterhead`, in m.
density_fresh : float, optional
    the density of freshwater (1000 kg/m3), or a different value if
    different units are used, or a different density reference is required.

Returns
-------
freshwaterhead : float or xr.DataArray of floats

imod.evaluate.facebudget
========================
Computes net face flow into a control volume, as defined by ``budgetzone``.

Returns a three dimensional DataArray with in- and outgoing flows for every
cell that is located on the edge of the control volume, thereby calculating
the flow through the control surface of the control volume.

Front, lower, and right arguments refer to iMOD face flow budgets, in cubic
meters per day. In terms of flow direction these are defined as:

* ``front``: positive with ``y`` (negative with row index)
* ``lower``: positive with ``layer`` (positive with layer index)
* ``right``: negative with ``x`` (negative with column index)

Only a single face flow has to be defined, for example, if only vertical
fluxes (``lower``) are to be considered.

Note that you generally should not define a budgetzone that is only one cell
wide. In that case, flow will both enter and leave the cell, resulting in a
net flow of zero (given there are no boundary conditions present).

The output of this function defines ingoing flow as positive, and outgoing
flow as negative. The output is a 3D array with net flows for every control
surface cell. You can select specific parts of the control surface, for
example only the east-facing side of the control volume. Please refer to the
examples.

Parameters
----------
budgetzone: xr.DataAray of floats
    Array defining zones. Non-zones should be with a ``np.nan`` value.
    Dimensions must be exactly ``("layer", "y", "x")``.
front: xr.DataArray of floats, optional
    Dimensions must be exactly ``("layer", "y", "x")`` or
    ``("time", "layer", "y", "x")``.
lower: xr.DataArray of floats, optional
    Dimensions must be exactly ``("layer", "y", "x")`` or
    ``("time", "layer", "y", "x")``.
right: xr.DataArray of floats, optional
    Dimensions must be exactly ``("layer", "y", "x")`` or
    ``("time", "layer", "y", "x")``.
netflow : bool, optional
    Whether to split flows by direction (front, lower, right).
    True: sum all flows. False: return individual directions.

Returns
-------
facebudget_front, facebudget_lower, face_budget_right : xr.DataArray of floats
    Only returned if `netflow` is False.
facebudget_net : xr.DataArray of floats
    Only returned if `netflow` is True.

Examples
--------
Load the face flows, and select the last time using index selection:

>>> import imod
>>> lower = imod.idf.open("bdgflf*.idf").isel(time=-1)
>>> right = imod.idf.open("bdgfrf*.idf").isel(time=-1)
>>> front = imod.idf.open("bdgfff*.idf").isel(time=-1)

Define the zone of interest, e.g. via rasterizing a shapefile:

>>> import geopandas as gpd
>>> gdf = gpd.read_file("zone_of_interest.shp")
>>> zone2D = imod.prepare.rasterize(gdf, like=lower.isel(layer=0))

Broadcast it to three dimensions:

>>> zone = xr.ones_like(flow) * zone2D

Compute net flow through the (control) surface of the budget zone:

>>> flow = imod.evaluate.facebudget(
>>>     budgetzone=zone, front=front, lower=lower, right=right
>>> )

Or evaluate only horizontal fluxes:

>>> flow = imod.evaluate.facebudget(
>>>     budgetzone=zone, front=front, right=right
>>> )

Extract the net flow, only on the right side of the zone, for example as
defined by x > 10000:

>>> netflow_right = flow.where(flow["x"] > 10_000.0).sum()

Extract the net flow for only the first four layers:

>>> netflow_layers = netflow_right.sel(layer=slice(1, 4)).sum()

Extract the net flow to the right of an arbitrary diagonal in ``x`` and
``y`` is simple using the equation for a straight line:

>>> netflow_right_of_diagonal = flow.where(
>>>    flow["y"] < (line_slope * flow["x"] + line_intercept)
>>> )

There are many ways to extract flows for a certain part of the zone of
interest. The most flexible one with regards to the ``x`` and ``y``
dimensions is by drawing a vector file, rasterizing it, and using it to
select with ``xarray.where()``.

To get the flows per direction, pass ``netflow=False``.

>>> flowfront, flowlower, flowright = imod.evaluate.facebudget(
>>>    budgetzone=zone, front=front, lower=lower, right=right, netflow=False
>>> )

imod.evaluate.flow_velocity
===========================
Compute flow velocities (m/d) from budgets (m3/d).

Parameters
----------
front: xr.DataArray of floats, optional
    Dimensions must be exactly ``("layer", "y", "x")``.
lower: xr.DataArray of floats, optional
    Dimensions must be exactly ``("layer", "y", "x")``.
right: xr.DataArray of floats, optional
    Dimensions must be exactly ``("layer", "y", "x")``.
top_bot: xr.Dataset of floats, containing 'top', 'bot' and optionally
    'dz' of layers.
    Dimensions must be exactly ``("layer", "y", "x")``.
porosity: float or xr.DataArray of floats, optional (default 0.3)
    If xr.DataArray, dimensions must be exactly ``("layer", "y", "x")``.

Returns
-------
vx, vy, vz: xr.DataArray of floats
    Velocity components in x, y, z direction.

imod.evaluate.interpolate_value_boundaries
==========================================
Function that returns all exceedance and non-exceedance boundaries for
a given threshold in a 3D values DataArray. Returned z-coordinates are
linearly interpolated between cell mids. As many boundaries are returned as are maximally
present in the 3D values DataArray. Function returns xr.DataArray of exceedance boundaries
and xr.DataArray of z-coordinates where values fall below the set treshold.

Parameters
----------
values : 3D xr.DataArray
    The datarray containing the values to search for boundaries. Dimensions ``layer``, ``y``, ``x``
z : 1D or 3D xr.DataArray
    Datarray containing z-coordinates of cell midpoints. Dimensions ``layer``, ``y``, ``x``. Should contain a dz coordinate.
threshold : float
    Value threshold

Returns
-------
xr.DataArray
    Z locations of successive exceedances of threshold from the top down. Dimensions ``boundary``, ``y``, ``x``
xr.DataArray
    Z locations of successive instances of falling below threshold from the top down. Dimensions ``boundary``, ``y``, ``x``

imod.evaluate.quiver_line
=========================
Obtain the u and v components for quiver plots for a line cross section
through a three-dimensional flux field. The u and v components are obtained
by first projecting the threedimensional flux components onto the provided
cross-section.

Parameters
----------
frf: `xarray.DataArray`
    Three- (or higher) dimensional dataarray of flow component along the rows (FLOW RIGHT FACE).
fff: `xarray.DataArray`
    Three- (or higher) dimensional dataarray of flow component along the columns (FLOW FRONT FACE).
flf: `xarray.DataArray`
    Three- (or higher) dimensional dataarray of flow component along the layers (FLOW LOWER FACE).
start: (2, ) array_like
    A latitude-longitude pair designating the start point of the cross
    section.
end: (2, ) array_like
    A latitude-longitude pair designating the end point of the cross
    section.

Returns
-------
u: `xarray.DataArray`
    The u component (x-component) of the flow projection on the cross-section between start and end coordinate,
    with new dimension "s" along the cross-section. The cellsizes along "s" are given in
    the "ds" coordinate.
v: `xarray.DataArray`
    The v component (y-component) of the flow projection on the cross-section between start and end coordinate,
    with new dimension "s" along the cross-section. The cellsizes along "s" are given in
    the "ds" coordinate.

Notes
-----
Use imod.evaluate.flow_velocity() first to obtain groundwater velocities
as input for this function. Velocity in x direction is, however, inverted and must
be re-inverted before using as input here.

imod.evaluate.quiver_linestring
===============================
Obtain the u and v components for quiver plots for a linestring cross section
through a three-dimensional flow field. The u and v components are obtained
by first projecting the threedimensional flow components onto the provided
cross-section.

Parameters
----------
frf: `xarray.DataArray`
    Three- (or higher) dimensional dataarray of flow component along the rows (FLOW RIGHT FACE).
fff: `xarray.DataArray`
    Three- (or higher) dimensional dataarray of flow component along the columns (FLOW FRONT FACE).
flf: `xarray.DataArray`
    Three- (or higher) dimensional dataarray of flow component along the layers (FLOW LOWER FACE).
linestring : shapely.geometry.LineString
    Shapely geometry designating the linestring along which to sample the
    cross section.

Returns
-------
u: `xarray.DataArray`
    The u component (x-component) of the flow projection on the cross-section between start and end coordinate,
    with new dimension "s" along the cross-section. The cellsizes along "s" are given in
    the "ds" coordinate.
v: `xarray.DataArray`
    The v component (y-component) of the flow projection on the cross-section between start and end coordinate,
    with new dimension "s" along the cross-section. The cellsizes along "s" are given in
    the "ds" coordinate.

Notes
-----
Use imod.evaluate.flow_velocity() first to obtain groundwater velocities
as input for this function. Velocity in x direction is, however, inverted and must
be re-inverted before using as input here.

imod.evaluate.streamfunction_line
=================================
Obtain the streamfunction for a line cross section through
a three-dimensional flow field. The streamfunction is obtained
by first projecting the horizontal flow components onto the provided
cross-section. The streamfunction can be contoured to visualize stream lines.
Stream lines are an efficient way to visualize groundwater flow.

Note, however, that the streamfunction is only defined in 2D, non-diverging,
steady-state flow without sources and sinks. These assumption are violated even
in a 2D model, but even more so in the approach followed here. Flow perpendicular
to the cross-section will not be visualized. It is up to the user to choose
cross-sections as perpendicular to the main flow direction as possible.

The 2D streamfunction and stream line visualization is based on work of Em. Prof. Olsthoorn.

Parameters
----------
frf: `xarray.DataArray`
    Three- (or higher) dimensional dataarray of flow component along the rows (FLOW RIGHT FACE).
fff: `xarray.DataArray`
    Three- (or higher) dimensional dataarray of flow component along the columns (FLOW FRONT FACE).
start: (2, ) array_like
    A latitude-longitude pair designating the start point of the cross
    section.
end: (2, ) array_like
    A latitude-longitude pair designating the end point of the cross
    section.

Returns
-------
`xarray.DataArray`
    The streamfunction projected on the cross-section between start and end coordinate,
    with new dimension "s" along the cross-section. The cellsizes along "s" are given in
    the "ds" coordinate.

imod.evaluate.streamfunction_linestring
=======================================
Obtain the streamfunction for a linestring cross section through
a three-dimensional flow field. The streamfunction is obtained
by first projecting the horizontal flow components onto the provided
cross-section. The streamfunction can be contoured to visualize stream lines.
Stream lines are an efficient way to visualize groundwater flow.

Note, however, that the streamfunction is only defined in 2D, non-diverging,
steady-state flow without sources and sinks. These assumption are violated even
in a 2D model, but even more so in the approach followed here. Flow perpendicular
to the cross-section will not be visualized. It is up to the user to choose
cross-sections as perpendicular to the main flow direction as possible.

The 2D streamfunction and stream line visualization is based on work of Em. Prof. Olsthoorn.

Parameters
----------
frf: `xarray.DataArray`
    Three- (or higher) dimensional dataarray of flow component along the rows (FLOW RIGHT FACE).
fff: `xarray.DataArray`
    Three- (or higher) dimensional dataarray of flow component along the columns (FLOW FRONT FACE).

linestring : shapely.geometry.LineString
    Shapely geometry designating the linestring along which to sample the
    cross section.

Returns
-------
`xarray.DataArray`
    The streamfunction projected on the cross-section defined by provided linestring,
    with new dimension "s" along the cross-section. The cellsizes along "s" are given in
    the "ds" coordinate.

imod.evaluate.intra_cell_boundary_conditions
============================================
Function to pre-check boundary-conditions against one another for large intra-cell fluxes.
ghb and river can function as source and sink, drn only as sink.

Parameters
----------
top_bot : xr.Dataset of floats
    'top_bot' should at least contain `top` and `bot` data_vars
porosity : float or xr.DataArray of floats, optional
    Effective porosity of model cells
riv : (dict or list of) imod.RiverPackage, optional
ghb : (dict or list of) imod.GeneralHeadBoundaryPackage, optional
drn : (dict or list of) imod.DrainagePackage, optional
drop_allnan : boolean, optional
    Whether source-sink combinations without overlap should be dropped from result (default True)

Returns
-------
dt_min: xr.DataArray of floats
    `dt_min` is the minimum calculated timestep over all combinations of boundary conditions
dt_all: xr.DataArray of floats
    `dt_all` is the calculated timestep for all combinations of boundary conditions

imod.evaluate.stability_constraint_advection
============================================
Computes advection stability constraint as applied in MT3D for adaptive
timestepping (Zheng & Wang, 1999 p54):

.. math:: \Delta t \leq \frac{R}{\frac{\left | v_{x} \right |}{\Delta x}+\frac{\left | v_{y} \right |}{\Delta y}+\frac{\left | v_{z} \right |}{\Delta z}}

This function can be used to select
which cells necessitate a small timestap, thereby slowing down calculations.

Front, lower, and right arguments refer to iMOD face flow budgets, in cubic
meters per day. In terms of flow direction these are defined as:

* ``front``: positive with ``y`` (negative with row index)
* ``lower``: positive with ``layer`` (positive with layer index)
* ``right``: negative with ``x`` (negative with column index)

Returns the minimum timestep that is required to satisfy this constraint.
The resulting dt xr.DataArray is the minimum timestep over all three directions,
dt_xyz is an xr.Dataset containing minimum timesteps for the three directions
separately.

Parameters
----------
front: xr.DataArray of floats, optional
    Dimensions must be exactly ``("layer", "y", "x")``.
lower: xr.DataArray of floats, optional
    Dimensions must be exactly ``("layer", "y", "x")``.
right: xr.DataArray of floats, optional
    Dimensions must be exactly ``("layer", "y", "x")``.
top_bot: xr.Dataset of floats, containing 'top', 'bot' and optionally
    'dz' of layers.
    Dimensions must be exactly ``("layer", "y", "x")``.
porosity: float or xr.DataArray of floats, optional (default 0.3)
    If xr.DataArray, dimensions must be exactly ``("layer", "y", "x")``.
R: Retardation factor, optional (default)
    Only when sorption is a factor.

Returns
-------
dt: xr.DataArray of floats
dt_xyz: xr.Dataset of floats

imod.evaluate.stability_constraint_wel
======================================
Computes sink/source stability constraint as applied in MT3D for adaptive
timestepping (Zheng & Wang, 1999 p54).

.. math:: \Delta t \leq \frac{R\theta }{\left | q_{s} \right |}

For the WEL package, a flux is known
beforehand, so we can evaluate beforehand if a flux assigned to a cell
will necessitate a small timestap, thereby slowing down calculations.

Returns a ipf DataFrame that includes a column for the specific discharge and
resulting minimum timestep.

Parameters
----------
wel: pd.DataFrame
    pd.DataFrame that defines a WEL package. Minimally includes
    x, y, layer and Q column.
top_bot: xr.Dataset of floats, containing 'top', 'bot' and optionally
    'dz' of layers.
    Dimensions must be exactly ``("layer", "y", "x")``.
porosity: float or xr.DataArray of floats, optional (default 0.3)
    If xr.DataArray, dimensions must be exactly ``("layer", "y", "x")``.
R: Retardation factor, optional (default)
    Only when sorption is a factor.

Returns
-------
wel: pd.DataFrame containing addition qs (specific discharge) and
    dt (minimum timestep) columns

