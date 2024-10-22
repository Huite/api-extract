imod.msw.GridData
=================
This contains the grid data of MetaSWAP.

This class is responsible for the file `area_svat.inp`

Parameters
----------
area: array of floats (xr.DataArray)
    Describes the area of SVAT units. This array must have a subunit coordinate
    to describe different landuses.
landuse: array of integers (xr.DataArray)
    Describes the landuse type of SVAT units.
    This array must have a subunit coordinate.
rootzone_depth: array of floats (xr.DataArray)
    Describes the rootzone depth of SVAT units.
    This array must have a subunit coordinate to describe different landuses.
surface_elevation: array of floats (xr.DataArray)
    Describes the surface elevation of SVAT units.
    This array must not have a subunit coordinate.
soil_physical_unit: array of integers (xr.DataArray)
    Describes the physical parameters of SVAT units.
    These parameters will be looked up in a table according to the given integers.
    This array must not have a subunit coordinate.
active: array of bools (xr.DataArray)
    Describes whether SVAT units are active or not.
    This array must not have a subunit coordinate.

imod.msw.GridData Class Members
===============================
   * imod.msw.GridData.generate_index_array
   * imod.msw.GridData.write
   * imod.msw.GridData.write_dataframe_fixed_width

imod.msw.GridData.generate_index_array
======================================
Generate index arrays to be used on other packages

imod.msw.GridData.write
=======================
Write MetaSWAP package to its corresponding fixed format file. This has
the `.inp` extension.

imod.msw.GridData.write_dataframe_fixed_width
=============================================
Write dataframe to fixed format file.

imod.msw.Infiltration
=====================
This contains the infiltration data.

This class is responsible for the file `infi_svat.inp`

Parameters
----------
infiltration_capacity: array of floats (xr.DataArray)
    Describes the infiltration capacity of SVAT units. This array must have
    a subunit coordinate to describe different land uses.
downward_resistance: array of floats (xr.DataArray)
    Describes the downward resisitance of SVAT units. Set to -9999.0 to make
    MetaSWAP ignore this resistance. This array must not have a subunit
    coordinate.
upward_resistance: array of floats (xr.DataArray)
    Describes the upward resistance of SVAT units. Set to -9999.0 to make
    MetaSWAP ignore this resistance. This array must not have a subunit
    coordinate.
bottom_resistance: array of floats (xr.DataArray)
    Describes the infiltration capacity of SVAT units. Set to -9999.0 to
    make MetaSWAP ignore this resistance. This array must not have a subunit
    coordinate.
extra_storage_coefficient: array of floats (xr.DataArray)
    Extra storage coefficient of phreatic layer. This array must not have a
    subunit coordinate.
active: array of bools (xr.DataArray)
    Describes whether SVAT units are active or not. This array must not have
    a subunit coordinate.

imod.msw.Infiltration Class Members
===================================
   * imod.msw.Infiltration.write
   * imod.msw.Infiltration.write_dataframe_fixed_width

imod.msw.Ponding
================
Set ponding related parameters for MetaSWAP. This class is responsible for
the svat2swnr_roff.inp file. Currently, we do not support ponds coupled to
SIMGRO's surface water module.

Parameters
----------
ponding_depth: array of floats (xr.DataArray)
    Ponding depth of the SVAT units in meters. If set to 0. water can freely
    flow over the soil surface. Runoff is disable by setting the ponding
    depth to 9999 m. Large values, e.g. 1000 m, should be avoided becauses
    this causes excess memory use. This array must have a subunit coordinate
    to describe different land uses.
runoff_resistance: array of floats (xr.DataArray)
    Runoff resistance of SVAT units in days. This array must have a subunit
    coordinate to describe different land uses.
runon_resistance: array of floats (xr.DataArray)
    Runon resistance of SVAT units in days. This array must have a subunit
    coordinate to describe different land uses.

imod.msw.Ponding Class Members
==============================
   * imod.msw.Ponding.write
   * imod.msw.Ponding.write_dataframe_fixed_width

imod.msw.ScalingFactors
=======================
This package allows you to do three things:
    1. Set scaling factors for some inputs in the soil physical database,
       namely the soil moisture content and the saturated hydraulic
       conductivity.
    2. Set a scaling factor for pressure head related parameters in the
       landuse class lookup table (LUSE_SVAT.INP).
    3. Set the depth of the perched watertable base.

This class is useful for sensitivity and uncertainty analyses, as well as
model calibration. Scaling factors are multiplied with their corresponding
parameters in the soil physical database.

Parameters
----------
scale_soil_moisture: array of floats (xr.DataArray)
    Scaling factor which adjusts the saturated soil moisture content, the
    residual soil moisture content, and the soil moisture content of
    macropores. This array must have a subunit coordinate to describe
    different landuses.
scale_hydraulic_conductivity: array of floats (xr.DataArray)
    Scaling factor which adjusts the (vertical) saturated hydraulic
    conductivity of the soil. This array must have a subunit coordinate to describe
    different landuses.
scale_pressure_head: array of floats (xr.DataArray)
    Scaling factor which adjusts the pressure head applied to the pressure
    parameters defined in LUSE_SVAT.INP. This array must have a subunit coordinate to describe
    different landuses.
depth_perched_water_table: array of floats (xr.DataArray)
    Sets the depth of the perched watertable base. If the groundwater depth
    exeeds this depth, the capillary rise is set to zero. This option has
    been included in the model on the request of a specific project (MIPWA),
    and is only sound for depths exceeding 2 meters. For more shallow
    presences of loam causing a perched watertable, it is advised to
    generate a new soil physical unit. This array must not have a subunit
    coordinate.

imod.msw.ScalingFactors Class Members
=====================================
   * imod.msw.ScalingFactors.write
   * imod.msw.ScalingFactors.write_dataframe_fixed_width

imod.msw.Sprinkling
===================
This contains the sprinkling capacities of links between SVAT units and
groundwater/surface water locations.

This class is responsible for the file `scap_svat.inp`

Parameters
----------
max_abstraction_groundwater: array of floats (xr.DataArray)
    Describes the maximum abstraction of groundwater to SVAT units in m3 per
    day. This array must not have a subunit coordinate.
max_abstraction_surfacewater: array of floats (xr.DataArray)
    Describes the maximum abstraction of surfacewater to SVAT units in m3
    per day. This array must not have a subunit coordinate.
well: WellDisStructured
    Describes the sprinkling of SVAT units coming groundwater.

imod.msw.Sprinkling Class Members
=================================
   * imod.msw.Sprinkling.write
   * imod.msw.Sprinkling.write_dataframe_fixed_width

imod.msw.IdfMapping
===================
Describes svat location in the IDF grid.

Note that MetaSWAP can only write equidistant grids.

imod.msw.IdfMapping Class Members
=================================
   * imod.msw.IdfMapping.write
   * imod.msw.IdfMapping.write_dataframe_fixed_width

imod.msw.TimeOutputControl
==========================
Specify the accumulation periods which will be used to write output. For
example, say the model computes on a daily timestep, but timesteps two days
apart are specified, the summed fluxes of each two days are written by
MetaSWAP.

Parameters
----------
time: xr.DataArray
    Timesteps at which to write output.

imod.msw.TimeOutputControl Class Members
========================================
   * imod.msw.TimeOutputControl.write
   * imod.msw.TimeOutputControl.write_dataframe_fixed_width

imod.msw.VariableOutputControl
==============================
Control which variables will be created as output. The variable names used
in this class provide a condensed water balance. You can use additional
keyword arguments to set more variables by using their specific name, e.g.
`vcr = True` for the water balance error. For all possibilities see the
SIMGRO Input and Output description.

All budgets will be written in m unit to in `.idf` files and to mm unit in
`.csv` files.

Parameters
----------
Pm: bool
    Write measured precipitation
Psgw: bool
    Write sprinkling precipitation, from groundwater
Pssw: bool
    Write sprinkling precipitation, from surface water
qrun: bool
    Write runon
qdr: bool
    Write net infiltration of surface water
qspgw: bool
    Groundwater extraction for sprinkling from layer
qmodf: bool
    Sum of all MODFLOW stresses on groundwater
ETact: bool
    Write total actual evapotranspiration, which is the sum of the
    sprinkling evaporation (Esp), interception evaporation (Eic), ponding
    evaporation (Epd) bare soil evaporation (Ebs), and actual transpiration
    (Tact).
**kwargs: bool
    Additional variables to let MetaSWAP write

imod.msw.VariableOutputControl Class Members
============================================
   * imod.msw.VariableOutputControl.write
   * imod.msw.VariableOutputControl.write_dataframe_fixed_width

imod.msw.InitialConditionsEquilibrium
=====================================
Use an equilibrium profile to initialize the model.

This class is responsible for the file `init_svat.inp`

imod.msw.InitialConditionsEquilibrium Class Members
===================================================
   * imod.msw.InitialConditionsEquilibrium.write
   * imod.msw.InitialConditionsEquilibrium.write_dataframe_fixed_width

imod.msw.InitialConditionsPercolation
=====================================
The precipitation intensity at the starting time (iybg, tdbg in
PARA_SIM.INP) is used for initializing the percolation flux in the profiles.
This type of initialization is normally done separately from the actual run,
using a specially prepared meteo-input file. After letting the model reach
near equilibrium by letting it run for a number of years, the saved state is
used for the initialization of subsequent runs.

This class is responsible for the file `init_svat.inp`

imod.msw.InitialConditionsPercolation Class Members
===================================================
   * imod.msw.InitialConditionsPercolation.write
   * imod.msw.InitialConditionsPercolation.write_dataframe_fixed_width

imod.msw.InitialConditionsRootzonePressureHead
==============================================
Use the pF-value of the root zone pressure head as initial condition.

This class is responsible for the file `init_svat.inp`

Parameters
----------
initial_pF: float
    Initial pF value to be used for all soil columns.

imod.msw.InitialConditionsRootzonePressureHead Class Members
============================================================
   * imod.msw.InitialConditionsRootzonePressureHead.write
   * imod.msw.InitialConditionsRootzonePressureHead.write_dataframe_fixed_width

imod.msw.InitialConditionsSavedState
====================================
Use saved state of a previous MetaSWAP run as initial condition.

This class is responsible for the file `init_svat.inp`

Parameters
----------
saved_state: Path or str
    Path to a previously saved state. This file will be copied to
    init_svat.inp.

imod.msw.InitialConditionsSavedState Class Members
==================================================
   * imod.msw.InitialConditionsSavedState.write
   * imod.msw.InitialConditionsSavedState.write_dataframe_fixed_width

imod.msw.InitialConditionsSavedState.write
==========================================
Write MetaSWAP package to its corresponding fixed format file. This has
the `.inp` extension.

imod.msw.LanduseOptions
=======================
Land use options. This object is responsible for luse_svat.inp

Parameters
----------
landuse_name: array of strings (xr.DataArray)
    Names of land use
vegetation_index: array of integers (xr.DataArray)
    Vegetation indices
jarvis_o2_stress: array of floats (xr.DataArray)
    Jarvis parameter for oxygen stress
jarvis_drought_stress: array of floats (xr.DataArray)
    Jarvis parameter for drought stress
feddes_p1: array of floats (xr.DataArray)
    p1 (m) in Feddes function for transpiration reduction
feddes_p2: array of floats (xr.DataArray)
    p2 (m) in Feddes function for transpiration reduction
feddes_p3h: array of floats (xr.DataArray)
    p3h (m) in Feddes function for transpiration reduction
feddes_p3l: array of floats (xr.DataArray)
    p3l (m) in Feddes function for transpiration reduction
feddes_p4: array of floats (xr.DataArray)
    p4 (m) in Feddes function for transpiration reduction
feddes_t3h: array of floats (xr.DataArray)
    t3h (mm/d) in Feddes function for transpiration reduction
feddes_t3l: array of floats (xr.DataArray)
    t3l (mm/d) in Feddes function for transpiration reduction
threshold_sprinkling: array of floats (xr.DataArray)
    If <0, pressure head (m) at which sprinkling begins. If >0 drought
    stress at which sprinkling begins.
fraction_evaporated_sprinkling: array of floats (xr.DataArray)
    Fraction evaporated sprinkling water
gift: array of floats (xr.DataArray)
    Gift (mm) during rotational period
gift_duration: array of floats (xr.DataArray)
    Gift duration (d)
rotational_period: array of floats (xr.DataArray)
    Rotational period (d)
start_sprinkling_season: array of floats (xr.DataArray)
    Day of year at which sprinkling season starts (d)
end_sprinkling_season: array of floats (xr.DataArray)
    Day of year at which sprinkling season ends (d)
interception_option: array of integers (xr.DataAray)
    Choose interception model. 0=Rutter, 1=Von Hoyningen. NOTE: option
    2=GASH, but this is not supported by MetaSWAP v8.1.0.3 and lower
interception_capacity_per_LAI: array of floats (xr.DataArray)
    Interception capacity (mm/LAI) will be set for both Rutter and Von
    Hoyningen.
interception_intercept: array of floats (xr.DataArray)
    Intercept of the interception evaporation curve. Pun unintended.

Notes
-----
No Penman-Monteith is supported in iMOD Python, so albedo, rsc, rsw, rsoil,
kdif, and kdir cannot be specified. (We might create a seperate object for
this if there is a demand for it.)

The GASH model (interception_option = 2) and salt stress parameters Maas &
Hoffman are not supported by MetaSWAP at the time of writing this class. So
these are not supported.

imod.msw.LanduseOptions Class Members
=====================================
   * imod.msw.LanduseOptions.write
   * imod.msw.LanduseOptions.write_dataframe_fixed_width

imod.msw.AnnualCropFactors
==========================
For each vegetation type specify a yearly trend in vegetation factors and
interception characteristics. These are used if WOFOST is not used.

This class is responsible for the file `fact_svat.inp`.

Parameters
----------
soil_cover: array of floats (xr.DataArray)
    Soil cover in m2/m2. Must have a "vegetation_index" and "day_of_year" a
    coordinates.
leaf_area_index: array of floats (xr.DataArray)
    Leaf area index in m2/m2. Must have a "vegetation_index" and
    "day_of_year" a coordinates.
interception_capacity: array of floats (xr.DataArray)
    Interception capacity in m3/m2. Must have a "vegetation_index" and
    "day_of_year" a coordinates.
vegetation_factor: array of floats (xr.DataArray)
    Vegetation factor. Must have a "vegetation_index" and "day_of_year" a
    coordinates.
interception_factor: array of floats (xr.DataArray)
    Interception evaporation factor. Must have a "vegetation_index" and
    "day_of_year" a coordinates.
bare_soil_factor: array of floats (xr.DataArray)
    Bare soil evaporation factor. Must have a "vegetation_index" and
    "day_of_year" a coordinates.
ponding_factor: array of floats (xr.DataArray)
    Ponding factor. Must have a "vegetation_index" and "day_of_year" a
    coordinates.

imod.msw.AnnualCropFactors Class Members
========================================
   * imod.msw.AnnualCropFactors.write
   * imod.msw.AnnualCropFactors.write_dataframe_fixed_width

imod.msw.MeteoGrid
==================
This contains the meteorological grid data. Grids are written to ESRI ASCII
files. The meteorological data requires a time coordinate. Next to a
MeteoGrid instance, instances of PrecipitationMapping and
EvapotranspirationMapping are required as well to specify meteorological
information to MetaSWAP.

This class is responsible for `mete_grid.inp`.

Parameters
----------
precipitation: array of floats (xr.DataArray)
    Contains the precipitation grids in mm/d. A time coordinate is required.
evapotranspiration: array of floats (xr.DataArray)
    Contains the evapotranspiration grids in mm/d. A time coordinate is
    required.

imod.msw.MeteoGrid Class Members
================================
   * imod.msw.MeteoGrid.check_string_lengths
   * imod.msw.MeteoGrid.write
   * imod.msw.MeteoGrid.write_dataframe_fixed_width
   * imod.msw.MeteoGrid.write_free_format_file

imod.msw.MeteoGrid.check_string_lengths
=======================================
Check if strings lengths do not exceed 256 characters.
With absolute paths this might be an issue.

imod.msw.MeteoGrid.write
========================
Write mete_grid.inp and accompanying ASCII grid files.

Parameters
----------
directory: str or Path
    directory to write file in.

imod.msw.MeteoGrid.write_free_format_file
=========================================
Write free format file. The mete_grid.inp file is free format.

imod.msw.EvapotranspirationMapping
==================================
This contains the data to connect evapotranspiration grid cells to MetaSWAP
svats. The evapotranspiration grid does not have to be equal to the metaswap
grid: connections between the evapotranspiration cells to svats will be
established using a nearest neighbour lookup.

This class is responsible for the file `svat2etrefgrid.inp`.

Parameters
----------
evapotransporation: array of floats (xr.DataArray)
    Describes the evapotransporation data. The extend of the grid must be
    larger than the MetaSvap grid. The data must also be coarser than the
    MetaSvap grid.

imod.msw.EvapotranspirationMapping Class Members
================================================
   * imod.msw.EvapotranspirationMapping.write
   * imod.msw.EvapotranspirationMapping.write_dataframe_fixed_width

imod.msw.PrecipitationMapping
=============================
This contains the data to connect precipitation grid cells to MetaSWAP
svats. The precipitation grid does not have to be equal to the metaswap
grid: connections between the precipitation cells to svats will be
established using a nearest neighbour lookup.

This class is responsible for the file `svat2precgrid.inp`.

Parameters
----------
precipitation: array of floats (xr.DataArray)
    Describes the precipitation data. The extend of the grid must be larger
    than the MetaSvap grid. The data must also be coarser than the MetaSvap
    grid.

imod.msw.PrecipitationMapping Class Members
===========================================
   * imod.msw.PrecipitationMapping.write
   * imod.msw.PrecipitationMapping.write_dataframe_fixed_width

imod.msw.CouplerMapping
=======================
This contains the data to connect MODFLOW 6 cells to MetaSWAP svats.

This class is responsible for the file `mod2svat.inp`. It also includes
connection to wells.

Parameters
----------
modflow_dis: StructuredDiscretization
    Modflow 6 structured discretization
well: WellDisStructured (optional)
    If given, this parameter describes sprinkling of SVAT units from MODFLOW
    cells.

imod.msw.CouplerMapping Class Members
=====================================
   * imod.msw.CouplerMapping.write
   * imod.msw.CouplerMapping.write_dataframe_fixed_width

imod.msw.MetaSwapModel
======================
Contains data and writes consistent model input files

Parameters
----------
unsaturated_database: Path-like or str
    Path to the MetaSWAP soil physical database folder.

imod.msw.MetaSwapModel Class Members
====================================
   * imod.msw.MetaSwapModel.clear
   * imod.msw.MetaSwapModel.get
   * imod.msw.MetaSwapModel.items
   * imod.msw.MetaSwapModel.keys
   * imod.msw.MetaSwapModel.pop
   * imod.msw.MetaSwapModel.popitem
   * imod.msw.MetaSwapModel.setdefault
   * imod.msw.MetaSwapModel.update
   * imod.msw.MetaSwapModel.values
   * imod.msw.MetaSwapModel.write

imod.msw.MetaSwapModel.clear
============================
D.clear() -> None.  Remove all items from D.

imod.msw.MetaSwapModel.get
==========================
D.get(k[,d]) -> D[k] if k in D, else d.  d defaults to None.

imod.msw.MetaSwapModel.items
============================
D.items() -> a set-like object providing a view on D's items

imod.msw.MetaSwapModel.keys
===========================
D.keys() -> a set-like object providing a view on D's keys

imod.msw.MetaSwapModel.pop
==========================
D.pop(k[,d]) -> v, remove specified key and return the corresponding value.
If key is not found, d is returned if given, otherwise KeyError is raised.

imod.msw.MetaSwapModel.popitem
==============================
D.popitem() -> (k, v), remove and return some (key, value) pair
as a 2-tuple; but raise KeyError if D is empty.

imod.msw.MetaSwapModel.setdefault
=================================
D.setdefault(k[,d]) -> D.get(k,d), also set D[k]=d if k not in D

imod.msw.MetaSwapModel.update
=============================
D.update([E, ]**F) -> None.  Update D from mapping/iterable E and F.
If E present and has a .keys() method, does:     for k in E: D[k] = E[k]
If E present and lacks .keys() method, does:     for (k, v) in E: D[k] = v
In either case, this is followed by: for k, v in F.items(): D[k] = v

imod.msw.MetaSwapModel.values
=============================
D.values() -> an object providing a view on D's values

imod.msw.MetaSwapModel.write
============================
Write packages and simulation settings (para_sim.inp).

Parameters
----------
directory: Path or str
    directory to write model in.

