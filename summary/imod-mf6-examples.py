"""
Circle
======

This example illustrates how to setup a very simple unstructured groundwater
model using the ``imod`` package and associated packages.

In overview, we'll set the following steps:

    * Create a triangular mesh for a disk geometry.
    * Create the xugrid UgridDataArrays containg the MODFLOW6 parameters.
    * Feed these arrays into the imod mf6 classes.
    * Write to modflow6 files.
    * Run the model.
    * Open the results back into UgridDataArrays.
    * Visualize the results.
"""

# sphinx_gallery_thumbnail_number = -1

# %%
# We'll start with the following imports:

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import xugrid as xu
from pandas import isnull

import imod

# %%
# Create a mesh
# -------------
#
# The first steps consists of generating a mesh. In this example, we'll use data
# included with iMOD Python for a circular mesh. Note that this is a `Ugrid2D
# object. <https://deltares.github.io/xugrid/api/xugrid.Ugrid2d.html>`_
# For more information on working with unstructured grids see the
# `Xugrid documentation <https://deltares.github.io/xugrid/index.html>`_

grid = imod.data.circle()

grid


# %%
# We can plot this object as follows:

fig, ax = plt.subplots()
xu.plot.line(grid, ax=ax)
ax.set_aspect(1)

# %%
# Create UgridDataArray
# ---------------------
#
# Now that we have defined the grid, we can start defining the model parameter
# data.
#
# Our goal here is to define a steady-state model with:
#
# * Uniform conductivities of 1.0 m/d;
# * Two layers of 5.0 m thick;
# * Uniform recharge of 0.001 m/d on the top layer;
# * Constant heads of 1.0 m along the exterior edges of the mesh.
#
# From these boundary conditions, we would expect circular mounding of the
# groundwater; with small flows in the center and larger flows as the recharge
# accumulates while the groundwater flows towards the exterior boundary.

nface = grid.n_face
nlayer = 2

idomain = xu.UgridDataArray(
    xr.DataArray(
        np.ones((nlayer, nface), dtype=np.int32),
        coords={"layer": [1, 2]},
        dims=["layer", grid.face_dimension],
    ),
    grid=grid,
)
icelltype = xu.full_like(idomain, 0)
k = xu.full_like(idomain, 1.0, dtype=float)
k33 = k.copy()
rch_rate = xu.full_like(idomain.sel(layer=1), 0.001, dtype=float)
bottom = idomain * xr.DataArray([5.0, 0.0], dims=["layer"])

# %%
# All the data above have been constants over the grid. For the constant head
# boundary, we'd like to only set values on the external border. We can
# `py:method:xugrid.UgridDataset.binary_dilation` to easily find these cells:

chd_location = xu.zeros_like(idomain.sel(layer=2), dtype=bool).ugrid.binary_dilation(
    border_value=True
)
constant_head = xu.full_like(idomain.sel(layer=2), 1.0, dtype=float).where(chd_location)

fig, ax = plt.subplots()
constant_head.ugrid.plot(ax=ax)
xu.plot.line(grid, ax=ax, color="black")
ax.set_aspect(1)

# %%
# Write the model
# ---------------
#
# The first step is to define an empty model, the parameters and boundary
# conditions are added in the form of the familiar MODFLOW packages.

gwf_model = imod.mf6.GroundwaterFlowModel()
gwf_model["disv"] = imod.mf6.VerticesDiscretization(
    top=10.0, bottom=bottom, idomain=idomain
)
gwf_model["chd"] = imod.mf6.ConstantHead(
    constant_head, print_input=True, print_flows=True, save_flows=True
)
gwf_model["ic"] = imod.mf6.InitialConditions(start=0.0)
gwf_model["npf"] = imod.mf6.NodePropertyFlow(
    icelltype=icelltype,
    k=k,
    k33=k33,
    save_flows=True,
)
gwf_model["sto"] = imod.mf6.SpecificStorage(
    specific_storage=1.0e-5,
    specific_yield=0.15,
    transient=False,
    convertible=0,
)
gwf_model["oc"] = imod.mf6.OutputControl(save_head="all", save_budget="all")
gwf_model["rch"] = imod.mf6.Recharge(rch_rate)

simulation = imod.mf6.Modflow6Simulation("circle")
simulation["GWF_1"] = gwf_model
simulation["solver"] = imod.mf6.Solution(
    modelnames=["GWF_1"],
    print_option="summary",
    outer_dvclose=1.0e-4,
    outer_maximum=500,
    under_relaxation=None,
    inner_dvclose=1.0e-4,
    inner_rclose=0.001,
    inner_maximum=100,
    linear_acceleration="cg",
    scaling_method=None,
    reordering_method=None,
    relaxation_factor=0.97,
)
simulation.create_time_discretization(additional_times=["2000-01-01", "2000-01-02"])

# %%
# We'll create a new directory in which we will write and run the model.

modeldir = imod.util.temporary_directory()
simulation.write(modeldir)

# %%
# Run the model
# -------------
#
# .. note::
#
#   The following lines assume the ``mf6`` executable is available on your PATH.
#   :ref:`The Modflow 6 examples introduction <mf6-introduction>` shortly
#   describes how to add it to yours.

simulation.run()

# %%
# Open the results
# ----------------
#
# First, we'll open the heads (.hds) file.

head = simulation.open_head()

head

# %%
# For a DISV MODFLOW6 model, the heads are returned as a UgridDataArray.  While
# all layers are timesteps are available, they are only loaded into memory as
# needed.
#
# We may also open the cell-by-cell flows (.cbc) file.

cbc = simulation.open_flow_budget()

print(cbc.keys())

# %%
# The flows are returned as a dictionary of UgridDataArrays. This dictionary
# contains all entries that are stored in the CBC file, but like for the heads
# file the data are only loaded into memory when needed.
#
# The horizontal flows are stored on the edges of the UgridDataArray topology.
# The other flows are generally stored on the faces; this includes the
# flow-lower-face.
#
# We'll create a dataset for the horizontal flows for further analysis.

cbc_grid = cbc["flow-horizontal-face-x"].grid
ds = xu.UgridDataset(grids=cbc_grid)
ds["u"] = cbc["flow-horizontal-face-x"]
ds["v"] = cbc["flow-horizontal-face-y"]

# %%
# Visualize the results
# ---------------------
#
# We can quickly and easily visualize the output with the plotting functions
# provided by xarray and xugrid. We'll add some some edge coordinates to the
# dataset so that they can be used to place the arrows in the quiver plot.

ds = ds.ugrid.assign_edge_coords()
fig, ax = plt.subplots()
head.isel(time=0, layer=0).compute().ugrid.plot(ax=ax)
ds.isel(time=0, layer=0).plot.quiver(
    x="mesh2d_edge_x", y="mesh2d_edge_y", u="u", v="v", color="white"
)
ax.set_aspect(1)

# %%
# As would be expected from our model input, we observe circular groundwater
# mounding and increasing flows as we move from the center to the exterior.

# %%
# Slice the model domain
# ----------------------
#
# We may also quickly setup a smaller model. We'll select half of the original
# domain. To set up the boundary conditions on the clipped edges you can provide
# a states_for_boundary dictionary. In this case we add the head values of the
# computed full domain simulation as the clipped boundary values

states_for_boundary = {
    "GWF_1": head.compute(),
}

half_simulation = simulation.clip_box(
    x_max=0.0, states_for_boundary=states_for_boundary
)

# %%
# Let's run the model, read the results, and visualize.

modeldir = imod.util.temporary_directory()
half_simulation.write(modeldir)
half_simulation.run()
head = half_simulation.open_head()

# %%
# Let's add constant head boundaries together and plot them

half_simulation_constant_head = half_simulation["GWF_1"]["chd"]["head"]

clipped_half_simulation_constant_head = (
    half_simulation["GWF_1"]["chd_clipped"]["head"].sel(layer=2).isel(time=0)
)

all_boundaries_constant_head = half_simulation_constant_head.where(
    ~isnull(half_simulation_constant_head), clipped_half_simulation_constant_head
)

# plot boundary conditions
fig, ax = plt.subplots()
all_boundaries_constant_head.ugrid.plot(ax=ax)
ax.set_aspect(1)

# %%
# plot computed heads
fig, ax = plt.subplots()
head.isel(time=0, layer=0).compute().ugrid.plot(ax=ax)
ax.set_aspect(1)
# %%


"""
Circle partitioned
==================

This example illustrates a circular model that is split into 3 submodels.
The split method returns a simulation object that can be run as is. In this
case the 3 submodels are roughly equal sized partitions that have the shape
of pie pieces.
"""

import matplotlib.pyplot as plt
from example_models import create_circle_simulation

import imod
from imod.mf6.multimodel.partition_generator import get_label_array

simulation = create_circle_simulation()
tmp_path = imod.util.temporary_directory()
simulation.write(tmp_path / "original", False)

idomain = simulation["GWF_1"]["disv"].dataset["idomain"]

number_partitions = 5
submodel_labels = get_label_array(simulation, number_partitions)

# Create a simulation that is split in subdomains according to the label array.
new_sim = simulation.split(submodel_labels)
# %%
# Write the simulation input files for the new simulation.
new_sim.write(tmp_path, False)

# run the split simulation
new_sim.run()
# %%
# Visualize the computed heads in the top layer.
fig, ax = plt.subplots()
head = new_sim.open_head()

head["head"].isel(layer=0, time=-1).ugrid.plot.contourf(ax=ax)
# %%
# Visualize the flow-horizontal-face-x componenty of the balances.
fig, ax = plt.subplots()
balances = new_sim.open_flow_budget()

balances["flow-horizontal-face-x"].isel(layer=0, time=-1).ugrid.plot()
pass

# %%


"""
Freshwater lens (circle)
========================

This example illustrates how to setup a very simple unstructured groundwater
transport model using the ``imod`` package and associated packages.

In overview, we'll set the following steps:

    * Setting up the flow model, just like in the circle.py example
    * set up the transport model
    * Run the simulation.
    * Visualize the results.
"""

# sphinx_gallery_thumbnail_number = -1

# %%
# We'll start with the following imports:

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
import xugrid as xu

import imod

# %%
# Parameters
# ----------

porosity = 0.10
max_concentration = 35.0
min_concentration = 0.0
max_density = 1025.0
min_density = 1000.0
k_value = 10.0

# %%
# Create a mesh
# -------------
#
# The first steps consists of generating a mesh. In this example, we'll use data
# included with iMOD Python for a circular mesh. Note that this is a `Ugrid2D
# object. <https://deltares.github.io/xugrid/api/xugrid.Ugrid2d.html>`_
# For more information on working with unstructured grids see the
# `Xugrid documentation <https://deltares.github.io/xugrid/index.html>`_

grid_triangles = imod.data.circle()

fig, ax = plt.subplots()
xu.plot.line(grid_triangles)
ax.set_aspect(1)

# %%
# However a triangular grid has the issue that the direction of the fluxes
# between cell centres is not perpendicular to the cell vertices. The default
# formulation of Modflow 6 does not account for this, which causes mass balance
# errors. The XT3D formulation is able to account for this, but the last version
# of Modflow 6 (6.3 at time of writing) does not support this in combination
# with the Buoyancy package and using XT3D comes with an extra, significant
# computational burden. It is therefore easier to use a voronoi grid, for which
# `Xugrid <https://deltares.github.io/xugrid/index.html>`_ has a very convenient
# method.

grid = grid_triangles.tesselate_centroidal_voronoi()

fig, ax = plt.subplots()
xu.plot.line(grid)
ax.set_aspect(1)

# %%
# Create arrays

nface = grid.n_face
nlayer = 15

layer = np.arange(nlayer, dtype=int) + 1

idomain = xu.UgridDataArray(
    xr.DataArray(
        np.ones((nlayer, nface), dtype=np.int32),
        coords={"layer": layer},
        dims=["layer", grid.face_dimension],
    ),
    grid=grid,
)
icelltype = xu.full_like(idomain, 0)
k = xu.full_like(idomain, k_value, dtype=float)
k33 = k.copy()

top = 0.0
bottom = xr.DataArray(top - (layer * 10.0), dims=["layer"])

# %%
# Recharge
# --------
#
# We need a recharge rate for the fluid and a recharge rate for the solute. The
# fluid recharge rate is volumetric and per unit area, so the unit is
# length/time. The solute recharge rate is the concentration of solute in the
# recharge, and has concentration units.

rch_rate = xu.full_like(idomain.sel(layer=1), 0.001, dtype=float)
rch_concentration = xu.full_like(rch_rate, min_concentration)
rch_concentration = rch_concentration.expand_dims(species=["salinity"])


# %%
# Unlike a recharge boundary, with a prescribed head boundary we don't know a
# priori whether water will flow in over the boundary or leave across the
# boundary. If water flows into the model over the boundary, it carries a
# prescribed solute concentration. If it leaves, it leaves with the
# concentration that was computed for the cell.
#
# In this example we set the prescribed head value to 0.0 and the external
# concentration to 35.0 as well. The boundary only operates on the top layer.

chd_location = xu.zeros_like(idomain.sel(layer=1), dtype=bool).ugrid.binary_dilation(
    border_value=True
)
constant_head = xu.full_like(idomain, 0.0, dtype=float).where(chd_location)
# Approximate face area
face_area = (1000.0 / 6) ** 2 * 0.5

conductance = xu.full_like(idomain, face_area * k_value, dtype=float).where(
    chd_location
)

constant_concentration = xu.full_like(constant_head, max_concentration).where(
    chd_location
)
constant_concentration = constant_concentration.expand_dims(species=["salinity"])


# %%
# Add flow model to simulation
# ----------------------------
#
# See the circle.py example for more information.

gwf_model = imod.mf6.GroundwaterFlowModel()
gwf_model["disv"] = imod.mf6.VerticesDiscretization(
    top=top, bottom=bottom, idomain=idomain
)
gwf_model["ghb"] = imod.mf6.GeneralHeadBoundary(
    constant_head,
    conductance=conductance,
    concentration=constant_concentration,
    print_input=True,
    print_flows=True,
    save_flows=True,
)
gwf_model["ic"] = imod.mf6.InitialConditions(start=0.0)
gwf_model["npf"] = imod.mf6.NodePropertyFlow(
    icelltype=icelltype,
    k=k,
    k33=k33,
    save_flows=True,
)
gwf_model["sto"] = imod.mf6.SpecificStorage(
    specific_storage=1.0e-5,
    specific_yield=0.15,
    transient=False,
    convertible=0,
)
gwf_model["oc"] = imod.mf6.OutputControl(save_head="last", save_budget="last")
gwf_model["rch"] = imod.mf6.Recharge(
    rch_rate, concentration=rch_concentration, print_flows=True, save_flows=True
)

simulation = imod.mf6.Modflow6Simulation("circle")
simulation["flow"] = gwf_model
simulation["flow_solver"] = imod.mf6.Solution(
    modelnames=["flow"],
    print_option="summary",
    outer_dvclose=1.0e-4,
    outer_maximum=500,
    under_relaxation=None,
    inner_dvclose=1.0e-4,
    inner_rclose=0.001,
    inner_maximum=100,
    linear_acceleration="bicgstab",
    scaling_method=None,
    reordering_method=None,
    relaxation_factor=0.97,
)

# %%
# Set the timesteps, we want output each year, so we specify stress periods
# which last 1 year. However, timesteps of 1 year yield unstable results, so we
# set ``n_timesteps`` to 10, which sets the amount of timesteps within a stress
# period.

simtimes = pd.date_range(start="2000-01-01", end="2030-01-01", freq="As")
simulation.create_time_discretization(additional_times=simtimes)
simulation["time_discretization"]["n_timesteps"] = 10

# %%
# Buoyancy
# --------
# Since we are solving a variable density problem, we need to add the buoyancy
# package. It will use the species "salinity" that we are simulating with a
# transport model defined below.

slope = (max_density - min_density) / (max_concentration - min_concentration)
gwf_model["buoyancy"] = imod.mf6.Buoyancy(
    reference_density=min_density,
    modelname=["transport"],
    reference_concentration=[min_concentration],
    density_concentration_slope=[slope],
    species=["salinity"],
)

# %%
# Add transport model to simulation
# ---------------------------------
#
# The transport model requires the flow field inside the domain computed by the
# flow model. It also needs the fluxes over the boundary. It uses the same
# discretization as the flow model. Here we create a transport model for
# salinity, derive sources and sinks based from the flow model, and tell it to
# use the same discretization.

transport_model = imod.mf6.GroundwaterTransportModel()
transport_model["ssm"] = imod.mf6.SourceSinkMixing.from_flow_model(
    gwf_model, "salinity"
)
transport_model["disv"] = gwf_model["disv"]

# %%
# Now we define some transport packages for simulating the physical processes
# of advection, mechanical dispersion, and molecular diffusion dispersion. This
# example is transient, and the volume available for storage is the porosity,
# in this case 0.10.

al = 0.001

transport_model["dsp"] = imod.mf6.Dispersion(
    diffusion_coefficient=1e-4,
    longitudinal_horizontal=al,
    transversal_horizontal1=al * 0.1,
    transversal_vertical=al * 0.01,
    xt3d_off=False,
    xt3d_rhs=False,
)
transport_model["adv"] = imod.mf6.AdvectionUpstream()
transport_model["mst"] = imod.mf6.MobileStorageTransfer(porosity)

# %%
# Define the maximum concentration as the initial conditions, also output
# options for the transport model, and assign the transport model to the
# simulation as well.

transport_model["ic"] = imod.mf6.InitialConditions(start=max_concentration)
transport_model["oc"] = imod.mf6.OutputControl(
    save_concentration="last", save_budget="last"
)

simulation["transport"] = transport_model
simulation["transport_solver"] = imod.mf6.Solution(
    modelnames=["transport"],
    print_option="summary",
    outer_dvclose=1.0e-4,
    outer_maximum=500,
    under_relaxation=None,
    inner_dvclose=1.0e-4,
    inner_rclose=0.001,
    inner_maximum=100,
    linear_acceleration="bicgstab",
    scaling_method=None,
    reordering_method=None,
    relaxation_factor=0.97,
)

# %%
# We'll create a new directory in which we will write and run the model.

modeldir = imod.util.temporary_directory()
simulation.write(modeldir, binary=False)

# %%
# Run the model
# -------------
#
# .. note::
#
#   The following lines assume the ``mf6`` executable is available on your
#   PATH. The examples introduction shortly describes how to add it to yours.

simulation.run()

# %%
# Open the results
# ----------------

sim_concentration = simulation.open_concentration().compute()
sim_head = simulation.open_head().compute()

# %%
# Assign coordinates to output
# ----------------------------
#
# The model output does not feature very useful coordinate values for ``time``
# and ``z``, therefore it is best to assign these to the datasets for more
# understandable plots.
#
# First we have to compute values for a z coordinate. The
interfaces = np.concatenate([[top], bottom.values])
z = (interfaces[:-1] + interfaces[1:]) / 2

z

# %%
# Assign these new coordinate values to the dataset
coords = {"time": simtimes[1:], "z": ("layer", z)}

sim_head = sim_head.assign_coords(**coords)
sim_concentration = sim_concentration.assign_coords(**coords)


# %%
# Visualize the results
# ---------------------
#
# We can quickly and easily visualize the output with the plotting functions
# provided by xarray and xugrid:

fig, ax = plt.subplots()
sim_head.isel(time=-1, layer=0).ugrid.plot(ax=ax)
ax.set_aspect(1)

# %%
# We can draw a crossection through the center by selecting y=0, for which we
# can plot the contours as follows:

fig, ax = plt.subplots()
sim_concentration.isel(time=-1).ugrid.sel(y=0).plot.contourf(
    ax=ax, x="mesh2d_x", y="z", cmap="RdYlBu_r"
)
# %%


"""
Different ways to regrid models
================================

This example focuses on regridding. It uses the TWRI model from modflow6 (`Harbaugh, 2005`).
More information about this model can be found in an example dedicated to building this model ( ex01_twri.py)

In overview, we'll set the following steps:
* we build a new grid, onto which we want to regrid the twri model
* we regrid the model in 3 different ways
** first by regridding the simulation itself. This automatically regrids the model and all the packages in the model using
default regridding methods for each field in each package
** second, we show how a model within a simulation can be regridded, also using default regridding methods
** third, we show how regridding can be done package per package. We illustrate that for one package, and we use a non-default
regridding method for the horizontal conductivity field

"""

# %%
# We'll start with the usual imports. As this is an simple (synthetic)
# structured model, we can make due with few packages.

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from example_models import create_twri_simulation

from imod.mf6.regrid import NodePropertyFlowRegridMethod
from imod.mf6.utilities.regrid import RegridderType, RegridderWeightsCache

# %%
# Now we create the twri simulation itself. It yields a simulation of a flow problem, with a grid of 3 layers and 15 cells in both x and y directions.
# To better illustrate the regridding, we replace the K field with a lognormal random K field. The original k-field is a constant per layer.
simulation = create_twri_simulation()

idomain = simulation["GWF_1"]["dis"]["idomain"]
heterogeneous_k = xr.zeros_like(idomain, dtype=np.double)
heterogeneous_k.values = np.random.lognormal(4.0, 1.0, heterogeneous_k.shape)
simulation["GWF_1"]["npf"]["k"] = heterogeneous_k

# %%
# Let's plot the k-field. This is going to be the input for the regridder, and the regridded output should somewhat resemble it.
fig, ax = plt.subplots()
heterogeneous_k.sel(layer=1).plot(y="y", yincrease=False, ax=ax)

# %%
# Now we create a new grid for this simulation. It has 3 layers,  45 rows and 20 columns.
# The length of the domain is slightly different from the input grid. That had a coordinate difference between the first and last cellcentre on the
# x axis and y axis of  15*5000 = 75000 on both axes, but the new grid that is  75020 on the x axis and  75015 on the y axis

nlay = 3
nrow = 21
ncol = 12
shape = (nlay, nrow, ncol)

dx = 6251
dy = -3572
xmin = 0.0
xmax = dx * ncol
ymin = 0.0
ymax = abs(dy) * nrow
dims = ("layer", "y", "x")

layer = np.array([1, 2, 3])
y = np.arange(ymax, ymin, dy) + 0.5 * dy
x = np.arange(xmin, xmax, dx) + 0.5 * dx
coords = {"layer": layer, "y": y, "x": x, "dx": dx, "dy": dy}
target_grid = xr.DataArray(np.ones(shape, dtype=int), coords=coords, dims=dims)

# %%
# A first way to regrid the twri model is to regrid the whole simulation object. This is the most straightforward method,
# and it uses default regridding methods for each input field. To see which ones are used, look at the _regrid_method
# class attribute of the relevant package. For example the _regrid_method attribute  of the NodePropertyFlow package
# specifies that field "k" uses an OVERLAP regridder in combination with the averaging function "geometric_mean".
new_simulation = simulation.regrid_like("regridded_twri", target_grid=target_grid)

# %%
# Let's plot the k-field. This is the regridded output, and it should should somewhat resemble the original k-field plotted earlier.
regridded_k_1 = new_simulation["GWF_1"]["npf"]["k"]
fig, ax = plt.subplots()
regridded_k_1.sel(layer=1).plot(y="y", yincrease=False, ax=ax)
# %%
# A second way to regrid  twri  is to regrid the groundwater flow model.

model = simulation["GWF_1"]
new_model = model.regrid_like(target_grid)

regridded_k_2 = new_model["npf"]["k"]
fig, ax = plt.subplots()
regridded_k_2.sel(layer=1).plot(y="y", yincrease=False, ax=ax)


# %% Finally, we can regrid package per package. This allows us to choose the
# regridding method as well. in this example we'll regrid the npf package
# manually and the rest of the packages using default methods.
#
# Note that we create a RegridderWeightsCache here. This will store the weights
# of the regridder. Using the same cache to regrid another package will lead to
# a performance increase if that package uses the same regridding method,
# because initializing a regridder is costly.

regridder_types = NodePropertyFlowRegridMethod(k=(RegridderType.CENTROIDLOCATOR,))
regrid_context = RegridderWeightsCache()
npf_regridded = model["npf"].regrid_like(
    target_grid=target_grid,
    regrid_context=regrid_context,
    regridder_types=regridder_types,
)
new_model["npf"] = npf_regridded


regridded_k_3 = new_model["npf"]["k"]
fig, ax = plt.subplots()
regridded_k_3.sel(layer=1).plot(y="y", yincrease=False, ax=ax)
# %%


"""
TWRI
====

This example has been converted from the `MODFLOW6 Example problems`_.  See the
`description`_ and the `notebook`_ which uses `FloPy`_ to setup the model.

This example is a modified version of the original MODFLOW example
("`Techniques of Water-Resources Investigation`_" (TWRI)) described in
(`McDonald & Harbaugh, 1988`_) and duplicated in (`Harbaugh & McDonald, 1996`_).
This problem is also is distributed with MODFLOW-2005 (`Harbaugh, 2005`_). The
problem has been modified from a quasi-3D problem, where confining beds are not
explicitly simulated, to an equivalent three-dimensional problem.

In overview, we'll set the following steps:

    * Create a structured grid for a rectangular geometry.
    * Create the xarray DataArrays containg the MODFLOW6 parameters.
    * Feed these arrays into the imod mf6 classes.
    * Write to modflow6 files.
    * Run the model.
    * Open the results back into DataArrays.
    * Visualize the results.

"""

# %%
# We'll start with the usual imports. As this is an simple (synthetic)
# structured model, we can make due with few packages.

import numpy as np
import xarray as xr

import imod

# %%
# Create grid coordinates
# -----------------------
#
# The first steps consist of setting up the grid -- first the number of layer,
# rows, and columns. Cell sizes are constant throughout the model.

nlay = 3
nrow = 15
ncol = 15
shape = (nlay, nrow, ncol)

dx = 5000.0
dy = -5000.0
xmin = 0.0
xmax = dx * ncol
ymin = 0.0
ymax = abs(dy) * nrow
dims = ("layer", "y", "x")

layer = np.array([1, 2, 3])
y = np.arange(ymax, ymin, dy) + 0.5 * dy
x = np.arange(xmin, xmax, dx) + 0.5 * dx
coords = {"layer": layer, "y": y, "x": x}

# %%
# Create DataArrays
# -----------------
#
# Now that we have the grid coordinates setup, we can start defining model
# parameters. The model is characterized by:
#
# * a constant head boundary on the left
# * a single drain in the center left of the model
# * uniform recharge on the top layer
# * a number of wells scattered throughout the model.

idomain = xr.DataArray(np.ones(shape, dtype=int), coords=coords, dims=dims)
bottom = xr.DataArray([-200.0, -300.0, -450.0], {"layer": layer}, ("layer",))

# Constant head
constant_head = xr.full_like(idomain, np.nan, dtype=float).sel(layer=[1, 2])
constant_head[..., 0] = 0.0

# Drainage
elevation = xr.full_like(idomain.sel(layer=1), np.nan, dtype=float)
conductance = xr.full_like(idomain.sel(layer=1), np.nan, dtype=float)
elevation[7, 1:10] = np.array([0.0, 0.0, 10.0, 20.0, 30.0, 50.0, 70.0, 90.0, 100.0])
conductance[7, 1:10] = 1.0

# Recharge
rch_rate = xr.full_like(idomain.sel(layer=1), 3.0e-8, dtype=float)

# Well
# fmt: off
wells_x = [52500.0, 27500.0, 57500.0, 37500.0, 47500.0, 57500.0, 67500.0, 37500.0,
           47500.0, 57500.0, 67500.0, 37500.0, 47500.0, 57500.0, 67500.0, ]
wells_y = [52500.0, 57500.0, 47500.0, 32500.0, 32500.0, 32500.0, 32500.0, 22500.0,
           22500.0, 22500.0, 22500.0, 12500.0, 12500.0, 12500.0, 12500.0, ]
screen_top = [-300.0, -200.0, -200.0, 200.0, 200.0, 200.0, 200.0, 200.0,
              200.0, 200.0, 200.0, 200.0, 200.0, 200.0, 200.0, ]
screen_bottom = [-450.0, -300.0, -300.0, -200.0, -200.0, -200.0, -200.0, -200.0,
                 -200.0, -200.0, -200.0, -200.0, -200.0, -200.0, -200.0, ]
rate_wel = [-5.0, -5.0, -5.0, -5.0, -5.0, -5.0, -5.0, -5.0,
            -5.0, -5.0, -5.0, -5.0, -5.0, -5.0, -5.0, ]
# fmt: on

# Node properties
icelltype = xr.DataArray([1, 0, 0], {"layer": layer}, ("layer",))
k = xr.DataArray([1.0e-3, 1.0e-4, 2.0e-4], {"layer": layer}, ("layer",))
k33 = xr.DataArray([2.0e-8, 2.0e-8, 2.0e-8], {"layer": layer}, ("layer",))

# %%
# Write the model
# ---------------
#
# The first step is to define an empty model, the parameters and boundary
# conditions are added in the form of the familiar MODFLOW packages.

gwf_model = imod.mf6.GroundwaterFlowModel()
gwf_model["dis"] = imod.mf6.StructuredDiscretization(
    top=200.0, bottom=bottom, idomain=idomain
)
gwf_model["chd"] = imod.mf6.ConstantHead(
    constant_head, print_input=True, print_flows=True, save_flows=True
)
gwf_model["drn"] = imod.mf6.Drainage(
    elevation=elevation,
    conductance=conductance,
    print_input=True,
    print_flows=True,
    save_flows=True,
)
gwf_model["ic"] = imod.mf6.InitialConditions(start=0.0)
gwf_model["npf"] = imod.mf6.NodePropertyFlow(
    icelltype=icelltype,
    k=k,
    k33=k33,
    variable_vertical_conductance=True,
    dewatered=True,
    perched=True,
    save_flows=True,
)
gwf_model["sto"] = imod.mf6.SpecificStorage(
    specific_storage=1.0e-5,
    specific_yield=0.15,
    transient=False,
    convertible=0,
)
gwf_model["oc"] = imod.mf6.OutputControl(save_head="all", save_budget="all")
gwf_model["rch"] = imod.mf6.Recharge(rch_rate)

gwf_model["wel"] = imod.mf6.Well(
    x=wells_x,
    y=wells_y,
    screen_top=screen_top,
    screen_bottom=screen_bottom,
    rate=rate_wel,
    minimum_k=0.0001,
)

# Attach it to a simulation
simulation = imod.mf6.Modflow6Simulation("ex01-twri")
simulation["GWF_1"] = gwf_model
# Define solver settings
simulation["solver"] = imod.mf6.Solution(
    modelnames=["GWF_1"],
    print_option="summary",
    outer_dvclose=1.0e-4,
    outer_maximum=500,
    under_relaxation=None,
    inner_dvclose=1.0e-4,
    inner_rclose=0.001,
    inner_maximum=100,
    linear_acceleration="cg",
    scaling_method=None,
    reordering_method=None,
    relaxation_factor=0.97,
)
# Collect time discretization
simulation.create_time_discretization(
    additional_times=["2000-01-01", "2000-01-02", "2000-01-03", "2000-01-04"]
)

# %%
# We'll create a new directory in which we will write and run the model.

modeldir = imod.util.temporary_directory()
simulation.write(modeldir)

# %%
# Run the model
# -------------
#
# .. note::
#
#   The following lines assume the ``mf6`` executable is available on your PATH.
#   :ref:`The Modflow 6 examples introduction <mf6-introduction>` shortly
#   describes how to add it to yours.

simulation.run()

# %%
# Open the results
# ----------------
#
# We'll open the heads (.hds) file.

head = imod.mf6.open_hds(
    modeldir / "GWF_1/GWF_1.hds",
    modeldir / "GWF_1/dis.dis.grb",
)

# %%
# Visualize the results
# ---------------------

head.isel(layer=0, time=0).plot.contourf()


# %%
# .. _MODFLOW6 example problems: https://github.com/MODFLOW-USGS/modflow6-examples
# .. _description: https://modflow6-examples.readthedocs.io/en/master/_examples/ex-gwf-twri.html
# .. _notebook: https://github.com/MODFLOW-USGS/modflow6-examples/tree/master/notebooks/ex-gwf-twri.ipynb
# .. _Techniques of Water-Resources Investigation: https://pubs.usgs.gov/twri/twri7-c1/
# .. _McDonald & Harbaugh, 1988: https://pubs.er.usgs.gov/publication/twri06A1
# .. _Harbaugh & McDonald, 1996: https://pubs.er.usgs.gov/publication/ofr96485
# .. _Harbaugh, 2005: https://pubs.er.usgs.gov/publication/tm6A16
# .. _FloPy: https://github.com/modflowpy/flopy

# %%


"""
1D Solute Transport Benchmarks
==============================

This example is taken from the MODFLOW6 Examples, number 35.

As explained there, the setup is a simple 1d homogeneous aquifer with a steady
state flow field of constant velocity. The benchmark consists of four transport
problems that are modeled using this flow field. Here we have modeled these
four transport problems as a single simulation with multiple species. In all
cases the initial concentration in the domain is zero, but water entering the
domain has a concentration of one:

* species a is transported with zero diffusion or dispersion and the
  concentration distribution should show a sharp front, but due to the
  numerical method we see some smearing, which is expected.
* species b has a sizeable dispersivity and hence shows more smearing than
  species a but the same centre of mass.
* Species c has linear sorption and therefore the concentration doesn't enter
  the domain as far as species a or species b, but the front of the solute
  plume has the same overall shape as for species a or species b.
* Species d has linear sorption and first order decay, and this changes the
  shape of the front of the solute plume.

"""

# %%

import numpy as np
import pandas as pd
import xarray as xr

import imod


def create_transport_model(flowmodel, speciesname, dispersivity, retardation, decay):
    """
    Function to create a transport model, as we intend to create four similar
    transport models.

    Parameters
    ----------
    flowmodel: GroundwaterFlowModel
    speciesname: str
    dispersivity: float
    retardation: float
    decay: float

    Returns
    -------
    transportmodel: GroundwaterTransportModel
    """

    rhobulk = 1150.0
    porosity = 0.25

    tpt_model = imod.mf6.GroundwaterTransportModel()
    tpt_model["ssm"] = imod.mf6.SourceSinkMixing.from_flow_model(
        flowmodel, speciesname, save_flows=True
    )
    tpt_model["adv"] = imod.mf6.AdvectionUpstream()
    tpt_model["dsp"] = imod.mf6.Dispersion(
        diffusion_coefficient=0.0,
        longitudinal_horizontal=dispersivity,
        transversal_horizontal1=0.0,
        xt3d_off=False,
        xt3d_rhs=False,
    )

    # Compute the sorption coefficient based on the desired retardation factor
    # and the bulk density. Because of this, the exact value of bulk density
    # does not matter for the solution.
    if retardation != 1.0:
        sorption = "linear"
        kd = (retardation - 1.0) * porosity / rhobulk
    else:
        sorption = None
        kd = 1.0

    tpt_model["mst"] = imod.mf6.MobileStorageTransfer(
        porosity=porosity,
        decay=decay,
        decay_sorbed=decay,
        bulk_density=rhobulk,
        distcoef=kd,
        first_order_decay=True,
        sorption=sorption,
    )

    tpt_model["ic"] = imod.mf6.InitialConditions(start=0.0)
    tpt_model["oc"] = imod.mf6.OutputControl(
        save_concentration="all", save_budget="last"
    )
    tpt_model["dis"] = flowmodel["dis"]
    return tpt_model


# %%
# Create the spatial discretization.

nlay = 1
nrow = 2
ncol = 101
dx = 10.0
xmin = 0.0
xmax = dx * ncol
layer = [1]
y = [1.5, 0.5]
x = np.arange(xmin, xmax, dx) + 0.5 * dx

grid_dims = ("layer", "y", "x")
grid_coords = {"layer": layer, "y": y, "x": x}
grid_shape = (nlay, nrow, ncol)
grid = xr.DataArray(np.ones(grid_shape, dtype=int), coords=grid_coords, dims=grid_dims)
bottom = xr.full_like(grid, -1.0, dtype=float)

gwf_model = imod.mf6.GroundwaterFlowModel()
gwf_model["ic"] = imod.mf6.InitialConditions(0.0)

# %%
# Create the input for a constant head boundary and its associated concentration.
constant_head = xr.full_like(grid, np.nan, dtype=float)
constant_head[..., 0] = 60.0
constant_head[..., 100] = 0.0

constant_conc = xr.full_like(grid, np.nan, dtype=float)
constant_conc[..., 0] = 1.0
constant_conc[..., 100] = 0.0
constant_conc = constant_conc.expand_dims(
    species=["species_a", "species_b", "species_c", "species_d"]
)

gwf_model["chd"] = imod.mf6.ConstantHead(constant_head, constant_conc)

# %%
# Add other flow packages.

gwf_model["npf"] = imod.mf6.NodePropertyFlow(
    icelltype=1,
    k=xr.full_like(grid, 1.0, dtype=float),
    variable_vertical_conductance=True,
    dewatered=True,
    perched=True,
)
gwf_model["dis"] = imod.mf6.StructuredDiscretization(
    top=0.0,
    bottom=bottom,
    idomain=grid,
)
gwf_model["oc"] = imod.mf6.OutputControl(save_head="all", save_budget="all")
gwf_model["sto"] = imod.mf6.SpecificStorage(
    specific_storage=1.0e-5,
    specific_yield=0.15,
    transient=False,
    convertible=0,
)
# %%
# Create the simulation.

simulation = imod.mf6.Modflow6Simulation("1d_tpt_benchmark")
simulation["flow"] = gwf_model

# %%
# Add four transport simulations, and setup the solver flow and transport.

simulation["tpt_a"] = create_transport_model(gwf_model, "species_a", 0.0, 1.0, 0.0)
simulation["tpt_b"] = create_transport_model(gwf_model, "species_b", 10.0, 1.0, 0.0)
simulation["tpt_c"] = create_transport_model(gwf_model, "species_c", 10.0, 5.0, 0.0)
simulation["tpt_d"] = create_transport_model(gwf_model, "species_d", 10.0, 5.0, 0.002)

simulation["flow_solver"] = imod.mf6.Solution(
    modelnames=["flow"],
    print_option="summary",
    outer_dvclose=1.0e-4,
    outer_maximum=500,
    under_relaxation=None,
    inner_dvclose=1.0e-4,
    inner_rclose=0.001,
    inner_maximum=100,
    linear_acceleration="bicgstab",
    scaling_method=None,
    reordering_method=None,
    relaxation_factor=0.97,
)
simulation["transport_solver"] = imod.mf6.Solution(
    modelnames=["tpt_a", "tpt_b", "tpt_c", "tpt_d"],
    print_option="summary",
    outer_dvclose=1.0e-4,
    outer_maximum=500,
    under_relaxation=None,
    inner_dvclose=1.0e-4,
    inner_rclose=0.001,
    inner_maximum=100,
    linear_acceleration="bicgstab",
    scaling_method=None,
    reordering_method=None,
    relaxation_factor=0.97,
)

duration = pd.to_timedelta("2000d")
start = pd.to_datetime("2000-01-01")
simulation.create_time_discretization(additional_times=[start, start + duration])
simulation["time_discretization"]["n_timesteps"] = 100

# %%
# Run the simulation.
modeldir = imod.util.temporary_directory()
simulation.write(modeldir, binary=False)
simulation.run()

# %%
# Open the concentration results and store them in a single DataArray.

concentration = simulation.open_concentration(species_ls=["a", "b", "c", "d"])
mass_budgets = simulation.open_transport_budget(species_ls=["a", "b", "c", "d"])

# %%
# Visualize the last concentration profiles of the model run for the different
# species.

concentration.isel(time=-1, y=0).plot(x="x", hue="species")

# %%


"""
Example models
==============

This source file contains functions that create a simulation that can be used in
examples that are not focused on building a simulation, but on doing something
with it (such as regridding).

"""

import numpy as np
import xarray as xr
import xugrid as xu

import imod


def create_twri_simulation() -> imod.mf6.Modflow6Simulation:
    """There is a separate example contained in :doc:`/examples/mf6/ex01_twri`
    that you should look at if you are interested in the model building. The
    TWRI model has 3 layers and contains wells, a drain and recharge.
    Geometrically it is rectangular with prescribed head on some of the
    boundaries. Conductivity is highly anisotropic but constant in each layer.
    Simulation is steady state.
    """

    nlay = 3
    nrow = 15
    ncol = 15
    shape = (nlay, nrow, ncol)

    dx = 5000.0
    dy = 5000.0
    xmin = 0.0
    xmax = dx * ncol
    ymin = 0.0
    ymax = dy * nrow
    dims = ("layer", "y", "x")

    layer = np.array([1, 2, 3])
    y = np.arange(ymax, ymin, -dy) - 0.5 * dy
    x = np.arange(xmin, xmax, dx) + 0.5 * dx
    coords = {"layer": layer, "y": y, "x": x, "dx": dx, "dy": -dy}

    # %%
    # Create DataArrays
    # -----------------
    #
    # Now that we have the grid coordinates setup, we can start defining model
    # parameters. The model is characterized by:
    #
    # * a constant head boundary on the left
    # * a single drain in the center left of the model
    # * uniform recharge on the top layer
    # * a number of wells scattered throughout the model.

    idomain = xr.DataArray(np.ones(shape, dtype=int), coords=coords, dims=dims)
    bottom = xr.DataArray([-200.0, -300.0, -450.0], {"layer": layer}, ("layer",))

    # Constant head
    constant_head = xr.full_like(idomain, np.nan, dtype=float).sel(layer=[1, 2])
    constant_head[..., 0] = 0.0

    # Drainage
    elevation = xr.full_like(idomain.sel(layer=1), np.nan, dtype=float)
    conductance = xr.full_like(idomain.sel(layer=1), np.nan, dtype=float)
    elevation[7, 1:10] = np.array([0.0, 0.0, 10.0, 20.0, 30.0, 50.0, 70.0, 90.0, 100.0])
    conductance[7, 1:10] = 1.0

    # Recharge
    rch_rate = xr.full_like(idomain.sel(layer=1), 3.0e-8, dtype=float)

    # Well
    screen_layer = [2, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

    # we set the screen top and bottoms such that each well falls in one layer and is long enough not to be filtered out
    perforation_length = 50
    delta_z = 0.1

    screen_bottom = bottom[screen_layer] + delta_z
    screen_top = screen_bottom + delta_z + perforation_length

    # we compute the x and y cooordinates of the wells based on the row and column indices presented in the original twri model
    well_y = (
        ymax
        - np.array(
            [
                5.0,
                4.0,
                6.0,
                9.0,
                9.0,
                9.0,
                9.0,
                11.0,
                11.0,
                11.0,
                11.0,
                13.0,
                13.0,
                13.0,
                13.0,
            ]
        )
        * abs(dy)
        + dy / 2
    )
    well_x = (
        np.array(
            [
                11.0,
                6.0,
                12.0,
                8.0,
                10.0,
                12.0,
                14.0,
                8.0,
                10.0,
                12.0,
                14.0,
                8.0,
                10.0,
                12.0,
                14.0,
            ]
        )
        * dx
        - dx / 2
    )
    well_rate = [-5.0] * 15

    # Node properties
    icelltype = xr.DataArray([1, 0, 0], {"layer": layer}, ("layer",))
    k = xr.DataArray([1.0e-3, 1.0e-4, 2.0e-4], {"layer": layer}, ("layer",))
    k33 = xr.DataArray([2.0e-8, 2.0e-8, 2.0e-8], {"layer": layer}, ("layer",))

    # %%
    # Build the model
    # ---------------
    #
    # The first step is to define an empty model, the parameters and boundary
    # conditions are added in the form of the familiar MODFLOW packages.

    gwf_model = imod.mf6.GroundwaterFlowModel()
    gwf_model["dis"] = imod.mf6.StructuredDiscretization(
        top=200.0, bottom=bottom, idomain=idomain
    )
    gwf_model["chd"] = imod.mf6.ConstantHead(
        constant_head, print_input=True, print_flows=True, save_flows=True
    )
    gwf_model["drn"] = imod.mf6.Drainage(
        elevation=elevation,
        conductance=conductance,
        print_input=True,
        print_flows=True,
        save_flows=True,
    )
    gwf_model["ic"] = imod.mf6.InitialConditions(start=0.0)
    gwf_model["npf"] = imod.mf6.NodePropertyFlow(
        icelltype=icelltype,
        k=k,
        k33=k33,
        variable_vertical_conductance=True,
        dewatered=True,
        perched=True,
        save_flows=True,
    )
    gwf_model["sto"] = imod.mf6.SpecificStorage(
        specific_storage=1.0e-5,
        specific_yield=0.15,
        transient=False,
        convertible=0,
    )
    gwf_model["oc"] = imod.mf6.OutputControl(save_head="all", save_budget="all")
    gwf_model["rch"] = imod.mf6.Recharge(rch_rate)
    gwf_model["wel"] = imod.mf6.Well(
        x=well_x,
        y=well_y,
        screen_top=screen_top,
        screen_bottom=screen_bottom,
        rate=well_rate,
        minimum_k=0.0001,
    )

    # %%
    # Attach it to a simulation
    # ---------------

    simulation = imod.mf6.Modflow6Simulation("ex01-twri")
    simulation["GWF_1"] = gwf_model
    # Define solver settings
    simulation["solver"] = imod.mf6.Solution(
        modelnames=["GWF_1"],
        print_option="summary",
        outer_dvclose=1.0e-4,
        outer_maximum=500,
        under_relaxation=None,
        inner_dvclose=1.0e-4,
        inner_rclose=0.001,
        inner_maximum=100,
        linear_acceleration="cg",
        scaling_method=None,
        reordering_method=None,
        relaxation_factor=0.97,
    )
    # Collect time discretization
    simulation.create_time_discretization(
        additional_times=["2000-01-01", "2000-01-02", "2000-01-03", "2000-01-04"]
    )
    return simulation


def create_circle_simulation():
    """
    There is a separate example contained in :doc:`/examples/mf6/circle`
    that you should look at if you are interested in the model building. The
    circle model uses an unstructured grid. It has 2 layers of constant
    thickness. In conductivity, it is isotropic, and constant in space.
    Boundary conditions include recharge and constant head. It is a transient
    simulation.
    """

    # %%
    # Create a mesh
    # -------------

    grid = imod.data.circle()

    nface = grid.n_face
    nlayer = 2

    idomain = xu.UgridDataArray(
        xr.DataArray(
            np.ones((nlayer, nface), dtype=np.int32),
            coords={"layer": [1, 2]},
            dims=["layer", grid.face_dimension],
        ),
        grid=grid,
    )

    # Create model and packages
    icelltype = xu.full_like(idomain, 0)
    k = xu.full_like(idomain, 1.0, dtype=float)
    k33 = k.copy()
    rch_rate = xu.full_like(idomain.sel(layer=1), 0.001, dtype=float)
    bottom = idomain * xr.DataArray([5.0, 0.0], dims=["layer"])

    chd_location = xu.zeros_like(
        idomain.sel(layer=2), dtype=bool
    ).ugrid.binary_dilation(border_value=True)
    constant_head = xu.full_like(idomain.sel(layer=2), 1.0, dtype=float).where(
        chd_location
    )

    gwf_model = imod.mf6.GroundwaterFlowModel()
    gwf_model["disv"] = imod.mf6.VerticesDiscretization(
        top=10.0, bottom=bottom, idomain=idomain
    )
    gwf_model["chd"] = imod.mf6.ConstantHead(
        constant_head, print_input=True, print_flows=True, save_flows=True
    )
    gwf_model["ic"] = imod.mf6.InitialConditions(start=0.0)
    gwf_model["npf"] = imod.mf6.NodePropertyFlow(
        icelltype=icelltype, k=k, k33=k33, save_flows=True, save_specific_discharge=True
    )
    gwf_model["sto"] = imod.mf6.SpecificStorage(
        specific_storage=1.0e-5,
        specific_yield=0.15,
        transient=False,
        convertible=0,
    )
    gwf_model["oc"] = imod.mf6.OutputControl(save_head="all", save_budget="all")
    gwf_model["rch"] = imod.mf6.Recharge(rch_rate)

    # Create simulation
    simulation = imod.mf6.Modflow6Simulation("circle")
    simulation["GWF_1"] = gwf_model
    simulation["solver"] = imod.mf6.Solution(
        modelnames=["GWF_1"],
        print_option="summary",
        outer_dvclose=1.0e-4,
        outer_maximum=500,
        under_relaxation=None,
        inner_dvclose=1.0e-4,
        inner_rclose=0.001,
        inner_maximum=100,
        linear_acceleration="cg",
        scaling_method=None,
        reordering_method=None,
        relaxation_factor=0.97,
    )
    simulation.create_time_discretization(additional_times=["2000-01-01", "2000-01-02"])

    return simulation


"""
Henry
=====

This example illustrates how to setup a variable density groundwater flow and
transport model using the ``imod`` package and associated packages.

In overview, we'll set the following steps:

    * Create a suitable 2d (x, z) grid.
    * Create a groundwater flow model, with variable density.
    * Create a solute transport model.
    * Combine these models into a single MODFLOW6 simulation.
    * Write to modflow6 files.
    * Run the model.
    * Open the results back into xarray DataArrays.
    * Visualize the results.

We are simulating the Henry problem, although not the original one but the one
outlined in the MODFLOW 6 manual (jupyter notebooks example 51). This is the
modified Henry problem with half the inflow rate.

The domain is a vertically oriented two dimensional rectangle, which is 2 m
long and 1 m high. Water flows in over the left boundary with a fixed rate,
which is represented by a Well package. The right boundary is in direct contact
with hydrostatic seawater with a density of 1025 kg m:sup:`-3`. This is
represented by a General Head Boundary package.
"""

# sphinx_gallery_thumbnail_number = -1

# %%
# We'll start with the usual imports. As this is a simple (synthetic)
# structured model, we can make due with few packages.

import numpy as np
import pandas as pd
import xarray as xr

import imod

# %%
# We'll start by defining the (vertical) rectangular domain and the physical
# parameters of the model.

nlay = 40
nrow = 1
ncol = 80
shape = (nlay, nrow, ncol)

total_flux = 5.7024  # m3/d
k = 864.0  # m/d
porosity = 0.35
max_concentration = 35.0
min_concentration = 0.0
max_density = 1025.0
min_density = 1000.0
diffusion_coefficient = 0.57024
longitudinal_horizontal = 0.1
transversal_horizontal1 = 0.01

# Time
start_date = pd.to_datetime("2020-01-01")
duration = pd.to_timedelta("0.5d")

# Domain size
xmax = 2.0
xmin = 0.0
dx = (xmax - xmin) / ncol
zmin = 0.0
zmax = 1.0
dz = (zmax - zmin) / nlay

x = np.arange(xmin, xmax, dx) + 0.5 * dx
y = np.array([0.5])
layer = np.arange(1, 41, 1)

dy = -1.0
coords = {"layer": layer, "y": y, "x": x, "dy": dy, "dx": dx}
dims = ("layer", "y", "x")
idomain = xr.DataArray(np.ones(shape, dtype=int), coords=coords, dims=dims)

top = 1.0
bottom = xr.DataArray(
    np.arange(zmin, zmax, dz)[::-1],
    {"layer": layer},
    ("layer",),
)

# %%
# Now make the flow model. We'll start with the non-boundary condition
# packages.

gwf_model = imod.mf6.GroundwaterFlowModel()
gwf_model["dis"] = imod.mf6.StructuredDiscretization(
    top=top, bottom=bottom, idomain=idomain
)
gwf_model["npf"] = imod.mf6.NodePropertyFlow(
    icelltype=0,
    k=k,
)
gwf_model["sto"] = imod.mf6.SpecificStorage(
    specific_storage=1.0e-4,
    specific_yield=0.15,
    transient=False,
    convertible=0,
)
gwf_model["ic"] = imod.mf6.InitialConditions(start=0.0)
gwf_model["oc"] = imod.mf6.OutputControl(save_head="last", save_budget="last")

# %%
# Now let's make the constant head boundary condition.

ghb_head = xr.ones_like(idomain, dtype=float)
ghb_head[:, :, :-1] = np.nan

ghb_conc = xr.full_like(idomain, max_concentration, dtype=float)
ghb_conc[:, :, :-1] = np.nan
ghb_conc = ghb_conc.expand_dims(species=["salinity"])

conductance = xr.full_like(idomain, 864.0 * 2.0, dtype=float)
conductance[:, :, :-1] = np.nan

gwf_model["right_boundary"] = imod.mf6.GeneralHeadBoundary(
    head=ghb_head,
    conductance=conductance,
    concentration=ghb_conc,
    concentration_boundary_type="AUX",
    print_input=True,
    print_flows=True,
    save_flows=True,
)

# %%
# ... and the constant flux condition.
from imod.prepare.layer import create_layered_top

screen_top = create_layered_top(bottom, top)

flux_concentration = xr.DataArray(
    data=np.full((1, nlay), min_concentration),
    dims=["species", "index"],
    coords={"species": ["salinity"], "index": layer},
)

gwf_model["left_boundary"] = imod.mf6.Well(
    x=np.full_like(layer, xmin, dtype=float),
    y=np.full_like(layer, y[0], dtype=float),
    screen_top=screen_top.values,
    screen_bottom=bottom.values,
    rate=np.full_like(layer, 0.5 * (total_flux / nlay), dtype=float),
    concentration=flux_concentration,
    concentration_boundary_type="AUX",
    minimum_thickness=0.002,
)

# %%
# Since we are solving a variable density problem, we need to add the buoyancy
# package. It will use the species "salinity" that we are simulating with a
# transport model defined below.

slope = (max_density - min_density) / (max_concentration - min_concentration)
gwf_model["buoyancy"] = imod.mf6.Buoyancy(
    reference_density=min_density,
    modelname=["transport"],
    reference_concentration=[min_concentration],
    density_concentration_slope=[slope],
    species=["salinity"],
)

# %%
# Now let's make the transport model. It contains the standard packages of
# storage, dispersion and advection as well as initial condiations and output
# control. Sinks and sources are automatically determined based on packages
# provided in the flow model.

gwt_model = imod.mf6.GroundwaterTransportModel()
gwt_model["ssm"] = imod.mf6.SourceSinkMixing.from_flow_model(gwf_model, "salinity")
gwt_model["adv"] = imod.mf6.AdvectionTVD()
gwt_model["dsp"] = imod.mf6.Dispersion(
    diffusion_coefficient=0.57024,
    longitudinal_horizontal=0.1,
    transversal_horizontal1=0.01,
    xt3d_off=False,
    xt3d_rhs=False,
)

gwt_model["mst"] = imod.mf6.MobileStorageTransfer(
    porosity=porosity,
)
gwt_model["ic"] = imod.mf6.InitialConditions(start=max_concentration)
gwt_model["oc"] = imod.mf6.OutputControl(save_concentration="last", save_budget="last")
gwt_model["dis"] = gwf_model["dis"]

# %%
# now let's define a simulation using the flow and transport models.

# Attach it to a simulation
simulation = imod.mf6.Modflow6Simulation("henry")

simulation["flow"] = gwf_model
simulation["transport"] = gwt_model

# %%
# Define solver settings. We need to define separate solutions for the flow and
# transport models. In this case, we'll use the same settings, but generally
# convergence settings should differ: the transport model has very different
# units from the flow model.

simulation["flow_solver"] = imod.mf6.Solution(
    modelnames=["flow"],
    print_option="summary",
    outer_dvclose=1.0e-6,
    outer_maximum=500,
    under_relaxation=None,
    inner_dvclose=1.0e-6,
    inner_rclose=1.0e-10,
    inner_maximum=100,
    linear_acceleration="bicgstab",
    scaling_method=None,
    reordering_method=None,
    relaxation_factor=0.97,
)
simulation["transport_solver"] = imod.mf6.Solution(
    modelnames=["transport"],
    print_option="summary",
    outer_dvclose=1.0e-6,
    outer_maximum=500,
    under_relaxation=None,
    inner_dvclose=1.0e-6,
    inner_rclose=1.0e-10,
    inner_maximum=100,
    linear_acceleration="bicgstab",
    scaling_method=None,
    reordering_method=None,
    relaxation_factor=0.97,
)
# Collect time discretization
times = [start_date, start_date + duration]
simulation.create_time_discretization(additional_times=times)

# %%
# Increase the number of time steps for the single stress period:
simulation["time_discretization"]["n_timesteps"] = 500

# %%
# We'll create a new directory in which we will write and run the model.
modeldir = imod.util.temporary_directory()
simulation.write(modeldir, binary=False)

# %%
# Run the model
# -------------
#
# This takes about 20 seconds.
#
# .. note::
#
#   The following lines assume the ``mf6`` executable is available on your PATH.
#   :ref:`The Modflow 6 examples introduction <mf6-introduction>` shortly
#   describes how to add it to yours.

simulation.run()

# %%
# Open the results
# ----------------
#
# We'll open the head and concentration files.

head = simulation.open_head()

conc = simulation.open_concentration()

# %%
# Visualize the results
# ---------------------

head.isel(y=0, time=-1).plot.contourf(yincrease=False)

# %%
# We can check the concentration to see that a fresh-saline interface has been
# formed:

conc.isel(y=0, time=-1).plot.contourf(yincrease=False, cmap="RdYlBu_r")

# %%


"""
Regional model
==============

This example shows a simplified script for building a groundwater model in the
northeast of the Netherlands. A primary feature of this area is an ice-pushed
ridge called the Hondsrug. This examples demonstrates modifying external data
for use in a MODFLOW6 model.

In overview, the model features:

    * Thirteen layers: seven aquifers and six aquitards;
    * A dense ditch network in the east;
    * Pipe drainage for agriculture;
    * Precipitation and evapotranspiration.

"""

# sphinx_gallery_thumbnail_number = -1

# %% Import packages
# We'll start with the usual imports, and an import from scipy.

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.ndimage
import xarray as xr

import imod

# %%
# Before starting to create the input data, we will create the groundwater
# :doc:`/api/generated/mf6/imod.mf6.GroundwaterFlowModel`.
# The data from all the model packages will be added to this variable.

gwf_model = imod.mf6.GroundwaterFlowModel()

# %%
# This package allows specifying a regular MODFLOW grid. This grid is assumed
# to be rectangular horizontally, but can be distorted vertically.
#
# Load data
# ---------
#
# We'll load the data from the examples that come with this package.

layermodel = imod.data.hondsrug_layermodel()

# Make sure that the idomain is provided as integers
idomain = layermodel["idomain"].astype(int)

# We only need to provide the data for the top as a 2D array. Modflow 6 will
# compare the top against the uppermost active bottom cell.
top = layermodel["top"].max(dim="layer")

bot = layermodel["bottom"]

top.plot.imshow()

# %%
# Adding information to the DIS package
# -------------------------------------
#
# The following step is to add the previously created discretization data to
# the gwf_model variable.  This is done using the function
# :doc:`/api/generated/mf6/imod.mf6.StructuredDiscretization`.
# The data to include is the top of the model domain, the bottom of the layers,
# and the idomain. All this information comes from the previously imported
# tifs (now converted to `xarray.DataArray
# <http://xarray.pydata.org/en/stable/generated/xarray.DataArray.html#xarray.DataArray>`_.

gwf_model["dis"] = imod.mf6.StructuredDiscretization(
    top=top, bottom=bot, idomain=idomain
)

# %%
# Node property flow package - NPF
# =================================
#
# This package contains the information related to the aquifer properties used to calculate
# hydraulic conductance. This package replaces the Layer Property Flow (LPF),
# Block-Centered Flow (BCF), and Upstream Weighting (UPW) packages from previous MODFLOW versions.
#
# Hydraulic conductivity
# ----------------------
#
k = layermodel["k"]

# %%
# icelltype
# ----------
#
# The cell type to be used in the model (confined or convertible) can be
# defined under ICELLTYPE, which is an input to the NPF package.  ICELLTYPE ==
# 0: *Confined cell* - Constant transmissivity ICELLTYPE != 0: *Convertible
# cell* - Transmissivity varies depending on the calculated head in the cell
# (based on the saturated cell thickness)
#
# In this example, all layers have a ICELLTYPE equal to 0, indicating confined
# cells.  This is defined in the following section.
#
# Adding information to the NPF package
# --------------------------------------
#
# The information for the NPF package is added to the gwf_model variable using
# :doc:`/api/generated/mf6/imod.mf6.NodePropertyFlow`.
# The information included is the icelltype value (equal to zero), the array
# for  the hydraulic conductivity (considered to be the same for the horizontal
# and vertical direction) and, optionally, the
# variable_vertical_conductance, dewatered, perched and save_flows options have
# been activated.  For more details about the meaning of these variables and
# other variables available to be used within this package, please refer to the
# :doc:`documentation </api/generated/mf6/imod.mf6.NodePropertyFlow>`.

gwf_model["npf"] = imod.mf6.NodePropertyFlow(
    icelltype=0,
    k=k,
    k33=k,
    variable_vertical_conductance=True,
    dewatered=True,
    perched=True,
    save_flows=True,
)

# %%
# Initial conditions package - IC
# ================================
#
# This package reads the starting heads for a simulation.
#
# Starting heads interpolation
# ----------------------------
#
# The starting heads to be used in this model are based on the interpolation of
# x-y head measurements, which were interpolated on a larger area.  This
# example was created in this example --insert-link-here--
#
# The heads were interpolated on a larger area, therefore these have to be
# clipped first

initial = imod.data.hondsrug_initial()
interpolated_head_larger = initial["head"]

xmin = 237_500.0
xmax = 250_000.0
ymin = 559_000.0
ymax = 564_000.0

interpolated_head = interpolated_head_larger.sel(
    x=slice(xmin, xmax), y=slice(ymax, ymin)
)

# Plotting the clipped interpolation
fig, ax = plt.subplots()
interpolated_head.plot.imshow(ax=ax)

# %%
# The final step is to assign the 2D heads interpolation to all the
# model layers (as a reference value) by using the xarray tool
# `xarray.full_like <http://xarray.pydata.org/en/stable/generated/xarray.full_like.html#xarray.full_like>`_.
# The 3d idomain array is used as reference for the geometry and then
# its original values are replaced by NaNs.
# This array is combined with the interpolated_head array using the xarray
# `DataArray.combine_first <http://xarray.pydata.org/en/stable/generated/xarray.DataArray.combine_first.html#xarray.DataArray.combine_first>`_
# option.
# The final result is an starting_heads xarray where all layers have the 2d interpolated_head information.

# Assign interpolated head values to all the model layers
like_3d = xr.full_like(idomain, np.nan, dtype=float)
starting_head = like_3d.combine_first(interpolated_head)
# Consequently ensure no data is specified in inactive cells:
starting_head = starting_head.where(idomain == 1)

starting_head

# %%
# Adding information to the IC package
# ------------------------------------
#
# The function for indicating the initial conditions is
# :doc:`/api/generated/mf6/imod.mf6.InitialConditions`.
# It is necessary to indicate the value(s) to be considered as the initial
# (starting) head of the simulation.
# In this case, this value is equal to the previously created starting_head array.

gwf_model["ic"] = imod.mf6.InitialConditions(starting_head)

# %%
# Constant head package - CHD
# ===========================
#
# This package allows to indicate if the head varies with time,
# if it is constant or if it is inactive.
#
# Constant head edge
# -------------------
#
# The previously interpolated starting_head array will be used to define
# the constant head value which will be used along the model boundaries.
# A function is defined to indicate the location of the outer edge
# (returning a boolean array).


def outer_edge(da):
    data = da.copy()
    from_edge = scipy.ndimage.binary_erosion(data)
    is_edge = (data == 1) & (from_edge == 0)
    return is_edge.astype(bool)


# %%
# For the next calculations, it is necessary to create a template array
# which can be used for assigning the corresponding geometry to other arrays.
# In this case, a 2d template is created based on the idomain layer information
# and filled with ones.

like_2d = xr.full_like(idomain.isel(layer=0), 1)
like_2d

# %%
# Using the previously created function and the 2d template,
# the outer edge is defined for this example.

edge = outer_edge(xr.full_like(like_2d.drop_vars("layer"), 1))

# %%
# Adding information to the CHD package
# --------------------------------------
#
# To add the information to the CHD package within the gwf_model variable, the
# :doc:`/api/generated/mf6/imod.mf6.ConstantHead`.
# function is used.
# The required information is the head array for this boundary condition.
# In this example, the starting_head array is selected where the idomain is > 0 (active)
# and it is located in the edge array.
#
# It is also possible (and optional) to indicate if the CHD information will be written
# to the listing file after it is read (print_input), if the constant head flow rates will
# be printed to the listing file for every stress period time step
# in which “BUDGET PRINT” is specified in Output Control (print_flows)
# and if the constant head flow terms will be written to the file
# specified with “BUDGET FILEOUT” in Output Control (save_flows).
# By default, these three options are set to False.

gwf_model["chd"] = imod.mf6.ConstantHead(
    starting_head.where((idomain > 0) & edge),
    print_input=False,
    print_flows=True,
    save_flows=True,
)

# %%
#
# Recharge
# ========
#
# This package is used to represent areally distributed recharge to the groundwater system.
# To calculate the recharge, the precipitation and evapotranspiration
# information from the KNMI website has been downloaded for the study area.
# This information is saved in netCDF files, which have been imported using the
# xarray function
# `xr.open_dataset <http://xarray.pydata.org/en/stable/generated/xarray.open_dataset.html#xarray.open_dataset>`_,
# slicing the area to the model's miminum and maximum dimensions.
#
# Note that the meteorological data has mm/d as unit and
# this has to be converted to m/d for Modflow 6.

xmin = 230_000.0
xmax = 257_000.0
ymin = 550_000.0
ymax = 567_000.0

meteorology = imod.data.hondsrug_meteorology()
pp = meteorology["precipitation"]
evt = meteorology["evapotranspiration"]

pp = pp.sel(x=slice(xmin, xmax), y=slice(ymax, ymin)) / 1000.0  # from mm/d to m/d
evt = evt.sel(x=slice(xmin, xmax), y=slice(ymax, ymin)) / 1000.0  # from mm/d to m/d

# %%
# Recharge - Steady state
# -----------------------
#
# For the steady state conditions of the model,
# the data from the period 2000 to 2009 was considered as reference.
# The initial information was sliced to this time period and averaged
# to obtain the a mean value grid. This process was done for both
# precipitation and evapotranspiration datasets.
#
# **Precipitation**
pp_ss = pp.sel(time=slice("2000-01-01", "2009-12-31"))
pp_ss_mean = pp_ss.mean(dim="time")

fig, ax = plt.subplots()
pp_ss_mean.plot(ax=ax)

# %%
# **Evapotranspiration**
evt_ss = evt.sel(time=slice("2000-01-01", "2009-12-31"))
evt_ss_mean = evt_ss.mean(dim="time")

fig, ax = plt.subplots()
evt_ss_mean.plot(ax=ax)

# %%
# For the recharge calculation, a first estimate
# is the difference between the precipitation and evapotranspiration values.

rch_ss = pp_ss_mean - evt_ss_mean

fig, ax = plt.subplots()
rch_ss.plot.imshow(ax=ax)

# %%
# Recharge - Transient
# --------------------
#
# The transient model will encompass the period from 2010 to 2015.
# The initial pp and evt datasets have been sliced to this time frame.

pp_trans = pp.sel(time=slice("2010-01-01", "2015-12-31"))
evt_trans = evt.sel(time=slice("2010-01-01", "2015-12-31"))

# %%
# As previously done, it is assumed that the recharge is equal
# to the difference between precipitation and evapotranspiration as a first estimate.
# Furthermore, the negative values found after doing this calculation have been
# replaced by zeros, as the recharge should not have a negative value.

rch_trans = pp_trans - evt_trans
rch_trans = rch_trans.where(rch_trans > 0, 0)  # check negative values

# %%
# The original information is on a daily step, so it is going to be
# resampled to a yearly step by using the xarray function
# `Dataset.resample <http://xarray.pydata.org/en/stable/generated/xarray.Dataset.resample.html#xarray.Dataset.resample>`_.

rch_trans_yr = rch_trans.resample(time="A", label="left").mean()
rch_trans_yr

# %%
# To create the final recharge for the transient simulation,
# the steady state information needs to be concatenated to the transient recharge data.
# The steady state simulation will be run for one second.
# This is achieved by using the numpy
# `Timedelta function <https://numpy.org/doc/stable/reference/arrays.datetime.html>`_,
# first creating a time delta of 1 second, which is assigned to the steady state recharge information.
# This dataset is then concatenated using the xarray function
# `xarray.concat <http://xarray.pydata.org/en/stable/generated/xarray.concat.html#xarray.concat>`_
# to the transient information and indicating that the dimension to join is "time".

starttime = "2009-12-31"

# Add first steady-state
timedelta = np.timedelta64(1, "s")  # 1 second duration for initial steady-state
starttime_steady = np.datetime64(starttime) - timedelta
rch_ss = rch_ss.assign_coords(time=starttime_steady)

rch_ss_trans = xr.concat([rch_ss, rch_trans_yr], dim="time")
rch_ss_trans

# %%
# The data obtained from KNMI has different grid dimensions
# than the one considered in this example. To fix this,
# imod-python includes the option
# :doc:`/api/generated/prepare/imod.prepare.Regridder`,
# which modifies the original grid dimensions to a different one.
# It is also possible to define the regridding method such as
# ``nearest``, ``multilinear``, ``mean``, among others.
# In this case, ``mean`` was selected and the 2d template (like_2d)
# was used as reference, as this is the geometry to be considered in the model.

rch_ss_trans = imod.prepare.Regridder(method="mean").regrid(rch_ss_trans, like_2d)
rch_ss_trans

# %%
# The previously created recharge array is a 2D array
# that needs to be assigned to a 3D array. This is done using the xarray
# `DataArray.where <http://xarray.pydata.org/en/stable/generated/xarray.DataArray.where.html#xarray.DataArray.where>`_
# option, where the recharge values are applied to the cells where the
# idomain value is larger than zero (that is, the active cells) and for the uppermost
# active cell (indicated by the minimum layer number).

rch_total = rch_ss_trans.where(
    idomain["layer"] == idomain["layer"].where(idomain > 0).min("layer")
)
rch_total

# %%
# Finally, transposing the array dimensions using
# `DataArray.transpose <http://xarray.pydata.org/en/stable/generated/xarray.DataArray.transpose.html#xarray.DataArray.transpose>`_
# so they are in the correct order.

rch_total = rch_total.transpose("time", "layer", "y", "x")
rch_total

fig, ax = plt.subplots()
rch_total.isel(layer=2, time=6).plot.imshow(ax=ax)

# %%
# Adding information to the RCH package
# --------------------------------------
#
# The information for the RCH package is added with the function
# :doc:`/api/generated/mf6/imod.mf6.Recharge`.
# It is required to insert the recharge flux rate, and it is optional
# to include the print_input, print_flows and save_flows information.

gwf_model["rch"] = imod.mf6.Recharge(rch_total)

# %%
# Drainage package - DRN
# =======================
#
# The drain package is used to simulate features that remove water from the aquifer,
# such as agricultural drains or springs.
# This occurs at a rate proportional to the head difference between the head in the
# aquifer and the drain elevation
# (the head in the aquifer has to be above that elevation).
# The conductance is the proportionality constant.
#
# Import drainage information
# ----------------------------

drainage = imod.data.hondsrug_drainage()

pipe_cond = drainage["conductance"]
pipe_elev = drainage["elevation"]

pipe_cond

# %%
# Adding information to the DRN package
# -------------------------------------
#
# To add the information to the DRN package within the gwf_model variable, the
# :doc:`/api/generated/mf6/imod.mf6.Drainage`.
# function is used. It is required to add the previously created arrays for
# the drain elevation and the drain conductance.
# It is optional to insert the information for
# ``print_input``, ``print_flows`` and ``save_flows``
# which are set to False by default.

gwf_model["drn-pipe"] = imod.mf6.Drainage(conductance=pipe_cond, elevation=pipe_elev)

# %%
# River package - RIV
# ===================
#
# This package simulates the effects of flow between
# surface-water features and groundwater systems.
#
# Import river information
# ------------------------

river = imod.data.hondsrug_river()
riv_cond = river["conductance"]
riv_stage = river["stage"]
riv_bot = river["bottom"]

# %%
# Adding information to the RIV package
# -------------------------------------
#
# The data is assigned to the gwf_model variable by using
# :doc:`/api/generated/mf6/imod.mf6.River`,
# based on the previously imported conductance, stage and bottom arrays.

gwf_model["riv"] = imod.mf6.River(
    conductance=riv_cond, stage=riv_stage, bottom_elevation=riv_bot
)

# %%
# Storage package - STO
# ======================
#
# When the STO Package is included in a model, storage changes
# will be calculated, and thus, the model will be transient.

ss = 0.0003
layer = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13])
sy = xr.DataArray(
    [0.16, 0.16, 0.16, 0.16, 0.15, 0.15, 0.15, 0.15, 0.14, 0.14, 0.14, 0.14, 0.14],
    {"layer": layer},
    ("layer",),
)
times_sto = np.array(
    [
        "2009-12-30T23:59:59.00",
        "2009-12-31T00:00:00.00",
        "2010-12-31T00:00:00.00",
        "2011-12-31T00:00:00.00",
        "2012-12-31T00:00:00.00",
        "2013-12-31T00:00:00.00",
        "2014-12-31T00:00:00.00",
    ],
    dtype="datetime64[ns]",
)

transient = xr.DataArray(
    [False, True, True, True, True, True, True], {"time": times_sto}, ("time",)
)

# %%
# Adding information to the STO package
# -------------------------------------
#
# The data is assigned to the gwf_model variable by using
# :doc:`/api/generated/mf6/imod.mf6.SpecificStorage`.
# It is necessary to indicate the values of specific storage,
# specific yield and if the layers are convertible.

gwf_model["sto"] = imod.mf6.SpecificStorage(
    specific_storage=ss,
    specific_yield=sy,
    transient=transient,
    convertible=0,
    save_flows=True,
)

# %%
# Output Control package - OC
# ===========================
#
# This package determines how and when heads are printed to the
# listing file and/or written to a separate binary output file
#
# Adding information to the OC package
# ------------------------------------
#
# The function
# :doc:`/api/generated/mf6/imod.mf6.OutputControl`
# is used to store the information for this package.
# It is possible to indicate if the heads and budget information is saved
# at the end of each stress period (``last``),
# for all timesteps a stress period (``all``),
# or at the start of a stress period (``first``)

gwf_model["oc"] = imod.mf6.OutputControl(save_head="last", save_budget="last")

# %%
# Model simulation
# ================
#
# In MODFLOW 6, the concept of "model" is that part of the program
# that solves a hydrologic process. MODFLOW 6 documentation supports
# one type of model: the GWF Model.
# It is possible within the MODFLOW 6 framewotk to solve multiple,
# tightly coupled, numerical models in a single system of equation,
# which may be multiple models of the same type or of different types.
#
# The previously created gwf_model variable now contains
# the information from all the variables.

gwf_model

# %%
# Attach the model information to a simulation
# --------------------------------------------
#
# The function
# :doc:`/api/generated/mf6/imod.mf6.Modflow6Simulation`
# allows to assign models to a simulation (in this case, the gwf_model).

simulation = imod.mf6.Modflow6Simulation("mf6-mipwa2-example")
simulation["GWF_1"] = gwf_model

# %%
# Solver settings
# ---------------
#
# The solver settings are indicated using
# :doc:`/api/generated/mf6/imod.mf6.Solution`.
# If the values are not indicated manually, the defaults values will be considered.

simulation["solver"] = imod.mf6.Solution(
    modelnames=["GWF_1"],
    print_option="summary",
    outer_dvclose=1.0e-4,
    outer_maximum=500,
    under_relaxation=None,
    inner_dvclose=1.0e-4,
    inner_rclose=0.001,
    inner_maximum=100,
    linear_acceleration="cg",
    scaling_method=None,
    reordering_method=None,
    relaxation_factor=0.97,
)

# %%
# Assign time discretization
# --------------------------
#
# The time discretization of this model is 6 years.

simulation.create_time_discretization(
    additional_times=["2009-12-30T23:59:59.000000000", "2015-12-31T00:00:00.000000000"]
)

# %%
# Run the model
# -------------
#
# .. note::
#
#   The following lines assume the ``mf6`` executable is available on your PATH.
#   :ref:`The Modflow 6 examples introduction <mf6-introduction>` shortly
#   describes how to add it to yours.

modeldir = imod.util.temporary_directory()
simulation.write(modeldir, binary=False)
simulation.run()

# %%
# Results visualization
# =====================
#
# The next section indicated how to visualize the model results.
#
# Import heads results
# --------------------
#
# The heads results are imported using
# :doc:`/api/generated/mf6/imod.mf6.open_hds`.
# on the background.

hds = simulation.open_head()

# %%
# We can plot the data of an individual layer as follows
fig, ax = plt.subplots()
hds.sel(layer=3).isel(time=3).plot(ax=ax)

# %%
# As you can see layer 3 has some missing cells in the west
# Whereas layer 4 only contains active cells in the
# eastern peatland area
fig, ax = plt.subplots()
hds.sel(layer=4).isel(time=3).plot(ax=ax)
# %%
# Layer 5 contains more data towards the west,
# but has no active cells in the centre.
fig, ax = plt.subplots()
hds.sel(layer=5).isel(time=3).plot(ax=ax)

# %%
# As you can see the data is individual layers
# have lots of inactive in different places.
#
# It is difficult for this model to get a good idea
# what is happening across the area based on 1 layer alone.
# Luckily xarray allows us to compute the mean across a selection
# of layers and plot this.
#
# By first selecting 3 layers with ``sel```,
# and then computing the mean across the layer dimension
# with ``mean(dim="layer")``.

fig, ax = plt.subplots()
hds.sel(layer=slice(3, 5)).mean(dim="layer").isel(time=3).plot(ax=ax)

# %%
# Assign dates to head
# --------------------
#
# MODFLOW6 has no concept of a calendar, so the output is not labelled only
# in terms of "time since start" in floating point numbers. For this model
# the time unit is days and we can assign a date coordinate as follows:

starttime = pd.to_datetime("2000-01-01")
timedelta = pd.to_timedelta(hds["time"], "D")
hds = hds.assign_coords(time=starttime + timedelta)

# %%
# Extract head at points
# ----------------------
#
# A typical operation is to extract simulated heads at point locations to
# compare them with measurements. In this example, we select the heads at
# two points:

x = [240_000.0, 244_000.0]
y = [560_000.0, 562_000.0]
selection = imod.select.points_values(hds, x=x, y=y)

# %%
# The result can be converted into a pandas dataframe for timeseries analysis,
# or written to a variety of tabular file formats.

dataframe = selection.to_dataframe().reset_index()
dataframe = dataframe.rename(columns={"index": "id"})
dataframe


"""
Partitioning a regional model
=============================

This example shows how a model can be partitioned into submodels. This will
allow parallelization when solving the model. The example used is the Hondsrug
model. It is partitioned into 3 rectangular parts. In the example we first run
the original, unpartitioned model. Then we partition the model and run the
resulting simulation. Finally we merge the head output of the submodels into a
head array for the whole grid. we print both the heads obtained without
partitioning, and the merged heads of the partitioned simulation, ' for
comparison.
"""

# %% Import packages
import matplotlib.pyplot as plt

import imod
from imod.mf6.multimodel.partition_generator import get_label_array

# %%
# Obtain the simulation, write it, run it, and plot some heads.
# There is a separate example contained in
# :doc:`hondsrug </examples/mf6/hondsrug>`
# that you should look at if you are interested in the model building
tmpdir = imod.util.temporary_directory()

gwf_simulation = imod.data.hondsrug_simulation(tmpdir / "hondsrug_saved")

# %%
# Write the model and run it (before partitioning, so we can compare if the
# results are similar).
original_modeldir = tmpdir / "original"

gwf_simulation.write(original_modeldir)
gwf_simulation.run()

# %%
# Plot the simulation results of the unpartitioned model.
hds_original = gwf_simulation.open_head()

fig, ax = plt.subplots()

hds_original.sel(layer=3).isel(time=6).plot(ax=ax)
ax.set_title("hondsrug original ")
# %%
# Now we partition the Hondsrug model
idomain = gwf_simulation["GWF"].domain
number_partitions = 16
submodel_labels = get_label_array(gwf_simulation, number_partitions)

# %%
# plot the partitioning array. It shows how the model will be partitioned.
fig, ax = plt.subplots()
submodel_labels.plot(ax=ax)
ax.set_title("hondsrug partitioning geometry")

split_simulation = gwf_simulation.split(submodel_labels)

# %%
# Now we  write and run the partitioned model
split_modeldir = tmpdir / "split"

split_simulation.write(split_modeldir)
split_simulation.run()


# %%
# Load and plot the simulation results. Also plot the differences with the original model
hds_split = split_simulation.open_head()["head"]
fig, ax = plt.subplots()
hds_split.sel(layer=3).isel(time=6).plot(ax=ax)
ax.set_title("hondsrug partitioned ")

diff = hds_split - hds_original
diff_for_plot = diff.max(dim=["time", "layer"])
fig, ax = plt.subplots()
diff_for_plot.plot(ax=ax)
ax.set_title("hondsrug diff ")

# %%


"""
Lake package example
====================

This is a synthetic example (using invented, not necesarily physical data) of how to use the
lake package api to generate models with lakes.

In overview, we'll set the following steps:

    * Create a structured grid for a rectangular geometry.
    * Create a constant head boundary
    * Create packages for  initial conditions, output control, storage, and node property flow
    * Create a lake package with a time-dependent rainfall
    * Write to modflow6 files.
    * Run the model.
    * Open the results back into DataArrays.
    * Visualize the results.

"""

# %%
# We'll start with the usual imports. As this is an simple (synthetic)
# structured model, we can make due with few packages.

import numpy as np
import xarray as xr

import imod
import imod.mf6.lak as lak

nlay = 3
nrow = 15
ncol = 15
shape = (nlay, nrow, ncol)

dx = 5000.0
dy = -5000.0
xmin = 0.0
xmax = dx * ncol
ymin = 0.0
ymax = abs(dy) * nrow
dims = ("layer", "y", "x")

layer = np.array([1, 2, 3])
y = np.arange(ymax, ymin, dy) + 0.5 * dy
x = np.arange(xmin, xmax, dx) + 0.5 * dx
coords = {"layer": layer, "y": y, "x": x}

idomain = xr.DataArray(np.ones(shape, dtype=int), coords=coords, dims=dims)

# %% Define lake location
lake_layer = 1
lake_x = x[4:7]
lake_y = y[4:7]
is_lake = xr.full_like(idomain, fill_value=False, dtype=bool)
is_lake.loc[{"layer": lake_layer, "x": lake_x, "y": lake_y}] = True
is_lake.sel(layer=1).plot.imshow()
# %% Specify lake data

VERTICAL = 1
connectionType = xr.where(is_lake, VERTICAL, np.nan)
bed_leak = xr.where(is_lake, 0.2, np.nan)
top_elevation = xr.where(is_lake, 0.4, np.nan)
bot_elevation = xr.where(is_lake, 0.1, np.nan)
connection_length = xr.where(is_lake, 0.5, np.nan)
connection_width = xr.where(is_lake, 0.6, np.nan)

times_rainfall = [
    np.datetime64("2000-01-01"),
    np.datetime64("2000-03-01"),
    np.datetime64("2000-05-01"),
]

rainfall = xr.DataArray(
    np.full((len(times_rainfall)), 0.001),
    coords={"time": times_rainfall},
    dims=["time"],
)


lake = lak.LakeData(
    10.0,
    "Nieuwkoopse_plas",
    connectionType,
    bed_leak,
    top_elevation,
    bot_elevation,
    connection_length,
    connection_width,
    None,
    None,
    rainfall,
    None,
    None,
    None,
    None,
    None,
)


# %%
# Create grid coordinates
# -----------------------
#
# The first steps consist of setting up the grid -- first the number of layer,
# rows, and columns. Cell sizes are constant throughout the model.


# %%
# We'll create a new directory in which we will write and run the model.
modeldir = imod.util.temporary_directory()

# %%
# Create DataArrays
# -----------------
#
# Now that we have the grid coordinates setup, we can start defining model
# parameters. The model is characterized by:
#
# * a constant head boundary on the left
# * a single drain in the center left of the model
# * uniform recharge on the top layer

bottom = xr.DataArray([-200.0, -300.0, -450.0], {"layer": layer}, ("layer",))

# Constant head
constant_head = xr.full_like(idomain, np.nan, dtype=float).sel(layer=[1, 2])
constant_head[..., 0] = 0.0

# Node properties
icelltype = xr.DataArray([1, 0, 0], {"layer": layer}, ("layer",))
k = xr.DataArray([1.0e-3, 1.0e-4, 2.0e-4], {"layer": layer}, ("layer",))
k33 = xr.DataArray([2.0e-8, 2.0e-8, 2.0e-8], {"layer": layer}, ("layer",))

gwf_model = imod.mf6.GroundwaterFlowModel()
gwf_model["dis"] = imod.mf6.StructuredDiscretization(
    top=200.0, bottom=bottom, idomain=idomain
)
gwf_model["chd"] = imod.mf6.ConstantHead(
    constant_head, print_input=True, print_flows=True, save_flows=True
)
gwf_model["ic"] = imod.mf6.InitialConditions(head=0.0)
gwf_model["npf"] = imod.mf6.NodePropertyFlow(
    icelltype=icelltype,
    k=k,
    k33=k33,
    variable_vertical_conductance=True,
    dewatered=True,
    perched=True,
    save_flows=True,
)
gwf_model["sto"] = imod.mf6.SpecificStorage(
    specific_storage=1.0e-5,
    specific_yield=0.15,
    transient=True,
    convertible=0,
)
gwf_model["oc"] = imod.mf6.OutputControl(save_head="all", save_budget="all")

# Attach it to a simulation
simulation = imod.mf6.Modflow6Simulation("ex01-twri")
simulation["GWF_1"] = gwf_model
# Define solver settings
simulation["solver"] = imod.mf6.Solution(
    modelnames=["GWF_1"],
    print_option="summary",
    outer_dvclose=1.0e-4,
    outer_maximum=500,
    under_relaxation=None,
    inner_dvclose=1.0e-4,
    inner_rclose=0.001,
    inner_maximum=100,
    linear_acceleration="cg",
    scaling_method=None,
    reordering_method=None,
    relaxation_factor=0.97,
)

gwf_model["lake"] = lak.Lake.from_lakes_and_outlets(
    [lake],
    print_input=True,
    print_stage=True,
    print_flows=True,
    save_flows=True,
    stagefile=modeldir / "GWF_1/stagefile.lak",
    budgetcsvfile=modeldir / "GWF_1/budgetcsvfile.lak",
    package_convergence_filename=modeldir / "GWF_1/convergence.lak",
)

# Collect time discretization
simulation.create_time_discretization(
    additional_times=["2000-01-01", "2000-01-02", "2000-01-03", "2013-06-04"]
)


simulation.write(modeldir)

# %%
# Run the model
# -------------
#
# .. note::
#
#   The following lines assume the ``mf6`` executable is available on your PATH.
#   :ref:`The Modflow 6 examples introduction <mf6-introduction>` shortly
#   describes how to add it to yours.

simulation.run()
# %%
# Open the results
# ----------------
#
# We'll open the heads (.hds) file.

head = imod.mf6.open_hds(
    modeldir / "GWF_1/GWF_1.hds",
    modeldir / "GWF_1/dis.dis.grb",
)

# %%
# Visualize the results
# ---------------------

head.isel(layer=0, time=4).plot.contourf()


"""
Transport 2d example
====================

The simulation shown here comes from the  1999 MT3DMS report, p 138:
Two-Dimensional Transport in a Uniform Flow of solute injected
continuously from a point source in a steady-state uniform flow field.

In this example, we build up the model, and the we run the model as is.
Next, we split the model in 4 partitions and run that as well.
Finally, we show that the difference in outcome for the partitioned and unpartitioned models
is small.

MT3DMS: A Modular Three-Dimensional
Multispecies Transport Model for Simulation
of Advection, Dispersion, and Chemical
Reactions of Contaminants in Groundwater
Systems; Documentation and User's Guide
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

import imod
from imod.mf6.multimodel.partition_generator import get_label_array
from imod.typing.grid import nan_like, zeros_like

# %%
# Set some grid dimensions
nlay = 1  # Number of layers
nrow = 31  # Number of rows
ncol = 46  # Number of columns
delr = 10.0  # Column width ($m$)
delc = 10.0  # Row width ($m$)
delz = 10.0  # Layer thickness ($m$)
shape = (nlay, nrow, ncol)
top = 10.0
dims = ("layer", "y", "x")

# %%
# Construct the "idomain" array, and then the discretization package which represents the model grid.
y = np.arange(delr * nrow, 0, -delr)
x = np.arange(0, delc * ncol, delc)
coords = {"layer": [1], "y": y, "x": x, "dx": delc, "dy": -delr}
idomain = xr.DataArray(np.ones(shape, dtype=int), coords=coords, dims=dims)

bottom = xr.DataArray([0.0], {"layer": [1]}, ("layer",))
gwf_model = imod.mf6.GroundwaterFlowModel()
gwf_model["dis"] = imod.mf6.StructuredDiscretization(
    top=10.0, bottom=bottom, idomain=idomain
)


# %%
# Construct the other flow packages. Flow is steady-state in this simulation,
# meaning specific storage is set to zero.
#
gwf_model["sto"] = imod.mf6.SpecificStorage(
    specific_storage=0.0,
    specific_yield=0.0,
    transient=False,
    convertible=0,
)
gwf_model["npf"] = imod.mf6.NodePropertyFlow(
    icelltype=zeros_like(idomain), k=1.0, save_flows=True, save_specific_discharge=True
)
gwf_model["oc"] = imod.mf6.OutputControl(save_head="all", save_budget="all")
gwf_model["ic"] = imod.mf6.InitialConditions(start=10.0)

# %%
# Set up the boundary conditions. We have: 2 constant head boundaries at
# the left and right, chosen so that the velocity is 1/3 m/day
# and a well that injects 1 m3 per day, with a concentration of 1000
Lx = 460
v = 1.0 / 3.0
prsity = 0.3
q = v * prsity
h1 = q * Lx
chd_field = nan_like(idomain)
chd_field.values[0, :, 0] = h1
chd_field.values[0, :, -1] = 0.0
chd_concentration = nan_like(idomain)
chd_concentration.values[0, :, 0] = 0.0
chd_concentration.values[0, :, -1] = 0.0
chd_concentration = chd_concentration.expand_dims(species=["Au"])


gwf_model["chd"] = imod.mf6.ConstantHead(
    chd_field,
    concentration=chd_concentration,
    print_input=True,
    print_flows=True,
    save_flows=True,
)
injection_concentration = xr.DataArray(
    [[1000.0]],
    coords={
        "species": ["Au"],
        "index": [0],
    },
    dims=("species", "index"),
)
gwf_model["wel"] = imod.mf6.Well(
    x=[150.0],
    y=[150.0],
    screen_top=[10.0],
    screen_bottom=[0.0],
    rate=[1.0],
    concentration=injection_concentration,
    concentration_boundary_type="aux",
)


# %%
# Now construct the transport simulation. The flow boundaries
# already have inflow concentration data associated, so the transport
# boundaries can be imported using the ssm package, and the rest of the
# transport model definition is straightforward.
tpt_model = imod.mf6.GroundwaterTransportModel()
tpt_model["ssm"] = imod.mf6.SourceSinkMixing.from_flow_model(
    gwf_model, species="Au", save_flows=True
)
tpt_model["adv"] = imod.mf6.AdvectionUpstream()
tpt_model["dsp"] = imod.mf6.Dispersion(
    diffusion_coefficient=0.0,
    longitudinal_horizontal=10.0,
    transversal_horizontal1=3.0,
    xt3d_off=False,
    xt3d_rhs=False,
)
tpt_model["mst"] = imod.mf6.MobileStorageTransfer(
    porosity=0.3,
)

tpt_model["ic"] = imod.mf6.InitialConditions(start=0.0)
tpt_model["oc"] = imod.mf6.OutputControl(save_concentration="all", save_budget="last")
tpt_model["dis"] = gwf_model["dis"]

# %%
# Create a simulation and add the flow and transport models to it.
# Then define some ims packages: 1 for every type of model.
# Finally create 365 time steps of 1 day each.
simulation = imod.mf6.Modflow6Simulation("ex01-twri")
simulation["GWF_1"] = gwf_model
simulation["TPT_1"] = tpt_model


simulation["flow_solver"] = imod.mf6.Solution(
    modelnames=["GWF_1"],
    print_option="summary",
    outer_dvclose=1.0e-4,
    outer_maximum=500,
    under_relaxation=None,
    inner_dvclose=1.0e-4,
    inner_rclose=0.001,
    inner_maximum=100,
    linear_acceleration="cg",
    scaling_method=None,
    reordering_method=None,
    relaxation_factor=0.97,
)
simulation["transport_solver"] = imod.mf6.Solution(
    modelnames=["TPT_1"],
    print_option="summary",
    outer_dvclose=1.0e-4,
    outer_maximum=500,
    under_relaxation=None,
    inner_dvclose=1.0e-4,
    inner_rclose=0.001,
    inner_maximum=100,
    linear_acceleration="bicgstab",
    scaling_method=None,
    reordering_method=None,
    relaxation_factor=0.97,
)
# Collect time discretization

duration = pd.to_timedelta("365d")
start = pd.to_datetime("2002-01-01")
simulation.create_time_discretization(additional_times=[start, start + duration])
simulation["time_discretization"]["n_timesteps"] = 365
modeldir = imod.util.temporary_directory()
simulation.write(modeldir, binary=False)

# %%
# To split the model in 4 partitions, we must create a label array.
# We use the utility function  ``get_label_array`` for that.

label_array = get_label_array(simulation, 4)
fig, ax = plt.subplots()
label_array.plot(ax=ax)

split_simulation = simulation.split(label_array)
# %%
# Run the unsplit model and load the simulation results.
simulation.run()
conc = simulation.open_concentration()

# %%
# Run the split model and load the simulation results.
split_modeldir = modeldir / "split"
split_simulation.write(split_modeldir, binary=False)
split_simulation.run()
split_conc = split_simulation.open_concentration()["concentration"]
fig, ax = plt.subplots()
split_conc.isel(layer=0, time=364).plot.contourf(ax=ax, levels=[0.1, 1, 10])

# %%
# Compute the difference between the split and unsplit simulation results for transport at the
# end of the simulation, and print them
diff = abs(conc - split_conc)
fig, ax = plt.subplots()
diff.isel(layer=0, time=364).plot.contourf(ax=ax)
