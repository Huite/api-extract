imod.visualize.cross_section
============================
Wraps matplotlib.pcolormesh to draw cross-sections, drawing cell boundaries
accurately. Aquitards can be plotted on top of the cross-section, by providing
a DataArray with the aquitard location for `aquitards`.

Parameters
----------
da : xr.DataArray
    Two dimensional DataArray containing data of the cross section. One
    dimension must be "layer", and the second dimension will be used as the
    x-axis for the cross-section.

    Coordinates "top" and "bottom" must be present, and must have at least the
    "layer" dimension (voxels) or both the "layer" and x-coordinate dimension.

    *Use imod.select.cross_section_line() or cross_section_linestring() to obtain
    the required DataArray.*
colors : list of str, or list of RGB tuples
    Matplotlib acceptable list of colors. Length N.
    Accepts both tuples of (R, G, B) and hexidecimal (e.g. "#7ec0ee").

    Looking for good colormaps? Try: http://colorbrewer2.org/
    Choose a colormap, and use the HEX JS array.
levels : listlike of floats or integers
    Boundaries between the legend colors/classes. Length: N - 1.
layers : boolean, optional
    Whether to draw lines separating the layers.
aquitards : xr.DataArray, optional
    Datarray containing data on location of aquitard layers.
kwargs_pcolormesh : dict
    Other optional keyword arguments for matplotlib.pcolormesh.
kwargs_colorbar : dict
    If optional key ``plot_colorbar`` is set to False, no colorbar is drawn. Defaults to True.
    Optional keyword argument ``whiten_triangles`` whitens respective colorbar triangle if
    data is not larger/smaller than legend_levels-range. Defaults to True.
    Other arguments are forwarded to fig.colorbar()
kwargs_aquitards: dict
    These arguments are forwarded to matplotlib.fill_between to draw the
    aquitards.
return_cmap_norm : boolean, optional
    Return the cmap and norm of the plot, default False
fig : matplotlib Figure instance, optional
    Figure to write plot to. If not supplied, a Figure instance is created
ax : matplotlib Axes instance, optional
    Axes to write plot to. If not supplied, an Axes instance is created

Returns
-------
fig : matplotlib.figure
ax : matplotlig.ax
if return_cmap_norm == True:
cmap : matplotlib.colors.ListedColormap
norm : matplotlib.colors.BoundaryNorm

Examples
--------

Basic cross section:

>>> imod.visualize.cross_section(da, colors, levels)

Aquitards can be styled in multiple ways. For a transparent grey overlay
(the default):

>>> kwargs_aquitards = {"alpha": 0.5, "facecolor": "grey"}
>>> imod.visualize.cross_section(da, colors, levels, aquitards=aquitards, kwargs_aquitards)

For a hatched overlay:

>>> kwargs_aquitards = {"hatch": "/", "edgecolor": "k"}
>>> imod.visualize.cross_section(da, colors, levels, aquitards=aquitards, kwargs_aquitards)

imod.visualize.plot_map
=======================
Plot raster on a map with optional vector overlays and basemap.

Parameters
----------
raster : xr.DataArray
    2D grid to plot.
colors : list of str, list of RGBA/RGBA tuples, colormap name (str), or
    LinearSegmentedColormap.

    If list, it should be a Matplotlib acceptable list of colors. Length N.
    Accepts both tuples of (R, G, B) and hexidecimal (e.g. `#7ec0ee`).

    If str, use an existing Matplotlib colormap. This function will
    autmatically add distinctive colors for pixels lower or high than the
    given min respectivly max level.

    If LinearSegmentedColormap, you can use something like
    `matplotlib.cm.get_cmap('jet')` as input. This function will not alter
    the colormap, so add under- and over-colors yourself.

    Looking for good colormaps? Try: http://colorbrewer2.org/ Choose a
    colormap, and use the HEX JS array.
levels : listlike of floats or integers
    Boundaries between the legend colors/classes. Length: N - 1.
overlays : list of dicts, optional
    Dicts contain geodataframe (key is "gdf"), and the keyword arguments for
    plotting the geodataframe.
basemap : bool or contextily._providers.TileProvider, optional
    When `True` or a `contextily._providers.TileProvider` object: plot a
    basemap as a background for the plot and make the raster translucent. If
    `basemap=True`, then `CartoDB.Positron` is used as the default provider.
    If not set explicitly through kwargs_basemap, plot_map() will try and
    infer the crs from the raster or overlays, or fall back to EPSG:28992
    (Amersfoort/RDnew).

    *Requires contextily*

kwargs_raster : dict of keyword arguments, optional
    These arguments are forwarded to ax.imshow()
kwargs_colorbar : dict of keyword arguments, optional
    These arguments are forwarded to fig.colorbar(). The key label can be
    used to label the colorbar. Key whiten_triangles can be set to False to
    alter the default behavior of coloring the min / max triangles of the
    colorbar white if the value is not present in the map.
kwargs_basemap : dict of keyword arguments, optional
    Except for "alpha", these arguments are forwarded to
    contextily.add_basemap(). Parameter "alpha" controls the transparency of
    raster.
figsize : tuple of two floats or integers, optional
    This is used in plt.subplots(figsize)
return_cbar : boolean, optional
    Return the matplotlib.Colorbar instance. Defaults to False.
fig : matplotlib.figure, optional
    If provided, figure to which to add the map
ax : matplot.ax, optional
    If provided, axis to which to add the map

Returns
-------
fig : matplotlib.figure ax : matplotlib.ax if return_cbar == True: cbar :
matplotlib.Colorbar

Examples
--------
Plot with an overlay:

>>> overlays = [{"gdf": geodataframe, "edgecolor": "black", "facecolor": "None"}]
>>> imod.visualize.plot_map(raster, colors, levels, overlays)

Label the colorbar:

>>> imod.visualize.plot_map(raster, colors, levels, kwargs_colorbar={"label":"Head aquifer (m)"})

Plot with a basemap:

>>> import contextily as ctx
>>> src = ctx.providers.Stamen.TonerLite
>>> imod.visualize.plot_map(raster, colors, levels, basemap=src, kwargs_basemap={"alpha":0.6})

imod.visualize.imshow_topview
=============================
Automatically colors by quantile.

Dumps PNGs into directory of choice.

imod.visualize.read_imod_legend
===============================
Parameters
----------
path : str
    Path to iMOD .leg file.

Returns
-------
colors : List of hex colors of length N.
levels : List of floats of length N-1. These are the boundaries between
    the legend colors/classes.

imod.visualize.quiver
=====================
Wraps matplotlib.quiver to draw quiver plots. Function can be used to draw
flow quivers on top of a cross-section.

Parameters
----------
u : xr.DataArray
    Two dimensional DataArray containing u component of quivers. One
    dimension must be "layer", and the second dimension will be used as the
    x-axis for the cross-section.

    Coordinates "top" and "bottom" must be present, and must have at least the
    "layer" dimension (voxels) or both the "layer" and x-coordinate dimension.

    *Use imod.evaluate.quiver_line() or quiver_linestring() to obtain
    the required DataArray.*
v : xr.DataArray
    Two dimensional DataArray containing v component of quivers. One
    dimension must be "layer", and the second dimension will be used as the
    x-axis for the cross-section.

    Coordinates "top" and "bottom" must be present, and must have at least the
    "layer" dimension (voxels) or both the "layer" and x-coordinate dimension.

    *Use imod.evaluate.quiver_line() or quiver_linestring() to obtain
    the required DataArray.*
ax : matplotlib Axes instance
    Axes to write plot to.
kwargs_quiver : dict
    Other optional keyword arguments for matplotlib.quiver.

Returns
-------
matplotlib.quiver.Quiver
    The drawn quivers.

Examples
--------

First: apply evaluate.quiver_line to get the u and v components of the quivers
from a three dimensional flow field. Assign top and bottom coordinates if these are
not already present in the flow field data arrays.

>>> u, v = imod.evaluate.quiver_line(right, front, lower, start, end)
>>> u.assign_coords(top=top, bottom=bottom)
>>> v.assign_coords(top=top, bottom=bottom)

The quivers can then be plotted over a cross section created by imod.visualize.cross_section():

>>> imod.visualize.quiver(u, v, ax)

Quivers can easily overwhelm your plot, so it is a good idea to 'thin out' some of the quivers:

>>> # Only plot quivers at every 5th cell and every 3rd layer
>>> thin = {"s": slice(0, None, 5), "layer": slice(0, None, 3)}
>>> imod.visualize.quiver(u.isel(**thin), v.isel(**thin), ax)

imod.visualize.streamfunction
=============================
Wraps matplotlib.contour to draw stream lines. Function can be used to draw
stream lines on top of a cross-section.

Parameters
----------
da : xr.DataArray
    Two dimensional DataArray containing data of the cross section. One
    dimension must be "layer", and the second dimension will be used as the
    x-axis for the cross-section.

    Coordinates "top" and "bottom" must be present, and must have at least the
    "layer" dimension (voxels) or both the "layer" and x-coordinate dimension.

    *Use imod.evaluate.streamfunction_line() or streamfunction_linestring() to obtain
    the required DataArray.*
ax : matplotlib Axes instance
    Axes to write plot to.
n_streamlines : int or array_like
    Determines the number and positions of the contour lines / regions.

    If an int n, use n data intervals; i.e. draw n+1 contour lines. The level heights are automatically chosen.

    If array-like, draw contour lines at the specified levels. The values must be in increasing order.
kwargs_contour : dict
    Other optional keyword arguments for matplotlib.contour.

Returns
-------
matplotlib.contour.QuadContourSet
    The drawn contour lines.

imod.visualize.waterbalance_barchart
====================================
Parameters
----------
df : pandas.DataFrame
    The dataframe containing the water balance data.
inflows : listlike of str
outflows : listlike of str
datecolumn : str, optional
format : str, optional,
ax : matplotlib.Axes, optional
unit : str, optional
colors : listlike of strings or tuples

Returns
-------
ax : matplotlib.Axes

Examples
--------

>>> fig, ax = plt.subplots()
>>> imod.visualize.waterbalance_barchart(
>>>    ax=ax,
>>>    df=df,
>>>    inflows=["Rainfall", "River upstream"],
>>>    outflows=["Evapotranspiration", "Discharge to Sea"],
>>>    datecolumn="Time",
>>>    format="%Y-%m-%d",
>>>    unit="m3/d",
>>>    colors=["#ca0020", "#f4a582", "#92c5de", "#0571b0"],
>>>    )
>>> fig.savefig("Waterbalance.png", dpi=300, bbox_inches="tight")

imod.visualize.grid_3d
======================
Constructs a 3D PyVista representation of a DataArray.
DataArrays should be two-dimensional or three-dimensional:

* 2D: dimensions should be ``{"y", "x"}``. E.g. a DEM.
* 3D: dimensions should be ``{"z", "y", "x"}``, for a voxel model.
* 3D: dimensions should be ``{"layer", "y", "x"}``, with coordinates
    ``"top"({"layer", "y", "x"})`` and ``"bottom"({"layer", "y", "x"})``.

Parameters
----------
da : xr.DataArray
vertical_exaggeration : float, default 30.0
exterior_only : bool, default True
    Whether or not to only draw the exterior. Greatly speeds up rendering,
    but it means that pyvista slices and filters produce "hollowed out"
    results.
exterior_depth : int, default 1
    How many cells to consider as exterior. In case of large jumps, holes
    can occur. By settings this argument to a higher value, more of the
    inner cells will be rendered, reducing the chances of gaps occurring.
return_index : bool, default False

Returns
-------
pyvista.UnstructuredGrid

Examples
--------

>>> grid = imod.visualize.grid_3d(da)

To plot the grid, call the ``.plot()`` method.

>>> grid.plot()

Use ``.assign_coords`` to assign tops and bottoms to layer models:

>>> top = imod.idf.open("top*.idf")
>>> bottom = imod.idf.open("bot*.idf")
>>> kd = imod.idf.open("kd*.idf")
>>> kd = kd.assign_coords(top=(("layer", "y", "x"), top))
>>> kd = kd.assign_coords(bottom=(("layer", "y", "x"), bottom))
>>> grid = imod.visualize.grid_3d(kd)
>>> grid.plot()

Refer to the PyVista documentation on how to customize plots:
https://docs.pyvista.org/index.html

imod.visualize.line_3d
======================
Returns the exterior line of a shapely polygon.

Parameters
----------
polygon : shapely.geometry.Polygon
z : float or xr.DataArray
    z-coordinate to assign to line. If DataArray, assigns z-coordinate
    based on xy locations in DataArray.

Returns
-------
pyvista.PolyData

imod.visualize.GridAnimation3D
==============================
Class to easily setup 3D animations for transient data. Use the
``imod.visualize.StaticGridAnimation3D`` when the location of the displayed
cells is constant over time: it will render much faster.

You can iteratively add or change settings to the plotter, until you're
satisfied. Call the ``.peek()`` method to take a look. When satisfied, call
``.output()`` to write to a file.


Parameters
----------
da : xr.DataArray
    The dataarray with transient data. Must contain a "time" dimension.
vertical_exaggeration : float, defaults to 30.0 mesh_kwargs : dict
    keyword arguments that are forwarded to the pyvista mesh representing
    "da". If "stitle" is given as one of the arguments, the special keyword
    "timestamp" can be used to render the plotted time as part of the title.
    See example. For a full list of kwargs supported, see the
    `plotter.add_mesh
    <https://docs.pyvista.org/api/plotting/_autosummary/pyvista.Plotter.add_mesh.html#pyvista.Plotter.add_mesh>`_
    method documentation.
plotter_kwargs : dict
    keyword arguments that are forwarded to the pyvista plotter. For a full
    list of of kwargs supported, see the `Plotter constructor
    <https://docs.pyvista.org/api/plotting/_autosummary/pyvista.Plotter.html>`_
    documention.

Examples
--------

Initialize the animation:

>>> animation = imod.visualize.GridAnimation3D(concentration, mesh_kwargs=dict(cmap="jet"))

Check what it looks like (if a window pops up: press "q" instead of the X to
return):

>>> animation.peek()

Change the camera position, add bounding box, and check the result:

>>> animation.plotter.camera_position = (2, 1, 0.5)
>>> animation.plotter.add_bounding_box()
>>> animation.peek()

When it looks good, write to a file:

>>> animation.write("example.mp4")

If you've made some changes that don't look good, call ``.reset()`` to start
over:

>>> animation.reset()

Note that ``.reset()`` is automatically called when the animation has
finished writing.

You can use "stitle" in mesh_kwargs in conjunction with the "timestamp"
keyword to print a formatted timestamp in the animation:

>>> animation = imod.visualize.GridAnimation3D(concentration, mesh_kwargs=dict(stitle="Concentration on {timestamp:%Y-%m-%d}"))

imod.visualize.GridAnimation3D Class Members
============================================
   * imod.visualize.GridAnimation3D.peek
   * imod.visualize.GridAnimation3D.reset
   * imod.visualize.GridAnimation3D.write

imod.visualize.GridAnimation3D.peek
===================================
Display the current state of the animation plotter.

imod.visualize.GridAnimation3D.reset
====================================
Reset the plotter to its base state.

imod.visualize.GridAnimation3D.write
====================================
Write the animation to a video or gif.

Resets the plotter when finished animating.

Parameters
----------
filename : str, pathlib.Path
    Filename to write the video to. Should be an .mp4 or .gif.
framerate : int, optional
    Frames per second. Not honoured for gif.

imod.visualize.StaticGridAnimation3D
====================================
Class to easily setup 3D animations for transient data;
Should only be used when the location of the displayed cells is constant
over time. It will render much faster than ``imod.visualize.GridAnimation3D``.

Refer to examples of ``imod.visualize.GridAnimation3D``.

imod.visualize.StaticGridAnimation3D Class Members
==================================================
   * imod.visualize.StaticGridAnimation3D.peek
   * imod.visualize.StaticGridAnimation3D.reset
   * imod.visualize.StaticGridAnimation3D.write

