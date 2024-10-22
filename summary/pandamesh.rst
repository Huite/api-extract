pandamesh.Preprocessor
======================
Utilities to clean-up geometries before meshing.

Each processing method return a new Preprocessor instance. This enables
method chaining.

Note: many methods require exhaustive checking of geometries. Processing
may be slow if your geometry contains thousands of geometries or if it has
excessive detail.

Parameters
----------
geometry: np.ndarray of shapely geometries of length N.
    Should be points, linestrings, or polygons.
values: np.ndarray of length N, optional
    Values associated with the geometries (e.g. cell sizes for meshing).
grid_size: float, optional
    Forwarded to shapely operations. If grid_size is nonzero, input
    coordinates will be snapped to a precision grid of that size and
    resulting coordinates will be snapped to that same grid for shapely
    operations such as ``intersection, difference, union``, etc. If 0, the
    operations will use double precision coordinates.

Examples
--------
This class allows for method chaining to flexibly combine pre-processing
operations:

>>> processed = (
>>>     Preprocessor(geometry=gdf.geometry)
>>>     .interpolate_lines_to_points(distance=5.0)
>>>     .snap_points(distance=2.0)
>>>     .unify_polygons()
>>>     .merge_polygons()
>>>     .clip_points(distance=5.0)
>>>     .to_geodataframe()
>>> )

An intermediate result can be checked and inspect by calling ``.to_geodataframe()``
at any time.

>>> check = Preprocessor(geometry=gdf.geometry).interpolate_lines_to_points(5.0).to_geodataframe()

pandamesh.Preprocessor Class Members
====================================
   * pandamesh.Preprocessor.clip_lines
   * pandamesh.Preprocessor.clip_points
   * pandamesh.Preprocessor.interpolate_lines_to_points
   * pandamesh.Preprocessor.merge_polygons
   * pandamesh.Preprocessor.snap_points
   * pandamesh.Preprocessor.to_geodataframe
   * pandamesh.Preprocessor.unify_lines
   * pandamesh.Preprocessor.unify_polygons

pandamesh.Preprocessor.clip_lines
=================================
Remove line segments that are outside or that are near polygon
segments.

Returns
-------
processed: pandamesh.Preprocessor

pandamesh.Preprocessor.clip_points
==================================
Remove points that are outside of a polygon or near line or
polygon segments.

Parameters
----------
distance: float, optional
    Minimum distance to a line or polygon segment.

Returns
-------
processed: pandamesh.Preprocessor

pandamesh.Preprocessor.interpolate_lines_to_points
==================================================
Convert lines into points.

This method adds vertices if needed, but it does not discard vertices
that are located close together: ``distance`` is the maximum distance,
not a minimum.

Parameters
----------
distance: float, optional
    Additional vertices will be added so that all line segments are no
    longer than this value. Must be greater than 0.
values_as_distance: bool, optional
    If true, ignores the value of ``distance`` and uses ``values``
    provided during initialization as the distance instead. Errors if
    no values have been provided.

Returns
-------
processed: pandamesh.Preprocessor

pandamesh.Preprocessor.merge_polygons
=====================================
Merge polygons with the same value for indexer (the same value if
``values`` was provided at initialization).

Returns
-------
processed: pandamesh.Preprocessor

pandamesh.Preprocessor.snap_points
==================================
Snap points together that are within tolerance of each other.

Will use Numba to accelerate the snapping if it is installed. This may
be significantly faster in case of snapping a large (>10 000) number of points.

Returns
-------
processed: pandamesh.Preprocessor

pandamesh.Preprocessor.to_geodataframe
======================================
Return the processed geometries as a ``geopandas.GeoDataFrame``.

Returns
-------
gdf: geopandas.GeoDataFrame
    Contains columns geometry and indexer. The indexer column can be
    used to (re-)associate the geometry with the original values.
    If ``values`` were provided at initialization, a values column will
    be added to the geodataframe.

pandamesh.Preprocessor.unify_lines
==================================
Ensure intersections between lines are present.

Returns
-------
processed: pandamesh.Preprocessor

pandamesh.Preprocessor.unify_polygons
=====================================
Resolve polygon overlaps and intersections.

In overview, this method takes the following steps:

1. collect all linear rings (exterior and interior boundaries), as well
   as the linestrings.
2. create a unary union of all the linework. This ensures intersections
   between lines are represented by a point on the lines.
3. polygonize the union linework. This creates a polygon for each ring that
   is encountered, including holes.
4. collect sampling points for the newly created polygons. Use these to
   locate in which polygon (or which hole!) the newly created polygon is
   located.
5. In case of overlapping polygons, the sampling point may be present
   in more than one of the original polygons. We choose the one with
   the lowest indexer value or the smallest value in case ``.values``
   was provided at initialization if ``first=True``; the highest
   indexer value or largest value is taken for ``first=False``.
6. re-associate with the original indexer and discard hole polygons.

Unify polygons may generate many neighboring sub-polygons with the same
indexer value. The can be merged with ``.merge_polygons``.

Parameters
----------
first: bool, optional
    Which value or index to assign in case of polygon overlap.

    In case ``values`` were provided at initialization:

    * ``first=True``: take the smallest value among the overlapping
      polygons.
    * ``first=False``: take the largest value among the
      overlapping polygons.

    If not provided:

    * ``first=True``: take the first overlapping polygon of the input
      geometry.
    * ``first=False``: take the last overlapping polygon of the input
      geometry.

    See the examples.

Returns
-------
processed: pandamesh.Preprocessor

Examples
--------
Resolve overlapping polygons, assigning the smallest cell size values
in case of overlap:

>>> processed = (
>>>     pandamesh.Preprocessor(geometry=gdf["geometry"], values=gdf["cellsize"])
>>>     .unify_polygons()
>>>     .merge_polygons(first=True)
>>>     .to_geodataframe()
>>> )

Assign the largest cell size value instead:

>>> processed = (
>>>     pandamesh.Preprocessor(geometry=gdf["geometry"], values=gdf["cellsize"])
>>>     .unify_polygons()
>>>     .merge_polygons(first=False)
>>>     .to_geodataframe()
>>> )

Alternatively, to control the result of the merging with values, we can
sort the ``gdf`` prior to processing.

>>> sorted_gdf = gdf.sort_values(["a", "b"])
>>> processor = pandamesh.Preprocessor(sorted_gdf["geometry"])

Afterwards, the returned indexer can be used to fetch the data associated
with the merged results.

>>> out = processor.to_geodataframe()
>>> processed = geopandas.GeoDataFrame(
>>>     data=sorted_gdf.iloc[out["indexer"]].loc[["a", "b"]],
>>>     geometry=out["geometry"],
>>> )

pandamesh.find_edge_intersections
=================================
Find all unresolved intersections between polygon boundaries, linestring,
and linearring edges.

A "resolved" intersection is one where the intersection of two lines is
represented by a vertex in both lines. Unresolved means: an intersection
which is not represented by an explicit vertex in the geometries.

Parameters
----------
geometry: gpd.GeoSeries
    Points, lines, polygons.

Returns
-------
intersections: gpd.GeoSeries
    Locations (points) of intersections.

pandamesh.find_proximate_perimeter_points
=========================================
Detect vertices in polygon perimeters (exteriors and interiors) that are
very close to each other.

Note that dangling edges can be detected through self-intersection: whether
a geometry is simple or not. However, some slivers will almost form a
dangling edge, where the sliver still have a very small thickness. This may
result in problems during mesh generation, as tiny triangles will be
required locally.

Note that sliver concavities are allowed: the vertex spacing **along** the
perimeter is not necessarily small.

Parameters
----------
geometry : geopandas.Geoseries
    Points, lines, polygons.
minimum_spacing : float, default is 1.0e-3.
    The minimum allowed distance between vertices, or the minimum width of
    slivers.

pandamesh.TriangleMesher
========================
Wrapper for the python bindings to Triangle. This class must be initialized
with a geopandas GeoDataFrame containing at least one polygon, and a column
named ``"cellsize"``.

Optionally, multiple polygons with different cell sizes can be included in
the geodataframe. These can be used to achieve local mesh refinement.

Linestrings and points may also be included. The segments of linestrings
will be directly forced into the triangulation. Points can also be forced
into the triangulation. The cell size values associated with these
geometries willl not be used.

Triangle cannot automatically resolve overlapping polygons, or points
located exactly on segments. During initialization, the geometries of
the geodataframe are checked:

    * Polygons should not have any overlap with each other.
    * Linestrings should not intersect each other, unless the intersection
      vertex is present in both.
    * Every linestring should be fully contained by a single polygon;
      a linestring may not intersect two or more polygons.
    * Linestrings and points should not "touch" / be located on
      polygon borders.
    * Holes in polygons are fully supported, but they must not contain
      any linestrings or points.

If such cases are detected, the initialization will error: use the
:class:`pandamesh.Preprocessor` to clean up geometries beforehand.

For more details on Triangle, see:
https://www.cs.cmu.edu/~quake/triangle.defs.html

Parameters
----------
gdf: gpd.GeoDataFrame
    GeoDataFrame containing the vector geometry. Must contain a "cellsize"
    column.
shift_origin: bool, optional, default is True.
    If True, temporarily shifts the coordinate system origin to the centroid
    of the geometry's bounding box during mesh generation. This helps mitigate
    floating-point precision issues. The resulting mesh vertices are
    automatically translated back to the original coordinate system.
intersecting_edges: str, optional, default is "error"
    String indicating how to report unresolved line segment intersections:

    * "ignore": skip check.
    * "warning": emit a warning.
    * "error": raise a ValueError.

minimum_perimeter_spacing: float, default is 1.0e-3.
    Errors if spacing of vertices on polygon perimeters is less or equal to
    minimum spacing. A distance of 0.0 indicates a dangling edge or a
    repeated vertex. Such features may cause a crash during mesh
    generation.

pandamesh.TriangleMesher Class Members
======================================
   * pandamesh.TriangleMesher.conforming_delaunay
   * pandamesh.TriangleMesher.consistency_check
   * pandamesh.TriangleMesher.delaunay_algorithm
   * pandamesh.TriangleMesher.generate
   * pandamesh.TriangleMesher.generate_geodataframe
   * pandamesh.TriangleMesher.generate_ugrid
   * pandamesh.TriangleMesher.maximum_steiner_points
   * pandamesh.TriangleMesher.minimum_angle
   * pandamesh.TriangleMesher.suppress_exact_arithmetic

pandamesh.TriangleMesher.conforming_delaunay
============================================
Conforming Delaunay: use this switch if you want all triangles in the
mesh to be Delaunay, and not just constrained Delaunay; or if you want
to ensure that all Voronoi vertices lie within the triangulation.

pandamesh.TriangleMesher.consistency_check
==========================================
Check the consistency of the final mesh. Uses exact arithmetic for
checking, even if ``suppress_exact_arithmetic`` is set to ``False``.
Useful if you suspect Triangle is buggy.

pandamesh.TriangleMesher.delaunay_algorithm
===========================================
Sets the Delaunay algorithm. Can be set to one of:
:py:class:`pandamesh.DelaunayAlgorithm`:

.. code::

    DIVIDE_AND_CONQUER = ""
    INCREMENTAL = "i"
    SWEEPLINE = "F"

pandamesh.TriangleMesher.generate
=================================
Generate a mesh of triangles.

Returns
-------
vertices: np.ndarray of floats with shape ``(n_vertex, 2)``
triangles: np.ndarray of integers with shape ``(n_triangle, 3)``

pandamesh.TriangleMesher.generate_geodataframe
==============================================
Generate a mesh and return it as a geopandas GeoDataFrame.

Returns
-------
mesh: geopandas.GeoDataFrame

pandamesh.TriangleMesher.generate_ugrid
=======================================
Generate a mesh and return it as an xugrid Ugrid2d.

Returns
-------
mesh: xugrid.Ugrid2d

pandamesh.TriangleMesher.maximum_steiner_points
===============================================
Specifies the maximum number of added Steiner points

See:
https://www.cs.cmu.edu/~quake/triangle.S.html

pandamesh.TriangleMesher.minimum_angle
======================================
Minimum allowed angle for any triangle in the mesh.

See:
https://www.cs.cmu.edu/~quake/triangle.q.html

pandamesh.TriangleMesher.suppress_exact_arithmetic
==================================================
Suppresses exact arithmetic.

See:
https://www.cs.cmu.edu/~quake/triangle.exact.html

pandamesh.DelaunayAlgorithm
===========================
The type of Delaunay algorithm for Triangle.

pandamesh.DelaunayAlgorithm Class Members
=========================================
   * pandamesh.DelaunayAlgorithm.DIVIDE_AND_CONQUER
   * pandamesh.DelaunayAlgorithm.INCREMENTAL
   * pandamesh.DelaunayAlgorithm.SWEEPLINE

pandamesh.DelaunayAlgorithm.DIVIDE_AND_CONQUER
==============================================
The type of Delaunay algorithm for Triangle.

pandamesh.DelaunayAlgorithm.INCREMENTAL
=======================================
The type of Delaunay algorithm for Triangle.

pandamesh.DelaunayAlgorithm.SWEEPLINE
=====================================
The type of Delaunay algorithm for Triangle.

pandamesh.GmshMesher
====================
Wrapper for the python bindings to Gmsh. This class must be initialized
with a geopandas GeoDataFrame containing at least one polygon, and a column
named ``"cellsize"``.

Optionally, multiple polygons with different cell sizes can be included in
the geodataframe. These can be used to achieve local mesh refinement.

Linestrings and points may also be included. The segments of linestrings
will be directly forced into the triangulation. Points can also be forced
into the triangulation. Unlike Triangle, the cell size values associated
with these geometries **will** be used.

Gmsh cannot automatically resolve overlapping polygons, or points
located exactly on segments. During initialization, the geometries of
the geodataframe are checked:

    * Polygons should not have any overlap with each other.
    * Linestrings should not intersect each other, unless the intersection
      vertex is present in both.
    * Every linestring should be fully contained by a single polygon;
      a linestring may not intersect two or more polygons.
    * Linestrings and points should not "touch" / be located on
      polygon borders.
    * Holes in polygons are fully supported, but they must not contain
      any linestrings or points.

If such cases are detected, the initialization will error: use the
:class:`pandamesh.Preprocessor` to clean up geometries beforehand.

For more details on Gmsh, see:
https://gmsh.info/doc/texinfo/gmsh.html

A helpful index can be found near the bottom:
https://gmsh.info/doc/texinfo/gmsh.html#Syntax-index

.. note::

    This meshers uses the Gmsh Python API, which is global. To avoid a
    situation where multiple GmshMeshers have been iniatilized and are
    mutating each other's (global) variables, the ``.finalize()`` method
    must be called before instantiating a new mesher.

Parameters
----------
gdf: gpd.GeoDataFrame
    GeoDataFrame containing the vector geometry. Must contain a "cellsize"
    column.
shift_origin: bool, optional, default is True.
    If True, temporarily shifts the coordinate system origin to the centroid
    of the geometry's bounding box during mesh generation. This helps mitigate
    floating-point precision issues. The resulting mesh vertices are
    automatically translated back to the original coordinate system.
intersecting_edges: str, optional, default is "error"
    String indicating how to report unresolved line segment intersections:

    * "ignore": skip check.
    * "warning": emit a warning.
    * "error": raise a ValueError.
minimum_perimeter_spacing: float, default is 1.0e-3.
    Errors if spacing of vertices on polygon perimeters is less or equal to
    minimum spacing. A distance of 0.0 indicates a dangling edge or a
    repeated vertex. Such features may cause a crash during mesh
    generation.
read_config_files: bool
    Gmsh initialization option: Read system Gmsh configuration files
    (gmshrc and gmsh-options).
interruptible: bool
    Gmsh initialization option.

pandamesh.GmshMesher Class Members
==================================
   * pandamesh.GmshMesher.add_matheval_distance_field
   * pandamesh.GmshMesher.add_structured_field
   * pandamesh.GmshMesher.add_structured_field_from_dataarray
   * pandamesh.GmshMesher.add_threshold_distance_field
   * pandamesh.GmshMesher.field_combination
   * pandamesh.GmshMesher.fields
   * pandamesh.GmshMesher.finalize_gmsh
   * pandamesh.GmshMesher.general_verbosity
   * pandamesh.GmshMesher.generate
   * pandamesh.GmshMesher.generate_geodataframe
   * pandamesh.GmshMesher.generate_ugrid
   * pandamesh.GmshMesher.mesh_algorithm
   * pandamesh.GmshMesher.mesh_size_extend_from_boundary
   * pandamesh.GmshMesher.mesh_size_from_curvature
   * pandamesh.GmshMesher.mesh_size_from_points
   * pandamesh.GmshMesher.recombine_all
   * pandamesh.GmshMesher.subdivision_algorithm
   * pandamesh.GmshMesher.write

pandamesh.GmshMesher.add_matheval_distance_field
================================================
Add a matheval distance field to the mesher.

The geometry of these fields are not forced into the mesh, but they are
used to specify zones of with cell sizes.

Uses the MathEval functionality in Gmsh, which relies on the SSCILIB
math expression evaluator.

https://gitlab.onelab.info/gmsh/gmsh/-/blob/master/contrib/MathEx/mathex.cpp

https://sscilib.sourceforge.net/

Parameters
----------
gdf: geopandas.GeoDataFrame
    Location of the features to measure distance to. Should contain
    ``spacing`` column to specify the spacing of interpolated vertices
    along linestrings and polygon boundaries, and a ``function`` column
    to specify the function to control cell size as a function of
    distance, ``"distance"`` must be present as an argument in the
    function string. Note that the distance must never evaluate to 0,
    since 0 sized cells are not allowed. See the examples.

Examples
--------
Generate a number of points:

>>> x = np.arange(0.0, 10.0)
>>> y = np.arange(0.0, 10.0)
>>> points = gpd.points_from_xy(x, y)
>>> field = gpd.GeoDataFrame(geometry=points)

Add spacing (dummy value for points) and a function:

>>> field["spacing"] = np.nan
>>> field["function"] = "max(distance^2, 1.0)"

Note that the ``max`` function is used to ensure that the cell size is
never smaller than 1.0

Apply it:

>>> mesher.add_matheval_distance_field(field)

The following mathematical operators are supported:

Basic Operators

- Arithmetic: ``+``, ``-``, ``*``, ``/``, ``%`` (modulo), ``^`` (power)
- Comparison: ``<``, ``>``

Mathematical Functions

- Absolute value: ``abs(x)``
- Square root: ``sqrt(x)``
- Exponential: ``exp(x)``
- Natural logarithm: ``log(x)``
- Base-10 logarithm: ``log10(x)``
- Power: ``pow(x,y)``

Statistical Functions

- Minimum: ``min(x, y, ...)``
- Maximum: ``max(x, y, ...)``
- Sum: ``sum(x, y, ...)``
- Average: ``med(x, y, ...)``

Trigonometric Functions

- Standard: ``sin(x)``, ``cos(x)``, ``tan(x)``
- Inverse: ``asin(x)``, ``acos(x)``, ``atan(x)``
- Hyperbolic: ``sinh(x)``, ``cosh(x)``, ``tanh(x)``

Rounding Functions

- ``floor(x)``, ``ceil(x)``, ``round(x)``, ``trunc(x)``

Constants

- Pi: ``pi``
- Euler's number: ``e``

pandamesh.GmshMesher.add_structured_field
=========================================
Add an equidistant structured field specifying cell sizes for mesh generation.

This method defines a grid of cell sizes that Gmsh will use as control points
for mesh generation. Gmsh interpolates between these points to determine the
desired cell size at any location within the grid.

Parameters
----------
cellsize: FloatArray with shape ``(n_y, n_x``)
    Specifies the cell size on a structured grid. The location of this grid
    is determined by ``xmin, ymin, dx, dy``.
xmin: float
    x-origin.
ymin: float
    y-origin.
dx: float
    Spacing along the x-axis.
dy: float
    Spacing along the y-axis.
outside_value: Union[float, None]
    Value outside of the window ``(xmin, xmax)`` and ``(ymin, ymax)``.
    Default value is None.

pandamesh.GmshMesher.add_structured_field_from_dataarray
========================================================
Add an equidistant structured field as an xarray DataArray specifying
cell sizes for mesh generation.

This method defines a grid of cell sizes that Gmsh will use as control
points for mesh generation. Gmsh interpolates between these points to
determine the desired cell size at any location within the grid.

Parameters
----------
da: xarray.DataArray
    Values are used as cell sizes. Must have dimensions `("y", "x")`.
outside_value: Union[float, None]
    Value outside of the window ``(xmin, xmax)`` and ``(ymin, ymax)``.
    Default value is None.

pandamesh.GmshMesher.add_threshold_distance_field
=================================================
Add a distance field to the mesher.

The of geometry of these fields are not forced into the mesh, but they
can be used to specify zones of with cell sizes.

Parameters
----------
gdf: geopandas.GeoDataFrame
    Location of the features to measure distance to. Should contain
    ``spacing`` column to specify the spacing of interpolated vertices
    along linestrings and polygon boundaries.

pandamesh.GmshMesher.field_combination
======================================
Controls how cell size fields are combined when they are found at the
same location. Can be set to one of
:py:class:`pandamesh.FieldCombination`:

.. code::

    MIN = "Min"
    MAX = "Max"

pandamesh.GmshMesher.fields
===========================
Read-only access to fields.

Use ``.clear_fields`` to remove fields from the mesher.

pandamesh.GmshMesher.finalize_gmsh
==================================
Finalize Gmsh.

pandamesh.GmshMesher.general_verbosity
======================================
Controls level of information printed. Can be set to one of
:py:class:`pandamesh.GeneralVerbosity`:

.. code::

    SILENT = 0
    ERRORS = 1
    WARNINGS = 2
    DIRECT = 3
    INFORMATION = 4
    STATUS = 5
    DEBUG = 99

pandamesh.GmshMesher.generate
=============================
Generate a mesh of triangles or quadrangles.

Parameters
----------
finalize: bool, default False
    Automatically finalize after generating.

Returns
-------
vertices: np.ndarray of floats with shape ``(n_vertex, 2)``
faces: np.ndarray of integers with shape ``(n_face, nmax_per_face)``
    ``nmax_per_face`` is 3 for exclusively triangles and 4 if
    quadrangles are included. A fill value of -1 is used as a last
    entry for triangles in that case.

pandamesh.GmshMesher.generate_geodataframe
==========================================
Generate a mesh and return it as a geopandas GeoDataFrame.

Parameters
----------
finalize: bool, default False
    Automatically finalize after generating.

Returns
-------
mesh: geopandas.GeoDataFrame

pandamesh.GmshMesher.generate_ugrid
===================================
Generate a mesh and return it as an xugrid Ugrid2d.

Parameters
----------
finalize: bool, default False
    Automatically finalize after generating.

Returns
-------
mesh: xugrid.Ugrid2d

pandamesh.GmshMesher.mesh_algorithm
===================================
Can be set to one of :py:class:`pandamesh.MeshAlgorithm`:

.. code::

    MESH_ADAPT = 1
    AUTOMATIC = 2
    INITIAL_MESH_ONLY = 3
    FRONTAL_DELAUNAY = 5
    BAMG = 7
    FRONTAL_DELAUNAY_FOR_QUADS = 8
    PACKING_OF_PARALLELLOGRAMS = 9
    QUASI_STRUCTURED_QUAD = 11

Each algorithm has its own advantages and disadvantages.

pandamesh.GmshMesher.mesh_size_extend_from_boundary
===================================================
Forces the mesh size to be extended from the boundary, or not, per
surface.

pandamesh.GmshMesher.mesh_size_from_curvature
=============================================
Automatically compute mesh element sizes from curvature, using the value as
the target number of elements per 2 * Pi radians.

pandamesh.GmshMesher.mesh_size_from_points
==========================================
Compute mesh element sizes from values given at geometry points.

pandamesh.GmshMesher.recombine_all
==================================
Apply recombination algorithm to all surfaces, ignoring per-surface
spec.

pandamesh.GmshMesher.subdivision_algorithm
==========================================
All meshes can be subdivided to generate fully quadrangular cells. Can
be set to one of :py:class:`pandamesh.SubdivisionAlgorithm`:

.. code::

    NONE = 0
    ALL_QUADRANGLES = 1
    BARYCENTRIC = 3

pandamesh.GmshMesher.write
==========================
Write a gmsh .msh file

Parameters
----------
path: Union[str, pathlib.Path

pandamesh.FieldCombination
==========================
Controls how cell size fields are combined in Gmsh when they are found at
the same location.

pandamesh.FieldCombination Class Members
========================================
   * pandamesh.FieldCombination.MAX
   * pandamesh.FieldCombination.MIN

pandamesh.FieldCombination.MAX
==============================
Controls how cell size fields are combined in Gmsh when they are found at
the same location.

pandamesh.FieldCombination.MIN
==============================
Controls how cell size fields are combined in Gmsh when they are found at
the same location.

pandamesh.GeneralVerbosity
==========================
Gmsh level of information printed.

pandamesh.GeneralVerbosity Class Members
========================================
   * pandamesh.GeneralVerbosity.DEBUG
   * pandamesh.GeneralVerbosity.DIRECT
   * pandamesh.GeneralVerbosity.ERRORS
   * pandamesh.GeneralVerbosity.INFORMATION
   * pandamesh.GeneralVerbosity.SILENT
   * pandamesh.GeneralVerbosity.STATUS
   * pandamesh.GeneralVerbosity.WARNINGS

pandamesh.GeneralVerbosity.DEBUG
================================
Gmsh level of information printed.

pandamesh.GeneralVerbosity.DIRECT
=================================
Gmsh level of information printed.

pandamesh.GeneralVerbosity.ERRORS
=================================
Gmsh level of information printed.

pandamesh.GeneralVerbosity.INFORMATION
======================================
Gmsh level of information printed.

pandamesh.GeneralVerbosity.SILENT
=================================
Gmsh level of information printed.

pandamesh.GeneralVerbosity.STATUS
=================================
Gmsh level of information printed.

pandamesh.GeneralVerbosity.WARNINGS
===================================
Gmsh level of information printed.

pandamesh.MeshAlgorithm
=======================
Gmsh meshing algorithm. Each algorithm has its own advantages and
disadvantages.

For all 2D unstructured algorithms a Delaunay mesh that contains all
the points of the 1D mesh is initially constructed using a
divide-and-conquer algorithm. Missing edges are recovered using edge
swaps. After this initial step several algorithms can be applied to
generate the final mesh:

* The MeshAdapt algorithm is based on local mesh modifications. This
  technique makes use of edge swaps, splits, and collapses: long edges
  are split, short edges are collapsed, and edges are swapped if a
  better geometrical configuration is obtained.
* The Delaunay algorithm is inspired by the work of the GAMMA team at
  INRIA. New points are inserted sequentially at the circumcenter of
  the element that has the largest adimensional circumradius. The mesh
  is then reconnected using an anisotropic Delaunay criterion.
* The Frontal-Delaunay algorithm is inspired by the work of S. Rebay.
* Other experimental algorithms with specific features are also
  available. In particular, Frontal-Delaunay for Quads is a variant of
  the Frontal-Delaunay algorithm aiming at generating right-angle
  triangles suitable for recombination; and BAMG allows to generate
  anisotropic triangulations.

For very complex curved surfaces the MeshAdapt algorithm is the most robust.
When high element quality is important, the Frontal-Delaunay algorithm should
be tried. For very large meshes of plane surfaces the Delaunay algorithm is
the fastest; it usually also handles complex mesh size fields better than the
Frontal-Delaunay. When the Delaunay or Frontal-Delaunay algorithms fail,
MeshAdapt is automatically triggered. The Automatic algorithm uses
Delaunay for plane surfaces and MeshAdapt for all other surfaces.

pandamesh.MeshAlgorithm Class Members
=====================================
   * pandamesh.MeshAlgorithm.AUTOMATIC
   * pandamesh.MeshAlgorithm.BAMG
   * pandamesh.MeshAlgorithm.FRONTAL_DELAUNAY
   * pandamesh.MeshAlgorithm.FRONTAL_DELAUNAY_FOR_QUADS
   * pandamesh.MeshAlgorithm.INITIAL_MESH_ONLY
   * pandamesh.MeshAlgorithm.MESH_ADAPT
   * pandamesh.MeshAlgorithm.PACKING_OF_PARALLELLOGRAMS
   * pandamesh.MeshAlgorithm.QUASI_STRUCTURED_QUAD

pandamesh.MeshAlgorithm.AUTOMATIC
=================================
Gmsh meshing algorithm. Each algorithm has its own advantages and
disadvantages.

For all 2D unstructured algorithms a Delaunay mesh that contains all
the points of the 1D mesh is initially constructed using a
divide-and-conquer algorithm. Missing edges are recovered using edge
swaps. After this initial step several algorithms can be applied to
generate the final mesh:

* The MeshAdapt algorithm is based on local mesh modifications. This
  technique makes use of edge swaps, splits, and collapses: long edges
  are split, short edges are collapsed, and edges are swapped if a
  better geometrical configuration is obtained.
* The Delaunay algorithm is inspired by the work of the GAMMA team at
  INRIA. New points are inserted sequentially at the circumcenter of
  the element that has the largest adimensional circumradius. The mesh
  is then reconnected using an anisotropic Delaunay criterion.
* The Frontal-Delaunay algorithm is inspired by the work of S. Rebay.
* Other experimental algorithms with specific features are also
  available. In particular, Frontal-Delaunay for Quads is a variant of
  the Frontal-Delaunay algorithm aiming at generating right-angle
  triangles suitable for recombination; and BAMG allows to generate
  anisotropic triangulations.

For very complex curved surfaces the MeshAdapt algorithm is the most robust.
When high element quality is important, the Frontal-Delaunay algorithm should
be tried. For very large meshes of plane surfaces the Delaunay algorithm is
the fastest; it usually also handles complex mesh size fields better than the
Frontal-Delaunay. When the Delaunay or Frontal-Delaunay algorithms fail,
MeshAdapt is automatically triggered. The Automatic algorithm uses
Delaunay for plane surfaces and MeshAdapt for all other surfaces.

pandamesh.MeshAlgorithm.BAMG
============================
Gmsh meshing algorithm. Each algorithm has its own advantages and
disadvantages.

For all 2D unstructured algorithms a Delaunay mesh that contains all
the points of the 1D mesh is initially constructed using a
divide-and-conquer algorithm. Missing edges are recovered using edge
swaps. After this initial step several algorithms can be applied to
generate the final mesh:

* The MeshAdapt algorithm is based on local mesh modifications. This
  technique makes use of edge swaps, splits, and collapses: long edges
  are split, short edges are collapsed, and edges are swapped if a
  better geometrical configuration is obtained.
* The Delaunay algorithm is inspired by the work of the GAMMA team at
  INRIA. New points are inserted sequentially at the circumcenter of
  the element that has the largest adimensional circumradius. The mesh
  is then reconnected using an anisotropic Delaunay criterion.
* The Frontal-Delaunay algorithm is inspired by the work of S. Rebay.
* Other experimental algorithms with specific features are also
  available. In particular, Frontal-Delaunay for Quads is a variant of
  the Frontal-Delaunay algorithm aiming at generating right-angle
  triangles suitable for recombination; and BAMG allows to generate
  anisotropic triangulations.

For very complex curved surfaces the MeshAdapt algorithm is the most robust.
When high element quality is important, the Frontal-Delaunay algorithm should
be tried. For very large meshes of plane surfaces the Delaunay algorithm is
the fastest; it usually also handles complex mesh size fields better than the
Frontal-Delaunay. When the Delaunay or Frontal-Delaunay algorithms fail,
MeshAdapt is automatically triggered. The Automatic algorithm uses
Delaunay for plane surfaces and MeshAdapt for all other surfaces.

pandamesh.MeshAlgorithm.FRONTAL_DELAUNAY
========================================
Gmsh meshing algorithm. Each algorithm has its own advantages and
disadvantages.

For all 2D unstructured algorithms a Delaunay mesh that contains all
the points of the 1D mesh is initially constructed using a
divide-and-conquer algorithm. Missing edges are recovered using edge
swaps. After this initial step several algorithms can be applied to
generate the final mesh:

* The MeshAdapt algorithm is based on local mesh modifications. This
  technique makes use of edge swaps, splits, and collapses: long edges
  are split, short edges are collapsed, and edges are swapped if a
  better geometrical configuration is obtained.
* The Delaunay algorithm is inspired by the work of the GAMMA team at
  INRIA. New points are inserted sequentially at the circumcenter of
  the element that has the largest adimensional circumradius. The mesh
  is then reconnected using an anisotropic Delaunay criterion.
* The Frontal-Delaunay algorithm is inspired by the work of S. Rebay.
* Other experimental algorithms with specific features are also
  available. In particular, Frontal-Delaunay for Quads is a variant of
  the Frontal-Delaunay algorithm aiming at generating right-angle
  triangles suitable for recombination; and BAMG allows to generate
  anisotropic triangulations.

For very complex curved surfaces the MeshAdapt algorithm is the most robust.
When high element quality is important, the Frontal-Delaunay algorithm should
be tried. For very large meshes of plane surfaces the Delaunay algorithm is
the fastest; it usually also handles complex mesh size fields better than the
Frontal-Delaunay. When the Delaunay or Frontal-Delaunay algorithms fail,
MeshAdapt is automatically triggered. The Automatic algorithm uses
Delaunay for plane surfaces and MeshAdapt for all other surfaces.

pandamesh.MeshAlgorithm.FRONTAL_DELAUNAY_FOR_QUADS
==================================================
Gmsh meshing algorithm. Each algorithm has its own advantages and
disadvantages.

For all 2D unstructured algorithms a Delaunay mesh that contains all
the points of the 1D mesh is initially constructed using a
divide-and-conquer algorithm. Missing edges are recovered using edge
swaps. After this initial step several algorithms can be applied to
generate the final mesh:

* The MeshAdapt algorithm is based on local mesh modifications. This
  technique makes use of edge swaps, splits, and collapses: long edges
  are split, short edges are collapsed, and edges are swapped if a
  better geometrical configuration is obtained.
* The Delaunay algorithm is inspired by the work of the GAMMA team at
  INRIA. New points are inserted sequentially at the circumcenter of
  the element that has the largest adimensional circumradius. The mesh
  is then reconnected using an anisotropic Delaunay criterion.
* The Frontal-Delaunay algorithm is inspired by the work of S. Rebay.
* Other experimental algorithms with specific features are also
  available. In particular, Frontal-Delaunay for Quads is a variant of
  the Frontal-Delaunay algorithm aiming at generating right-angle
  triangles suitable for recombination; and BAMG allows to generate
  anisotropic triangulations.

For very complex curved surfaces the MeshAdapt algorithm is the most robust.
When high element quality is important, the Frontal-Delaunay algorithm should
be tried. For very large meshes of plane surfaces the Delaunay algorithm is
the fastest; it usually also handles complex mesh size fields better than the
Frontal-Delaunay. When the Delaunay or Frontal-Delaunay algorithms fail,
MeshAdapt is automatically triggered. The Automatic algorithm uses
Delaunay for plane surfaces and MeshAdapt for all other surfaces.

pandamesh.MeshAlgorithm.INITIAL_MESH_ONLY
=========================================
Gmsh meshing algorithm. Each algorithm has its own advantages and
disadvantages.

For all 2D unstructured algorithms a Delaunay mesh that contains all
the points of the 1D mesh is initially constructed using a
divide-and-conquer algorithm. Missing edges are recovered using edge
swaps. After this initial step several algorithms can be applied to
generate the final mesh:

* The MeshAdapt algorithm is based on local mesh modifications. This
  technique makes use of edge swaps, splits, and collapses: long edges
  are split, short edges are collapsed, and edges are swapped if a
  better geometrical configuration is obtained.
* The Delaunay algorithm is inspired by the work of the GAMMA team at
  INRIA. New points are inserted sequentially at the circumcenter of
  the element that has the largest adimensional circumradius. The mesh
  is then reconnected using an anisotropic Delaunay criterion.
* The Frontal-Delaunay algorithm is inspired by the work of S. Rebay.
* Other experimental algorithms with specific features are also
  available. In particular, Frontal-Delaunay for Quads is a variant of
  the Frontal-Delaunay algorithm aiming at generating right-angle
  triangles suitable for recombination; and BAMG allows to generate
  anisotropic triangulations.

For very complex curved surfaces the MeshAdapt algorithm is the most robust.
When high element quality is important, the Frontal-Delaunay algorithm should
be tried. For very large meshes of plane surfaces the Delaunay algorithm is
the fastest; it usually also handles complex mesh size fields better than the
Frontal-Delaunay. When the Delaunay or Frontal-Delaunay algorithms fail,
MeshAdapt is automatically triggered. The Automatic algorithm uses
Delaunay for plane surfaces and MeshAdapt for all other surfaces.

pandamesh.MeshAlgorithm.MESH_ADAPT
==================================
Gmsh meshing algorithm. Each algorithm has its own advantages and
disadvantages.

For all 2D unstructured algorithms a Delaunay mesh that contains all
the points of the 1D mesh is initially constructed using a
divide-and-conquer algorithm. Missing edges are recovered using edge
swaps. After this initial step several algorithms can be applied to
generate the final mesh:

* The MeshAdapt algorithm is based on local mesh modifications. This
  technique makes use of edge swaps, splits, and collapses: long edges
  are split, short edges are collapsed, and edges are swapped if a
  better geometrical configuration is obtained.
* The Delaunay algorithm is inspired by the work of the GAMMA team at
  INRIA. New points are inserted sequentially at the circumcenter of
  the element that has the largest adimensional circumradius. The mesh
  is then reconnected using an anisotropic Delaunay criterion.
* The Frontal-Delaunay algorithm is inspired by the work of S. Rebay.
* Other experimental algorithms with specific features are also
  available. In particular, Frontal-Delaunay for Quads is a variant of
  the Frontal-Delaunay algorithm aiming at generating right-angle
  triangles suitable for recombination; and BAMG allows to generate
  anisotropic triangulations.

For very complex curved surfaces the MeshAdapt algorithm is the most robust.
When high element quality is important, the Frontal-Delaunay algorithm should
be tried. For very large meshes of plane surfaces the Delaunay algorithm is
the fastest; it usually also handles complex mesh size fields better than the
Frontal-Delaunay. When the Delaunay or Frontal-Delaunay algorithms fail,
MeshAdapt is automatically triggered. The Automatic algorithm uses
Delaunay for plane surfaces and MeshAdapt for all other surfaces.

pandamesh.MeshAlgorithm.PACKING_OF_PARALLELLOGRAMS
==================================================
Gmsh meshing algorithm. Each algorithm has its own advantages and
disadvantages.

For all 2D unstructured algorithms a Delaunay mesh that contains all
the points of the 1D mesh is initially constructed using a
divide-and-conquer algorithm. Missing edges are recovered using edge
swaps. After this initial step several algorithms can be applied to
generate the final mesh:

* The MeshAdapt algorithm is based on local mesh modifications. This
  technique makes use of edge swaps, splits, and collapses: long edges
  are split, short edges are collapsed, and edges are swapped if a
  better geometrical configuration is obtained.
* The Delaunay algorithm is inspired by the work of the GAMMA team at
  INRIA. New points are inserted sequentially at the circumcenter of
  the element that has the largest adimensional circumradius. The mesh
  is then reconnected using an anisotropic Delaunay criterion.
* The Frontal-Delaunay algorithm is inspired by the work of S. Rebay.
* Other experimental algorithms with specific features are also
  available. In particular, Frontal-Delaunay for Quads is a variant of
  the Frontal-Delaunay algorithm aiming at generating right-angle
  triangles suitable for recombination; and BAMG allows to generate
  anisotropic triangulations.

For very complex curved surfaces the MeshAdapt algorithm is the most robust.
When high element quality is important, the Frontal-Delaunay algorithm should
be tried. For very large meshes of plane surfaces the Delaunay algorithm is
the fastest; it usually also handles complex mesh size fields better than the
Frontal-Delaunay. When the Delaunay or Frontal-Delaunay algorithms fail,
MeshAdapt is automatically triggered. The Automatic algorithm uses
Delaunay for plane surfaces and MeshAdapt for all other surfaces.

pandamesh.MeshAlgorithm.QUASI_STRUCTURED_QUAD
=============================================
Gmsh meshing algorithm. Each algorithm has its own advantages and
disadvantages.

For all 2D unstructured algorithms a Delaunay mesh that contains all
the points of the 1D mesh is initially constructed using a
divide-and-conquer algorithm. Missing edges are recovered using edge
swaps. After this initial step several algorithms can be applied to
generate the final mesh:

* The MeshAdapt algorithm is based on local mesh modifications. This
  technique makes use of edge swaps, splits, and collapses: long edges
  are split, short edges are collapsed, and edges are swapped if a
  better geometrical configuration is obtained.
* The Delaunay algorithm is inspired by the work of the GAMMA team at
  INRIA. New points are inserted sequentially at the circumcenter of
  the element that has the largest adimensional circumradius. The mesh
  is then reconnected using an anisotropic Delaunay criterion.
* The Frontal-Delaunay algorithm is inspired by the work of S. Rebay.
* Other experimental algorithms with specific features are also
  available. In particular, Frontal-Delaunay for Quads is a variant of
  the Frontal-Delaunay algorithm aiming at generating right-angle
  triangles suitable for recombination; and BAMG allows to generate
  anisotropic triangulations.

For very complex curved surfaces the MeshAdapt algorithm is the most robust.
When high element quality is important, the Frontal-Delaunay algorithm should
be tried. For very large meshes of plane surfaces the Delaunay algorithm is
the fastest; it usually also handles complex mesh size fields better than the
Frontal-Delaunay. When the Delaunay or Frontal-Delaunay algorithms fail,
MeshAdapt is automatically triggered. The Automatic algorithm uses
Delaunay for plane surfaces and MeshAdapt for all other surfaces.

pandamesh.SubdivisionAlgorithm
==============================
Controls how Gmsh recombines triangles to form quads.

The default recombination algorithm might leave some triangles in the mesh,
if recombining all the triangles leads to badly shaped quads. In such
cases, to generate full-quad meshes, you can either subdivide the resulting
hybrid mesh (ALL_QUADRANGLES), or use the full-quad recombination
algorithm, which will automatically perform a coarser mesh followed by
recombination, smoothing and subdivision.

pandamesh.SubdivisionAlgorithm Class Members
============================================
   * pandamesh.SubdivisionAlgorithm.ALL_QUADRANGLES
   * pandamesh.SubdivisionAlgorithm.BARYCENTRIC
   * pandamesh.SubdivisionAlgorithm.NONE

pandamesh.SubdivisionAlgorithm.ALL_QUADRANGLES
==============================================
Controls how Gmsh recombines triangles to form quads.

The default recombination algorithm might leave some triangles in the mesh,
if recombining all the triangles leads to badly shaped quads. In such
cases, to generate full-quad meshes, you can either subdivide the resulting
hybrid mesh (ALL_QUADRANGLES), or use the full-quad recombination
algorithm, which will automatically perform a coarser mesh followed by
recombination, smoothing and subdivision.

pandamesh.SubdivisionAlgorithm.BARYCENTRIC
==========================================
Controls how Gmsh recombines triangles to form quads.

The default recombination algorithm might leave some triangles in the mesh,
if recombining all the triangles leads to badly shaped quads. In such
cases, to generate full-quad meshes, you can either subdivide the resulting
hybrid mesh (ALL_QUADRANGLES), or use the full-quad recombination
algorithm, which will automatically perform a coarser mesh followed by
recombination, smoothing and subdivision.

pandamesh.SubdivisionAlgorithm.NONE
===================================
Controls how Gmsh recombines triangles to form quads.

The default recombination algorithm might leave some triangles in the mesh,
if recombining all the triangles leads to badly shaped quads. In such
cases, to generate full-quad meshes, you can either subdivide the resulting
hybrid mesh (ALL_QUADRANGLES), or use the full-quad recombination
algorithm, which will automatically perform a coarser mesh followed by
recombination, smoothing and subdivision.

