# Meshing in Octave

Tutorial from [link](http://wiki.octave.org/Geometry_package#Relation_to_matGeom).

## Installation

* Dependencies: `gmsh` and `liboctave-dev`

* In `octave --no-gui` prompt, install the packages
```
pkg install -forge geometry msh fpl
```
If all the dependencies are needed, do:
```
pkg install -forge matgeom geometry splines msh
```

### Q: do we need to install `octave-geometry` using `apt-get install`?

* Define a polygon using
```
P = [10 10; 40 10; 40 30; 20 25; 15 30];
```

#### Plotting a polygon
* Load the geom package
```
pkg load geometry
```
*  Plot the polygon
```
drawPolygon(P, '-o')
```

### Creating `.geo` files from polygons
* Generate the `.geo` file from the polygon using `data2geo()` function from `geometry` package:
```
meshsize = 2;
data2geo (P, meshsize, "output", "tempfile.geo")
```
Some kind of suggestion to take the mesh size is given by (what is that?)
```
meshsize = sqrt (mean (sumsq (diff (P, 1, 1), 2)))/2;
```
### Generating a mesh data structure from a `.geo` file
* Load the meshing library
```
pkg load msh
```
* Generate the mesh using (without the extension `.geo`. Necessary?)
```
T = msh2m_gmsh("tempfile")
```
See [link](https://octave.sourceforge.io/msh/overview.html) for more help and other types of mesh.

The output data `T` has a _MESH_ data structure, which is compatible with that of Matlab PDE toolbox. The data contains 3 information (see [link](https://octave.sourceforge.io/msh/function/msh2m_structured_mesh.html)):

 1. `T.p`: coordinate of the nodes. It is a matrix with size 2 times number of mesh points.
------- | -------------
1st row | x-coordinates of the points.
2nd row | y-coordinates of the points. 
------- | -------------
 2. `T.e`: boundary-edges (edges that lie on the boundary of the region). It is a matrix with size 7 times number of mesh side edges.
------- | -------------
1st row | index of the first vertex of the side edge.
2nd row | index of the second vertex of the side edge.
3rd row | set to 0, present for compatibility with MatLab PDE-tool.
4th row | set to 0, present for compatibility with MatLab PDE-tool.
5th row | index of the geometrical border containing the side edge. (This index comes from the boundary provided in the `.geo` file)
6th row | index of the geometrical surface to the right of side edge.
7th row | index of the geometrical surface to the left of the side edge. 
------- | -------------
 3. `T.t`: triangular elements. It is a matrix with size 4 times number of mesh elements.
------- | -------------
1st row | index of the first vertex of the element.
2nd row | index of the second vertex of the element.
3rd row | index of the third vertex of the element.
4th row | index of the geometrical surface containing the element. 
------- | -------------

**Note:** Here, nodes are the vertices of the triangles, elements are the triangles, and the edges are lines between vertices. Sides are boundaries.

**Note:** the order of the edges and the elements seem to preserve the right-handed orientation, i.e. nodes of an element are listed in the counterclockwise direction (provided the polygon data in the `.geo` file was provided in that direction).
	

### Plotting the mesh data structure
* Load the FEM plotting library `fpl`
```
pkg load fpl
```
* Plot the mesh (see [link](https://octave.sourceforge.io/fpl/function/pdemesh.html))
```
pdemesh (T.p, T.e, T.t);
```
* `pdemesh` plots in 3D. View it in 2D using 
```
view(2)
```
**Note:** `T.e` is actually a dummy input. It is only added to make the syntax consistent with Matlab's PDE toolbox function `pdemesh()`.

### Boundary (sides):
* Total number of boundaries = `max(T.e(5,:))` (verify)
* **Nodes:** Even though it is easy to compute, nodes that lie on certain borders (e.g. border indices 1, 2, and 4) of the region can be returned with 
```
msh2m_nodes_on_sides(T,[1 2 4])
```
* **Edges:** All edges on the boundary
```
msh2m_topological_properties(T, "sides")
```
* **Elements:** Element indices that share an edge with the boundary of the region
```
msh2m_topological_properties(T, "boundary")
```
The first row of this output contains the element indices. The second row of this output contains the boundary index of the corresponding edge.
**Note:** There will be repeated entries for elements which shares more than one boundary indices.

### Geometric properties of the mesh data structure
* Get the coordinates of the centroid of the triangular elements
```
T_centroid = msh2m_geometrical_properties(mesh, "bar")
```
* Extract multiple geometric information from the _MESH_ data structure (see [link](https://octave.sourceforge.io/msh/function/msh2m_geometrical_properties.html)) by comma-separating the strings
```
[T_centroid, T_area] = msh2m_geometrical_properties(mesh, "bar", "area")
```

### Topological properties of the mesh data structure
* Get the element indices of neighbors of each element (`NaN` for missing values)
```
msh2m_topological_properties(T, "n")
```

### Extra

#### SVG to polygon
* Get polygon from a `.svg` file using
```
octavesvg = svg ("octave.svg").normalize();
ids       = octavesvg.pathid();
P         = octavesvg.path2polygon (ids{1}, 12)(1:end-1,:);
P         = bsxfun (@minus, P, centroid (P));
```
* If some complex polygon has too many nodes, simplify them using Douglas-Peucker algorithm using
```
P  = simplifypolygon(P, 'tol', 1e-3);
```

#### VTK output format
* `fpl` allows output in VTK format [link](https://octave.sourceforge.io/fpl/overview.html). Use `savevtk`, `savevtkvector` etc.

#### Modifying the mesh
```
T_J = msh2m_jiggle_mesh (T, 100)
```
returns a new mesh by considering the edges are springs and solves a PDE to get the quasi static equilibrium, where the second argument is the number of steps for which the PDE is simulated.
See also: `msh2m_equalize_mesh()`, `msh2m_displacement_smoothing()` and `mshm_refine()` (which requires `dolfin` package, which is part of `fem-fenics` package (verify))

# Peridynamic simulations

The goal is to simulate multi-particle interaction.

Steps

* Generate a mesh, choose the nodes, compute the nodal volume with `script.m`
* Compute a neighborhood array that contains the indices of the neighbors of each of the nodes using `gen_nbdlist.m`
* Simulate the dynamics with the peridynamic force function

### Code
* `genmesh()` (Octave) produces a mesh where the centroid of the triangular element is a node and the nodal volume is the area of the triangle.
* `genmesh2()` (Octave) produces a mesh where each vertex of the triangle is a node and the nodal volume is the sum of 1/3 of each of the neighboring triangular elements.

##### Running the code
1. In Octave, run `full_proj`, close the program after various matrices are saved to the disk (e.g. `Pos.m`  etc). The mesh is generated in this step and the nodal volumes are computed (`Vol.mat`).
1. In Matlab, run `full_proj_matlab` to load the information and simulate the dynamics.


# Multithreading
* See the number of threads allowed using `maxNumCompThreads`. Set the number using `maxNumCompThreads(4)`. Check the effect using a long computation like this:
```
X=ones(10000);
sqrt(X).^sqrt(X).^sqrt(X).^sqrt(X).^sqrt(X).^sqrt(X)
```

# Physics with numerical schemes

* When using velocity Verlet scheme, if the body is under a constant gravity, we need to specify that gravity in initial acceleration (`uolddotdot`) as well as the external force density (`extforce`) as `extforce = zeros(totalnodes, 2) + [0 -10] * rho`. To see that we need to specify this is both places, we can consider the ODE _u''(t) = a_ whose solution is _u(t) = 1/2 * a * t^2_. This is the case when the total internal force is zero. So, after `n` iterations, the displacement should be _1/2*a * (n *dt)^2_ where _dt_ is the timestep size. Follow the steps of the scheme and make sure that that is the case.
