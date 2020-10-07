import numpy as np
import pygmsh

# geom = pygmsh.built_in.Geometry()
# poly = geom.add_polygon([
        # [0.1, 0],
        # [0.2, 0],
        # [0.2, 0.04],
        # [0.1, 0.04]
        # ],
        # lcar = 0.05
# )
geom = pygmsh.built_in.Geometry()
poly = geom.add_polygon([
    [ 0.0,  0.5, 0.0],
    [-0.1,  0.1, 0.0],
    [-0.5,  0.0, 0.0],
    [-0.1, -0.1, 0.0],
    [ 0.0, -0.5, 0.0],
    [ 0.1, -0.1, 0.0],
    [ 0.5,  0.0, 0.0],
    [ 0.1,  0.1, 0.0]
    ],
    lcar=0.05
)

