function drawmesh(N, dx, dy, nx, ny)
[pos, A] = coord(N, dx, dy, nx, ny);
gplot(A, pos)
