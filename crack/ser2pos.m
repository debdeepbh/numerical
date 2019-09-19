function A = ser2pos(y,dx,dy,nx, ny)
    B = ser2ind(y,nx, ny);
    A = ind2pos(B(1),B(2), dx, dy);
