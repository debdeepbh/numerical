function A = ser2pos(y)
    B = ser2ind(y);
    A = ind2pos(B(1),B(2));