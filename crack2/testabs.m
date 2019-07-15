function testabs(Nbd)

ss = size(Nbd);
vect = (1:ss(1))';

u_r = ser2pos1(Nbd);
r = ser2pos1(vect);

diff = (u_r - r).*(~~Nbd);
max(max(abs(diff)))

