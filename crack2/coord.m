function [pos, A] = coord(Nbd)
global length;
global width;
global dx;
global dy;

nx = length/dx +1;
ny = width/dy+1;
nodes=nx*ny;
pos = zeros(nodes,2);
A = zeros(nodes,nodes);
for i=1:nodes
        pos(i,:) = ser2pos(i);
        nbdsers = find(Nbd(i,:));
        for k=nbdsers
            nbd = Nbd(i,k);
            A(i,nbd) = 1;
        end
end
        
   