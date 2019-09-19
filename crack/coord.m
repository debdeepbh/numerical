function [pos, A] = coord(Nbd, dx, dy, nx, ny)

% updated to depend only on ser2pos1, ser2pos2

nodes=nx*ny;
pos = zeros(nodes,2);
A = zeros(nodes,nodes);
for i=1:nodes
        %pos(i,:) = ser2pos(i, dx, dy, nx, ny);
        pos(i,1) = ser2pos1(i, dx, nx);
        pos(i,2) = ser2pos2(i, dy, nx);
        nbdsers = find(Nbd(i,:));
        for k=nbdsers
            nbd = Nbd(i,k);
            A(i,nbd) = 1;
        end
end
        
   
