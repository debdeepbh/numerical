ss = size(ref);

%diffabs = zeros(ss);
%diffabs2 = zeros(ss);
%diffabs3 = zeros(ss);
%for i=1:ss(1)
%	for j = 1:ss(2)
%		if ref(i,j) > 0
%		diffabs(i,j) = abs(ser2pos1(ref(i,j)) - ser2pos1(i));
%		diffabs2(i,j) = abs(ser2pos2(ref(i,j)) - ser2pos2(i));
%	end
%	end
%end
%
%max(max(diffabs)) > 0.02
%max(max(diffabs2))> 0.02

t1 = ser2pos1(ref);
rowv = ser2pos1((1:ss(1))');

max(max(abs((t1 - rowv).*(~~ref))))

%badp = find(abs((t1 - rowv).*(~~ref)) > 0.002);

%max(t1(badp) - t1(badp)) > 0 

%for i=1:ss(1)
%	for j=1:ss(2)
%		if ref(i,j) > 0
%			diffabs3(i,j) = abs(t1(i,j) - ser2pos1(i));
%	else
%		diffabs3(i,j) = 0;
%		end
%	end
%end
%max(max(diffabs3))

%for i = 1:ss(1)
%	for j=1:ss(2)
%		realpos(i,j) = ser2pos1(ref(i,j));
%	end
%end
%
%fakepos = ser2pos1(ref);
%
%sum(realpos ~= fakepos)

