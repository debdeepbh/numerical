function intrsct=intersects(A1,A2,A3,A4) %Ai = [x y], the coordinate
% does NOT work for disjoint segments of same line
A = [A1;A2;A3;A4];
x= A(:,1);
y= A(:,2);

% if parallel, make sure distance is not zero(y-intercept)
 D1 = A2-A1;
 D2 = A4-A3;

    
        
    
    
%     slope1 = inf
% else
%     slope1 = D(2)/D(1);
%     slope2 = D(3)
    
dt1=det([1,1,1;x(1),x(2),x(3);y(1),y(2),y(3)])*det([1,1,1;x(1),x(2),x(4);y(1),y(2),y(4)]);
dt2=det([1,1,1;x(1),x(3),x(4);y(1),y(3),y(4)])*det([1,1,1;x(2),x(3),x(4);y(2),y(3),y(4)]);

if(dt1<=0 & dt2<=0)
  intrsct=1;         %If lines intesect
  if (all(y == y(1)))   %check if disjoint segment of a horizontal line intersects
      %disp('y values are same');
    xmin = min(x);
    if numel((find(x == xmin))) == 1
        i = find(x == xmin);
        if i == 1
            j = 2;
            if((x(j) < x(3)) &&( x(j) <x(4)))  
                intrsct = 0;
            end
        elseif i == 2
            j = 1;
            if((x(j) < x(3)) &&( x(j) <x(4)))  
                intrsct = 0;
            end
        elseif i == 3
            j = 4;
            if((x(j) < x(1)) &&( x(j) <x(2)))  
                intrsct = 0;
            end
        elseif i ==4
            j = 3;
            if((x(j) < x(1)) &&( x(j) <x(2)))  
                intrsct = 0;
            end
        end
    else
        intrsct = 1;
    end
end

else
 intrsct=0;
end