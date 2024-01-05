function [X,Y] = gethexind(n,s)
% generate the xy coordinates of the centers of the hexagonal grid elements
% n is the dimension of the array, i.e. the number of grid elements is n^2
% coordinates are output as row vectors of length n^2
if nargin < 2
    s = 1;
end;
X=[];
Y=[];
xinit = 0;
yinit = 0; 
xvec=0:2*s:2*(n-1)*s;  % vector of x-centers along a row
for i=1:n % for each row
    if mod(i,2)==1  % for odd numbered rows
        X=[X,xinit+xvec];
    else
        X=[X,xinit+s+xvec]; % for even numbered rows
    end;
    Y = [Y,2*(i-1)*s*ones(1,n)];
end

