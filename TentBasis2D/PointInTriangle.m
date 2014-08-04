function flag = PointInTriangle(P, A,B,C) 
% This script will decide whether point P is inside the triangle defined by
% point A, B, C or not (edge included). 
% The Input P can also be a Nx2 vector, each row is a point
%
% This method is sometimes referred as the "Barycentric Technique"
% for more info. see the following link
% http://en.wikipedia.org/wiki/Barycentric_coordinates_%28mathematics%29
% or the one I used
% http://www.blackpawn.com/texts/pointinpoly/default.html

% Implemented by Yuwei Cui, Oct 30. 2012
A = A(:)';
B = B(:)';
C = C(:)';
Npt = size(P,1);
% Compute vectors        
v0 = C - A;
v1 = B - A;
v2 = P - ones(Npt,1)*A;

% Compute dot products
dot00 = v0*v0';
dot01 = v0*v1';
dot11 = v1*v1';
dot02 = v2*v0';
dot12 = v2*v1';

% Compute barycentric coordinates
invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
u = (dot11 * dot02 - dot01 * dot12) * invDenom;
v = (dot00 * dot12 - dot01 * dot02) * invDenom;

% Check if point is in triangle
flag = (u >= 0) & (v >= 0) & (u + v <= 1);