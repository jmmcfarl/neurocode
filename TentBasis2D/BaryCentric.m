function lambda = BaryCentric(P,A,B,C)
% Calculate Barycenteric Coordinate of point P with respect to
% Triangle ABC
% Yuwei Cui, Oct 30, 2012
Npt = size(P,1);
lambda = zeros(Npt,3);

x = [A(1) B(1) C(1)];
y = [A(2) B(2) C(2)];

M = ((y(2)-y(3))*(x(1)-x(3))+(x(3)-x(2))*(y(1)-y(3)));
lambda(:,1) = ((y(2)-y(3))*(P(:,1)-x(3))+(x(3)-x(2))*(P(:,2)-y(3)))/M;
lambda(:,2) = ((y(3)-y(1))*(P(:,1)-x(3))+(x(1)-x(3))*(P(:,2)-y(3)))/M;

% lambda(:,1) = ((B(2)-C(2))*(P(:,1)-C(1))+(C(1)-B(1))*(P(:,2)-C(2)))./...
%     ((B(2)-C(2))*(A(1)-C(1))+(C(1)-B(1))*(A(2)-C(2)));
% lambda(:,2) = ((C(2)-A(2))*(P(:,1)-C(1))+(A(1)-C(1))*(P(:,2)-C(2)))./...
%     ((B(2)-C(2))*(A(1)-C(1))+(C(1)-B(1))*(A(2)-C(2)));

lambda(:,3) = ones(Npt,1)-lambda(:,1)-lambda(:,2);