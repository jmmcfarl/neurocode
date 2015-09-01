function h=gaussplot(mu, Sigma, ncontours,color)
% GAUSSDRAW Plot some contours of constant probability density of a 2D
% Gaussian centered o
% n the mean
% h=gaussdraw(mu, Sigma, ncontours)

if ~exist('ncontours'), ncontours = -1; end

% Modelled afer demgauss from netlab.

ngrid = 50;
sx = 1*Sigma(1,1);
sy = 1*Sigma(2,2);
Xmin = mu(1)-3*sx; Xmax = mu(1)+3*sx;
Ymin = mu(2)-3*sy; Ymax = mu(2)+3*sy;
Xvals = linspace(Xmin, Xmax, ngrid);
Yvals = linspace(Ymin, Ymax, ngrid);
[X1, X2] = meshgrid(Xvals, Yvals);
XX = [X1(:), X2(:)];
% the i'th row of XX is the (x,y) coord of the i'th point in the raster
% scan of the grid

probs = gauss(mu, Sigma, XX);
probs = reshape(probs, ngrid, ngrid);

if ncontours==-1
  [C,h]=contour(fliplr(X1), flipud(X2), probs, 'Color',color,'linewidth',2);
else
  [C,h]=contour(fliplr(X1), flipud(X2), probs, ncontours, 'Color',color,'linewidth',2);
end