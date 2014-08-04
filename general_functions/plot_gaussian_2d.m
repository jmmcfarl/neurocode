function h = plot_gaussian_2d(gmean,gsigma,sig_widths,color,linewidth)

if nargin < 3 || isempty(color)
    color = 'b';
end
if nargin < 4 || isempty(linewidth)
    linewidth = 1;
end
    
n_pts = 50;

tt=linspace(0,2*pi,n_pts)';
x = cos(tt); y=sin(tt);
ap = [x(:) y(:)]';
[v,d]=eig(gsigma); 

for ii = 1:length(sig_widths)
    cur_d = sig_widths(ii) * sqrt(d); % convert variance to sdwidth*sd
    bp = (v*cur_d*ap) + repmat(gmean, 1, size(ap,2));
    h(ii) = plot(bp(1,:), bp(2,:),'color',color,'linewidth',linewidth);
end