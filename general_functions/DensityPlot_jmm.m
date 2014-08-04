function [handles, details] = DensityPlot_jmm(x,y,varargin)
% handles = DensityPlot(x,y,varargin)
% plots a density plot of the distribution of x,y
%DensityPlot(..., 'nbins', N) bins data into NxN bins (default 50);
%DensityPlot(..., 'exact') N bins = range of data (useful for ints)
%DensityPlot(..., 'ynormal') sets y axis to go in  normal direction


%get rid of nanvalues
nan_vals = find(isnan(x) | isnan(y));
x(nan_vals) = []; y(nan_vals) = [];

yrange = minmax(y);
xrange = minmax(x);
setydir = 0;
nbins = [50 50];
sx = 0; sy = 0;
sqrt_scale = 0;
plot_on = 1;
handles = [];
details = [];
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'exact',4)
        nbins(1) = diff(xrange)+1;
        nbins(2) = diff(yrange)+1;
    elseif strncmpi(varargin{j},'nbins',4)
        j = j+1;
        nbins = varargin{j};
    elseif strncmpi(varargin{j},'sd',2)
        j = j+1;
        sx = varargin{j}(1);
        sy = varargin{j}(2);
    elseif strncmpi(varargin{j},'ynormal',5)
        setydir = 1;
    elseif strncmpi(varargin{j},'xrange',5)
        j = j+1;
        xrange = varargin{j};
    elseif strncmpi(varargin{j},'yrange',5)
        j = j+1;
        yrange = varargin{j};
    elseif strncmpi(varargin{j},'sqrtsc',6)
        sqrt_scale = 1; 
    elseif strncmpi(varargin{j},'noplot',6)
        plot_on = 0;
    end
    j = j + 1;
end

%For gaussian Kernel
[gx,gy] = meshgrid(-10:10,-10:10);
if sx == 0
    sx= std(x)/10;
end
if sy == 0
    sy= std(y)/10;
end
G = exp(-(gx).^2/sx - (gy).^2/sy);
G = G./sum(G(:));


dx = diff(xrange)./nbins(1);
dy = diff(yrange)./nbins(2);
z = zeros(nbins(1)+1,nbins(2)+1);
% ignore spikes where both are set to 0
xtra(1) = diff(xrange)./(nbins(1)-1);
xtra(2) = diff(yrange)./(nbins(2)-1);
xi = 1+floor(nbins(1) * (x-xrange(1))/diff(xrange));
yi = 1+floor(nbins(2) * (y-yrange(1))/diff(yrange));

bad_pts = find(xi < 1 | xi > nbins(2) | yi < 1 | yi > nbins(1) | isnan(xi) | isnan(yi));
xi(bad_pts) = [];
yi(bad_pts) = [];

%something like histc might make this a lot faster...
for j = 1:length(xi)
    z(yi(j),xi(j)) = z(yi(j),xi(j))+1;
end
if sx+sy > 0
    z = conv2(z,G,'same');
end
if sqrt_scale == 1
    z = sqrt(z);
end
if plot_on
handles(1) = imagesc([xrange(1) xrange(2)+xtra(1)],[yrange(1) yrange(2)+xtra(2)],z);
if setydir == 1
    set(gca,'ydir','normal');
end
end
[details.x,details.y] = meshgrid(xrange(1):dx:xrange(2),yrange(1):dy:yrange(2));
details.z = z;
