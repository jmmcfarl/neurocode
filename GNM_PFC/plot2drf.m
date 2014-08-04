function RF = plot2drf(k,sdim,zrange,to_transpose)
% USAGE: RF = plot2drf(k,flen)
% Plots the RF and gives a matrix containing it back

if (nargin==2)
    zrange=[-0.5,0.5];
end;

flen = length(k)/sdim;
RF   = reshape(k,flen,sdim);
if to_transpose
    RF = RF';
end
colormap(gray(256));
imagesc(RF,zrange);
if ~to_transpose
    set(gca,'XTick',[])
    set(gca,'YTick',[])
end
end

