function [ output_args ] = plot2drfmat(rf,varargin)
% USAGE: [ output_args ] = plot2drfmat(rf,varargin)
%  plot matrix, eventually with zlims
colormap(gray); 
imagesc(rf,varargin{:}); 
% pcolor(rf); shading interp
set(gca,'XTick',[])
set(gca,'YTick',[])
set(gca,'YDir','normal'); 
end

