function [x] = downsample2x(x)
    x = convn(x,[.25,.5,.25],'same');
    x = x(1:2:end-1);
end