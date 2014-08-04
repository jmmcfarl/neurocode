function smooth_mat = gauss_smooth2d(mat,sigma)

temp_axis = -2*sigma:2*sigma;
[X,Y] = meshgrid(temp_axis,temp_axis);
gaussfun = exp(-1/(2*sigma^2)*(X.^2+Y.^2));
gaussfun = gaussfun/sum(gaussfun(:));

smooth_mat = conv2(mat,gaussfun,'same');