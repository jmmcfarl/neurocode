function grating = gen_square_grating(imSize,lambda,theta,phase)


X = 1:imSize;                           % X is a vector from 1 to imageSize
X0 = (X / imSize) - .5;                 % rescale X -> -.5 to .5
freq = imSize/lambda;     

[Xm Ym] = meshgrid(X0, X0);             % 2D matrices

Xt = Xm * cos(theta);                % compute proportion of Xm for given orientation
Yt = Ym * sin(theta);                % compute proportion of Ym for given orientation
XYt = Xt + Yt;                      % sum X and Y components
XYf = XYt * freq * 2*pi;                % convert to radians and scale by frequency
grating = sin( XYf + phase);                   % make 2D sinewave
