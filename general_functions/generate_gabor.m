function Gfun = generate_gabor(XX,YY,x0,y0,theta,sigma,phase,lambda,gamma)

if nargin < 7 || isempty(phase)
    phase = rand*2*pi;
end
if nargin < 8 || isempty(lambda)
   lambda = sigma*2; 
end
if nargin < 9 || isempty(gamma)
    gamma = 1;
end
    
xp = (XX-x0)*cos(theta)+(YY-y0)*sin(theta);
yp = -(XX-x0)*sin(theta)+(YY-y0)*cos(theta);

Gfun = exp(-(xp.^2/2/sigma^2 + gamma*yp.^2/2/sigma^2)) .* cos(2*pi*xp/lambda+phase);
Gfun = Gfun - mean(Gfun(:)); %remove DC component

end