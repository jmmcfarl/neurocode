function Gfun = generate_gabor_sequence(XX,YY,x0,y0,theta,sigma,phase,lambda)

n_gabors = length(x0);

if nargin < 7 || isempty(phase)
    phase = rand(n_gabors,1)*2*pi;
end
if nargin < 8 || isempty(lambda)
   lambda = sigma*2; 
end
% if nargin < 9 || isempty(gamma)
%     gamma = ones(n_gabors,1);
% end
    
XX = bsxfun(@minus,XX,reshape(x0,[1 1 n_gabors]));
YY = bsxfun(@minus,YY,reshape(y0,[1 1 n_gabors]));

xp = bsxfun(@times,XX,reshape(cos(theta),[1 1 n_gabors])) + ...
    bsxfun(@times,YY,reshape(sin(theta),[1 1 n_gabors]));
yp = bsxfun(@times,-XX,reshape(sin(theta),[1 1 n_gabors])) + ...
    bsxfun(@times,YY,reshape(cos(theta),[1 1 n_gabors]));

Gfun = exp(-(bsxfun(@rdivide,xp.^2,reshape(2*sigma.^2,[1 1 n_gabors])) + ...
    bsxfun(@rdivide,yp.^2,reshape(2*sigma.^2,[1 1 n_gabors]))));
Gfun = Gfun.*cos(bsxfun(@plus,2*pi*bsxfun(@rdivide,xp,reshape(lambda,[1 1 n_gabors])),reshape(phase,[1 1 n_gabors]))); 

end