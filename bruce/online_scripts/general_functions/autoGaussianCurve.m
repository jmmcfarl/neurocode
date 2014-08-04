function [results] = autoGaussianCurve(xi,zi)

%% Adapted from Mineault's autoGaussianCurve function

    sz = size(zi);
    
    %Verify that the grid is regular
    if any(any(abs(diff(xi,2,2)) >=1e3*eps))
        error('xi is not a regular grid');
    end
    
    if any(size(zi)~=size(xi))
        error('xi and zi are not the same size');
    end
    
    if numel(zi) ~= length(zi)
        error('zi must be 1-dimensional');
    end
    
    xi = xi(:);
    boundx = [min(xi),max(xi)];
    
    %Find a minimum sigma based on number of elements, range of x and y
    rgx = diff(boundx);
    minsigma = rgx/sz(2)/5;
    maxsigma = rgx*.3;
   
    sigmas = exp(log(minsigma):.3:log(maxsigma));
    
    rgx = -ceil(sz(2)/2)+1:sz(2)/2;
    
    res = zeros(length(sigmas),5);
    
    %Run through all the different values for sigma
    for ii = 1:length(sigmas)
        thefiltx = exp(-rgx.^2/2/sigmas(ii));
        %convolve zi with gaussian filters and find the maximum response
        %(matched filtering)
        zi2 = reflectconv(zi,thefiltx);
        if opts.positive
            [~,pos] = max(zi2(:));
        else
            [~,pos] = max(abs(zi2(:)));
        end
        x0 = xi(pos);
        
        %Determine the residual error for the optimal x for this sigma
        G = exp(-((xi-x0).^2)/2/sigmas(ii)^2);
        X = [G,ones(length(G),1)];
        ps = X\zi(:);
        res(ii,:) = [sum((zi(:) - X*ps).^2),ps(:)',x0,sigmas(ii)];
    end
    
    %Find sigma with corresponding least error
    [~,optsigma] = min(res(:,1));
    params0 = res(optsigma,2:end)';
    
    varnames = {'a','b','x0','sigmax'};
%     lb = [-Inf,-Inf,xi(1),minsigma /1.01]';
     lb = [-Inf,0,xi(1),minsigma /1.01]';
   ub = [ Inf, Inf,xi(end),maxsigma + .01]';
    if opts.positive
        lb(1) = 0;
    end
    results = doFinalOptimization(@pointgaussian,xi(:),zi(:),params0,lb,ub,true(4,1),varnames,false(4,1),opts);
    results.G = reshape(results.G,size(xi));
end

%Convolution with reflected edge handling
function A = reflectconv(A,f)
    A = conv(A,f,'same');
end

function [a] = pointgaussian(x,xdata)
    a = exp(-.5*(xdata-x(1)).^2*(1./x(2).^2));
end