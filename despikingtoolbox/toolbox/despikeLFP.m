function [results] = despikeLFP(y,S,Bs,g,opts)
    %function [results] = despikeLFP(y,S,Bs,g,opts)
    %
    % Removes spikes from wideband signals to obtain a local field potential
    % estimated by the methods of Zanos, Mineault & Pack (2010). 
    %
    % The generative model for the waveform is:
    %  y = w + sum(conv2(S,Bs*phi,'same'),2) + mu + epsilon
    %
    % Meaning:
    % voltage = lfp + sum of spike contributions + offset + noise
    %
    % And: epsilon ~ N(0,sigma^2*I)
    %            w ~ N(0,gamma^2*C'(ifft(g)))
    %
    % Here C is a function that maps a vector onto a circulant matrix,
    % such that:
    %
    % C'(a)b = ifft(fft(a).*fft(b));
    %
    % The function estimates the parameters through MAP and determines the
    % hyperparameters sigma and gamma through evidence optimization. It
    % returns the waveforms phi, the offset mu, and the despiked wideband 
    % signal z defined as:
    %
    %           z = y - sum(conv2(S,Bs*phi,'same'),2) - mu
    %
    % To obtain a despiked LFP, filter z with some low-pass zero phase
    % filter.
    %
    % Parameters:
    %
    % y : the wideband signal to despike as a vector of size nx1
    %     n must be even
    % S : An n x m matrix of zeros and 1 that specifies spike timings. m is 
    %     the number of neurons emitting spikes. The element S(i,j) is 1 if 
    %     a spike from the j'th neuron begins at the i'th time sample
    % Bs: A basis matrix of size q x k that specifies the subspace that the spike
    %     waveforms can take. q is the length of a spike in samples. 
    % g : a nx1 method that embodies the prior assumptions on w
    %
    % opts: a struct with options
    %     .hptol  - tolerance for minimum log change in hyperparameters to 
    %               stop iterating
    %     .cgtol  - tolerance to use in CG iterations, default 1e-6
    %     .auxtol - tolerance for an auxillary equation used to determine 
    %               convergence, default 1e-4
    %     .mltol  - tolerance for minimum increase in log marginal
    %               likelihood to continue iterating, default 1
    %     .maxiter - Max number of iterations, default 10
    %     .maxcgiter - Max number of iterations for conjugate gradients
    %                  (default 7/8*size(Bs,2))
    %     .displaylevel - either 0 (no output)
    %                            1 (some output, default)
    %                            2 (full debug info)
    %     .sigma - the sigma to start with, default 1/3*std(y)
    %     .gamma - the gamma to start with, default 1/10
    % 
    % Returns a struct with keys:
    %
    % z: the despiked wideband signal
    % mu: the offset
    % phi:  the estimated spike waveforms, a k x m matrix
    % sigma, gamma : hyperparameters
    % opts: the options given to the method
    %
    % History:
    %
    % 09/06/2010: v1.0
    %
    % Created by Patrick Mineault
    %
    %See also despikeLFPbyChunks
    
    thestart = tic;
    
    %Sanitize input
    if mod(length(y),2) ~= 0
        error('despikeLFP:oddy','The length of y must be even');
    end
    
    if size(y,2) ~= 1
        error('despikeLFP:wrongsizey','y is of an incorrect size. This method requires y to be of size n x 1, where n is even');
    end
    
    if size(Bs,1) > 500
        warning('despikeLFP:largeBs','The height of Bs is larger than 500. This could take a while.');
    end
    
    Bs = full(Bs);
    
    %Orthogonalize Bs
    [U,S1,V] = svd(Bs);
    S2 = S1;
    S2(S1~=0) = 1;
    Bs = U*S2*V';
    
    %With phi expressed in this basis, phi0 can be expressed in the
    %original basis as phi0 = V*diag(1./diag(S1))*V'*phi
    
    if any(size(g) ~= size(y))
        error('despikeLFP:wrongg','g is of an incorrect size. The correct size is (%dx1)',length(y));
    end
    fftgam = g;
    clear g;
        
    %Defaults for opts
    defaults.hptol = 1e-2;
    defaults.cgtol = 1e-6;
    defaults.auxtol = 1e-4;
    defaults.maxiter = 10;
    defaults.maxcgiter = ceil(7/8*size(Bs,2));
    defaults.displaylevel = 1;
    defaults.mltol = 1;
    defaults.nbasis = 10;
    defaults.minfreq = .5;
    stdy = std(y);
    defaults.sigma = 1/3*stdy;
    defaults.gamma = 1/10/sqrt(max(fftgam))*stdy;
    
    
    
    if nargin < 5
        sgspec = 0;
        opts = defaults;
    else
        sgspec = isfield(opts,'sigma') && isfield(opts,'gamma');
        opts = mergeoptions(opts,defaults);
    end
    
    if opts.displaylevel > 1
        fminlevel = 'Iter';
    else
        fminlevel = 'Off';
    end    

    ffty = fft(y);
    psdz = ffty.*conj(ffty);

    
    if sgspec
        if opts.displaylevel > 1
            fprintf('   Using opts.sigma, opts.gamma\n');
        end
        sigma = opts.sigma;
        gamma = opts.gamma;
        lml0 = marginallikelihood(log([sigma;gamma]),fftgam,psdz);
    else
        %Find a good start value for gamma
        if opts.displaylevel > 1
            fprintf('   Initializing gamma, sigma\n');
        end
        %Optimize sigma, gamma through marginal likelihood
        ps = log([opts.sigma;opts.gamma]);
        fminuncopts = optimset('TolFun',1e-9,'Display',fminlevel,'GradObj','on','typicalx',ps,'DerivativeCheck','off','Hessian','On','MaxIter',5);
        [ps lml0] = fminunc(@(x) marginallikelihood(x,fftgam,psdz),ps,fminuncopts);
        sigma = exp(ps(1));
        gamma = exp(ps(2));
    end
    
    clear psdz;
    
    % Initialize the rest of the parameters
    phi = zeros(size(Bs,2),size(S,2));
    mu = 0;
    
    %Dirac delta at the 1st position
    dirac1 = sparse(1,1,length(y),length(y),1);
    
    fftS = fft(full(S));
    S = sparse(S);
    fftss = zeros(size(fftS));
    fftz = ffty - sum(fftss,2) - mu*dirac1;
    
    iter = 0;
    if opts.displaylevel > 0 
        fprintf('Iteration %2d, sigma = %.3e, gamma = %.3e, marginal likelihood = %10.2f\n',iter,sigma,gamma,-lml0);
    end
    
    %Start the ball rolling
    %The condition for stopping the outer loop is that sigma and gamma are
    %stable
    sigmaold = sigma*1e4;
    gammaold = gamma*1e4;
    
    lmlold = lml0 + 1e3;
    lml = lml0;
    iter = 1;
        
    while (abs(log(sigma)-log(sigmaold)) > opts.hptol ...
        || abs(log(gamma)-log(gammaold)) > opts.hptol) ...
       && -(lml - lmlold) > opts.mltol ...
       && iter <= opts.maxiter
        lmlold = lml;
        
        phiold    = phi;
        gammaold  = gamma;
        sigmaold  = sigma;
        
        %The filtering kernels
        ffth = sigma^2./(fftgam*gamma^2+sigma^2);
        %This is only used for conditioning the phi equation better
        iffth = ifft(1+fftgam*gamma^2/(1e-12+sigma^2));
        iffth = [iffth((size(Bs,1)+1):-1:2);iffth(1:(size(Bs,1)+1))];
        
        %Do a loop on the relative error measured by the aux equation,
        %decreasing the tolerance of conjugate graidents if the error is
        %too large
        CGTOLERANCESTEP = .3;        
        maxerror = 1;
        opts.cgtol = opts.cgtol/CGTOLERANCESTEP;
        
        while maxerror > opts.auxtol
            opts.cgtol = opts.cgtol*CGTOLERANCESTEP;
            
            if opts.displaylevel > 1
                fprintf('   Computing phi\n');
            end
            
            %Compute phi
            for ii = 1:size(S,2)
                fftu = ffty - sum(fftss,2) + fftss(:,ii);
                phi(:,ii) = solvephi(fftu,fftS(:,ii),Bs,ffth,iffth,opts.maxcgiter,opts.cgtol,opts.displaylevel>1);
            end
            clear fftu;

            if opts.displaylevel > 1
                fprintf('   Computing w\n');
            end

            %Compute mu and w
            clear fftss;
            fftss = fftS.*fft(Bs*phi,length(y));
            %The first element of an FFT vector is the sum, hence
            %sum(fftss(1,:))/length(y)
            mu = mean(y) - sum(fftss(1,:))/length(y);
            fftz = ffty - sum(fftss,2) - mu*dirac1;
            %fftw = (1-ffth).*fftz;
            
            %Check for lack of convergence using an auxillary equation
            %fftwe = fftw - (fftw
            wp = zeros(size(phi));
            for ii = 1:size(S,2)
                wp(:,ii) = 1/sum(S(:,ii))*Bs'*sts(S(:,ii),ifft((fftz - fftz(1)*dirac1).*ffth),size(Bs,1));
            end
            maxerror = max(max(abs(Bs*wp)))/stdy;
            
            clear fftw;
            
            if maxerror > opts.auxtol && opts.displaylevel > 0
                opts.maxcgiter = ceil(opts.maxcgiter/sqrt(1-CGTOLERANCESTEP));
                warning('despikeLFP:auxeq','Auxilliary residual too large (%.5f > %.5f); Setting cgtol to %.2e and raising maxiter to %d. Seeing several of these errors (>5) in a row indicates bad conditioning or some bug. Getting one or two is normal.\n',maxerror,opts.auxtol,opts.cgtol/10,opts.maxcgiter);
            end
        end        
        
        %Now refit hyperparameters
        psdz = fftz.*conj(fftz);
        
        ps = [sigma;gamma];
        ps = log(ps);

        if opts.displaylevel > 1
            fprintf('   Refitting gamma, sigma\n');
        end

        %Optimize the marginal likelihood of the model using fminunc
        fminopts = optimset('TolFun',1e-9,'Display',fminlevel,'GradObj','on','typicalx',ps,'DerivativeCheck','Off','Hessian','On','MaxIter',3);
        [ps lml] = fminunc(@(x) marginallikelihood(x,fftgam,psdz),ps,fminopts);

        clear psdz;
        
        sigma = exp(ps(1));
        gamma = exp(ps(2));
        
        if opts.displaylevel > 0
            fprintf('Iteration %2d, sigma = %.3e, gamma = %.3e, marginal likelihood = %10.2f\n',iter,sigma,gamma,-lml);
        end
        
        iter = iter + 1;
        
        if lml > lml0 + opts.mltol
            %Some sort of unforeseen runaway condition happened. 
            %Regression of death
            warning('despikeLFP:regressionofdeath',...
                    ['The marginal likelihood decreased during the last iteration. This should never ' ...
                    'happen. Reverting to the last good values for z, phi etc. which may still be usable. ' ...
                    'If this continues to happen create a minimal example showing the bug by saving a .mat file ' ...
                    'with the necessary variables and a .m file calling this function with the arguments used, ' ...
                    'and send to the author, after first sending an email saying that you encountered this bug ' ...
                    'and that you will be sending an attachment in the next mail. If the file is > 10 MB, upload it ' ...
                    'to a file sharing service and send the link instead.']);
                
            phi    = phiold;
            gamma  = gammaold;
            sigma  = sigmaold;
            
            fftss = fftS.*fft(Bs*phi,length(y));
            mu = mean(y) - fftss(1)/length(y);
            fftz = ffty - sum(fftss,2) - mu*dirac1;
            
            break;
        end
    end
    
    if opts.displaylevel > 0
        fprintf('Done! Optimization took %.3f s\n',toc(thestart));
    end
    
    results.z = ifft(fftz);
    phi0 = V*diag(1./diag(S1))*V'*phi;
    results.phi = phi0;
    results.gamma = gamma;
    results.sigma = sigma;
    results.mu = mu;
    results.opts = opts;
end

%Merge a struct with default options
function [opts] = mergeoptions(opts, defaults)
    names = fieldnames(defaults);
    for ii = 1:length(names)
        if ~isfield(opts,names{ii})
            if isnan(defaults.(names{ii}))
                error('despikeLFP:requiredOption','Required option not set: %s',defaults.(names{ii}));
            end
            opts.(names{ii}) = defaults.(names{ii});
        end
    end
    
    %Check for unrecognized options
    names = fieldnames(opts);
    for ii = 1:length(names)
        if ~isfield(defaults,names{ii})
            warning('despikeLFP:unknownOption','Unrecognized option %s',names{ii});
        end
    end
end

%Compute the marginal likelihood of the model
%params:
%     sigma = exp(params(1));
%     gamma = exp(params(2));
%g: As in the prior of the model, Gamma = gamma^2*C'(F^-1(g))
%psdz: The PSD of the variable z defined in the Appendix
%      e.g. conj(fftz).*fftz
function [ml grad H] = marginallikelihood(params,g,psdz)
    h = psdz;
    sigma = exp(params(1));
    gamma = exp(params(2));
    
    n = length(g);
    
    %{
    %The obvious way
    fftw = (sigma^-2./(sigma^-2+1./(gamma^2*fftgam)).*fftz);
    wmap = ifft(fftw);
    werror = sum((z-wmap).^2);
    hipassw = ifft(fftw./fftgam);
    wpriorerror = wmap'*hipassw;
    %}
    
    %The faster way
    divider1 = 1./(sigma^2+gamma^2*g);
    divider2 = divider1.^2;
    divider3 = divider1.*divider2;
    divider4 = divider2.^2;
    
    hg = h.*g;
    sh2 = sum(h.*divider2);
    shg2 = sum(hg.*divider2);
    
    E1 = 1/n*sigma^2*sh2;
    E2 = 1/n*gamma^2*shg2;
    
    ml = 1/2*sum(log(sigma^2+gamma^2*g));
    ml = ml + 1/2*E1 + 1/2*E2;
    
    %{
    fftw = (sigma^-2./(sigma^-2+1./(gamma^2*g)).*fftz);
    wmap = ifft(fftw);
    werror = sum((ifft(fftz)-wmap).^2);
    hipassw = ifft(fftw./g);
    wpriorerror = wmap'*hipassw;
    %}
    
    if nargout > 1
        grad = zeros(2,1);
        
        %This is over-optimized. See commented out text at bottom for readable version
        s1 = sum(divider1);
        sg1 = sum(g.*divider1);
        sh3 = sum(h.*divider3);
        shg3 = sum(hg.*divider3);
        hgg3 = hg.*g.*divider3;
        shgg3 = sum(hgg3);
        

        dE1dp1 =  4/n*sigma^2*(sh2-sigma^2*sh3);
        dE2dp1 = -4/n*sigma^2*gamma^2*shg3;
        dE1dp2 = -4/n*sigma^2*gamma^2*shg3;
        dE2dp2 =  4/n*gamma^2*(shg2-gamma^2*shgg3);
        
        grad(1) = sigma^2*s1    - E1 + 1/2*dE1dp1 + 1/2*dE2dp1;
        grad(2) = gamma^2*sg1   - E2 + 1/2*dE1dp2 + 1/2*dE2dp2;
        
        if nargout > 2
            
            hgg4 = hgg3.*divider1;
            shgg4 = sum(hgg4);
            shg4 = sum(hg.*divider4);
            sh4 = sum(h.*divider4);
            shggg4 = sum(hgg4.*g);
            s2 = sum(divider2);
            sg2 = sum(g.*divider2);
            sgg2 = sum(g.*g.*divider2);

            
            dE1dp1dp1 =  4/n*sigma^2        *(4*sh2 - 10*sh3*sigma^2 + 6*sigma^4*sh4);
            dE1dp1dp2 =  4/n*sigma^2*gamma^2*(-4*shg3 + 6*sigma^2*shg4); %
            dE1dp2dp2 =  4/n*sigma^2*gamma^2*(-2*shg3 + 6*gamma^2*shgg4); %
            dE2dp1dp1 =  4/n*sigma^2*gamma^2*(-2*shg3 + 6*sigma^2*shg4); %
            dE2dp1dp2 =  4/n*sigma^2*gamma^2*(-4*shg3 + 6*gamma^2*shgg4); %
            dE2dp2dp2 =  4/n*gamma^2        *(4*shg2 - 10*shgg3*gamma^2 + 6*gamma^4*shggg4);
        
            H11 =  2*sigma^2*(s1 - sigma^2*s2)    + 2*E1 - 2*dE1dp1 + 1/2*(dE1dp1dp1+dE2dp1dp1);
            H21 = -2*gamma^2*sigma^2*sg2 - (dE1dp2+dE2dp1)          + 1/2*(dE1dp1dp2+dE2dp1dp2);
            H22 =  2*gamma^2*(sg1 - gamma^2*sgg2) + 2*E2 - 2*dE2dp2 + 1/2*(dE1dp2dp2+dE2dp2dp2);

            H = [H11,H21;H21,H22];
        end

        
        %{
        dE1dp1 =  4/n*sigma^4*        sum(h.*   (1-sigma^2   *divider1).*divider2);
        dE2dp1 = -4/n*sigma^2*gamma^4*sum(h.*g.*divider3);
        dE1dp2 = -4/n*sigma^4*gamma^2*sum(h.*g.*divider3);
        dE2dp2 =  4/n*gamma^4*        sum(h.*g.*(1-gamma^2*g.*divider1).*divider2);
        
        dE1dp1dp1 =  4/n*sigma^4        *sum( 4*h.*   divider2-10*h*      sigma^2.*divider3+6*h*sigma^4.*divider4);
        dE1dp1dp2 =  4/n*sigma^4*gamma^2*sum(-4*h.*g.*divider3+ 6*h.*g*   sigma^2.*divider4);
        dE1dp2dp2 =  4/n*sigma^4*gamma^2*sum(-2*h.*g.*divider3+ 6*h.*g.*g*gamma^2.*divider4);
        dE2dp1dp1 =  4/n*sigma^2*gamma^4*sum(-2*h.*g.*divider3+ 6*h.*g*   sigma^2.*divider4);
        dE2dp1dp2 =  4/n*sigma^2*gamma^4*sum(-4*h.*g.*divider3+ 6*h.*g.*g*gamma^2.*divider4);
        dE2dp2dp2 =  4/n*gamma^4        *sum( 4*h.*g.*divider2-10*h.*g.*g*gamma^2.*divider3+6*gamma^4*h.*g.*g.*g.*divider4);
        
        grad(1) = sigma^2*sum(divider1)    - 1/sigma^2*E1 + 1/(2*sigma^2)*dE1dp1 + 1/(2*gamma^2)*dE2dp1;
        grad(2) = gamma^2*sum(divider1.*g) - 1/gamma^2*E2 + 1/(2*sigma^2)*dE1dp2 + 1/(2*gamma^2)*dE2dp2;
        
        H11 =  2*sigma^2*sum(   divider1.*(1-sigma^2   *divider1)) + 2/sigma^2*E1 - 2/sigma^2*dE1dp1 + 1/2*(1/sigma^2*dE1dp1dp1+1/gamma^2*dE2dp1dp1);
        H21 = -2*gamma^2*sigma^2*sum(g.*divider2) - (1/sigma^2*dE1dp2+1/gamma^2*dE2dp1)                + 1/2*(1/sigma^2*dE1dp1dp2+1/gamma^2*dE2dp1dp2);
        H22 =  2*gamma^2*sum(g.*divider1.*(1-gamma^2*g.*divider1)) + 2/gamma^2*E2 - 2/gamma^2*dE2dp2 + 1/2*(1/sigma^2*dE1dp2dp2+1/gamma^2*dE2dp2dp2);
        
        H = [H11,H21;H21,H22];
        %}
    end
end

%Solve for phi
function phi = solvephi(fftu,ffts,Bs,ffth,iffth,maxcgiter,cgtol,debug)
    %Solve for phi. See paper for equations
    thekern = ifft(conj(ffts).*ffts.*ffth);
    p = ffts(1);
    d = sum(Bs,1);
    symmkern = [thekern((size(Bs,1)+1):-1:2);thekern(1:(size(Bs,1)+1))];
    lhs = Bs'*conv2(Bs,symmkern,'same') - 1/length(fftu)*p^2*ffth(1)*(d'*d);
    
    rhs = ifft(fftu.*ffth.*conj(ffts));
    rhs = rhs - mean(rhs);
    rhs = rhs(1:size(Bs,1));
    rhs = Bs'*rhs;
    
    %iffth is the inverse filter for better conditioning
    precond = @(x) conv2(x,iffth,'same');
    
    %Solve using conjugate gradients instead of \. This way is slightly
    %faster and allows for the tolerance to be controlled
    [phi, flag, relres, iter] = pcg(lhs,rhs,cgtol,maxcgiter,precond);
    
    if flag ~= 0
        if flag == 1
            warning('despikeLFPfirstorder:PCGwarning','pcg (Conjugate gradient descent function) did not converge to desired relative residual (%.1e) within %d iterations. Actual residual %.2e. Raise opts.maxcgiter',cgtol,iter,relres);
        else
            error('despikeLFPfirstorder:PCGerror','pcg (Conjugate gradient descent function) returned with flag %d; type help pcg to interpret flag',flag);
        end
    elseif debug
        fprintf('        pcg converged in %d iterations, with rel residual %.2e\n', iter, relres);
    end
end

%Compute the spike-triggered sum of the stimulus a for a spike train s
%n is the length of the window of the sum
function [thests] = sts(s,a,n)
    a = [a;a(1:n)];
    spktimes = find(s);
    if length(spktimes) ==1
        thests = a(spktimes + (0:(n-1)));
    else
        snipidx = bsxfun(@plus,spktimes,0:(n-1));
        thesnips = a(snipidx);
        thests = sum(thesnips,1)';
    end
end