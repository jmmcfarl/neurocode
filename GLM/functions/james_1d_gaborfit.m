function [gabor_fit,gabor_fitvals] = james_1d_gaborfit(z, sdim, hyperparams, init_guess)


max_lambda = 500;
min_lambda = 2;

options.Display = 'off';
options.Algorithm = 'interior-point';
x = (1:sdim)';

%%
if length(unique(z)) < 2 %check if the z values are well-behaved
    disp('Function not well-behaved, aborting fit!');
    gabor_fitvals = zeros(size(z));
        gabor_fit.b = nan;
        gabor_fit.x0 = nan;
        gabor_fit.lambda = nan;
        gabor_fit.psi = nan;
        gabor_fit.sigma = nan;
        gabor_fit.a = nan;
        gabor_fit.ssq_err = nan;
        gabor_fit.ssq = nan;
        gabor_fit.fit_type = 'fail';   
else    
    if nargin < 4
        %% create initial parameters for the gaussian fit
        znmat = conv(z,epanechikov(3),'same');
        zgauss = znmat(:);
        zn = abs(zgauss)/sum(abs(zgauss)); %treat z as a distribution to estimate first and second moments
        % x0 = x'*zn;
        % y0 = y'*zn;
        [~,mloc] = max(zn);
        x0 = x(mloc);        
        sig = sqrt((x'-x0).^2*zn);
        K0 = [max(z) x0 sig median(z)];        
        %fit 1d gaussian model
        A = []; b = []; Aeq = []; beq = [];
        lb = [-1*max(abs(z)) (1-sdim/2) 1 min(z)];
        ub = [1*max(abs(z)) (sdim+sdim/2) sdim max(z)];
        if isempty(hyperparams)
            [gxfit,gauss_ssq] = fmincon(@(K) sum((z-gauss_1dfun(x,K)).^2),K0,A,b,Aeq,beq,lb,ub,[],options);
        else
            [gxfit,gauss_ssq] = fmincon(@(K) gauss_posterior(x,z,K,hyperparams),K0,A,b,Aeq,beq,lb,ub,[],options);
        end
        
    else %if initial parameters are specified
        if strcmp(init_guess.fit_type,'gauss')
            K0 = [init_guess.b init_guess.x0 init_guess.sigma init_guess.a];
            %fit 1d gaussian model
            A = []; b = []; Aeq = []; beq = [];
            lb = [-1*max(abs(z)) (1-sdim/2) 1 min(z)];
            ub = [1*max(abs(z)) (sdim+sdim/2) sdim max(z)];
            if isempty(hyperparams)
                [gxfit,gauss_ssq] = fmincon(@(K) sum((z-gauss_1dfun(x,K)).^2),K0,A,b,Aeq,beq,lb,ub,[],options);
            else
                [gxfit,gauss_ssq] = fmincon(@(K) gauss_posterior(x,z,K,hyperparams),K0,A,b,Aeq,beq,lb,ub,[],options);
            end            
        else strcmp(init_guess.fit_type,'gabor')
            gauss_ssq = inf;
        end
    end
    
    %% Now initialize parameters of gabor fit
    if nargin < 4 %if no initial guess specified
        b0 = abs(gxfit(1)); %initial estimate of scale
                
        %estimate fundamental non-zero frequency
        zm   = z-mean(z);
        D1   = fftshift(fft(zm)); %compute FFT
        rabs = abs(D1); %amplitude spectrum
        tfreqs   = [(sdim/2):-1:-((sdim/2)-1)]'/sdim; %spatial frequency axis (nyquist to minus nyquist)        
        [xmax,imax] = extrema(rabs); %find local maxima
        if ~isempty(xmax)
            f0 = tfreqs(imax(find(xmax > xmax(1)*0.1 & tfreqs(imax)~=0,1,'first'))); 
            if ~isempty(f0)
                lambda0 = 1/f0;
            else
                lambda0 = max_lambda;
            end
        else
            lambda0 = max_lambda;
        end
            
        x0 = [2*gxfit(1) gxfit(2) lambda0 0 gxfit(3) median(z)];        
        % now fit 2d gabor function
        A = [];
        b = [];
        Aeq = [];
        beq = [];
        lb = [0 (1-sdim/2) min_lambda 0 0.2 min(z)];
        ub = [max(abs(z)) (sdim+sdim/2) max_lambda 2*pi 25 max(z)];
        
        %now find best starting guess for psi by a uniform scan
        cur_x0 = x0;
        psi_range = linspace(0,2*pi,4);
        psi_fit = zeros(size(psi_range));
        %     fprintf('Trying Psi Initialization (of %d): ',length(psi_range));
        for i = 1:length(psi_range)
            %         fprintf(' %d ',i);
            cur_x0(6) = psi_range(i);
            if isempty(hyperparams)
                [cur_xfit(i,:),cur_ssq(i)] = fmincon(@(K) sum((z-gabor_1dfun(x,K)).^2),cur_x0,A,b,Aeq,beq,lb,ub,[],options);
            else
                [cur_xfit(i,:),cur_ssq(i)] = fmincon(@(K) gabor_posterior(x,z,K,hyperparams),cur_x0,A,b,Aeq,beq,lb,ub,[],options);
            end
        end
        %     fprintf('\n');
        [~,mind] = min(cur_ssq);
        xfit = cur_xfit(mind,:);
        ssq = cur_ssq(mind);
        
    else %if there is a specified initial guess
        if strcmp(init_guess.fit_type,'gabor') %fit a gabor using these params, only if the specified function type is gabor
            x0(1) = init_guess.b;
            x0(2) = init_guess.x0;
            x0(3) = init_guess.lambda;
            x0(4) = init_guess.psi;
            x0(5) = init_guess.sigma;
            x0(6) = init_guess.a;
            
            %% now fit 2d gabor function
            A = [];
            b = [];
            Aeq = [];
            beq = [];
            lb = [0 (1-sdim/2) min_lambda 0 0.2 min(z)];
            ub = [max(abs(z)) (sdim+sdim/2) max_lambda 2*pi 25 max(z)];
            
            %now find best starting guess for psi by a uniform scan
            cur_x0 = x0;
            psi_range = linspace(0,2*pi,4);
            psi_fit = zeros(size(psi_range));
            for i = 1:length(psi_range)
                cur_x0(4) = psi_range(i);
                if isempty(hyperparams)
                    [cur_xfit(i,:),cur_ssq(i)] = fmincon(@(K) sum((z-gabor_1dfun(x,K)).^2),cur_x0,A,b,Aeq,beq,lb,ub,[],options);
                else
                    [cur_xfit(i,:),cur_ssq(i)] = fmincon(@(K) gabor_posterior(x,z,K,hyperparams),cur_x0,A,b,Aeq,beq,lb,ub,[],options);
                end
            end
            [~,mind] = min(cur_ssq);
            xfit = cur_xfit(mind,:);
            ssq = cur_ssq(mind);
        else
            ssq = Inf;
        end
    end
    
    %create output structure
    if ssq > gauss_ssq        
        gabor_fit.b = gxfit(1);
        gabor_fit.x0 = gxfit(2);
        gabor_fit.lambda = nan;
        gabor_fit.psi = nan;
        gabor_fit.sigma = gxfit(5);
        gabor_fit.a = gxfit(6);
        gabor_fit.ssq_err = gauss_ssq;
        gabor_fit.ssq = sum(z.^2);
        gabor_fit.fit_type = 'gauss';
         if nargout > 1
            gabor_fitvals = gauss_2dfun(x,y,gxfit); %final gabor fit
        end
   else
        gabor_fit.b = xfit(1);
        gabor_fit.x0 = xfit(2);
        gabor_fit.lambda = xfit(3);
        gabor_fit.psi = xfit(4);
        gabor_fit.sigma = xfit(5);
        gabor_fit.a = xfit(6);
        gabor_fit.ssq_err = ssq;
        gabor_fit.ssq = sum(z.^2);
        gabor_fit.fit_type = 'gabor';
        if nargout > 1
            gabor_fitvals = gabor_1dfun(x,xfit); %final gabor fit
        end      
    end
end
end

%% helpr function, gaussian nLP
function nLP = gauss_posterior(x,z,K,hyperparams)
pred = gauss_1dfun(x,K);
nLL = sum((pred-z).^2)/hyperparams.noise_var;
nLP = nLL + -log(gampdf(K(3),hyperparams.sigma_shape,hyperparams.sigma_scale));
end

%% helpr function, gabor nLP
function nLP = gabor_posterior(x,z,K,hyperparams)
pred = gabor_1dfun(x,K);
nLL = sum((pred-z).^2)/hyperparams.noise_var;
sigma = K(5); 
nLP = nLL - log(gampdf(sigma,hyperparams.sigma_shape,hyperparams.sigma_scale));
end

%% helper function, 2d gaussian
function z = gauss_1dfun(x,K)
%assume no zero-offset, could easily add
b = K(1); %scale
x0 = K(2); %x-center
sig = K(3); %variance in x-dir
a = K(4); %offset
z = b*exp(-(x-x0).^2/(2*sig^2))+a;
end

%% helper function, defines Gabor function
function z = gabor_1dfun(x,K)
%x-direction defined as direction of spatial oscillation, must also be
%oriented with a principle axis of the gaussian envelope...  Could
%generalize beyond this.
b = K(1); %scale
x0 = K(2); %x-center
lambda = K(3); %spatial wavelength (in x-direction)
psi = K(4); %phase
sigma = K(5); %variance of gaussian, in x-direction
a = K(6);

xp = (x-x0);
z = b*exp(-(xp.^2/2/sigma^2)) .* (cos(2*pi*xp/lambda+psi))+a;
end

function [efilt] = epanechikov(len)
%T2's epanechikov function
refilt = 1- linspace(-1,1,len+2).^2;
efilt  = refilt(2:(len+1));
efilt  = efilt/sum(efilt);
end
