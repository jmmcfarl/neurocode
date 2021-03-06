function [gnm,old_scales,new_scales,init_grad,fin_grad] = fitGNM_filter_scales(gnm,X,spkbs,disp_type)

%INPUTS:
% gnm: model structure
% X: X-matrix
% spkbs: vector of spike bin indices
% disp_type: display type: 'tots','full','min','none'
% its_per_round: number of iterations to fit model before checking for user input to continue

% optTol = 1e-4;
% progTol = 1e-6;
nmods = length(gnm.mods);
options.Display = disp_type;

%stimulus parameters
[stimlen klen] = size(X);
kern_t = klen/gnm.stim_params.fsdim;

sdim = gnm.stim_params.sdim;
fsdim = gnm.stim_params.fsdim;
flen = klen/fsdim;

%compute binned spike vector
Robs = zeros(1,stimlen);
ftable = tabulate(spkbs);
Robs(ftable(:,1)) = ftable(:,2);


%compute initial fit parameters
initial_params = [];
for m = 1:nmods
    cur_kern = gnm.mods(m).k';
    old_scales(m) = norm(cur_kern);
    gnm.mods(m).k = gnm.mods(m).k/old_scales(m);
    initial_params = [initial_params old_scales(m)]; %add coefs to initial param vector
    
    NLx = gnm.mods(m).nlx;
    NL = gnm.mods(m).nly;
    if strcmp(gnm.mods(m).nltype,'uncon')
        %compute derivative of non-linearity
        fpr = zeros(1,length(NLx)-1);
        for n = 1:length(fpr)
            fpr(n) = (NL(n+1)-NL(n))/(NLx(n+1)-NLx(n));
        end
        fprimes{m} = fpr;
    else
        fprimes{m} = [];
    end
end
initial_params(end+1) = gnm.const; %add constant offset term to params
[init_LL,init_grad] = fitGNM_filter_scales_internal(initial_params', Robs, X, gnm,fprimes);
[params LL exitflag] = minFunc( @(K) fitGNM_filter_scales_internal(K, Robs, X, gnm,fprimes), initial_params', options );
[fin_LL,fin_grad] = fitGNM_filter_scales_internal(params, Robs, X, gnm,fprimes);

% Reassign variables
spk_theta.const = params(end);
for n = 1:nmods
    new_scales(n) = params(n);
    cur_kern = gnm.mods(n).k;
    gnm.mods(n).k = cur_kern(:)*new_scales(n);
end

%%

disp('Model Fitting Complete');


