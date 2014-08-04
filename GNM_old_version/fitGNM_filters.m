function [gnm] = fitGNM_filters(gnm,X,spkbs,disp_type,its_per_round,optTol,progTol,silent)

%INPUTS:
% gnm: model structure
% X: X-matrix
% spkbs: vector of spike bin indices
% disp_type: display type: 'tots','full','min','none'
% its_per_round: number of iterations to fit model before checking for user input to continue

if nargin < 6
    optTol = 1e-4;
end
if nargin < 7
    progTol = 1e-6;
end
if nargin < 8
    silent = 1;
end
nmods = length(gnm.mods);

if nargin < 5 || isempty(its_per_round)
    its_per_round = 5000; %default
    check_user_in = 0;
else
    check_user_in = 1;
end

[stimlen,klen]   = size(X); %stimulus dimensions

sdim = gnm.stim_params.sdim;
fsdim = gnm.stim_params.fsdim;
flen = klen/fsdim;

% fprintf('INITIAL MODEL\n');
[ll0, ll0p] = getLL_GNM(gnm,X,spkbs,disp_type);
gnm.LP = ll0p; gnm.LL = ll0;
gnm.LP_seq = [gnm.LP_seq ll0p];
gnm.LL_seq = [gnm.LL_seq ll0];

%compute binned spike vector
if strcmp(gnm.spk_nl,'gauss')
    Robs = spkbs;
    if size(X) ~= length(spkbs)
        error('Input data vector not spike bins for Gaussian likelihood')
    end
else
    Robs = zeros(1,stimlen);
    ftable = tabulate(spkbs);
    Robs(ftable(:,1)) = ftable(:,2);
end

% silent = 1;
eflag = 0;
n_opt_rounds = 0;
while eflag == 0
    n_opt_rounds = n_opt_rounds + 1;
%     fprintf('Optimization round %d\n',n_opt_rounds);
    [gnm eflag] = fitGNM_filters_helper(gnm,X,Robs,its_per_round,'lbgs',silent,optTol,progTol);
    if check_user_in == 1
        R = input('Continue? (y/n)','s');
        if R == 'n'
            break
        end
    else
        break
    end
    eflag = 0;
end

%%
%if any subunits have all coefs to 0, move the 0 to the subunit weight
for imod = 1:nmods
    if all(gnm.mods(imod).k == 0)
%         gnm.mods(imod).w = 0;
    end
end

[ll0, ll0p] = getLL_GNM(gnm,X,spkbs,disp_type);
gnm.LP_seq(end) = ll0p;
gnm.LL_seq = [gnm.LL_seq ll0];

gnm.LL = ll0; gnm.LP = ll0p;

for imod =1:length(gnm.mods);
    if max(abs(gnm.mods(imod).k)) > 0
        cur_inputs = X*gnm.mods(imod).k;
        [foy,fox] = ksdensity(cur_inputs);
        gnm.mods(imod).fox = fox;
        gnm.mods(imod).foy = foy;
        
        if ~strcmp(gnm.spk_nl,'gauss')
        [foys,foxs] = ksdensity(cur_inputs(spkbs));
        gnm.mods(imod).foxs = foxs;
        gnm.mods(imod).foys = foys;
        end
    else
        gnm.mods(imod).fox = linspace(-3.1,3.1,10);
        gnm.mods(imod).foy = zeros(1,10);
        gnm.mods(imod).foxs = gnm.mods(imod).fox;
        gnm.mods(imod).foys = gnm.mods(imod).foy;
    end
end

% disp('Model Fitting Complete');


