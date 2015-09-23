clear all
close all

N = 1e5;
d = 50; %number of spatial dims
sf = 0.1; %spatial freq
bandwidth = 1/2; %spatial freq bandwidth
env_sigma = 1/sf*bandwidth; %spatial envelope SD

xax = 1:d; xax = xax- mean(xax);

%make gabor filter
gabor_filt = exp(-xax.^2/(2*env_sigma^2)).*exp(sqrt(-1)*2*pi*xax*sf);
stim_filts = [real(gabor_filt(:)) imag(gabor_filt(:))]; %quad pair

%%
dds = 67; %bar density
X{1} = generate_RLS_stim(N,d,dds,1); %make RLS stimulus (1d white noise)
filt_outs = X{1}*stim_filts; %compute output of each filter
filt_outs = zscore(filt_outs);

%% apply upstream NLs
gamma = 1.5; %exponent
theta = 0.5; %offset

%make an energy-type model out of pairs of rectified power filters
pos_filts_out = (filt_outs - theta).^gamma;
pos_filts_out(filt_outs < theta) = 0;

neg_filts_out = (-filt_outs - theta).^gamma;
neg_filts_out(filt_outs > -theta) = 0;

%% compute overall model output
sp_beta = 1; %spiking NL beta
sp_theta = 0; %spiking NL theta

gen_signal = sum(pos_filts_out,2) + sum(neg_filts_out,2); %compute generating signal
gen_signal = zscore(gen_signal); %zscore normalize

rate = log(1 + exp(sp_beta*gen_signal + sp_theta));
Robs = poissrnd(rate);
%% add in a column vector of ones for estimating filter offsets
X{1} = [X{1} ones(N,1)];
d = d + 1; 

%%
n_subs = 4;
STIM_PARAMS = NIM.create_stim_params(d);
init_filts_mat = randn(d,n_subs)/d; %initialize fitler coefs with gaussian noise
init_filts = cell(n_subs,1);
for ii = 1:n_subs
    init_filts{ii} = init_filts_mat(:,ii);
    init_filts{ii}(end) = 0; %make initial offset terms 0
end

optim_params.optTol = 1e-5; %[minFunc] termination tolerance on first order optimality (max(abs(grad))
optim_params.progTol = 1e-8; %[minFunc] termination tolerance on function/parameter values
mod = NIM(STIM_PARAMS,'rectpow',ones(1,n_subs),'init_filts',init_filts);
mod = mod.fit_filters(Robs,X,[],'optim_params',optim_params);

% for ii = 1:n_subs
%    mod.subunits(ii).NLparam_con = [1 Inf]; 
% end
%%
% n_steps = 25;
% optim_params.MaxIter = 500;
% for ii = 1:n_steps
% mod = mod.fit_upstreamNLs(Robs,X,[],'optim_params',optim_params,'silent',1);
% mod = mod.fit_filters(Robs,X,[],'optim_params',optim_params,'silent',1);
% [mod.subunits(1).NLparams]
% end
