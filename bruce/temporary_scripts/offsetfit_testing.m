clear all
close all

N = 1e5;
d = 50;
sf = 0.1; %spatial freq
bandwidth = 1/2; %spatial freq bandwidth

env_sigma = 1/sf*bandwidth; %spatial envelope SD

xax = 1:d; xax = xax- mean(xax);

%make gabor filter
gabor_filt = exp(-xax.^2/(2*env_sigma^2)).*exp(sqrt(-1)*2*pi*xax*sf);
stim_filts = [real(gabor_filt(:)) imag(gabor_filt(:))];

%%
dds = 67;
X{1} = generate_RLS_stim(N,d,dds,1);
filt_outs = X{1}*stim_filts;
filt_outs = zscore(filt_outs);

%%
gamma = 1.5;
theta = 0.5;

pos_filts_out = (filt_outs - theta).^gamma;
pos_filts_out(filt_outs < theta) = 0;

neg_filts_out = (-filt_outs - theta).^gamma;
neg_filts_out(filt_outs > -theta) = 0;

%%
sp_beta = 1;
sp_theta = 0;

gen_signal = sum(pos_filts_out,2) + sum(neg_filts_out,2);
gen_signal = zscore(gen_signal);
% gen_signal = sum(

rate = log(1 + exp(sp_beta*gen_signal + sp_theta));
Robs = poissrnd(rate);
%%
X{1} = [X{1} ones(N,1)];
d = d + 1;
%%
n_subs = 4;
STIM_PARAMS = NIM.create_stim_params(d);
% init_filts_mat = randn(d,n_subs)/d; %initialize fitler coefs with gaussian noise

% init_filts = cell(n_subs,1);
% for ii = 1:n_subs
% %     init_filts{ii} = init_filts_mat(:,ii);
% %     init_filts{ii}(end) = 1;
% end
% for ii = 1:n_subs
%     init_filts{ii} = cat(1,init_filts{ii},0);
% end
mod = NIM(STIM_PARAMS,'rectpow',[ones(1,n_subs)]);
mod = mod.fit_filters(Robs,X);

%%
n_steps = 50;
optim_params.MaxIter = 500;
for ii = 1:n_steps
mod = mod.fit_upstreamNLs(Robs,X,[],'optim_params',optim_params,'silent',1);
mod = mod.fit_filters(Robs,X,[],'optim_params',optim_params,'silent',1);
[mod.subunits(1).NLparams]
end
