clear all
close all
addpath('~/NIMclass');

N = 1e4;
d = 30; %number of spatial dims
sf = 0.1; %spatial freq
bandwidth = 1/2; %spatial freq bandwidth
env_sigma = 1/sf*bandwidth; %spatial envelope SD

xax = 1:d; xax = xax- mean(xax);

%make gabor filter
gabor_filt = exp(-xax.^2/(2*env_sigma^2)).*exp(sqrt(-1)*2*pi*xax*sf);
stim_filts = [real(gabor_filt(:)) imag(gabor_filt(:))]; %quad pair
gabor_filt_ps = exp(-xax.^2/(2*env_sigma^2)).*exp(sqrt(-1)*(2*pi*xax*sf + pi/4));
stim_filts_ps = [real(gabor_filt(:)) imag(gabor_filt(:))]; %quad pair

%%
dds = 67; %bar density
X{1} = generate_RLS_stim(N,d,dds,1); %make RLS stimulus (1d white noise)
filt_outs = X{1}*stim_filts; %compute output of each filter
filt_outs = zscore(filt_outs);
filt_outs_ps = X{1}*stim_filts_ps; %compute output of each filter
filt_outs_ps = zscore(filt_outs_ps);

% apply upstream NLs
gamma = 2; %exponent
theta = 0; %offset

%make an energy-type model out of pairs of rectified power filters
pos_filts_out = (filt_outs - theta).^gamma;
pos_filts_out(filt_outs < theta) = 0;

neg_filts_out = (-filt_outs - theta).^gamma;
neg_filts_out(filt_outs > -theta) = 0;
% neg_filts_out = (-filt_outs_ps - theta).^gamma;
% neg_filts_out(filt_outs_ps > -theta) = 0;

% compute overall model output
sp_beta = 1; %spiking NL beta
sp_theta = 0; %spiking NL theta

gen_signal = sum(pos_filts_out,2) + sum(neg_filts_out,2); %compute generating signal
% gen_signal = sum(pos_filts_out,2); %compute generating signal
% gen_signal = sum(pos_filts_out,2); %compute generating signal
gen_signal = zscore(gen_signal); %zscore normalize

rate = log(1 + exp(sp_beta*gen_signal + sp_theta));
Robs = poissrnd(rate);

%%
n_subs = 2;
STIM_PARAMS = NIM.create_stim_params(d);
clear optim_params
optim_params.optTol = 1e-6; %[minFunc] termination tolerance on first order optimality (max(abs(grad))
optim_params.progTol = 1e-10; %[minFunc] termination tolerance on function/parameter values
optim_params.silent = 1;
mod = NIM(STIM_PARAMS,'rectpow',[1 1 1 1],'d2t',10);
mod = mod.fit_filters(Robs,X,[],'optim_params',optim_params,'fit_offsets',0,'silent',0);
% mod = NIM(STIM_PARAMS,'quad',[1 1],'d2t',0);
% mod = mod.fit_filters(Robs,X,[],'optim_params',optim_params,'fit_offsets',0,'silent',0);

%%
n_subs = length(mod.subunits);
[filt_SE,hessMat] = mod.compute_filter_SEs(Robs,X);
inv_hess = pinv(hessMat); %pseudo-inverse
allK_SE = sqrt(diag(inv_hess)); %take sqrt of diagonal component as SE
if any(~isreal(allK_SE))
    fprintf('complex\n\n\n\n\n');
end
close all
f1 = figure(); hold on
cmap = jet(n_subs);
% cmap(2,:) = [1 0 0];
for ii = 1:n_subs
    filt_range = (ii-1)*d + (1:d);
    errorbar(1:d,mod.subunits(ii).filtK,1.96*filt_SE{ii},'color',cmap(ii,:))
end
for ii = 1:n_subs
    plot(1:d,mod.subunits(ii).filtK,'k','linewidth',3); 
end

%%
mod_weights = [mod.subunits(:).weight];
[~,prate,mod_internals] = mod.eval_model(Robs,X);
G_plus_theta = mod_internals.G + mod.spkNL.theta;

F_fd_out = mod.apply_spkNL_deriv(G_plus_theta);
F_sd_out = mod.apply_spkNL_deriv2(G_plus_theta);

hess_mat = zeros(n_subs*d);
for ii = 1:n_subs
    irange = (ii-1)*d + (1:d);
    for jj = ii:n_subs
        jrange = (jj-1)*d + (1:d);
        
        filt_fd_ii = mod.subunits(ii).apply_NL_deriv(mod_internals.gint(:,ii));
        filt_fd_jj = mod.subunits(jj).apply_NL_deriv(mod_internals.gint(:,jj));
        filt_sd_ii = mod.subunits(ii).apply_NL_deriv2(mod_internals.gint(:,ii));
        filt_sd_jj = mod.subunits(jj).apply_NL_deriv2(mod_internals.gint(:,jj));

        if ii == jj
            interior = F_sd_out.*(filt_fd_ii).^2 + F_fd_out.*filt_sd_ii;
        else
            interior = mod_weights(ii)*mod_weights(jj)*F_sd_out.*filt_fd_ii.*filt_fd_jj;
        end
        
        dr2 = bsxfun(@times,X{1},interior);
        dr2 = bsxfun(@times,dr2,(Robs./prate - 1));
        
        dr_ii = bsxfun(@times,X{1},F_fd_out.*filt_fd_ii);
        dr_jj = bsxfun(@times,X{1},F_fd_out.*filt_fd_jj);
        dr_ii = bsxfun(@times,dr_ii,Robs./prate.^2);
        
        cur_hess = dr2'*X{1} + dr_ii'*dr_jj;
        
        hessMat(irange,jrange) = cur_hess;
        if ii ~= jj
           hessMat(jrange,irange) = cur_hess; 
        end
    end
end
inv_hess = pinv(hessMat);
filt_SEs = sqrt(diag(inv_hess))*1.96;


%%
alpha = 1;
beta = 1;
mod_weights = [mod.subunits(:).weight];
[~,prate,mod_internals] = mod.eval_model(Robs,X);
G_plus_theta = mod_internals.G + mod.spkNL.theta;

F_fd_out = alpha*beta*exp(beta*G_plus_theta)./(1 + exp(beta*G_plus_theta));
F_fd_out(beta*G_plus_theta > 50) = alpha*beta;

F_sd_out = alpha*beta^2*exp(beta*G_plus_theta)./(1 + exp(beta*G_plus_theta));
F_sd_out = F_sd_out.*(1 - exp(beta*G_plus_theta)./(1 + exp(beta*G_plus_theta)));
F_sd_out(beta*G_plus_theta > 50) = 0;

hess_mat = zeros(n_subs*d);
for ii = 1:n_subs
    irange = (ii-1)*d + (1:d);
    for jj = ii:n_subs
        jrange = (jj-1)*d + (1:d);
        
        filt_fd_ii = 2*mod_internals.gint(:,ii);
        filt_fd_jj = 2*mod_internals.gint(:,jj);
        filt_sd_ii = 2*ones(N,1);
        filt_sd_jj = 2*ones(N,1);
        
        if ii == jj
            interior = F_sd_out.*(filt_fd_ii).^2 + F_fd_out.*filt_sd_ii;
        else
            interior = mod_weights(ii)*mod_weights(jj)*F_sd_out.*filt_fd_ii.*filt_fd_jj;
        end
        
        dr2 = bsxfun(@times,X{1},interior);
        dr2 = bsxfun(@times,dr2,(Robs./prate - 1));
        
        dr_ii = bsxfun(@times,X{1},F_fd_out.*filt_fd_ii);
        dr_jj = bsxfun(@times,X{1},F_fd_out.*filt_fd_jj);
        dr_ii = bsxfun(@times,dr_ii,Robs./prate.^2);
        
        cur_hess = dr2'*X{1} + dr_ii'*dr_jj;
        
        hessMat(irange,jrange) = cur_hess;
        if ii ~= jj
           hessMat(jrange,irange) = cur_hess; 
        end
    end
end

% filt_fd1 = 2*mod_internals.gint(:,1);
% filt_fd2 = 2*mod_internals.gint(:,2);
% interior = F_sd_out.*filt_fd1.*filt_fd2;
% dr2 = bsxfun(@times,X{1},interior);
% dr2 = bsxfun(@times,dr2,(Robs./prate - 1));
% dr_f1 = bsxfun(@times,X{1},F_fd_out.*filt_fd1);
% dr_f2 = bsxfun(@times,X{1},F_fd_out.*filt_fd2);
% dr_times = bsxfun(@times,dr_f1,Robs./prate.^2);
% 
% hessMat(d + (1:d),(1:d)) = dr2'*X{1} + dr_f2'*dr_times;
% hessMat((1:d),d + (1:d)) = dr2'*X{1} + dr_f2'*dr_times;



% inv_hess = inv(hessMat);
inv_hess = pinv(hessMat);
filt_SEs = sqrt(diag(inv_hess))*1.96;
%%
f1 = figure(); hold on
cmap = jet(n_subs);
for ii = 1:n_subs
    filt_range = (ii-1)*d + (1:d);
    errorbar(1:d,mod.subunits(ii).filtK,filt_SEs(filt_range),'color',cmap(ii,:))
end
for ii = 1:n_subs
    plot(1:d,mod.subunits(ii).filtK,'k','linewidth',3); 
end

%%
mod = mod.fit_filters(Robs,X,[],'optim_params',optim_params,'fit_offsets',1);

%%
optim_params.optimizer = 'fminsearch';
optim_params.MaxIter = 500;
optim_params.silent = 0;
% for ii = 1:n_subs
%    mod.subunits(ii).NLoffset = -mod.subunits(ii).NLoffset; 
% end
mod = mod.fit_upstreamNLs(Robs,X,[],'optim_params',optim_params);
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
%     init_filts{ii}(end) = 0; %make initial offset terms 0
end

clear optim_params
optim_params.optTol = 1e-5; %[minFunc] termination tolerance on first order optimality (max(abs(grad))
optim_params.progTol = 1e-8; %[minFunc] termination tolerance on function/parameter values
% optim_params.silent = 0;
mod = NIM(STIM_PARAMS,'rectpow',ones(1,n_subs),'init_filts',init_filts);
mod = mod.fit_filters(Robs,X,[],'optim_params',optim_params,'silent',0);

% for ii = 1:n_subs
%    mod.subunits(ii).NLparam_con = [1 Inf]; 
% end
%%
silent = 1;
n_steps = 25;
optim_params.MaxIter = 500;
for ii = 1:n_steps
mod = mod.fit_upstreamNLs(Robs,X,[],'optim_params',optim_params,'silent',silent);
mod = mod.fit_filters(Robs,X,[],'optim_params',optim_params,'silent',silent);
[mod.subunits(1).NLparams]
end
