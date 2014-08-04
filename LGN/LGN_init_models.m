clear all
cd ~/Data/LGN
load('LGNrecsNAT.mat')
load('NATstim.mat');
load('stimres.mat');
%%
expt_num = 1;
block_num = 1;
[NT,SDIM] = size(NATstim);
slen = sqrt(SDIM);

n_sus = size(LGNrecsNAT{expt_num}.SUspks,1);
Robs_SU = nan(NT,n_sus);

for ii = 1:n_sus
   cur_spk_times = LGNrecsNAT{expt_num}.SUspks{ii,block_num};
   cur_spk_bins = floor(cur_spk_times/stimres);
   cur_spk_binned = hist(cur_spk_bins,1:NT);
   Robs_SU(:,ii) = cur_spk_binned;
end

n_mus = length(LGNrecsNAT{expt_num}.MUspks{1});
Robs_MU = nan(NT,n_mus);

for ii = 1:n_mus
   cur_spk_times = LGNrecsNAT{expt_num}.MUspks{block_num}{ii};
   cur_spk_bins = floor(cur_spk_times/stimres);
   cur_spk_binned = hist(cur_spk_bins,1:NT);
   Robs_MU(:,ii) = cur_spk_binned;
end

Robs_tot = [Robs_SU Robs_MU];
n_cells = n_sus + n_mus;
%%
flen = 10;
stim_params = NMMcreate_stim_params([flen slen slen],stimres);
norm_stim = bsxfun(@minus,NATstim,mean(NATstim));
norm_stim = bsxfun(@rdivide,norm_stim,std(norm_stim));
X = create_time_embedding(norm_stim,stim_params);
%%
base_lambda_d2XT = 1000;
base_lambda_L1 = 25;

n_squared_filts = 0;
mod_signs = ones(1,n_squared_filts+1);
% mod_signs(2) = -1;
% NL_types = [{'lin'} repmat({'quad'},1,n_squared_filts)];
NL_types = {'lin'};
init_d2XT = [500*ones(n_squared_filts+1,1)];
init_reg_params = NMMcreate_reg_params('lambda_d2XT',init_d2XT);
init_Xtargs = [ones(n_squared_filts+1,1)];
silent = 1;


for tcell = 1:n_cells;
    fprintf('Fitting model for Cell %d of %d\n',tcell,n_cells);
    gqm1 = NMMinitialize_model(stim_params,mod_signs,NL_types,init_reg_params,init_Xtargs);
    gqm1 = NMMfit_filters(gqm1,Robs_tot(:,tcell),X,[],[],silent);
    [LL, penLL, pred_rate, G, gint] = NMMmodel_eval(gqm1,Robs_tot(:,tcell),X);
    gqm1 = NMMadjust_regularization(gqm1,find(init_Xtargs==1),'lambda_d2XT',base_lambda_d2XT./var(gint)');
    gqm1 = NMMadjust_regularization(gqm1,find(init_Xtargs==1),'lambda_L1',base_lambda_L1./std(gint)');
    gqm1 = NMMfit_filters(gqm1,Robs_tot(:,tcell),X,[],[],silent);
    init_mod_fit(tcell) = gqm1;
     
    init_mod_fits_withspkNL(tcell) = NMMfit_logexp_spkNL(init_mod_fit(tcell),Robs_tot(:,tcell),X);
    [LL,~,~,~,~,~,nullLL] = NMMmodel_eval(init_mod_fits_withspkNL(tcell),Robs_tot(:,tcell),X);
    mod_nullLL(tcell) = nullLL;
    init_mod_LLimp(tcell) = LL-nullLL;
    
end

%%
save init_mods init_mod* stim_params init_reg_params base_* mod_*

%%
% sp_dx = 0.0497*2;
sp_dx = 0.1;
max_shift = 5;
dshift = 1;

x_shifts = -max_shift:dshift:max_shift;
y_shifts = -max_shift:dshift:max_shift;

[Xsh,Ysh] = meshgrid(x_shifts,y_shifts);
[Xsh_inds,Ysh_inds] = meshgrid(1:length(x_shifts),1:length(y_shifts));
SH = [Xsh(:) Ysh(:)];
SH_inds = [Xsh_inds(:) Ysh_inds(:)];
n_shifts = size(SH,1);
zero_frame = find(SH(:,1) == 0 & SH(:,2) == 0);

% temp_dx = [-2*max_shift:dshift:2*max_shift];
% shift_dx_edges = [temp_dx*sp_dx-dshift*sp_dx/2 temp_dx(end)*sp_dx+dshift*sp_dx/2];
% ovshift_dx_edges = [shifts*sp_dx-dshift*sp_dx/2 shifts(end)*sp_dx+dshift*sp_dx/2];

%generate shift matrices. Must be applied to the stimulus (not the filters)
It = speye(flen);
Ix = speye(slen);
Iy = speye(slen);
xshift_mat = cell(length(x_shifts),1);
for xx = 1:length(x_shifts)
    tempx = spdiags( ones(slen,1), -x_shifts(xx), slen, slen);
    xshift_mat{xx} = kron(Iy, kron(tempx, It));
end
yshift_mat = cell(length(y_shifts),1);
for yy = 1:length(y_shifts)
    tempy = spdiags( ones(slen,1), -y_shifts(yy), slen, slen);
    yshift_mat{yy} = kron(tempy, kron(Ix, It));
end
shift_mat = cell(n_shifts,1);
for xx = 1:n_shifts
    shift_mat{xx} = xshift_mat{SH_inds(xx,1)}*yshift_mat{SH_inds(xx,2)};
end

drift_dsf = 5;

%overall prior on shifts
eps_prior_sigma = 0.25; %0.125 start
leps_prior = -sum((SH*sp_dx).^2,2)/(2*eps_prior_sigma^2);
leps_prior = bsxfun(@minus,leps_prior,logsumexp(leps_prior)); %normalize
lA_tflip = repmat(leps_prior',n_shifts,1);

cdist = squareform(pdist(SH*sp_dx));
deps_sigma = 0.04; %.01
lA = -cdist.^2/(2*deps_sigma^2);
lA = bsxfun(@minus,lA,logsumexp(lA,2)); %normalize

%% PREPROCESS MODEL COMPONENTS
klen = flen*slen^2;
filt_bank = zeros(n_cells,klen,n_squared_filts+1);
mod_spkNL_params = nan(n_cells,3);
for ss = 1:n_cells
    cur_Xtargs = [init_mod_fits_withspkNL(ss).mods(:).Xtarget];
    cur_k = [init_mod_fits_withspkNL(ss).mods(cur_Xtargs == 1).filtK];
    n_used_filts = size(cur_k,2);
    filt_bank(ss,:,1:n_used_filts) = cur_k;
    mod_spkNL_params(ss,:) = init_mod_fits_withspkNL(ss).spk_NL_params;
end
filt_bank = permute(filt_bank,[2 1 3]);


%%
loo = 10;
tr_cells = 1:n_cells;
tr_cells(loo) = [];
% tr_cells = 1:26;
frame_LLs = nan(NT,n_shifts);
for xx = 1:n_shifts
    fprintf('Shift %d of %d\n',xx,n_shifts);
    cur_stim_shift = X*shift_mat{xx};
    
    %outputs of stimulus models at current X-matrix shift
    gfuns = ones(NT,n_cells);
    gfuns = bsxfun(@times,gfuns,mod_spkNL_params(:,1)');
    gfuns = gfuns + cur_stim_shift*squeeze(filt_bank(:,:,1));
    for ff = 2:(n_squared_filts+1)
        gfuns = gfuns + mod_signs(ff)*(cur_stim_shift*squeeze(filt_bank(:,:,ff))).^2;
    end
    
%     %add contributions from extra lin kernels
%     gfuns = gfuns + block_out;
    
    %incorporate beta
    gfuns = bsxfun(@times,gfuns,mod_spkNL_params(:,2)');
    
    %handle numerical overflow with log(1+exp)
    too_large = gfuns > 50;
    pred_rate = log(1+exp(gfuns));
    pred_rate(too_large) = gfuns(too_large);
    
    %incorporate alpha
    pred_rate = bsxfun(@times,pred_rate,mod_spkNL_params(:,3)');
    
    %enforce min predicted rate
    pred_rate(pred_rate < 1e-50) = 1e-50;
    
    frame_LLs(:,xx) = squeeze(nansum(Robs_tot(:,tr_cells).*log(pred_rate(:,tr_cells)) - ...
        pred_rate(:,tr_cells),2));
end

%%
cur_NT = NT;
use_prior = zeros(NT,1);
use_prior(1) = 1;

nt_pts = ceil(cur_NT/drift_dsf);
tset_inds = 1+floor((0:(cur_NT-1))/drift_dsf);
talpha=zeros(nt_pts,n_shifts);
tbeta = zeros(nt_pts,n_shifts);

tpt_loc = ceil(1:drift_dsf:nt_pts*drift_dsf);
tpt_loc(end) = cur_NT;

cur_LL_set = frame_LLs(1:cur_NT,:);
if mod(cur_NT,drift_dsf) ~= 0
    dangling_pts = nt_pts*drift_dsf-cur_NT;
    cur_LL_set = cat(1,cur_LL_set,zeros(dangling_pts,n_shifts));
end
cur_LL_set = reshape(cur_LL_set,[drift_dsf nt_pts n_shifts]);
cur_LL_set = squeeze(sum(cur_LL_set,1));

lalpha=zeros(nt_pts,n_shifts);
lbeta = zeros(nt_pts,n_shifts);
%compute rescaled forward messages
lalpha(1,:) = leps_prior' + frame_LLs(1,:);
for t=2:nt_pts
    if mod(t,100) == 0
        fprintf('Forward time %d of %d\n',t,nt_pts);
    end
    if use_prior(t)
        cur_lA = lA_tflip;
    else
        cur_lA = lA;
    end
    lalpha(t,:) = logmulexp(lalpha(t-1,:),cur_lA) + cur_LL_set(t,:);
end

%compute rescaled backward messages
lbeta(nt_pts,:)=log(ones(1,n_shifts));
for t=nt_pts-1:-1:1
    if mod(t,100) == 0
        fprintf('Backward time %d of %d\n',t,nt_pts);
    end
    if use_prior(t+1)
        cur_lA = lA_tflip;
    else
        cur_lA = lA;
    end
    lf1 = lbeta(t+1,:) + cur_LL_set(t+1,:);
    lbeta(t,:) = logmulexp(lf1,cur_lA');
end
lgamma= lalpha + lbeta;
lgamma = bsxfun(@minus,lgamma,logsumexp(lgamma,2));

gamma = exp(lgamma);

int_gamma = interp1(tpt_loc,gamma,1:cur_NT);

clear drift_post*
drift_post_mean(:,1) = sum(bsxfun(@times,gamma,SH(:,1)'),2);
drift_post_mean(:,2) = sum(bsxfun(@times,gamma,SH(:,2)'),2);

cur_xdiff = bsxfun(@minus,drift_post_mean(:,1)',SH(:,1)).^2';
drift_post_std(:,1) = sqrt(sum(cur_xdiff.*gamma,2));
cur_ydiff = bsxfun(@minus,drift_post_mean(:,2)',SH(:,2)).^2';
drift_post_std(:,2) = sqrt(sum(cur_ydiff.*gamma,2));
clear cur_xdiff cur_ydiff

% ninds = find(~isnan(drift_post_mean(:,1)));
% drift_post_mean = interp1(ninds,drift_post_mean(ninds,:),1:NT);
% drift_post_std = interp1(ninds,drift_post_std(ninds,:),1:NT);
% drift_post_mean(isnan(drift_post_mean(nn,1)),:) = 0;

%%
[~,drift_post_map] = max(int_gamma,[],2); 
X_shifted = X;
for i=1:NT
    X_shifted(i,:) = X(i,:)*shift_mat{drift_post_map(i)};
end
cur_X{1} = X_shifted;

%%
tcell = 1;

base_mod = init_mod_fit(tcell);
base_mod = NMMfit_filters(base_mod,Robs_tot(:,tcell),cur_X,[],[],silent);
base_mod_withspkNL = NMMfit_logexp_spkNL(base_mod,Robs_tot(:,tcell),cur_X);
[LL,~,~,~,~,~,nullLL] = NMMmodel_eval(base_mod_withspkNL,Robs_tot(:,tcell),cur_X);


