clear all; 
%%%%%%% INITIAL PROCESSING OF Tim Data
dt   = 19.9920031987/1000; 
sdim = 64; 

stimfnames = {'ns1-64p-50h-2m-mn125-ctNAT','ns1-64p-50h-2m-mnNAT-ct035', 'ns1-64p-50h-2m-mnNAT-ctNAT', 'ns2-64p-50h-2m-mn125-ctNAT', 'ns2-64p-50h-2m-mnNAT-ct035', 'ns2-64p-50h-2m-mnNAT-ctNAT', 'ns3-64p-50h-2m-mnNAT-ctNAT'}
stimfnames79  = cellfun(@(x)x(1:26),stimfnames,'UniformOutput',0); 
nstims = length(stimfnames79); 
cellfun(@disp,stimfnames79)

datdir = '/media/NTlab_data1/Data/blanche/rec_76/';
cd(datdir)

cd /media/NTlab_data1/Data/blanche/stimuli/
rstimparts= cell(nstims,1); 
for istim =1:nstims; 
	tstimfile = stimfnames79{istim};
	fid = fopen(tstimfile,'r');
    rawdat = fread(fid,inf,'uchar'); 
    fclose(fid);	
	rawdat = 2*((rawdat/255)-0.5); 
	rstimparts{istim} = rawdat; 
end; 

%% reshape, normalize & downsample the stimulus [-1,1]
dsr=0.5; 
stimps76  = cellfun(@(x)reshape(x,[sdim sdim 6000]),rstimparts,'UniformOutput',0); 
dstimps76 = cellfun(@(x)imresize(x,dsr),stimps76,'UniformOutput',0); 
dstimps76  = cellfun(@(x)reshape(x,[(sdim*dsr)^2 6000]),dstimps76,'UniformOutput',0); 

ulen = 32;
[XX,YY] = meshgrid(1:sdim*dsr);
uset = find(abs(XX) <= ulen & abs(YY) <= ulen);
dstimps76  = cellfun(@(x) x(uset,:)',dstimps76,'UniformOutput',0); 

%% Read timing data

fid  = fopen([datdir,'/stim_data/76_-_track_7c_tracking.din'],'r');
frt  = fread(fid,inf,'uint64'); fclose(fid); ntot=length(frt); 

rfrt = reshape(frt,2,ntot/2)'; rfrt(1:10,:)
ts0  = frt(1,1)
tims      = (rfrt(:,1)-ts0)/1000; 
fonsets   = reshape(tims,4,size(tims,1)/4)'; 

seqstarts = [0; fonsets( (1:(nstims-1))*6000+1 , 1)]/1000;


%% read the spikes
spikefiles = dir([datdir,'/spk_data/*.spk']); 
ncells = length(spikefiles);
rspks76 = cell(ncells,1); 
% Read spike times
for ispfile = 1:ncells
  fid = fopen([datdir,'/spk_data/',spikefiles(ispfile).name],'r');
  rawts = fread(fid,inf,'uint64'); fclose(fid);
  ts = (rawts - ts0)/1000;
  rspks76{ispfile} = ts(ts > 0)/1000;
end


%% total protocol len is nstim*stimlen + (n-1)*breaklen (5sec)
n_seqs = length(stimfnames);
seqlen = 6000;  
NT = seqlen*n_seqs;

Robs_mat = nan(NT,ncells);
for icell=1:ncells;
	tspks = rspks76{icell}; 
	for iep =1:n_seqs; 
		tstart = seqstarts(iep); tend=tstart+seqlen; 
        cur_rel_spike_times = tspks( tspks>tstart & tspks < tend) - tstart;
        uset = find(cur_rel_spike_times >= 0 & cur_rel_spike_times <= seqlen*dt);
        cur_spk_bins = floor(cur_rel_spike_times(uset)/dt);
        cur_binned_spikes = hist(cur_spk_bins,1:seqlen);
        cur_block_inds = (iep-1)*seqlen + (1:seqlen);
        Robs_mat(cur_block_inds,icell) = cur_binned_spikes;
	end; 
end; 

%%
sdim_ds = ulen;
flen = 5;
stim_params = NMMcreate_stim_params([flen sdim_ds sdim_ds],dt);

selstim   = [];
cur_means = [];
cur_vars = [];
X = [];
% stimmat = [];
block_vec = [];
for icond=1:n_seqs
    cur_stim = dstimps76{icond};
    cur_stim = bsxfun(@minus,cur_stim,mean(cur_stim));
    cur_means = [cur_means; mean(cur_stim(:))];
    cur_vars = [cur_vars; var(cur_stim(:))];
%     stimmat = [stimmat; reshape(cur_stim,[size(cur_stim,1) sdim_ds sdim_ds])];
    X = [X; create_time_embedding(reshape(cur_stim,[size(cur_stim,1) sdim_ds sdim_ds]),stim_params)];
    block_vec = [block_vec; ones(seqlen,1)*icond];
end

Xblock = zeros(NT,n_seqs);
for i = 1:n_seqs
    cur_set = find(block_vec==i);
    Xblock(cur_set,i) = 1;
end

cur_X{1} = X;
cur_X{2} = Xblock;

%%

base_lambda_d2XT = 500;
base_lambda_L1 = 10;
stim_params(2) = NMMcreate_stim_params([n_seqs],dt);

n_squared_filts = 1;
mod_signs = ones(1,n_squared_filts+2);
NL_types = [{'lin'} repmat({'quad'},1,n_squared_filts) {'lin'}];
init_d2XT = [500*ones(n_squared_filts+1,1); 0];
init_reg_params = NMMcreate_reg_params('lambda_d2XT',init_d2XT);
init_Xtargs = [ones(n_squared_filts+1,1); 2];
silent = 1;

for tcell = 1:ncells;
    fprintf('Fitting model for Cell %d of %d\n',tcell,ncells);
    gqm1 = NMMinitialize_model(stim_params,mod_signs,NL_types,init_reg_params,init_Xtargs);
    gqm1 = NMMfit_filters(gqm1,Robs_mat(:,tcell),cur_X,[],[],silent);
    [LL, penLL, pred_rate, G, gint] = NMMmodel_eval(gqm1,Robs_mat(:,tcell),cur_X);
    gqm1 = NMMadjust_regularization(gqm1,find(init_Xtargs==1),'lambda_d2XT',base_lambda_d2XT./var(gint)');
    gqm1 = NMMadjust_regularization(gqm1,find(init_Xtargs==1),'lambda_L1',base_lambda_L1./std(gint)');
    gqm1 = NMMfit_filters(gqm1,Robs_mat(:,tcell),cur_X,[],[],silent);
    init_mod_fit(tcell) = gqm1;
     
    init_mod_fits_withspkNL(tcell) = NMMfit_logexp_spkNL(init_mod_fit(tcell),Robs_mat(:,tcell),cur_X);
    [LL,~,~,~,~,~,nullLL] = NMMmodel_eval(init_mod_fits_withspkNL(tcell),Robs_mat(:,tcell),cur_X);
    mod_nullLL(tcell) = nullLL;
    init_mod_LLimp(tcell) = LL-nullLL;
    
end

%%
cd ~/Analysis/blanche/rec76/
save init_mods init_mod* stim_params init_reg_params base_* mod_*

%%
sp_dx = 0.0497*2;
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
Ix = speye(sdim_ds);
Iy = speye(sdim_ds);
xshift_mat = cell(length(x_shifts),1);
for xx = 1:length(x_shifts)
    tempx = spdiags( ones(sdim_ds,1), -x_shifts(xx), sdim_ds, sdim_ds);
    xshift_mat{xx} = kron(Iy, kron(tempx, It));
end
yshift_mat = cell(length(y_shifts),1);
for yy = 1:length(y_shifts)
    tempy = spdiags( ones(sdim_ds,1), -y_shifts(yy), sdim_ds, sdim_ds);
    yshift_mat{yy} = kron(tempy, kron(Ix, It));
end
shift_mat = cell(n_shifts,1);
for xx = 1:n_shifts
    shift_mat{xx} = xshift_mat{SH_inds(xx,1)}*yshift_mat{SH_inds(xx,2)};
end

drift_dsf = 10;

%overall prior on shifts
eps_prior_sigma = 0.25; %0.125 start
leps_prior = -sum((SH*sp_dx).^2,2)/(2*eps_prior_sigma^2);
leps_prior = bsxfun(@minus,leps_prior,logsumexp(leps_prior)); %normalize
lA_tflip = repmat(leps_prior',n_shifts,1);

cdist = squareform(pdist(SH*sp_dx));
deps_sigma = 0.05; %.01
lA = -cdist.^2/(2*deps_sigma^2);
lA = bsxfun(@minus,lA,logsumexp(lA,2)); %normalize

%% PREPROCESS MODEL COMPONENTS
klen = flen*sdim_ds^2;
filt_bank = zeros(ncells,klen,n_squared_filts+1);
lin_kerns = nan(ncells,n_seqs);
mod_spkNL_params = nan(ncells,3);
for ss = 1:ncells
    cur_Xtargs = [init_mod_fits_withspkNL(ss).mods(:).Xtarget];
    cur_k = [init_mod_fits_withspkNL(ss).mods(cur_Xtargs == 1).filtK];
    n_used_filts = size(cur_k,2);
    filt_bank(ss,:,1:n_used_filts) = cur_k;
    mod_spkNL_params(ss,:) = init_mod_fits_withspkNL(ss).spk_NL_params;
    lin_kerns(ss,:) = init_mod_fits_withspkNL(ss).mods(cur_Xtargs == 2).filtK;
end
filt_bank = permute(filt_bank,[2 1 3]);

%indicator predictions
block_out = Xblock*lin_kerns';

%%
loo = 13;
tr_cells = 1:ncells;
tr_cells(loo) = [];
frame_LLs = nan(NT,n_shifts);
for xx = 1:n_shifts
    fprintf('Shift %d of %d\n',xx,n_shifts);
    cur_stim_shift = X*shift_mat{xx};
    
    %outputs of stimulus models at current X-matrix shift
    gfuns = ones(NT,ncells);
    gfuns = bsxfun(@times,gfuns,mod_spkNL_params(:,1)');
    gfuns = gfuns + cur_stim_shift*squeeze(filt_bank(:,:,1));
    for ff = 2:(n_squared_filts+1)
        gfuns = gfuns + mod_signs(ff)*(cur_stim_shift*squeeze(filt_bank(:,:,ff))).^2;
    end
    
    %add contributions from extra lin kernels
    gfuns = gfuns + block_out;
    
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
    
    frame_LLs(:,xx) = squeeze(nansum(Robs_mat(:,tr_cells).*log(pred_rate(:,tr_cells)) - ...
        pred_rate(:,tr_cells),2));
end

%%
cur_NT = NT;
use_prior = zeros(NT,1);
use_prior((0:n_seqs-1)*seqlen+1) = 1;

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

ninds = find(~isnan(drift_post_mean(:,1)));
drift_post_mean = interp1(ninds,drift_post_mean(ninds,:),1:NT);
drift_post_std = interp1(ninds,drift_post_std(ninds,:),1:NT);
drift_post_mean(isnan(drift_post_mean(nn,1)),:) = 0;

%%
[~,drift_post_map] = max(gamma,[],2); 
X_shifted = X;
for i=1:NT
    X_shifted(i,:) = X(i,:)*shift_mat{drift_post_map(i)};
end
cur_X{1} = X_shifted;

%%
tcell = 13;

base_mod = init_mod_fit(tcell);
base_mod = NMMfit_filters(base_mod,Robs_mat(:,tcell),cur_X,[],[],silent);
base_mod_withspkNL = NMMfit_logexp_spkNL(base_mod,Robs_mat(:,tcell),cur_X);
[LL,~,~,~,~,~,nullLL] = NMMmodel_eval(base_mod_withspkNL,Robs_mat(:,tcell),cur_X);

