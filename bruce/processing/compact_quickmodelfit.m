clear all

global Expt_name bar_ori use_LOOXV monk_name rec_type

Expt_name = 'M009';
monk_name = 'jbe';
use_LOOXV = 1; %[0 no LOOXV; 1 SU LOOXV; 2 all LOOXV]
bar_ori = 0; %bar orientation to use (only for UA recs)

Expt_num = str2num(Expt_name(2:end));

data_dir = ['~/Data/bruce/' Expt_name];
if ~exist(data_dir,'dir')
    data_dir = ['/media/NTlab_data3/Data/bruce/' Expt_name];
end
if ~exist(data_dir,'dir');
    error('Couldnt find data directory');
end
Edata_file = strcat(data_dir,'/',monk_name,Expt_name,'Expts');
load(Edata_file);

%is this a laminar probe or utah array rec?
if strcmp(Expts{1}.Header.DataType,'GridData 96')
    rec_type = 'UA';
elseif strcmp(Expts{1}.Header.DataType,'Spike2')
    rec_type = 'LP';
end

data_name = sprintf('%s/packaged_data_ori%d',data_dir,bar_ori);
fprintf('Loading %s\n',data_name);
load(data_name);

%%
if strcmp(rec_type,'LP')
    switch Expt_num
        case 266
            cor_ori = 80;
        case 270
            cor_ori = 60;
        case 275
            cor_ori = 135;
        case 277
            cor_ori = 70;
        case 281
            cor_ori = 140;
        case 287
            cor_ori = 90;
        case 289
            cor_ori = 160;
        case 294
            cor_ori = 40;
        case 296
            cor_ori = 45;
        case 297
            cor_ori = [0 90];
        case 309
            cor_ori = 120;
        case 5
            cor_ori = 50;
        case 9
            cor_ori = 0;
    end
else
    cor_ori = [0 90];
end
if ~ismember(bar_ori,cor_ori)
    error('this bar ori is not allowed');
end

cd(data_dir);

if strcmp(rec_type,'LP')
    n_probes = 24;
elseif strcmp(rec_type,'UA')
    n_probes = 96;
end

load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data

%%
if strcmp(rec_type,'LP')
    use_nPix = 32;
elseif strcmp(rec_type,'UA')
    use_nPix = 16;
end

if ~isnan(params.rpt_seeds)
    xv_type = 'rpt';
    xv_frac = nan;
else
    xv_type = 'uni';
    xv_frac = 0;
end

flen = 12;

%these recs have larger bar widths
if ismember(Expt_num,[287 289 294])
    use_nPix = 15;
elseif ismember(Expt_num,[296 297])
    use_nPix = 22;
end

stim_fs = 100; %in Hz
dt = 0.01;
Fr = 1;

full_nPix=36;
switch Expt_num
    case 270
        full_nPix=32;
    case  287
        full_nPix = 22;
    case 289
        full_nPix = 22;
    case 294
        full_nPix = 20;
end

if full_nPix ~= params.full_nPix
    fprintf('Using full_nPix in params struct\n');
    full_nPix = params.full_nPix;
end
if use_nPix > full_nPix
    fprintf('Using npix == full_nPix\n');
    use_nPix = full_nPix;
end

%exclude data at beginning and end of each trial
trial_dur = 4;

stim_params = NIMcreate_stim_params([flen full_nPix],dt);

%%
base_sp_dx = mode(expt_data.expt_dw);
if length(unique(expt_data.expt_dw)) > 1
    fprintf('Warning, multiple different dws detected, using %.3f\n',base_sp_dx);
end
sp_dx = base_sp_dx/params.scale_fac; %model dx in deg

%%
all_stim_mat = decompressTernNoise(stimComp);

%%
all_Xmat = create_time_embedding(all_stim_mat,stim_params);

%% select submatrix with central pixels

%pick out the predictors (of the time-embedded stimulus-matrices)
%corresponding to the central pixels we're using in the models
[Xinds,~] = meshgrid(1:full_nPix,1:flen);
buffer_pix = floor((full_nPix - use_nPix)/2);
cur_use_pix = (1:use_nPix) + buffer_pix;
use_kInds = find(ismember(Xinds(:),cur_use_pix)); %use this set of dimensions in the Xmatrix

all_Xmat = all_Xmat(:,use_kInds);
%% BIN SPIKES FOR MU AND SU
all_binned_mua = spikes_int82double(spike_data.binned_mua);
all_binned_sua = spikes_int82double(spike_data.binned_sua);
Clust_data = spike_data.Clust_data;
su_probes = Clust_data.SU_probes;
SU_numbers = Clust_data.SU_numbers;

%%
NT = length(used_inds);
fullNT = size(all_binned_mua,1);
n_trials = length(time_data.trial_flip_ids);
n_blocks = length(expt_data.used_blocks);

all_t_axis = time_data.t_axis;
trial_start_inds = [1+time_data.trial_flip_inds];
trial_end_inds = [time_data.trial_flip_inds(2:end); fullNT];
all_trialvec = nan(fullNT,1);
for ii = 1:n_trials
    all_trialvec(trial_start_inds(ii):trial_end_inds(ii)) = time_data.trial_flip_ids(ii);
end

block_start_inds = [1+time_data.block_flip_inds];
block_end_inds = [time_data.block_flip_inds(2:end); fullNT];
all_blockvec = nan(fullNT,1);
for ii = 1:n_blocks
    all_blockvec(block_start_inds(ii):block_end_inds(ii)) = time_data.block_flip_ids(ii);
end
Xblock = zeros(fullNT,n_blocks);
for i = 1:n_blocks
    cur_set = find(all_blockvec==i);
    Xblock(cur_set,i) = 1;
end

%% Create set of TR and XV trials
use_trials = unique(all_trialvec(used_inds));
nuse_trials = length(use_trials);

if strcmp(xv_type,'rpt')
    xv_trials = find(ismember([trial_data(:).se],params.rpt_seeds));
    n_xv_trials = length(xv_trials);
else
    n_xv_trials = round(xv_frac*nuse_trials);
    xv_trials = randperm(nuse_trials);
    xv_trials(n_xv_trials+1:end) = [];
    xv_trials = use_trials(xv_trials);
end
tr_trials = setdiff(use_trials,xv_trials);
n_tr_trials = length(tr_trials);
fprintf('Initializing models with %d training trials and %d xval trials\n',n_tr_trials,n_xv_trials);

tr_inds = find(ismember(all_trialvec(used_inds),tr_trials));
xv_inds = find(ismember(all_trialvec(used_inds),xv_trials));

%%
% make Robs_mat
tot_sus = size(all_binned_sua,2);
Robs_mat = nan(length(used_inds),n_probes + tot_sus);
for ss = 1:size(Robs_mat,2)
    if ss > n_probes
        Robs_mat(:,ss) = all_binned_sua(used_inds,ss-n_probes);
    else
        Robs_mat(:,ss) = all_binned_mua(used_inds,ss);
    end
end

%%
tr_X{1} = abs(all_Xmat(used_inds,:));
tr_X{2} = all_Xmat(used_inds,:);

%%
init_stim_params(1) = NMMcreate_stim_params([flen use_nPix],dt);
init_stim_params(2) = NMMcreate_stim_params([flen use_nPix],dt);
silent = 1;
base_lambda_d2XT = 5;
base_lambda_L1 = 0;
init_optim_p.optTol = 1e-5; init_optim_p.progTol = 1e-9;

n_stim_filts = 2;
mod_signs = ones(1,n_stim_filts);
% NL_types = [{'lin'} repmat({'quad'},1,n_stim_filts-1)];
NL_types = [{'lin' 'lin'}];
init_d2XT = [ones(n_stim_filts,1)];
init_L2 = [zeros(n_stim_filts,1);];
init_reg_params = NMMcreate_reg_params('lambda_d2XT',init_d2XT);
init_Xtargs = [1; 2];

for ss = 1:size(Robs_mat,2);
    fprintf('Fitting model for MU %d of %d\n',ss,size(Robs_mat,2));
    Robs = Robs_mat(:,ss);
    uinds = find(~isnan(Robs));
    
    if ~isempty(uinds)
    gqm1 = NMMinitialize_model(init_stim_params,mod_signs,NL_types,init_reg_params,init_Xtargs);
    gqm1 = NMMfit_filters(gqm1,Robs,tr_X,[],uinds,silent,init_optim_p);
    
    [~, ~, ~, ~, gint] = NMMeval_model(gqm1,Robs,tr_X,[],uinds);
    gqm1 = NMMadjust_regularization(gqm1,find(init_Xtargs==1),'lambda_d2XT',base_lambda_d2XT./var(gint)');
    gqm1 = NMMadjust_regularization(gqm1,find(init_Xtargs==1),'lambda_L1',base_lambda_L1./std(gint)');
    gqm1 = NMMfit_filters(gqm1,Robs,tr_X,[],uinds,silent);
    
    all_mod_fits(ss) = gqm1;
    end
%     [LL, penLL, pred_rate, G, gint, fgint, nullLL] = NMMmodel_eval(all_mod_fits(ss),Robs,tr_X);
%     all_mod_LLimp(ss) = (LL-nullLL)/log(2);
    
end

%%
for ss = 1:size(Robs_mat,2)
    if ~isempty(all_mod_fits(ss).mods)
    all_filts(ss,:,:) = reshape(all_mod_fits(ss).mods(1).filtK,flen,use_nPix);
    end
end
temp_profiles = squeeze(std(all_filts,[],3));
space_profiles = squeeze(std(all_filts,[],2));
[~,best_lags] = max(temp_profiles,[],2);
[~,best_pos] = max(space_profiles,[],2);
for ss = 1:size(Robs_mat,2)
    best_temp_profiles(ss,:) = squeeze(all_filts(ss,:,best_pos(ss)));
    best_space_profiles(ss,:) = squeeze(all_filts(ss,best_lags(ss),:));
end
figure
imagesc(0:dt*1e3:(flen-1)*dt*1e3,1:size(Robs_mat,2),best_temp_profiles);
xlabel('Time (ms)','fontsize',12);
ylabel('Probe','fontsize',12);

figure
imagesc(1:use_nPix,1:size(Robs_mat,2),best_space_profiles);
xlabel('Bar position','fontsize',12);
ylabel('Probe','fontsize',12);