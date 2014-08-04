clear all
% close all

dir_prefix = '~';
Expt_name = 'M266';
Expt_num = str2num(Expt_name(2:end));
%% SET ANALYSIS PARAMETERS
stim_fs = 100; %in Hz
stim_dt = 0.01;
usfac = 2;
tbspace = 1;
dt = stim_dt/usfac;

min_trial_dur = 0.5;
beg_buffer = 0.3;

flen = round(12*usfac/tbspace);
use_nPix = 30;
% full_nPix = 36;
stim_params = NMMcreate_stim_params([flen use_nPix],dt,usfac,tbspace);
Fr = 1;

n_probes = 24;

%% LOOP OVER EXPTS

data_dir = [dir_prefix '/Data/bruce/' Expt_name];
cd(data_dir);

if Expt_name(1) == 'G'
    load(sprintf('jbe%sExpts.mat',Expt_name));
elseif Expt_name(1) == 'M'
    load(sprintf('lem%sExpts.mat',Expt_name));
end

save_dir = [dir_prefix '/Analysis/bruce/' Expt_name '/sac_mod'];
if ~exist(save_dir,'dir')
    system(['mkdir ' save_dir]);
end

trial_dur = 4;

load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data

%%
include_expts = {'rls.Fa', 'rls.FaXimi','rls.FaXFaXFs'};
for ii = 1:length(Expts)
    if ~isempty(Expts{ii})
        expt_names{ii} = Expts{ii}.Header.expname;
        expt_dds(ii) = Expts{ii}.Stimvals.dd;
        expt_bar_ori(ii) = Expts{ii}.Stimvals.or;
        expt_sac_dir(ii) = mod(Expts{ii}.Stimvals.Fa,180);
        expt_Fr(ii) = Expts{ii}.Stimvals.Fr;
        expt_imback(ii) = isfield(Expts{ii}.Trials,'imi');
        included_type(ii) = any(strcmp(expt_names{ii},include_expts));
    end
end
expt_has_ds(isnan(expt_has_ds)) = 0;
expt_has_ds(expt_has_ds == -1) = 0;
expt_binoc(isnan(expt_binoc)) = 0;
cur_block_set = find(included_type & ~expt_binoc' & expt_Fr == 1);

if strcmp(Expt_name,'M270')
    cur_block_set(cur_block_set == 5) = [];
    cur_block_set(cur_block_set == 19) = [];
end

sim_sac_expts = find(~expt_has_ds(cur_block_set));
imback_gs_expts = find(expt_has_ds(cur_block_set) & expt_imback(cur_block_set)');
grayback_gs_expts = find(expt_has_ds(cur_block_set) & ~expt_imback(cur_block_set)');

%%
if exist('./fin_aclust_data.mat','file')
    load('./fin_aclust_data.mat');
    [n_aclust_expts,n_aclust_probes] = size(autoclust);
else
    disp('No fin_aclust_data found.');
    autoclust = [];
    n_aclust_expts = 0; n_aclust_probes = 0;
end
aclust_probenums = [autoclust(cur_block_set(1),:).probenum];
aclust_CellList = aclust_CellList(cur_block_set,:);
su_probes = find(any(aclust_CellList > 0));
su_ids = unique(aclust_CellList(aclust_CellList > 0));

best_corresp_su_probe = nan(length(su_ids),1);
for ii = 1:length(su_ids)
    cnts = sum(aclust_CellList == su_ids(ii));
    [~,best_corresp_su_probe(ii)] = max(cnts);
end
%% COMPUTE TRIAL DATA
cd(data_dir);

fprintf('Computing prep data\n');
trial_cnt = 0;
cur_toffset = 0;

all_stim_times = [];
% all_Xmat = [];
all_stim_mat = [];
all_t_axis = [];
all_t_bin_edges = [];
all_tsince_start = [];
all_blockvec = [];
all_trialvec = [];
all_trial_Se = [];
all_trial_blk = [];
all_trial_start_times = [];
all_trial_end_times = [];
all_bin_edge_pts = [];
all_spk_times = cell(n_probes,1);
all_clust_ids = cell(length(su_probes),1);
trial_toffset = zeros(length(cur_block_set),1);
for ee = 1:length(cur_block_set);
    if ismember(ee,grayback_gs_expts)
        fprintf('Block %d of %d; grayback GS, ori:%d\n',ee,length(cur_block_set),expt_bar_ori(cur_block_set(ee)));
    elseif ismember(ee,imback_gs_expts)
        fprintf('Block %d of %d; imback GS, ori:%d\n',ee,length(cur_block_set),expt_bar_ori(cur_block_set(ee)));
    elseif ismember(ee,sim_sac_expts)
        fprintf('Block %d of %d; SimSac, ori:%d\n',ee,length(cur_block_set),expt_bar_ori(cur_block_set(ee)));
    else
        fprintf('Block %d of %d;  UNMATCHED EXPT TYPE\n',ee,length(cur_block_set));
    end
    cur_block = cur_block_set(ee);
    fname = sprintf('Expt%dClusterTimesDetails.mat',cur_block);
    load(fname);
    for cc = 1:n_probes
        all_spk_times{cc} = cat(1,all_spk_times{cc},ClusterDetails{cc}.t' + cur_toffset);
        cur_ind = find(su_probes == cc);
        if ~isempty(cur_ind)
            aclust_probe = find(aclust_probenums == cc);
            all_clust_ids{cur_ind} = cat(1,all_clust_ids{cur_ind},autoclust(cur_block_set(ee),aclust_probe).idx(:));
        end
    end
    
    trial_start_times = [Expts{cur_block}.Trials(:).TrialStart]/1e4;
    trial_end_times = [Expts{cur_block}.Trials(:).TrueEnd]/1e4;
    trial_durs = [Expts{cur_block}.Trials(:).dur]/1e4;
    trial_ids = [Expts{cur_block}.Trials(:).id];
    
    [un_ids,id_inds] = unique(trial_ids);
    rpt_trials = false;
    if length(un_ids) < length(trial_ids)
        rpt_trials = true;
        fprintf('Warning, repeat trial inds detected!\n');
        use_trials = [];
    else
        use_trials = find(trial_durs >= min_trial_dur);
    end
    %         trial_durs = trial_durs(id_inds);
    %         trial_start_times = trial_start_times(id_inds);
    %         trial_end_times = trial_end_times(id_inds);
    
    
    all_trial_start_times = cat(1,all_trial_start_times,trial_start_times(use_trials)' + cur_toffset);
    all_trial_end_times = cat(1,all_trial_end_times,trial_end_times(use_trials)' + cur_toffset);
    
    trial_Se = [Expts{cur_block}.Trials(:).se];
    trial_Se = trial_Se(id_inds);
    all_trial_Se = cat(1,all_trial_Se,trial_Se(use_trials)');
    all_trial_blk = cat(1,all_trial_blk,ones(length(use_trials),1)*ee);
    
    fname = sprintf('%s/stims/Expt%d_stim',data_dir,cur_block);
    load(fname);
    buffer_pix = floor((expt_npix(cur_block) - use_nPix)/2);
    cur_use_pix = (1:use_nPix) + buffer_pix;
    
    n_trials = length(use_trials);
    for tt = 1:n_trials
        cur_stim_times = Expts{cur_block}.Trials(use_trials(tt)).Start'/1e4;
        n_frames = size(left_stim_mats{use_trials(tt)},1);
        if length(cur_stim_times) == 1
            cur_stim_times = (cur_stim_times:dt*Fr:(cur_stim_times + (n_frames-1)*dt*Fr))';
%             cur_t_edges = [cur_stim_times; cur_stim_times(end) + dt*Fr];
        end
        cur_t_edges = (trial_start_times(use_trials(tt)):dt:trial_end_times(use_trials(tt)))';
        if length(cur_t_edges) > trial_dur/dt + 1
            cur_t_edges(round(trial_dur/dt+2):end) = [];
        end
        
        cur_t_axis = 0.5*cur_t_edges(1:end-1) + 0.5*cur_t_edges(2:end);
        temp(ee,tt) = length(cur_t_axis);
        cur_tsince_start = cur_t_axis - trial_start_times(use_trials(tt));
        
        if ~any(isnan(left_stim_mats{use_trials(tt)}(:))) && n_frames > min_trial_dur/dt
            use_frames = min([length(cur_stim_times),n_frames,trial_dur/stim_dt round(length(cur_t_axis)/usfac)]);
            temp2(ee,tt) = use_frames;
            cur_stim_mat = double(left_stim_mats{use_trials(tt)}(1:use_frames,cur_use_pix));
            
%             cur_Xmat = create_time_embedding(cur_stim_mat,stim_params);
            
            all_stim_mat = [all_stim_mat; cur_stim_mat];   
%             all_Xmat = [all_Xmat; cur_Xmat];
            all_stim_times = [all_stim_times; cur_stim_times' + cur_toffset];
            all_t_axis = [all_t_axis; cur_t_axis + cur_toffset];
            all_t_bin_edges = [all_t_bin_edges; cur_t_edges + cur_toffset];
            all_tsince_start = [all_tsince_start; cur_tsince_start];
            all_blockvec = [all_blockvec; ones(size(cur_t_axis))*ee];
            all_trialvec = [all_trialvec; ones(size(cur_t_axis))*(tt + trial_cnt)];
            all_bin_edge_pts = [all_bin_edge_pts; length(all_t_bin_edges)];
            
        end
    end
    trial_cnt = trial_cnt + n_trials;
    trial_toffset(ee) = all_t_bin_edges(end);
    cur_toffset = trial_toffset(ee);
end

%%
% mahal_thresh = su_data.mah_thresh;
all_binned_mua = nan(length(all_t_axis),n_probes);
%for only-MU probes
for cc = 1:n_probes
    [cur_spkhist,cur_bin_inds] = histc(all_spk_times{cc},all_t_bin_edges);
    cur_spkhist(all_bin_edge_pts) = [];
    all_binned_mua(:,cc) = cur_spkhist;
end

%for SU probes
fprintf('Using %d SUs\n',length(su_ids));
all_binned_sua = nan(length(all_t_axis),length(su_ids));
su_used_blocks = false(length(cur_block_set),length(su_ids));
probe_with_su = false(length(cur_block_set),n_probes);
for ss = 1:length(su_ids)
    all_su_spk_times{ss} = [];
    cur_su_probes = find(any(aclust_CellList(:,su_probes) == su_ids(ss)));
    for pp = 1:length(cur_su_probes)
        all_su_inds = find(all_clust_ids{cur_su_probes(pp)} == 1);
        cur_suahist = histc(all_spk_times{su_probes(cur_su_probes(pp))}(all_su_inds),all_t_bin_edges);
        cur_suahist(all_bin_edge_pts) = [];
        
        spk_block_inds = round(interp1(all_t_axis,all_blockvec,all_spk_times{su_probes(cur_su_probes(pp))}));
        for ee = 1:length(cur_block_set);
            cur_block_inds = find(all_blockvec == ee);
            %if SU is isolated in this block
            if aclust_CellList(ee,su_probes(cur_su_probes(pp))) == su_ids(ss)
                su_used_blocks(ee,ss) = true;
                probe_with_su(ee,su_probes(cur_su_probes(pp))) = true;
                all_binned_sua(cur_block_inds,ss) = cur_suahist(cur_block_inds);
                cur_set = all_su_inds(spk_block_inds(all_su_inds) == ee);
                all_su_spk_times{ss} = cat(1,all_su_spk_times{ss},all_spk_times{su_probes(cur_su_probes(pp))}(cur_set));
            end
        end
    end
end

%%
Fs = 1000;
niqf = Fs/2;
[bb,aa] = butter(2,[1 80]/niqf);

full_lfps = [];
% full_phasegrams = [];
% full_ampgrams = [];
cur_toffset = 0;
for ee = 1:length(cur_block_set);
    fprintf('Expt %d of %d\n',ee,length(cur_block_set));
    fname = sprintf('lemM%dA.%d.lfp.mat',Expt_num,cur_block_set(ee));
    load(fname);
    
    Fs = 1/LFP.Header.CRsamplerate;
    
    n_trials(ee) = length(LFP.Trials);
    %     lfp_trial_starts = [LFP.Trials(:).Start]/1e4;
    lfp_trial_starts = [LFP.Trials(:).ftime]/1e4;
    lfp_trial_ends = [LFP.Trials(:).End]/1e4;
    expt_lfp_t_axis = [];
    expt_lfps = [];
    for tt = 1:n_trials(ee)
        
        cur_npts = size(LFP.Trials(tt).LFP,1);
        cur_t_end(tt) = lfp_trial_starts(tt)+(cur_npts-1)/Fs;
        cur_t_axis = (lfp_trial_starts(tt):1/Fs:cur_t_end(tt)) + cur_toffset;
        
        if ~isempty(expt_lfp_t_axis)
            cur_sp = find(cur_t_axis > max(expt_lfp_t_axis),1,'first');
        else
            cur_sp = 1;
        end
        cur_t_axis = cur_t_axis(cur_sp:end);
        cur_LFP = [LFP.Trials(tt).LFP];
        cur_LFP = cur_LFP(cur_sp:end,:);
        cur_LFP = filtfilt(bb,aa,cur_LFP);
        
        expt_lfp_t_axis = [expt_lfp_t_axis; cur_t_axis(:)];
        expt_lfps = [expt_lfps; cur_LFP];
    end
    
    cur_set = find(all_blockvec==ee);
    interp_lfps = interp1(expt_lfp_t_axis,expt_lfps,all_t_axis(cur_set));
    
    full_lfps = [full_lfps; interp_lfps];

cur_toffset = trial_toffset(ee);
end

%% DEFINE DATA USED FOR ANALYSIS
used_inds = find(all_tsince_start >= beg_buffer);
used_inds(isnan(full_lfps(used_inds,1))) = [];

all_Xmat = create_time_embedding(all_stim_mat,stim_params);

n_blocks = length(cur_block_set);
Xblock = zeros(length(all_t_axis),n_blocks);
for i = 1:n_blocks
    cur_set = find(all_blockvec==i);
    Xblock(cur_set,i) = 1;
end

%% LOCATE REPEAT TRIALS AND SET THEM ASIDE FOR XV
% rpt_seed = 1e4;
% repeat_trials = find(all_trial_Se == rpt_seed);
% n_rpts = length(repeat_trials);
% tr_inds = used_inds(~ismember(all_trialvec(used_inds),repeat_trials));
% xv_inds = used_inds(ismember(all_trialvec(used_inds),repeat_trials));

%% MODELING PARAMS
init_stim_params = NMMcreate_stim_params([flen use_nPix],dt,usfac,tbspace);
init_stim_params(2) = NMMcreate_stim_params([flen use_nPix],dt,usfac,tbspace);

base_lambda_d2XT = 30;
base_lambda_L1 = 0;

mod_signs = [1 1];
NL_types = [{'lin','lin'}];
init_d2XT = [1; 1];
init_reg_params = NMMcreate_reg_params('lambda_d2XT',init_d2XT);
init_Xtargs = [1 2];
silent = 1;

cur_X{1} = all_Xmat(used_inds,:);
cur_X{2} = abs(all_Xmat(used_inds,:));

optim_p.optTol = 1e-5;
optim_p.progTol = 1e-15;
optim_p.maxIter = 50;
for cc = 1:24;
    fprintf('Fitting model for LFP %d of %d\n',cc,24);
    nim = NMMinitialize_model( init_stim_params, mod_signs, NL_types, init_reg_params, init_Xtargs,[], 'linear');
    nim = NMMfit_filters(nim,full_lfps(used_inds,cc),cur_X,[],[],silent,optim_p);
    
    [LL, penLL, pred_rate, G, gint] = NMMmodel_eval(nim,full_lfps(used_inds,cc),cur_X);
    nim = NMMadjust_regularization(nim,1:2,'lambda_d2XT',base_lambda_d2XT./var(gint)');
    nim = NMMadjust_regularization(nim,1:2,'lambda_L1',base_lambda_L1./std(gint)');
    lfp_nim(cc) = NMMfit_filters(nim,full_lfps(used_inds,cc),cur_X,[],[],silent);
end

%%
for cc = 1:24
    filts = [lfp_nim(cc).mods(:).filtK];
    filts = reshape(filts,[flen use_nPix 2]);
    spatial_profiles(cc,:,:) = squeeze(std(filts));
    temporal_profiles(cc,:,:) = squeeze(std(filts,[],2));
end
norm_temporal_profiles = bsxfun(@rdivide,temporal_profiles,max(temporal_profiles,[],2));
norm_spatial_profiles = bsxfun(@rdivide,spatial_profiles,max(spatial_profiles,[],2));
%%
cd(save_dir);
t_ax = (0:(flen-1))*dt;
dw = Expts{cur_block_set(1)}.Stimvals.dw;
pix_ax = ((1:use_nPix) - use_nPix/2)*dw

subplot(2,2,1)
imagesc(t_ax,1:24,squeeze(temporal_profiles(:,:,1)));
xlabel('Time lag (s)','fontsize',14);
ylabel('Probe number','fontsize',14);
title('Phase-dependent','fontsize',14);
subplot(2,2,3)
imagesc(t_ax,1:24,squeeze(temporal_profiles(:,:,2)));
xlabel('Time lag (s)','fontsize',14);
ylabel('Probe number','fontsize',14);
title('Phase-independent','fontsize',14);
subplot(2,2,2)
imagesc(t_ax,1:24,squeeze(norm_temporal_profiles(:,:,1)));
xlabel('Time lag (s)','fontsize',14);
ylabel('Probe number','fontsize',14);
title('Normalized Phase-independent','fontsize',14);
subplot(2,2,4)
imagesc(t_ax,1:24,squeeze(norm_temporal_profiles(:,:,2)));
xlabel('Time lag (s)','fontsize',14);
ylabel('Probe number','fontsize',14);
title('Normalized Phase-independent','fontsize',14);
fillPage(gcf,'papersize',[8 7]);
fname = 'LFP_stim_temporal_profiles';
print(fname,'-dpdf');

figure
subplot(2,2,1)
imagesc(pix_ax,1:24,squeeze(spatial_profiles(:,:,1)));
xlabel('Relative position (deg)','fontsize',14);
ylabel('Probe number','fontsize',14);
title('Phase-dependent','fontsize',14);
subplot(2,2,3)
imagesc(pix_ax,1:24,squeeze(spatial_profiles(:,:,2)));
xlabel('Relative position (deg)','fontsize',14);
ylabel('Probe number','fontsize',14);
title('Phase-independent','fontsize',14);
subplot(2,2,2)
imagesc(pix_ax,1:24,squeeze(norm_spatial_profiles(:,:,1)));
xlabel('Relative position (deg)','fontsize',14);
ylabel('Probe number','fontsize',14);
title('Normalized Phase-independent','fontsize',14);
subplot(2,2,4)
imagesc(pix_ax,1:24,squeeze(norm_spatial_profiles(:,:,2)));
xlabel('Relative position (deg)','fontsize',14);
ylabel('Probe number','fontsize',14);
title('Normalized Phase-independent','fontsize',14);
fillPage(gcf,'papersize',[8 7]);
fname = 'LFP_stim_spatial_profiles';
print(fname,'-dpdf');

