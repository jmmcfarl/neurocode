clear all
% close all
addpath('~/James_scripts/bruce/eye_tracking/');

Expt_num = 86;
Expt_name = sprintf('G0%d',Expt_num);
data_dir = ['~/Data/bruce/' Expt_name];
cd(data_dir);

load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct
load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data

anal_dir = ['~/Analysis/bruce/' Expt_name];
if ~exist(anal_dir,'dir')
    system(['mkdir ' anal_dir]);
end


%%
stim_fs = 100; %in Hz
dt = 0.01;
flen = 10;
use_nPix = 24;
stim_params = NIMcreate_stim_params([flen use_nPix],dt);
Fr = 1;
bar_ori = 90;
dds = 12;
min_trial_dur = 0.75;

beg_buffer = 0.2;
end_buffer = 0.05;
trial_dur = 4;


%%
if strcmp(Expt_name,'G093')
    include_expts = {'rls.FaXwi','rls.FaXwiXimi'};
else
    include_expts = {'rls.Fa', 'rls.FaXimi','rls.FaRC'};
end
for ii = 1:length(Expts)
    expt_names{ii} = Expts{ii}.Header.expname;
    expt_dds(ii) = Expts{ii}.Stimvals.dd;
    expt_bar_ori(ii) = Expts{ii}.Stimvals.or;
    expt_sac_dir(ii) = mod(Expts{ii}.Stimvals.Fa,180);
    expt_Fr(ii) = Expts{ii}.Stimvals.Fr;
    expt_imback(ii) = isfield(Expts{ii}.Trials,'imi');
    included_type(ii) = any(strcmp(expt_names{ii},include_expts));
end
expt_has_ds(isnan(expt_has_ds)) = 0;
expt_has_ds(expt_has_ds == -1) = 0;
expt_binoc(isnan(expt_binoc)) = 0;
cur_block_set = find(included_type & expt_binoc' & expt_Fr == 1 & expt_bar_ori == bar_ori & expt_sac_dir == bar_ori & (expt_dds == 12 | expt_dds == 67));

if exist('./fin_aclust_data.mat','file')
    load('./fin_aclust_data.mat');
    [n_aclust_expts,n_aclust_probes] = size(autoclust);
else
    disp('No fin_aclust_data found.');
    autoclust = [];
    n_aclust_expts = 0; n_aclust_probes = 0;
end


if strcmp(Expt_name,'G087')
    cur_block_set(cur_block_set == 15) = []; %only 6 trials and causes problems
end

if strcmp(Expt_name,'G087')
    cur_block_set(cur_block_set == 15) = []; %only 6 trials and causes problems
end
if strcmp(Expt_name,'G093')
    cur_block_set(cur_block_set ==  28) = []; %only 6 trials and causes problems
end

sim_sac_expts = find(~expt_has_ds(cur_block_set));
imback_gs_expts = find(expt_has_ds(cur_block_set) & expt_imback(cur_block_set)');
grayback_gs_expts = find(expt_has_ds(cur_block_set) & ~expt_imback(cur_block_set)');

cluster_dir = ['~/Analysis/bruce/' Expt_name '/clustering'];
n_probes = 96;
%% COMPUTE TRIAL DATA
cd(data_dir);

fprintf('Computing prep data\n');
trial_cnt = 0;

all_stim_times = [];
% all_left_Xmat = [];
% all_right_Xmat = [];
all_left_stim_mat = [];
all_right_stim_mat = [];
all_t_axis = [];
all_t_bin_edges = [];
all_tsince_start = [];
all_blockvec = [];
all_trialvec = [];
all_trial_wi = [];
all_trial_start_times = [];
all_trial_end_times = [];
all_bin_edge_pts = [];
all_spk_times = cell(96,1);
all_spk_inds = cell(96,1);
all_clust_ids = cell(96,1);
cur_toffset = 0;
cur_spkind_offset = 0;
for ee = 1:length(cur_block_set);
    if ismember(ee,grayback_gs_expts)
        fprintf('Expt %d Block %d of %d; grayback GS, ori:%d\n',Expt_num,ee,length(cur_block_set),expt_bar_ori(cur_block_set(ee)));
    elseif ismember(ee,imback_gs_expts)
        fprintf('Expt %d Block %d of %d; imback GS, ori:%d\n',Expt_num,ee,length(cur_block_set),expt_bar_ori(cur_block_set(ee)));
    elseif ismember(ee,sim_sac_expts)
        fprintf('Expt %d Block %d of %d; SimSac, ori:%d\n',Expt_num,ee,length(cur_block_set),expt_bar_ori(cur_block_set(ee)));
    else
        fprintf('Expt %d Block %d of %d;  UNMATCHED EXPT TYPE\n',Expt_num,ee,length(cur_block_set));
    end
    cur_block = cur_block_set(ee);
    
    fname = [cluster_dir sprintf('/Block%d_Clusters.mat',cur_block)];
    load(fname,'Clusters');
    for cc = 1:n_probes
        all_spk_times{cc} = cat(1,all_spk_times{cc},Clusters{cc}.times + cur_toffset);
        all_spk_inds{cc} = cat(1,all_spk_inds{cc},Clusters{cc}.spk_inds + cur_spkind_offset);
        all_clust_ids{cc} = cat(1,all_clust_ids{cc},Clusters{cc}.spike_clusts);
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
    
    all_trial_start_times = cat(1,all_trial_start_times,trial_start_times(use_trials)');
    all_trial_end_times = cat(1,all_trial_end_times,trial_end_times(use_trials)');
    
    if strcmp(Expt_name,'G093')
        trial_wi = [Expts{cur_block}.Trials(:).wi];
        trial_wi = trial_wi(id_inds);
        all_trial_wi = cat(1,all_trial_wi,trial_wi(use_trials)');
    end
    
    fname = sprintf('%s/stims/Expt%d_stim',data_dir,cur_block);
    load(fname);
    buffer_pix = floor((expt_npix(cur_block) - use_nPix)/2);
    cur_use_pix = (1:use_nPix) + buffer_pix;
    
    n_trials = length(use_trials);
    for tt = 1:n_trials
        cur_stim_times = Expts{cur_block}.Trials(use_trials(tt)).Start/1e4;
        n_frames = size(left_stim_mats{use_trials(tt)},1);
        n_frames_r = size(right_stim_mats{use_trials(tt)},1);
        if n_frames_r ~= n_frames
            error('Left-right mismatch!');
        end
        if n_frames > 0
            if length(cur_stim_times) == 1
                cur_stim_times = (cur_stim_times:dt*Fr:(cur_stim_times + (n_frames-1)*dt*Fr))';
                cur_stim_times(cur_stim_times > trial_end_times(use_trials(tt))) = [];
            end
        end
        cur_t_edges = [cur_stim_times; cur_stim_times(end) + dt*Fr];
        cur_t_axis = 0.5*cur_t_edges(1:end-1) + 0.5*cur_t_edges(2:end);
        
        cur_tsince_start = cur_t_axis - trial_start_times(use_trials(tt));
        
        if ~any(isnan(left_stim_mats{use_trials(tt)}(:))) && n_frames > min_trial_dur/dt
            use_frames = min(length(cur_stim_times),n_frames);
            %             cur_stim_mat = double(left_stim_mats{use_trials(tt)}(1:use_frames,cur_use_pix));
            cur_lstim_mat = double(left_stim_mats{use_trials(tt)}(1:use_frames,cur_use_pix));
            cur_rstim_mat = double(right_stim_mats{use_trials(tt)}(1:use_frames,cur_use_pix));
            %             left_Xmat = create_time_embedding(cur_lstim_mat,stim_params);
            %             right_Xmat = create_time_embedding(cur_rstim_mat,stim_params);
            
            if ~isempty(all_stim_times)
                if any(cur_stim_times < all_stim_times(end))
                    fprintf('Warn trial %d\n',tt);
                end
            end
            all_stim_times = [all_stim_times; cur_stim_times];
            all_t_axis = [all_t_axis; cur_t_axis];
            all_left_stim_mat = [all_left_stim_mat; cur_lstim_mat];
            all_right_stim_mat = [all_right_stim_mat; cur_rstim_mat];
            all_t_bin_edges = [all_t_bin_edges; cur_t_edges];
            all_tsince_start = [all_tsince_start; cur_tsince_start];
            all_blockvec = [all_blockvec; ones(size(cur_t_axis))*ee];
            all_trialvec = [all_trialvec; ones(size(cur_t_axis))*(tt + trial_cnt)];
            all_bin_edge_pts = [all_bin_edge_pts; length(all_t_bin_edges)];
        end
    end
    trial_cnt = trial_cnt + n_trials;
end

%% BIN SPIKES FOR MU AND SU
rec_type = 'UA';
clust_params.n_probes = n_probes;
if strcmp(rec_type,'LP')
    clust_params.exclude_adjacent = true;
else
    clust_params.exclude_adjacent = false;
end
[all_binned_mua,all_binned_sua,Clust_data] = ...
    get_binned_spikes(cluster_dir,all_spk_times,all_clust_ids,all_spk_inds,...
    all_t_axis,all_t_bin_edges,all_bin_edge_pts,cur_block_set,all_blockvec,clust_params);
SU_probes = Clust_data.SU_probes;
SU_numbers = Clust_data.SU_numbers;

%% DEFINE DATA USED FOR ANALYSIS
used_inds = find(all_tsince_start >= beg_buffer & (trial_dur-all_tsince_start) >= end_buffer);
if strcmp(Expt_name,'G093')
    un_wi_vals = unique(all_trial_wi);
    use_wi_trials = find(all_trial_wi == un_wi_vals(2));
    used_inds = used_inds(ismember(all_trialvec(used_inds),use_wi_trials));
end
NT = length(used_inds);

%% create and normalize stim mats
dd12_set = find(expt_dds(cur_block_set) == 12);
dd67_set = find(expt_dds(cur_block_set) == 67);
dd67_set(end) = [];

used_inds_12 = used_inds(ismember(all_blockvec(used_inds),dd12_set));
used_inds_67 = used_inds(ismember(all_blockvec(used_inds),dd67_set));
used_inds_67((length(used_inds_12)+1):end) = [];
%%
all_left_Xmat = create_time_embedding(all_left_stim_mat,stim_params);
all_right_Xmat = create_time_embedding(all_right_stim_mat,stim_params);



%% SUA STA/STC ANALYSIS
nneg = 5; npos = 5;
stc_thresh = -5e-3;

max_squared_filts = 2;
full_Xmat = [all_left_Xmat all_right_Xmat];
nmm_stim_params = NMMcreate_stim_params([flen 2*use_nPix],dt);
nmm_stim_params(2) = NMMcreate_stim_params([length(cur_block_set),1],dt);

%create L2 mat for left eye
L2_params = create_L2_params([],[1 flen*use_nPix],[flen use_nPix],2,3,[Inf 0]);
L2_mat = generate_L2_mat(L2_params,2*flen*use_nPix);
%add L2 mat for sine term
L2_params = create_L2_params([],[flen*use_nPix+1 2*flen*use_nPix],[flen use_nPix],2,3,[Inf 0]);
L2_mat = L2_mat + generate_L2_mat(L2_params,2*flen*use_nPix);

lambda_d2XT = 4000;
lambda_L1 = 100;
sq_sc_fac = 4;
block_L2 = 1;
sub_samp_fac = 5;
silent = 1;
optim_p.optTol = 1e-5; optim_p.progTol = 1e-9;

for ss = 1:length(SU_numbers)
    
    fprintf('Computing STC for SU %d of %d\n',ss,length(SU_numbers));
    
    Robs = all_binned_sua(used_inds_12,ss);
    cur_used_inds = used_inds_12(~isnan(Robs));
    if ~isempty(cur_used_inds)
        Robs = all_binned_sua(cur_used_inds,ss);
        avg_rate = mean(Robs);
        
        spikebins = convert_to_spikebins(Robs);
        spike_cond_stim = full_Xmat(cur_used_inds(spikebins),:);
        sta      = mean(spike_cond_stim) - mean(full_Xmat(cur_used_inds,:));
        sta = sta/norm(sta);
        proj_mat = sta'/(sta*sta')*sta;
        stim_proj = full_Xmat(cur_used_inds,:) - full_Xmat(cur_used_inds,:)*proj_mat;
        % stim_proj = stim_emb;
        stvcv = cov(stim_proj(spikebins,:));  utvcv = cov(stim_proj);
        [evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);
        stcs  = evecs(:,[nneg:-1:1,length(evs)-npos+1:end]); stcs  = stcs(:,end:-1:1);
        
        sua_data(ss).avg_rates = avg_rate;
        sua_data(ss).var_rates = nanvar(Robs);
        sua_data(ss).sta = sta;
        sua_data(ss).stcs = stcs;
        sua_data(ss).evals = diag(evals);
        cur_evec_diff = diff(flipud(sua_data(ss).evals));
        sua_data(ss).npos_stc = find(cur_evec_diff(1:250) < stc_thresh,1,'last');
        sua_data(ss).nneg_stc = find(cur_evec_diff(end:-1:250) < stc_thresh,1,'last');
        if isempty(sua_data(ss).npos_stc); sua_data(ss).npos_stc = 0; end;
        if isempty(sua_data(ss).nneg_stc); sua_data(ss).nneg_stc = 0; end;
        sua_data(ss).stc_use = true;
        sua_data(ss).avg_rate = mean(Robs);
        sua_data(ss).nspks = sum(Robs);
        
        %%
    else
        sua_data(ss).stc_use = false;
    end
end

%%
for ss = 1:length(SU_numbers)
    
    fprintf('Computing STC for SU %d of %d\n',ss,length(SU_numbers));
    
    Robs = all_binned_sua(used_inds_67,ss);
    cur_used_inds = used_inds_67(~isnan(Robs));
    if ~isempty(cur_used_inds)
        Robs = all_binned_sua(cur_used_inds,ss);
        avg_rate = mean(Robs);
        
        spikebins = convert_to_spikebins(Robs);
        spike_cond_stim = full_Xmat(cur_used_inds(spikebins),:);
        sta      = mean(spike_cond_stim) - mean(full_Xmat(cur_used_inds,:));
        sta = sta/norm(sta);
        proj_mat = sta'/(sta*sta')*sta;
        stim_proj = full_Xmat(cur_used_inds,:) - full_Xmat(cur_used_inds,:)*proj_mat;
        % stim_proj = stim_emb;
        stvcv = cov(stim_proj(spikebins,:));  utvcv = cov(stim_proj);
        [evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);
        stcs  = evecs(:,[nneg:-1:1,length(evs)-npos+1:end]); stcs  = stcs(:,end:-1:1);
        
        sua_data67(ss).avg_rates = avg_rate;
         sua_data67(ss).var_rates = nanvar(Robs);
       sua_data67(ss).sta = sta;
        sua_data67(ss).stcs = stcs;
        sua_data67(ss).evals = diag(evals);
        cur_evec_diff = diff(flipud(sua_data67(ss).evals));
        sua_data67(ss).npos_stc = find(cur_evec_diff(1:250) < stc_thresh,1,'last');
        sua_data67(ss).nneg_stc = find(cur_evec_diff(end:-1:250) < stc_thresh,1,'last');
        if isempty(sua_data67(ss).npos_stc); sua_data67(ss).npos_stc = 0; end;
        if isempty(sua_data67(ss).nneg_stc); sua_data67(ss).nneg_stc = 0; end;
        sua_data67(ss).stc_use = true;
        sua_data67(ss).avg_rate = mean(Robs);
        sua_data67(ss).nspks = sum(Robs);
        
        %%
    else
        sua_data67(ss).stc_use = false;
    end
end

%% MODEL FITS FOR MUA
stc_thresh = -1e-3;
for ss = 1:96
    
    fprintf('Computing STC for MU %d of %d\n',ss,96);
    su_probe_ind = find(su_probes == ss);
    if ~isempty(su_probe_ind)
        cur_used_blocks = find(su_used_blocks(:,su_probe_ind) == 0); %blocks when NO SU
        cur_poss_inds = find(ismember(all_blockvec,cur_used_blocks));
    else
        cur_used_blocks = 1:length(cur_block_set); %blocks when NO SU
        cur_poss_inds = find(ismember(all_blockvec,cur_used_blocks));
    end
    cur_used_inds = used_inds(ismember(used_inds,cur_poss_inds));
    cur_NT = length(cur_used_inds);
    
    if ~isempty(cur_used_inds)
        Robs = all_binned_spikes(cur_used_inds,ss);
        avg_rate = mean(Robs);
        
        spikebins = convert_to_spikebins(Robs);
        spike_cond_stim = full_Xmat(cur_used_inds(spikebins),:);
        sta      = mean(spike_cond_stim) - mean(full_Xmat(cur_used_inds,:));
        sta = sta/norm(sta);
        proj_mat = sta'/(sta*sta')*sta;
        stim_proj = full_Xmat(cur_used_inds,:) - full_Xmat(cur_used_inds,:)*proj_mat;
        % stim_proj = stim_emb;
        stvcv = cov(stim_proj(spikebins,:));  utvcv = cov(stim_proj);
        [evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);
        stcs  = evecs(:,[nneg:-1:1,length(evs)-npos+1:end]); stcs  = stcs(:,end:-1:1);
        
        mua_data(ss).block_avg_rates = block_mean_rates(cur_used_blocks,ss);
        mua_data(ss).sta = sta;
        mua_data(ss).stcs = stcs;
        mua_data(ss).evals = diag(evals);
        cur_evec_diff = diff(flipud(mua_data(ss).evals));
        mua_data(ss).npos_stc = find(cur_evec_diff(1:250) < stc_thresh,1,'last');
        mua_data(ss).nneg_stc = find(cur_evec_diff(end:-1:250) < stc_thresh,1,'last');
        if isempty(mua_data(ss).npos_stc); mua_data(ss).npos_stc = 0; end;
        if isempty(mua_data(ss).nneg_stc); mua_data(ss).nneg_stc = 0; end;
        mua_data(ss).stc_use = true;
        mua_data(ss).avg_rate = mean(Robs);
        mua_data(ss).nspks = sum(Robs);
        cur_n_grayback_blocks = sum(ismember(cur_used_blocks,grayback_gs_expts));
        cur_n_imback_blocks = sum(ismember(cur_used_blocks,imback_gs_expts));
        mua_data(ss).nbocks = [cur_n_grayback_blocks cur_n_imback_blocks];
        
        %%
        n_squared_filts = min(mua_data(ss).npos_stc,max_squared_filts);
        fprintf('Fitting NIM with %d exc squared\n',n_squared_filts);
        mod_signs = [1 ones(1,n_squared_filts) 1];
        NL_types = [{'lin'} repmat({'quad'},1,n_squared_filts) {'lin'}];
        Xtargets = [ones(1,n_squared_filts+1) 2];
        d2XT = [lambda_d2XT repmat(sq_sc_fac*sqrt(lambda_d2XT),1,n_squared_filts) 0];
        L2 = [zeros(1,n_squared_filts+1) block_L2];
        reg_params = NMMcreate_reg_params('lambda_custom',d2XT,'lambda_L2',L2);
        init_filts = cell(length(mod_signs),1); init_filts{end} = zeros(size(Xblock,2),1);
        fit0 = NMMinitialize_model(nmm_stim_params,mod_signs,NL_types,reg_params,Xtargets,init_filts); %initialize NIM
        sub_sample = randperm(cur_NT);
        sub_sample = sub_sample(1:round(cur_NT/sub_samp_fac));
        cur_Xmat{1} = full_Xmat(cur_used_inds(sub_sample),:); cur_Xmat{2} = Xblock(cur_used_inds(sub_sample),:);
        fit0 = NMMfit_filters(fit0,Robs(sub_sample),cur_Xmat,[],[],silent,[],L2_mat); %fit stimulus filters
        
        fit0 = NMMadjust_regularization(fit0,1:(1+n_squared_filts),'lambda_L1',[lambda_L1 sqrt(lambda_L1) sqrt(lambda_L1)]);
        cur_Xmat{1} = full_Xmat(cur_used_inds,:); cur_Xmat{2} = Xblock(cur_used_inds,:);
        fit0 = NMMfit_filters(fit0,Robs,cur_Xmat,[],[],silent,[],L2_mat); %fit stimulus filters
        [LL, penLL, pred_rate, G, gint, fgint, nullLL] = NMMmodel_eval( fit0, Robs, cur_Xmat);
        mua_data(ss).nimFit = fit0;
        mua_data(ss).nullLL = fit0;
        
    else
        mua_data(ss).stc_use = false;
    end
end


%%
cd(anal_dir)
modfit_used_blocks = cur_used_blocks;
fprintf('Saving data to %s\n',anal_name);
if exist('mua_data','var')
    save(anal_name,'sua_data','mua_data','modfit_used_blocks');
else
    fprintf('Saving SU data only\n');
    save(anal_name,'sua_data','modfit_used_blocks');
end
%% VIEW MODEL FITS COMPARED TO STA/STC
f1 = figure();
f2 = figure();
for ss = 1:length(su_probes)
    if sua_data(ss).stc_use
        fprintf('SU %d of %d\n',ss,length(su_probes));
        fprintf('%d blocks, %d spikes\n',sua_data(ss).nbocks(1),sua_data(ss).nspks);
        figure(f1); clf
        cur_sta = sua_data(ss).sta; ca = max(abs(cur_sta));
        subplot(3,3,1)
        imagesc(reshape(cur_sta,flen,2*use_nPix)); caxis([-ca ca]);
        cur_stcs = sua_data(ss).stcs; ca = max(abs(cur_stcs(:)));
        for ii = 1:3
            subplot(3,3,3+ii)
            imagesc(reshape(cur_stcs(:,ii),flen,2*use_nPix)); caxis([-ca ca]);
        end
        for ii = 1:3
            subplot(3,3,6+ii)
            imagesc(reshape(cur_stcs(:,end-3+ii),flen,2*use_nPix)); caxis([-ca ca]);
        end
        
        kmat = [sua_data(ss).nimFit.mods(1:end-1).filtK];
        figure(f2); clf
        subplot(2,2,1)
        imagesc(reshape(kmat(:,1),flen,2*use_nPix));
        ca = max(abs(kmat(:,1))); caxis([-ca ca]);
        for ii = 1:(size(kmat,2)-1)
            subplot(2,2,2+ii)
            imagesc(reshape(kmat(:,ii+1),flen,2*use_nPix));
            ca = max(abs(kmat(:,ii+1))); caxis([-ca ca]);
        end
    end
    pause
    
end
%% FOR MUA
f1 = figure();
f2 = figure();
for ss = 1:96
    if mua_data(ss).stc_use
        fprintf('MU %d of %d\n',ss,96);
        
        figure(f1); clf
        cur_sta = mua_data(ss).sta; ca = max(abs(cur_sta));
        subplot(3,3,1)
        imagesc(reshape(cur_sta,flen,2*use_nPix)); caxis([-ca ca]);
        cur_stcs = mua_data(ss).stcs; ca = max(abs(cur_stcs(:)));
        for ii = 1:3
            subplot(3,3,3+ii)
            imagesc(reshape(cur_stcs(:,ii),flen,2*use_nPix)); caxis([-ca ca]);
        end
        for ii = 1:3
            subplot(3,3,6+ii)
            imagesc(reshape(cur_stcs(:,end-3+ii),flen,2*use_nPix)); caxis([-ca ca]);
        end
        
        kmat = [mua_data(ss).nimFit.mods(1:end-1).filtK];
        figure(f2); clf
        subplot(2,2,1)
        imagesc(reshape(kmat(:,1),flen,2*use_nPix));
        ca = max(abs(kmat(:,1))); caxis([-ca ca]);
        for ii = 1:(size(kmat,2)-1)
            subplot(2,2,2+ii)
            imagesc(reshape(kmat(:,ii+1),flen,2*use_nPix));
            ca = max(abs(kmat(:,ii+1))); caxis([-ca ca]);
        end
        pause
    end
end

