clear all
% close all

dir_prefix = '~';
% dir_prefix = '/Volumes/james';

Expt_num = 86;

%% SET ANALYSIS PARAMETERS
stim_fs = 100; %in Hz
dt = 0.01;

min_trial_dur = 0.5;
beg_buffer = 0.3;

backlag = round(0.25/dt);
forwardlag = round(0.55/dt);

sua_sm_sig = (0.01/dt);
mua_sm_sig = (0.005/dt);

flen = 10;
use_nPix = 24;
full_nPix = 36;
stim_params = NIMcreate_stim_params([flen full_nPix],dt);
Fr = 1;
trial_dur = 4;
bar_ori = 90;

eye_data_name = 'binoc_eyecorr_vbar_eyeonly';

%% LOAD OVERALL SU DATA
load ~/Analysis/bruce/summary_analysis/su_data.mat


clear expt* included_type

Expt_name = sprintf('G0%d',Expt_num);
data_dir = [dir_prefix '/Data/bruce/' Expt_name];
anal_dir = ['~/Analysis/bruce/' Expt_name];
addpath([dir_prefix '/James_scripts/bruce/G081/'])

cd(data_dir);

load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct
if ~strcmp(Expt_name,'G081')
    load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data
end

save_dir = [dir_prefix '/Analysis/bruce/' Expt_name '/sac_mod'];
if ~exist(save_dir,'dir')
    system(['mkdir ' save_dir]);
end

if exist('./fin_aclust_data.mat','file')
    load('./fin_aclust_data.mat');
    [n_aclust_expts,n_aclust_probes] = size(autoclust);
else
    disp('No fin_aclust_data found.');
    autoclust = [];
    n_aclust_expts = 0; n_aclust_probes = 0;
end

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
% cur_block_set = find(included_type & ~expt_binoc' & expt_Fr == 1);
cur_block_set = find(included_type & expt_binoc' & expt_Fr == 1 & expt_bar_ori == bar_ori & expt_sac_dir == bar_ori);

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

%% load overall su data
cur_expt_id = find(su_data.expt_nums == Expt_num);
su_probes = find(su_data.is_su(cur_expt_id,:));
mua_probes = setdiff(1:96,su_probes); %probes with ONLY MU
aclust_probenums = [autoclust(cur_block_set(1),:).probenum];
autoclust = autoclust(:,ismember(aclust_probenums,su_probes));

%% COMPUTE TRIAL DATA
cd(data_dir);

fprintf('Computing prep data\n');
trial_cnt = 0;

all_stim_times = [];
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
all_clust_ids = cell(length(su_probes),1);
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
    fname = sprintf('Expt%dClusterTimesDetails.mat',cur_block);
    load(fname);
    for cc = 1:96
        all_spk_times{cc} = cat(1,all_spk_times{cc},ClusterDetails{cc}.t');
        cur_ind = find(su_probes == cc);
        if ~isempty(cur_ind)
            all_clust_ids{cur_ind} = cat(1,all_clust_ids{cur_ind},autoclust(cur_block_set(ee),cur_ind).idx(:));
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
    
    
    all_trial_start_times = cat(1,all_trial_start_times,trial_start_times(use_trials)');
    all_trial_end_times = cat(1,all_trial_end_times,trial_end_times(use_trials)');
    
    if strcmp(Expt_name,'G093')
        trial_wi = [Expts{cur_block}.Trials(:).wi];
        trial_wi = trial_wi(id_inds);
        all_trial_wi = cat(1,all_trial_wi,trial_wi(use_trials)');
    end
    
    
    fname = sprintf('%s/stims/Expt%d_stim',data_dir,cur_block);
    load(fname);
    buffer_pix = floor((expt_npix(cur_block) - full_nPix)/2);
    cur_use_pix = (1:full_nPix) + buffer_pix;
    
    n_trials = length(use_trials);
    for tt = 1:n_trials
        cur_stim_times = Expts{cur_block}.Trials(use_trials(tt)).Start'/1e4;
        n_frames = size(left_stim_mats{use_trials(tt)},1);
        if n_frames > 0
            if length(cur_stim_times) == 1
                cur_stim_times = (cur_stim_times:dt*Fr:(cur_stim_times + (n_frames-1)*dt*Fr))';
                cur_stim_times(cur_stim_times > trial_end_times(use_trials(tt))) = [];
                cur_t_edges = [cur_stim_times; cur_stim_times(end) + dt*Fr];
            end
        end
        %             cur_t_edges = (trial_start_times(use_trials(tt)):dt:trial_end_times(use_trials(tt)))';
        cur_t_axis = 0.5*cur_t_edges(1:end-1) + 0.5*cur_t_edges(2:end);
        
        cur_tsince_start = cur_t_axis - trial_start_times(use_trials(tt));
        
        if ~any(isnan(left_stim_mats{use_trials(tt)}(:))) && n_frames > min_trial_dur/dt
            use_frames = min(length(cur_stim_times),n_frames);
            cur_lstim_mat = double(left_stim_mats{use_trials(tt)}(1:use_frames,cur_use_pix));
            cur_rstim_mat = double(right_stim_mats{use_trials(tt)}(1:use_frames,cur_use_pix));
            
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
mahal_thresh = su_data.mah_thresh;
all_binned_spikes = nan(length(all_t_axis),96);
su_used_blocks = false(length(cur_block_set),length(su_probes));
%for only-MU probes
for cc = 1:96
    if ~ismember(cc,su_probes) %if probe doesn't have an SU
        [cur_spkhist,cur_bin_inds] = histc(all_spk_times{cc},all_t_bin_edges);
        cur_spkhist(all_bin_edge_pts) = [];
        all_binned_spikes(:,cc) = cur_spkhist;
    end
end
%for SU probes
for ss = 1:length(su_probes)
    %for MUA
    cur_muahist = histc(all_spk_times{su_probes(ss)},all_t_bin_edges);
    cur_muahist(all_bin_edge_pts) = [];
    
    all_su_inds = find(all_clust_ids{ss} == 1);
    cur_suahist = histc(all_spk_times{su_probes(ss)}(all_su_inds),all_t_bin_edges);
    cur_suahist(all_bin_edge_pts) = [];
    
    spk_block_inds = round(interp1(all_t_axis,all_blockvec,all_spk_times{su_probes(ss)}));
    for ee = 1:length(cur_block_set);
        cur_block_inds = find(all_blockvec == ee);
        %if SU is isolated in this block
        if autoclust(cur_block_set(ee),ss).mahal_d > mahal_thresh || autoclust(cur_block_set(ee),ss).man_code == 4
            su_used_blocks(ee,ss) = true;
            all_binned_spikes(cur_block_inds,su_probes(ss)) = cur_suahist(cur_block_inds);
        else %otherwise use MUA
            all_binned_spikes(cur_block_inds,su_probes(ss)) = cur_muahist(cur_block_inds);
        end
    end
end
%% DEFINE DATA USED FOR ANALYSIS
used_inds = find(all_tsince_start >= beg_buffer);
if strcmp(Expt_name,'G093')
    un_wi_vals = unique(all_trial_wi);
    use_wi_trials = find(all_trial_wi == un_wi_vals(2) | all_trial_wi == un_wi_vals(3));
    used_inds = used_inds(ismember(all_trialvec(used_inds),use_wi_trials));
end

%% PARSE EYE-TRACKING DATA (IGNORING RIGHT EYE BECAUSE OF NOISE)

emfile = ['jbe' Expt_name '.em.mat'];
load(emfile);

all_eye_vals = [];
all_eye_speed = [];
all_eye_ts = [];
all_eye_blockvec = [];
eye_smooth = 3;
for ee = 1:length(cur_block_set);
    fprintf('Loading ET data for expt %d, block %d of %d\n',Expt_num,ee,length(cur_block_set));
    cur_set = find(all_blockvec==ee);
    if ~isempty(cur_set)
        [eye_vals_interp,eye_ts_interp,eye_insample] = get_eye_tracking_data_new(Expt,all_t_axis(cur_set([1 end])));
        
        eye_dt = Expt.Header.CRrates(1);
        eye_fs = 1/eye_dt;
        lEyeXY = eye_vals_interp(:,1:2);
        rEyeXY = eye_vals_interp(:,3:4);
        
        %slight smoothing before computing speed
        sm_avg_eyepos = lEyeXY; eye_vel = lEyeXY; %initialization
        sm_avg_eyepos(:,1) = smooth(lEyeXY(:,1),eye_smooth);
        sm_avg_eyepos(:,2) = smooth(lEyeXY(:,2),eye_smooth);
        eye_vel(:,1) = [0; diff(sm_avg_eyepos(:,1))]/eye_dt;
        eye_vel(:,2) = [0; diff(sm_avg_eyepos(:,2))]/eye_dt;
        
        eye_speed = sqrt(eye_vel(:,1).^2+eye_vel(:,2).^2);
        
        all_eye_vals = [all_eye_vals; lEyeXY rEyeXY];
        all_eye_speed = [all_eye_speed; eye_speed];
        all_eye_ts = [all_eye_ts; eye_ts_interp'];
        all_eye_blockvec = [all_eye_blockvec; ee*ones(size(eye_speed))];
    end
end

back_pts = 1 + find(diff(all_eye_ts) <= 0);
double_samples = [];
for i = 1:length(back_pts)
    next_forward = find(all_eye_ts > all_eye_ts(back_pts(i)-1),1,'first');
    double_samples = [double_samples back_pts(i):next_forward];
end
all_eye_ts(double_samples) = [];
all_eye_speed(double_samples) = [];
all_eye_vals(double_samples,:) = [];

interp_eye_vals = interp1(all_eye_ts,all_eye_vals,all_t_axis);
interp_eye_speed = interp1(all_eye_ts,all_eye_speed,all_t_axis);

%% DETECT SACCADES
sac_thresh = 10;
peak_sig = [0; diff(sign(diff(all_eye_speed))); 0];
saccade_inds = find(peak_sig == -2 & all_eye_speed > sac_thresh);

peri_thresh = 3; %threshold eye speed for defining saccade boundary inds
thresh_cross_up = 1 + find(all_eye_speed(1:end-1) < peri_thresh & all_eye_speed(2:end) >= peri_thresh);
thresh_cross_down = 1 + find(all_eye_speed(1:end-1) >= peri_thresh & all_eye_speed(2:end) < peri_thresh);
sac_start_inds = nan(size(saccade_inds));
sac_stop_inds = nan(size(saccade_inds));
for ii = 1:length(saccade_inds)
    next_tc = find(thresh_cross_down > saccade_inds(ii),1,'first');
    if ~isempty(next_tc)
        sac_stop_inds(ii) = thresh_cross_down(next_tc);
    end
    prev_tc = find(thresh_cross_up < saccade_inds(ii),1,'last');
    if ~isempty(prev_tc)
        sac_start_inds(ii) = thresh_cross_up(prev_tc);
    end
    
end

%get rid of double-peaks
min_isi = 0.05; max_isi = Inf;
isis = [Inf; diff(sac_start_inds)]/eye_fs;
bad_isis = (isis < min_isi | isis > max_isi);
bad_sacs = find(isnan(sac_stop_inds) | isnan(sac_start_inds) | bad_isis);
saccade_inds(bad_sacs) = []; isis(bad_sacs) = []; sac_start_inds(bad_sacs) = []; sac_stop_inds(bad_sacs) = [];

saccade_times = all_eye_ts(saccade_inds);
sac_start_times = all_eye_ts(sac_start_inds);
sac_stop_times = all_eye_ts(sac_stop_inds);
sac_durs = sac_stop_times - sac_start_times;

sac_dbuff = round(0.05/eye_dt);
sac_pre_pos = nan(length(saccade_inds),4);
sac_post_pos = nan(length(saccade_inds),4);
for ii = 1:length(saccade_inds)
    pre_inds = (sac_start_inds(ii) - sac_dbuff):sac_start_inds(ii);
    pre_inds(pre_inds < 1) = 1;
    sac_pre_pos(ii,:) = median(all_eye_vals(pre_inds,:),1);
    post_inds = sac_stop_inds(ii):(sac_stop_inds(ii) + sac_dbuff);
    post_inds(post_inds > length(all_eye_ts)) = length(all_eye_ts);
    sac_post_pos(ii,:) = median(all_eye_vals(post_inds,:),1);
end

%use only left-eye signal here
sac_delta_pos = sac_post_pos(:,1:2) - sac_pre_pos(:,1:2);
sac_amps = sqrt(sum(sac_delta_pos.^2,2));
sac_dirs = atan2(sac_delta_pos(:,2),sac_delta_pos(:,1));

temp = ones(length(saccade_times),1);
saccades = struct('peak_time',mat2cell(saccade_times,temp),'start_time',mat2cell(sac_start_times,temp),...
    'stop_time',mat2cell(sac_stop_times,temp),'isi',mat2cell(isis,temp),...
    'duration',mat2cell(sac_durs,temp),'amplitude',mat2cell(sac_amps,temp),'direction',mat2cell(sac_dirs,temp),...
    'pre_Lx',mat2cell(sac_pre_pos(:,1),temp),'post_Lx',mat2cell(sac_post_pos(:,1),temp),...
    'pre_Ly',mat2cell(sac_pre_pos(:,2),temp),'post_Ly',mat2cell(sac_post_pos(:,2),temp),...
    'pre_Rx',mat2cell(sac_pre_pos(:,3),temp),'post_Rx',mat2cell(sac_post_pos(:,3),temp),...
    'pre_Ry',mat2cell(sac_pre_pos(:,4),temp),'post_Ry',mat2cell(sac_post_pos(:,4),temp));

%% BLOCK PREDICTOR
Xblock = zeros(length(all_stim_times),length(cur_block_set));
for i = 1:length(cur_block_set)
    cur_set = find(all_blockvec==i);
    Xblock(cur_set,i) = 1;
end
if Expt_num == 93
    Xwi = zeros(length(all_stim_times),2);
    cur_set = find(all_trial_wi == un_wi_vals(2));
    Xwi(cur_set,1) = 1;
    cur_set = find(all_trial_wi == un_wi_vals(3));
    Xwi(cur_set,2) = 1;
    
    Xblock = [Xblock Xwi];
end

[Xinds,Tinds] = meshgrid(1:full_nPix,1:flen);
buffer_pix = floor((full_nPix - use_nPix)/2);
cur_use_pix = (1:use_nPix) + buffer_pix;
use_kInds = find(ismember(Xinds(:),cur_use_pix));

%% LOAD INFERRED EYE POSITIONS
NT = length(used_inds);
cd(anal_dir)
load(eye_data_name);

sp_dx = 0.05;
inferred_eyepos = inferred_eyepos_cor/sp_dx;
inferred_eyepos = round(interp1(inferred_eyepos_taxis,inferred_eyepos,all_t_axis(used_inds)));

interp_eyepos_cor = zeros(length(all_stim_times),1);
interp_eyepos_cor(used_inds) = inferred_eyepos;
interp_eyepos_cor(isnan(interp_eyepos_cor)) = 0;

all_left_stim_mat_cor = all_left_stim_mat;
all_right_stim_mat_cor = all_right_stim_mat;
for ii = 1:NT
    all_left_stim_mat_cor(used_inds(ii),:) = shift_matrix_Nd(all_left_stim_mat_cor(used_inds(ii),:),...
        -interp_eyepos_cor(used_inds(ii)),2);
    all_right_stim_mat_cor(used_inds(ii),:) = shift_matrix_Nd(all_right_stim_mat_cor(used_inds(ii),:),...
        -interp_eyepos_cor(used_inds(ii)),2);
end
all_left_Xmat_cor = create_time_embedding(all_left_stim_mat_cor,stim_params);
all_right_Xmat_cor = create_time_embedding(all_right_stim_mat_cor,stim_params);
all_left_Xmat_cor = all_left_Xmat_cor(:,use_kInds);
all_right_Xmat_cor = all_right_Xmat_cor(:,use_kInds);
full_Xmat_cor = [all_left_Xmat_cor all_right_Xmat_cor];

all_left_Xmat = create_time_embedding(all_left_stim_mat,stim_params);
all_right_Xmat = create_time_embedding(all_right_stim_mat,stim_params);
all_left_Xmat = all_left_Xmat(:,use_kInds);
all_right_Xmat = all_right_Xmat(:,use_kInds);
full_Xmat = [all_left_Xmat all_right_Xmat];

%% SUA STA/STC ANALYSIS

dense_block_set = find(expt_dds(cur_block_set) == 67);
sparse_block_set = find(expt_dds(cur_block_set) == 12);
dense_inds = used_inds(ismember(all_blockvec(used_inds),dense_block_set));
sparse_inds = used_inds(ismember(all_blockvec(used_inds),sparse_block_set));

dense_std = mean(std(full_Xmat_cor(dense_inds,:)));
sparse_std = mean(std(full_Xmat_cor(sparse_inds,:)));

full_Xmat_cor(dense_inds,:) = full_Xmat_cor(dense_inds,:)/dense_std;
full_Xmat_cor(sparse_inds,:) = full_Xmat_cor(sparse_inds,:)/sparse_std;

stc_thresh = -5e-3;
nneg = 10;npos = 10;
for ss = 1:length(su_probes)
    %%
    fprintf('Computing STC for SU %d of %d\n',ss,length(su_probes));
    
    cur_used_blocks = find(su_used_blocks(:,ss)); %blocks when SU
    cur_poss_inds = find(ismember(all_blockvec,cur_used_blocks));
    cur_used_inds = used_inds(ismember(used_inds,cur_poss_inds));
    cur_NT = length(cur_used_inds);
    
    cur_used_sparse_inds = cur_used_inds(ismember(cur_used_inds,sparse_inds));
    cur_used_dense_inds = cur_used_inds(ismember(cur_used_inds,dense_inds));
    fprintf('%d sparse and %d dense\n',length(cur_used_sparse_inds),length(cur_used_dense_inds));
    if length(cur_used_dense_inds) > length(cur_used_sparse_inds)
        cur_used_dense_inds((length(cur_used_sparse_inds)+1):end) = [];
    else
        cur_used_sparse_inds((length(cur_used_dense_inds)+1):end) = [];
    end
    
    if ~isempty(cur_used_sparse_inds)
      
        Robs = all_binned_spikes(cur_used_sparse_inds,ss);
        avg_rate = mean(Robs);
        
%         spikebins = convert_to_spikebins(Robs);
%         spike_cond_stim = full_Xmat(cur_used_sparse_inds(spikebins),:);
%         sta      = mean(spike_cond_stim) - mean(full_Xmat(cur_used_sparse_inds,:));
%         sta = sta/norm(sta);
%         proj_mat = sta'/(sta*sta')*sta;
%         stim_proj = full_Xmat(cur_used_sparse_inds,:) - full_Xmat(cur_used_sparse_inds,:)*proj_mat;
%         % stim_proj = stim_emb;
%         stvcv = cov(stim_proj(spikebins,:));  utvcv = cov(stim_proj);
%         [evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);
%         stcs  = evecs(:,[nneg:-1:1,length(evs)-npos+1:end]); stcs  = stcs(:,end:-1:1);
%         
%         sua_data(ss).sparse_sta_unc = sta;
%         sua_data(ss).sparse_stc_unc = stcs;
%         sua_data(ss).sparse_evals_unc = diag(evals);
 
        spikebins = convert_to_spikebins(Robs);
        spike_cond_stim = full_Xmat_cor(cur_used_sparse_inds(spikebins),:);
        sta      = mean(spike_cond_stim) - mean(full_Xmat_cor(cur_used_sparse_inds,:));
        sta_norm = norm(sta);
        sta = sta/norm(sta);
        proj_mat = sta'/(sta*sta')*sta;
        stim_proj = full_Xmat_cor(cur_used_sparse_inds,:) - full_Xmat_cor(cur_used_sparse_inds,:)*proj_mat;
        % stim_proj = stim_emb;
        stvcv = cov(stim_proj(spikebins,:));  utvcv = cov(stim_proj);
        [evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);
        stcs  = evecs(:,[nneg:-1:1,length(evs)-npos+1:end]); stcs  = stcs(:,end:-1:1);
        
        sua_data(ss).sparse_sta = sta;
        sua_data(ss).sparse_stc = stcs;
        sua_data(ss).sparse_evals = diag(evals);
        sua_data(ss).sparse_sta_norm = sta_norm;
%         sua_data(ss).sparse_maxeval = sua_data(ss).sparse_evals(end);
        

%         cur_NT = length(cur_used_sparse_inds);
%         half_set = randperm(cur_NT);
%         half_set = half_set(1:round(cur_NT/10));
%         
%         Robs = all_binned_spikes(cur_used_sparse_inds(half_set),ss);
%         spikebins = convert_to_spikebins(Robs);
%         spike_cond_stim = full_Xmat_cor(cur_used_sparse_inds(half_set(spikebins)),:);
%         sta      = mean(spike_cond_stim) - mean(full_Xmat_cor(cur_used_sparse_inds(half_set),:));
%         sta_norm = norm(sta);
%         sta = sta/norm(sta);
% %         proj_mat = sta'/(sta*sta')*sta;
% %         stim_proj = full_Xmat_cor(cur_used_sparse_inds(half_set),:) - full_Xmat_cor(cur_used_sparse_inds(half_set),:)*proj_mat;
% %         % stim_proj = stim_emb;
% %         stvcv = cov(stim_proj(spikebins,:));  utvcv = cov(stim_proj);
% %         [evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);
% %         stcs  = evecs(:,[nneg:-1:1,length(evs)-npos+1:end]); stcs  = stcs(:,end:-1:1);
%         sua_data(ss).sparse_sta_half = sta;
%         sua_data(ss).sparse_stc_half = stcs;
%         sua_data(ss).sparse_evals_half = diag(evals);
%          sua_data(ss).sparse_sta_norm_half = sta_norm;
% %         sua_data(ss).sparse_maxeval_half = sua_data(ss).sparse_evals_half(end);
       
    end
    
    if ~isempty(cur_used_dense_inds)
        Robs = all_binned_spikes(cur_used_dense_inds,ss);
        avg_rate = mean(Robs);
        
%         spikebins = convert_to_spikebins(Robs);
%         spike_cond_stim = full_Xmat(cur_used_dense_inds(spikebins),:);
%         sta      = mean(spike_cond_stim) - mean(full_Xmat(cur_used_dense_inds,:));
%         sta = sta/norm(sta);
%         proj_mat = sta'/(sta*sta')*sta;
%         stim_proj = full_Xmat(cur_used_dense_inds,:) - full_Xmat(cur_used_dense_inds,:)*proj_mat;
%         % stim_proj = stim_emb;
%         stvcv = cov(stim_proj(spikebins,:));  utvcv = cov(stim_proj);
%         [evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);
%         stcs  = evecs(:,[nneg:-1:1,length(evs)-npos+1:end]); stcs  = stcs(:,end:-1:1);
%         
%         sua_data(ss).dense_sta_unc = sta;
%         sua_data(ss).dense_stc_unc = stcs;
%         sua_data(ss).dense_evals_unc = diag(evals);
%  
                spikebins = convert_to_spikebins(Robs);
        spike_cond_stim = full_Xmat_cor(cur_used_dense_inds(spikebins),:);
        sta      = mean(spike_cond_stim) - mean(full_Xmat_cor(cur_used_dense_inds,:));
        sta_norm = norm(sta);
        sta = sta/norm(sta);
        proj_mat = sta'/(sta*sta')*sta;
        stim_proj = full_Xmat_cor(cur_used_dense_inds,:) - full_Xmat_cor(cur_used_dense_inds,:)*proj_mat;
        % stim_proj = stim_emb;
        stvcv = cov(stim_proj(spikebins,:));  utvcv = cov(stim_proj);
        [evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);
        stcs  = evecs(:,[nneg:-1:1,length(evs)-npos+1:end]); stcs  = stcs(:,end:-1:1);
        
        sua_data(ss).dense_sta = sta;
        sua_data(ss).dense_stc = stcs;
        sua_data(ss).dense_evals = diag(evals);
         sua_data(ss).dense_sta_norm = sta_norm;
%         sua_data(ss).dense_maxeval = sua_data(ss).dense_evals(end);

%         cur_NT = length(cur_used_dense_inds);
%         half_set = randperm(cur_NT);
%         half_set = half_set(1:round(cur_NT/10));
%         
%         Robs = all_binned_spikes(cur_used_dense_inds(half_set),ss);
%         spikebins = convert_to_spikebins(Robs);
%         spike_cond_stim = full_Xmat_cor(cur_used_dense_inds(half_set(spikebins)),:);
%         sta      = mean(spike_cond_stim) - mean(full_Xmat_cor(cur_used_dense_inds(half_set),:));
%         sta_norm = norm(sta);
%         sta = sta/norm(sta);
%         proj_mat = sta'/(sta*sta')*sta;
%         stim_proj = full_Xmat_cor(cur_used_dense_inds(half_set),:) - full_Xmat_cor(cur_used_dense_inds(half_set),:)*proj_mat;
%         % stim_proj = stim_emb;
%         stvcv = cov(stim_proj(spikebins,:));  utvcv = cov(stim_proj);
%         [evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);
%         stcs  = evecs(:,[nneg:-1:1,length(evs)-npos+1:end]); stcs  = stcs(:,end:-1:1);
%         sua_data(ss).dense_sta_half = sta;
%         sua_data(ss).dense_stc_half = stcs;
%         sua_data(ss).dense_evals_half = diag(evals);
%           sua_data(ss).dense_sta_norm_half = sta_norm;
%         sua_data(ss).dense_maxeval_half = sua_data(ss).dense_evals_half(end);
   end
    
end

%% SUA STA/STC ANALYSIS FOR MUA

% dense_block_set = find(expt_dds(cur_block_set) == 67);
% sparse_block_set = find(expt_dds(cur_block_set) == 12);
% dense_inds = used_inds(ismember(all_blockvec(used_inds),dense_block_set));
% sparse_inds = used_inds(ismember(all_blockvec(used_inds),sparse_block_set));

stc_thresh = -5e-3;
nneg = 10;npos = 10;
n_squared_filts = 2;

nmm_stim_params = NMMcreate_stim_params([flen 2*use_nPix],dt);
nmm_stim_params(2) = NMMcreate_stim_params([length(cur_block_set),1],dt);

%create L2 mat for left eye
L2_params = create_L2_params([],[1 flen*use_nPix],[flen use_nPix],2,3,[Inf 0]);
L2_mat = generate_L2_mat(L2_params,2*flen*use_nPix);
%add L2 mat for sine term
L2_params = create_L2_params([],[flen*use_nPix+1 2*flen*use_nPix],[flen use_nPix],2,3,[Inf 0]);
L2_mat = L2_mat + generate_L2_mat(L2_params,2*flen*use_nPix);

lambda_d2XT = 200;
lambda_L1 = 10;
sq_sc_fac = 3;
block_L2 = 1;
sub_samp_fac = 10;
silent = 1;
optim_p.optTol = 1e-5; optim_p.progTol = 1e-9;

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
    
    cur_used_sparse_inds = cur_used_inds(ismember(cur_used_inds,sparse_inds));
    cur_used_dense_inds = cur_used_inds(ismember(cur_used_inds,dense_inds));
    fprintf('%d sparse and %d dense\n',length(cur_used_sparse_inds),length(cur_used_dense_inds));
    if length(cur_used_dense_inds) > length(cur_used_sparse_inds)
        cur_used_dense_inds((length(cur_used_sparse_inds)+1):end) = [];
    else
        cur_used_sparse_inds((length(cur_used_dense_inds)+1):end) = [];
    end
    
    if ~isempty(cur_used_sparse_inds)
        Robs = all_binned_spikes(cur_used_sparse_inds,ss);
        avg_rate = mean(Robs);
        
        spikebins = convert_to_spikebins(Robs);
        spike_cond_stim = full_Xmat_cor(cur_used_sparse_inds(spikebins),:);
        sta      = mean(spike_cond_stim) - mean(full_Xmat_cor(cur_used_sparse_inds,:));
        sta_norm = norm(sta);
        sta = sta/norm(sta);
        proj_mat = sta'/(sta*sta')*sta;
        stim_proj = full_Xmat_cor(cur_used_sparse_inds,:) - full_Xmat_cor(cur_used_sparse_inds,:)*proj_mat;
        % stim_proj = stim_emb;
        stvcv = cov(stim_proj(spikebins,:));  utvcv = cov(stim_proj);
        [evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);
        stcs  = evecs(:,[nneg:-1:1,length(evs)-npos+1:end]); stcs  = stcs(:,end:-1:1);
        outer_evals = diag(evals);
        outer_evals = outer_evals([end]);
        
%         mod_signs = [1 ones(1,n_squared_filts) 1];
%         NL_types = [{'lin'} repmat({'quad'},1,n_squared_filts) {'lin'}];
%         Xtargets = [ones(1,n_squared_filts+1) 2];
%         d2XT = [lambda_d2XT repmat(sqrt(lambda_d2XT),1,n_squared_filts) 0];
%         L2 = [zeros(1,n_squared_filts+1) block_L2];
%         reg_params = NMMcreate_reg_params('lambda_custom',d2XT,'lambda_L2',L2);
%         init_filts = cell(length(mod_signs),1); init_filts{end} = zeros(size(Xblock,2),1);
%         fit0 = NMMinitialize_model(nmm_stim_params,mod_signs,NL_types,reg_params,Xtargets,init_filts); %initialize NIM
%         cur_NT = length(cur_used_sparse_inds);
%         sub_sample = randperm(cur_NT);
%         sub_sample = sub_sample(1:round(cur_NT/sub_samp_fac));
%         cur_Xmat{1} = full_Xmat_cor(cur_used_sparse_inds(sub_sample),:); cur_Xmat{2} = Xblock(cur_used_sparse_inds(sub_sample),:);
%         fit0 = NMMfit_filters(fit0,Robs(sub_sample),cur_Xmat,[],[],silent,[],L2_mat); %fit stimulus filters
%         
%         fit0 = NMMadjust_regularization(fit0,1:(1+n_squared_filts),'lambda_L1',[lambda_L1 sqrt(lambda_L1) sqrt(lambda_L1)]);
%         cur_Xmat{1} = full_Xmat_cor(cur_used_sparse_inds,:); cur_Xmat{2} = Xblock(cur_used_sparse_inds,:);
%         fit0 = NMMfit_filters(fit0,Robs,cur_Xmat,[],[],silent,[],L2_mat); %fit stimulus filters
        
        mua_data(ss).sparse_sta = sta;
        mua_data(ss).sparse_stc = stcs;
        mua_data(ss).sparse_evals = diag(evals);
        mua_data(ss).sparse_sta_norm = sta_norm;
        mua_data(ss).sparse_outer_evals = outer_evals;
%         mua_data(ss).sparse_nimFit = fit0;
    end
    
    if ~isempty(cur_used_dense_inds)
        Robs = all_binned_spikes(cur_used_dense_inds,ss);
        avg_rate = mean(Robs);
        
        spikebins = convert_to_spikebins(Robs);
        spike_cond_stim = full_Xmat_cor(cur_used_dense_inds(spikebins),:);
        sta      = mean(spike_cond_stim) - mean(full_Xmat_cor(cur_used_dense_inds,:));
        sta_norm = norm(sta);
        sta = sta/norm(sta);
        proj_mat = sta'/(sta*sta')*sta;
        stim_proj = full_Xmat_cor(cur_used_dense_inds,:) - full_Xmat_cor(cur_used_dense_inds,:)*proj_mat;
        % stim_proj = stim_emb;
        stvcv = cov(stim_proj(spikebins,:));  utvcv = cov(stim_proj);
        [evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);
        stcs  = evecs(:,[nneg:-1:1,length(evs)-npos+1:end]); stcs  = stcs(:,end:-1:1);
                outer_evals = diag(evals);
        outer_evals = outer_evals([end]);

        
%         mod_signs = [1 ones(1,n_squared_filts) 1];
%         NL_types = [{'lin'} repmat({'quad'},1,n_squared_filts) {'lin'}];
%         Xtargets = [ones(1,n_squared_filts+1) 2];
%         d2XT = [lambda_d2XT repmat(sqrt(lambda_d2XT),1,n_squared_filts) 0];
%         L2 = [zeros(1,n_squared_filts+1) block_L2];
%         reg_params = NMMcreate_reg_params('lambda_custom',d2XT,'lambda_L2',L2);
%         init_filts = cell(length(mod_signs),1); init_filts{end} = zeros(size(Xblock,2),1);
%         fit0 = NMMinitialize_model(nmm_stim_params,mod_signs,NL_types,reg_params,Xtargets,init_filts); %initialize NIM
%         cur_NT = length(cur_used_dense_inds);
%         sub_sample = randperm(cur_NT);
%         sub_sample = sub_sample(1:round(cur_NT/sub_samp_fac));
%         cur_Xmat{1} = full_Xmat_cor(cur_used_dense_inds(sub_sample),:); cur_Xmat{2} = Xblock(cur_used_dense_inds(sub_sample),:);
%         fit0 = NMMfit_filters(fit0,Robs(sub_sample),cur_Xmat,[],[],silent,[],L2_mat); %fit stimulus filters
%         
%         fit0 = NMMadjust_regularization(fit0,1:(1+n_squared_filts),'lambda_L1',[lambda_L1 sqrt(lambda_L1) sqrt(lambda_L1)]);
%         cur_Xmat{1} = full_Xmat_cor(cur_used_dense_inds,:); cur_Xmat{2} = Xblock(cur_used_dense_inds,:);
%         fit0 = NMMfit_filters(fit0,Robs,cur_Xmat,[],[],silent,[],L2_mat); %fit stimulus filters
        
        mua_data(ss).dense_sta = sta;
        mua_data(ss).dense_stc = stcs;
        mua_data(ss).dense_evals = diag(evals);
                mua_data(ss).dense_sta_norm = sta_norm;
        mua_data(ss).dense_outer_evals = outer_evals;

%         mua_data(ss).dense_nimFit = fit0;
    end

end

%%
close all
f1 = figure('Name','sparse');
f2 = figure('Name','dense');
% f3 = figure('Name','half sparse');
% f4 = figure('Name','half dense');
for ss = 1:length(su_probes)
    
    %%
    if ~isempty(sua_data(ss).sparse_sta)
        fprintf('SU %d of %d\n',ss,length(su_probes));
        
        figure(f1); clf
        cur_sta = sua_data(ss).sparse_sta; ca = max(abs(cur_sta));
        subplot(3,3,1)
        imagesc(reshape(cur_sta,flen,2*use_nPix)); caxis([-ca ca]);
          yl = ylim();
        line([use_nPix use_nPix],yl,'color','k');
        cur_sta = abs_sua(ss).sparse_sta; ca = max(abs(cur_sta));
        subplot(3,3,2)
        imagesc(reshape(cur_sta,flen,2*use_nPix)); caxis([-ca ca]);
          yl = ylim();
        line([use_nPix use_nPix],yl,'color','k');
       cur_stcs = sua_data(ss).sparse_stc; ca = max(abs(cur_stcs(:)));
        for ii = 1:3
            subplot(3,3,3+ii)
            imagesc(reshape(cur_stcs(:,ii),flen,2*use_nPix)); caxis([-ca ca]);
          yl = ylim();
        line([use_nPix use_nPix],yl,'color','k');
       end
        for ii = 1:3
            subplot(3,3,6+ii)
            imagesc(reshape(cur_stcs(:,end-nneg+ii),flen,2*use_nPix)); caxis([-ca ca]);
          yl = ylim();
        line([use_nPix use_nPix],yl,'color','k');
       end
        
        figure(f2); clf
        cur_sta = sua_data(ss).dense_sta; ca = max(abs(cur_sta));
        subplot(3,3,1)
        imagesc(reshape(cur_sta,flen,2*use_nPix)); caxis([-ca ca]);
          yl = ylim();
        line([use_nPix use_nPix],yl,'color','k');
        cur_sta = abs_sua(ss).dense_sta; ca = max(abs(cur_sta));
        subplot(3,3,2)
        imagesc(reshape(cur_sta,flen,2*use_nPix)); caxis([-ca ca]);
          yl = ylim();
        line([use_nPix use_nPix],yl,'color','k');
       cur_stcs = sua_data(ss).dense_stc; ca = max(abs(cur_stcs(:)));
        for ii = 1:3
            subplot(3,3,3+ii)
            imagesc(reshape(cur_stcs(:,ii),flen,2*use_nPix)); caxis([-ca ca]);
          yl = ylim();
        line([use_nPix use_nPix],yl,'color','k');
       end
        for ii = 1:3
            subplot(3,3,6+ii)
            imagesc(reshape(cur_stcs(:,end-nneg+ii),flen,2*use_nPix)); caxis([-ca ca]);
         yl = ylim();
        line([use_nPix use_nPix],yl,'color','k');
        end
        
%          figure(f3); clf
%         cur_sta = sua_data(ss).sparse_sta_half; ca = max(abs(cur_sta));
%         subplot(3,3,1)
%         imagesc(reshape(cur_sta,flen,2*use_nPix)); caxis([-ca ca]);
%           yl = ylim();
%         line([use_nPix use_nPix],yl,'color','k');
%        cur_stcs = sua_data(ss).sparse_stc_half; ca = max(abs(cur_stcs(:)));
%         for ii = 1:3
%             subplot(3,3,3+ii)
%             imagesc(reshape(cur_stcs(:,ii),flen,2*use_nPix)); caxis([-ca ca]);
%           yl = ylim();
%         line([use_nPix use_nPix],yl,'color','k');
%        end
%         for ii = 1:3
%             subplot(3,3,6+ii)
%             imagesc(reshape(cur_stcs(:,end-nneg+ii),flen,2*use_nPix)); caxis([-ca ca]);
%           yl = ylim();
%         line([use_nPix use_nPix],yl,'color','k');
%        end
%         
%         figure(f4); clf
%         cur_sta = sua_data(ss).dense_sta_half; ca = max(abs(cur_sta));
%         subplot(3,3,1)
%         imagesc(reshape(cur_sta,flen,2*use_nPix)); caxis([-ca ca]);
%           yl = ylim();
%         line([use_nPix use_nPix],yl,'color','k');
%        cur_stcs = sua_data(ss).dense_stc_half; ca = max(abs(cur_stcs(:)));
%         for ii = 1:3
%             subplot(3,3,3+ii)
%             imagesc(reshape(cur_stcs(:,ii),flen,2*use_nPix)); caxis([-ca ca]);
%           yl = ylim();
%         line([use_nPix use_nPix],yl,'color','k');
%        end
%         for ii = 1:3
%             subplot(3,3,6+ii)
%             imagesc(reshape(cur_stcs(:,end-nneg+ii),flen,2*use_nPix)); caxis([-ca ca]);
%          yl = ylim();
%         line([use_nPix use_nPix],yl,'color','k');
%         end
       pause
    end    
end
    
%%
close all
f1 = figure('Name','sparse');
f2 = figure('Name','dense');
for ss = 1:96
    if ~isempty(mua_data(ss).sparse_sta)
        fprintf('MU %d of %d\n',ss,96);
        
        figure(f1); clf
        cur_sta = mua_data(ss).sparse_sta; ca = max(abs(cur_sta));
        subplot(3,3,1)
        imagesc(reshape(cur_sta,flen,2*use_nPix)); caxis([-ca ca]);
        cur_stcs = mua_data(ss).sparse_stc; ca = max(abs(cur_stcs(:)));
        yl = ylim();
        line([use_nPix use_nPix],yl,'color','k');
        for ii = 1:3
            subplot(3,3,3+ii)
            imagesc(reshape(cur_stcs(:,ii),flen,2*use_nPix)); caxis([-ca ca]);
         yl = ylim();
        line([use_nPix use_nPix],yl,'color','k');
       end
        for ii = 1:3
            subplot(3,3,6+ii)
            imagesc(reshape(cur_stcs(:,end-nneg+ii),flen,2*use_nPix)); caxis([-ca ca]);
         yl = ylim();
        line([use_nPix use_nPix],yl,'color','k');
       end
        
        figure(f2); clf
        cur_sta = mua_data(ss).dense_sta; ca = max(abs(cur_sta));
        subplot(3,3,1)
        imagesc(reshape(cur_sta,flen,2*use_nPix)); caxis([-ca ca]);
        yl = ylim();
        line([use_nPix use_nPix],yl,'color','k');
        cur_stcs = mua_data(ss).dense_stc; ca = max(abs(cur_stcs(:)));
        for ii = 1:3
            subplot(3,3,3+ii)
            imagesc(reshape(cur_stcs(:,ii),flen,2*use_nPix)); caxis([-ca ca]);
        yl = ylim();
        line([use_nPix use_nPix],yl,'color','k');
        end
        for ii = 1:3
            subplot(3,3,6+ii)
            imagesc(reshape(cur_stcs(:,end-nneg+ii),flen,2*use_nPix)); caxis([-ca ca]);
         yl = ylim();
        line([use_nPix use_nPix],yl,'color','k');
       end
        
        pause
        
    end
    
end
