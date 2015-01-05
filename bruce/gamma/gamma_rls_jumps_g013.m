clear all
% close all
addpath('~/James_scripts/bruce/eye_tracking/');

%run on 86 [0, 90];
Expt_num = 103;
Expt_name = sprintf('G%.3d',Expt_num);

use_LOOXV = 0; %[0 no LOOXV; 1 SU LOOXV; 2 all LOOXV]
bar_ori = 0; %bar orientation to use (only for UA recs)

recompute_init_mods = 0; %use existing initial models?
use_measured_pos = 0; %1 for init with coils, 2 for init with trial-sub coils, 0 for init to perfect fixation
use_sac_kerns = 1; %use sac-modulation kernels


data_dir = ['~/Data/bruce/' Expt_name];
cd(data_dir);

Expt_num = str2num(Expt_name(2:end));
if Expt_name(1) == 'M'
    rec_type = 'LP';
elseif Expt_name(1) == 'G'
    rec_type = 'UA';
end

if strcmp(rec_type,'LP')
    switch Expt_num
        case 266
            bar_ori = 80;
        case 270
            bar_ori = 60;
        case 275
            bar_ori = 135;
        case 277
            bar_ori = 70;
        case 281
            bar_ori = 140;
        case 287
            bar_ori = 90;
        case 289
            bar_ori = 160;
        case 294
            bar_ori = 40;
    end
end

if Expt_num >= 280
    data_dir = ['/media/NTlab_data3/Data/bruce/' Expt_name];
else
    data_dir = ['~/Data/bruce/' Expt_name];
end

cd(data_dir);

if strcmp(rec_type,'LP')
    if Expt_num >= 275
        rpt_seed = 1001; %M275 M277 M281
    else
        rpt_seed = 1e4; %m270 and 266
    end
    load(sprintf('lem%sExpts.mat',Expt_name)); %load in Expts struct
    good_coils = [1 1]; %which coils are usable
    use_coils = [1 1]; %[L R] Use info from coils?
    n_probes = 24;
    use_nPix = 32;
elseif strcmp(rec_type,'UA')
    rpt_seed = nan;
    load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct
    good_coils = [1 0]; %which coils are usable
    use_coils = [0 0]; %[L R] Use info from coils?
    n_probes = 96;
    use_nPix = 16;
end

load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data

anal_dir = ['~/Analysis/bruce/' Expt_name '/ET_final_imp/'];
if ~exist(anal_dir,'dir')
    system(['mkdir ' anal_dir]);
end
cluster_dir = ['~/Analysis/bruce/' Expt_name '/clustering'];

mod_data_name = 'full_eyetrack_initmods';
anal_name = 'full_eyetrack';

if bar_ori == 90 & strcmp(rec_type,'UA')
    mod_data_name = [mod_data_name '_vbars'];
    anal_name = [anal_name '_vbars'];
end

%dont fit stim models using these blocks
ignore_blocks = [];
switch Expt_num
    case 270
        ignore_blocks = [5 19];
    case 289
        ignore_blocks = [27]; %27 is off somehow
    case 294
        ignore_blocks = [37 38 39]; %37-39 have slightly different dw used in these blocks
    case 86
        ignore_blocks = [16 17 28 30];
    case 87
        ignore_blocks = [15];
    case 93
        ignore_blocks = [28];
end

%problem with M270 where vd was wrong, need this correction factor to get
%correct units
if Expt_num==270
    scale_fac = 1.72;
else
    scale_fac = 1;
end

%if using coil info
if any(use_coils > 0)
    anal_name = [anal_name '_Cprior'];
end
%if using coil initialization
if use_measured_pos == 1
    mod_data_name = [mod_data_name '_Cinit'];
    anal_name = [anal_name '_Cinit'];
end
%if using coil init with trial-sub
if use_measured_pos == 2
    mod_data_name = [mod_data_name '_CPinit'];
    anal_name = [anal_name '_CPinit'];
end

%%

min_trial_dur = 0.75;

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

%exclude data at beginning and end of each trial
beg_buffer = 0.15;
end_buffer = 0.05;
trial_dur = 4;

%%
include_expts = {'rls.sl'};
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
cur_block_set = find(included_type & ~expt_binoc' & expt_Fr == 1 & expt_bar_ori == bar_ori);


cur_block_set(ismember(cur_block_set,ignore_blocks)) = [];

n_blocks = length(cur_block_set);

%% COMPUTE TRIAL DATA
cd(data_dir);

fprintf('Computing prep data\n');
trial_cnt = 0;

all_stim_times = [];
% all_Xmat = [];
all_stim_mat = [];
all_t_axis = [];
all_t_bin_edges = [];
all_tsince_start = [];
all_blockvec = [];
all_trialvec = [];
all_trial_sl = [];
all_trial_se = [];
all_trial_blocknums = [];
all_trial_start_times = [];
all_trial_end_times = [];
all_bin_edge_pts = [];
all_spk_times = cell(n_probes,1);
all_spk_inds = cell(n_probes,1);
all_clust_ids = cell(n_probes,1);

trial_toffset = zeros(length(cur_block_set),1);
cur_spkind_offset = 0;
cur_toffset = 0;
for ee = 1:n_blocks;
    fprintf('Expt %d Block %d of %d;  UNMATCHED EXPT TYPE\n',Expt_num,ee,n_blocks);
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
    all_trial_blocknums = cat(1,all_trial_blocknums,ones(length(use_trials),1)*ee);
    
    trial_sl = [Expts{cur_block}.Trials(:).sl];
    trial_sl = trial_sl(id_inds);
    all_trial_sl = cat(1,all_trial_sl,trial_sl(use_trials)');

        trial_se = [Expts{cur_block}.Trials(:).se];
    trial_se = trial_se(id_inds);
    all_trial_se = cat(1,all_trial_se,trial_se(use_trials)');

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
        cur_t_axis = 0.5*cur_t_edges(1:end-1) + 0.5*cur_t_edges(2:end);
        
        cur_tsince_start = cur_t_axis - trial_start_times(use_trials(tt));
        
        if ~any(isnan(left_stim_mats{use_trials(tt)}(:))) && n_frames > min_trial_dur/dt
            use_frames = min(length(cur_stim_times),n_frames);
            cur_stim_mat = double(left_stim_mats{use_trials(tt)}(1:use_frames,cur_use_pix));
            %             bar_Xmat = create_time_embedding(cur_stim_mat,stim_params);
            
            if ~isempty(all_stim_times)
                if any(cur_stim_times < all_stim_times(end))
                    fprintf('Warn trial %d\n',tt);
                end
            end
            all_stim_times = [all_stim_times; cur_stim_times];
            all_t_axis = [all_t_axis; cur_t_axis];
            %             all_Xmat = [all_Xmat; bar_Xmat];
            all_stim_mat = [all_stim_mat; cur_stim_mat];
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
clust_params.n_probes = n_probes;
if strcmp(rec_type,'LP')
    clust_params.exclude_adjacent = true;
else
    clust_params.exclude_adjacent = false;
end
[all_binned_mua,all_binned_sua,Clust_data] = ...
    get_binned_spikes(cluster_dir,all_spk_times,all_clust_ids,all_spk_inds,...
    all_t_axis,all_t_bin_edges,all_bin_edge_pts,cur_block_set,all_blockvec,clust_params);
su_probes = Clust_data.SU_probes;
SU_numbers = Clust_data.SU_numbers;


%% DEFINE DATA USED FOR ANALYSIS
used_trials = find(all_trial_sl == 50);
used_inds = find(all_tsince_start >= beg_buffer & (trial_dur-all_tsince_start) >= end_buffer & ...
    ismember(all_trialvec,used_trials));
NT = length(used_inds);

%% CREATE EVENT PREDICTORS FOR REAL AND SIM SACCADES (POOLING MICRO AND MACRO SACS)
Xblock = zeros(length(all_stim_times),n_blocks);
for i = 1:n_blocks
    cur_set = find(all_blockvec==i);
    Xblock(cur_set,i) = 1;
end

%% PROCESS EYE TRACKING DATA
% [all_eye_vals,all_eye_ts,all_eye_speed,et_params] = process_ET_data(all_t_axis,all_blockvec,cur_block_set,Expt_name,trial_toffset);
[all_eye_vals,all_eye_ts,all_eye_speed,et_params] = process_ET_data_v2(all_t_axis,all_blockvec,cur_block_set,Expt_name,trial_toffset,good_coils);
interp_eye_speed = interp1(all_eye_ts,all_eye_speed,all_t_axis);

%compute corrected eye data in bar-oriented frame
[corrected_eye_vals,corrected_eye_vals_interp]  = get_corrected_ET_data(Expts(cur_block_set),all_eye_vals,all_eye_ts,...
    all_t_axis,all_blockvec,expt_bar_ori(cur_block_set),used_inds);

[saccades,et_params] = detect_saccades_v2(corrected_eye_vals,all_eye_vals,all_eye_speed,all_eye_ts,et_params);

is_blink = detect_blinks(all_eye_ts,all_eye_vals,saccades,et_params);

[saccades,is_blink] = merge_blinks(saccades,is_blink);

sac_start_times = [saccades(:).start_time];
sac_stop_times = [saccades(:).stop_time];
interp_sac_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),sac_start_times));
interp_sac_start_inds(isnan(interp_sac_start_inds)) = 1;
interp_sac_stop_inds = round(interp1(all_t_axis,1:length(all_t_axis),sac_stop_times));
interp_sac_stop_inds(isnan(interp_sac_stop_inds)) = 1;
sac_error = abs(sac_start_times - all_t_axis(interp_sac_start_inds)');
bad_sacs = find(isnan(interp_sac_start_inds) | sac_error > dt);
sac_start_times(bad_sacs) = [];
sac_stop_times(bad_sacs) = [];
saccades(bad_sacs) = [];
is_blink(bad_sacs) = [];
interp_sac_start_inds(bad_sacs) = [];
interp_sac_stop_inds(bad_sacs) = [];

%% CREATE SACCADE PREDICTOR MATS
saccade_start_inds = find(ismember(used_inds,interp_sac_start_inds));
used_saccade_set = find(ismember(interp_sac_start_inds,used_inds));
%nearest index in the used data set of the saccade stop time
saccade_stop_inds = round(interp1(used_inds,1:length(used_inds),interp_sac_stop_inds(used_saccade_set)))';

saccade_stop_inds(isnan(saccade_stop_inds)) = length(used_inds);
used_is_blink = is_blink(used_saccade_set);

sac_amps = [saccades(:).amplitude];
is_micro = sac_amps(used_saccade_set) < 1;
big_sacs = find(~is_micro & ~used_is_blink');
micro_sacs = find(is_micro & ~used_is_blink');
sac_durs = [saccades(:).duration];

saccade_trial_inds = all_trialvec(used_inds(saccade_start_inds));

%% LOAD LFP DATA
cd(data_dir)
if strcmp(rec_type,'LP')
    Fs = 1000;
    dsf = 5;
    Fsd = Fs/dsf;
    niqf = Fs/2;
    new_niqf = Fsd/2;
    [bb,aa] = butter(2,0.8*new_niqf/niqf);
    use_lfps = 1:n_probes;
elseif strcmp(rec_type,'UA')
    Fs = 400.0032;
    dsf = 2;
    Fsd = Fs/dsf;
    niqf = Fs/2;
    new_niqf = Fsd/2;
    [bb,aa] = butter(2,0.8*new_niqf/niqf);
    use_lfps = 1:4:n_probes;
    %     use_lfps = [28 SU_probes];
end

full_lfps = [];
full_lfp_taxis = [];
cur_toffset = 0;
ublock_set = 1:length(cur_block_set);
if Expt_num == 275
    bad_blocks = 9;
    ublock_set(bad_blocks) = [];
end
for ee = ublock_set
    
    fprintf('Loading LFPs, Expt %d of %d\n',ee,length(cur_block_set));
    
    if strcmp(rec_type,'LP')
        fname = sprintf('lemM%dA.%d.lfp.mat',Expt_num,cur_block_set(ee));
        load(fname);
        
        tlens = arrayfun(@(X) length(X.ftime),LFP.Trials);
        bad_trials = find(tlens == 0);
        LFP.Trials(bad_trials) = [];
        n_trials(ee) = length(LFP.Trials);
        lfp_trial_starts = [LFP.Trials(:).ftime]/1e4;
        lfp_trial_ends = [LFP.Trials(:).End]/1e4;
        expt_lfp_t_axis = [];
        expt_lfps = [];
        for tt = 1:n_trials(ee)
            %         tt
            cur_npts = size(LFP.Trials(tt).LFP,1);
            cur_t_end(tt) = lfp_trial_starts(tt)+(cur_npts-1)/Fs;
            cur_t_axis = (lfp_trial_starts(tt):1/Fs:cur_t_end(tt)) + cur_toffset;
            
            if ~isempty(expt_lfp_t_axis)
                cur_sp = find(cur_t_axis > max(expt_lfp_t_axis),1,'first');
            else
                cur_sp = 1;
            end
            cur_t_axis = cur_t_axis(cur_sp:end);
            if length(cur_t_axis) > 50 & size(LFP.Trials(tt).LFP,2) == n_probes
                cur_LFP = double([LFP.Trials(tt).LFP]);
                cur_LFP = cur_LFP(cur_sp:end,use_lfps);
                if dsf > 1
                    cur_LFP = filtfilt(bb,aa,cur_LFP);
                    cur_LFP = downsample(cur_LFP,dsf);
                    cur_t_axis = downsample(cur_t_axis,dsf);
                end
                if size(cur_LFP,2) == length(use_lfps)
                    expt_lfp_t_axis = [expt_lfp_t_axis; cur_t_axis(:)];
                    expt_lfps = [expt_lfps; cur_LFP];
                end
            end
        end
    else
        lfp_fname = sprintf('Expt%d_LFP.mat',cur_block_set(ee));
        load(lfp_fname);
        
        cur_lfps = bsxfun(@times,double(lfp_mat(:,use_lfps)),lfp_int2V(use_lfps)');
        if dsf > 1
            cur_lfps = filtfilt(bb,aa,cur_lfps);
            expt_lfps = downsample(cur_lfps,dsf);
            expt_lfp_t_axis = downsample(lfp_t_ax',dsf);
        end
    end
    
    cur_uset = find(all_blockvec == ee);
    if ~isempty(cur_uset)
        uinds = find(expt_lfp_t_axis >= all_t_axis(cur_uset(1)) & expt_lfp_t_axis <= all_t_axis(cur_uset(end)));
        full_lfps = cat(1,full_lfps,expt_lfps(uinds,:));
        full_lfp_taxis = cat(1,full_lfp_taxis,expt_lfp_t_axis(uinds));
    end
    cur_toffset = trial_toffset(ee);
end

%%
interp_lfps = interp1(full_lfp_taxis,full_lfps,all_t_axis_us);
% interp_lfps(ismember(all_blockvec,bad_blocks),:) = nan;
interp_lfps = interp_lfps(used_inds_us,:);
interp_lfps = nanzscore(interp_lfps);
