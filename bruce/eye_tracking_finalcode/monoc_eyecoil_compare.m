clear all
close all
addpath('~/James_scripts/bruce/eye_tracking/');
dir_prefix = '/home/james';
% dir_prefix = '/Volumes/james';

fig_dir = [dir_prefix '/Analysis/bruce/ET_final/'];
single_bar_expts = [232 235 239];

%%
% Expt_name ='G086';
% Expt_name ='M266';
Expt_name ='M275';
Expt_num = str2num(Expt_name(2:end));

% bar_ori = 0;
% bar_ori = 0;
% bar_ori = 135;
bar_ori = 135;
% bar_ori = 60;

use_sac_kerns = 1;
use_coils = [0 0]; %[L R]

flen = 12;
use_nPix = 16;
min_trial_dur = 0.75;
spatial_usfac = 4;
old_usfac = 2;

data_dir = [dir_prefix '/Data/bruce/' Expt_name];
cd(data_dir);

if Expt_name(1) == 'G'
    load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct
else
    load(sprintf('lem%sExpts.mat',Expt_name)); %load in Expts struct
end
if ~ismember(Expt_num,single_bar_expts)
    load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data
else
    load ./random_bar_eyedata_ftime.mat bar_expts
    load ./bar_params.mat
end

anal_dir = [dir_prefix '/Analysis/bruce/' Expt_name '/ET_final/'];
if ~exist(anal_dir,'dir')
    system(['mkdir ' anal_dir]);
end

ignore_blocks = [];
%dont fit stim models using these blocks
if Expt_num == 86
    ignore_blocks = [16 17 28 30]; %G086
elseif Expt_num == 87
    ignore_blocks = [15];
elseif Expt_num == 93
    ignore_blocks = [28];
end
%dont fit stim models using these blocks
if Expt_num == 270
    ignore_blocks = [5 19];
elseif Expt_num == 275
    ignore_blocks = 15;
end
%dont fit stim models using these blocks
if Expt_num == 235
    ignore_blocks = [51]; %G086
elseif Expt_num == 239
    ignore_blocks = [40];
end

if Expt_num==270
    scale_fac = 1.72;
else
    scale_fac = 1;
end

%%
stim_fs = 100; %in Hz
dt = 0.01;
full_nPix = 36;
Fr = 1;

beg_buffer = 0.2;
end_buffer = 0.05;
if ~ismember(Expt_num,single_bar_expts)
    trial_dur = 4;
else
    trial_dur = 2;
end

use_right_eye = false;

n_use_blocks = Inf;

use_nPix_us = use_nPix*spatial_usfac;
klen_us = use_nPix_us*flen;

sac_backlag = round(0.05/dt);
sac_forlag = round(0.3/dt);
sac_bincents = -sac_backlag:sac_forlag;
n_sac_bins = length(sac_bincents);

% if ~ismember(Expt_num,single_bar_expts)
%     sp_dx = 0.0565/spatial_usfac/scale_fac;
% else
%     sp_dx = 0.125/4;
% end

%%
if strcmp(Expt_name,'G093')
    include_expts = {'rls.FaXwi','rls.FaXwiXimi'};
else
    include_expts = {'rls.Fa', 'rls.FaXimi'};
end
if Expt_name(1) == 'M';
    include_expts = {'rls.Fa', 'rls.FaXimi','rls.FaXFaXFs','rls.AllSac','rls.imiXFa'};
end
if ~ismember(Expt_num,single_bar_expts)
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
    if Expt_num==275
        expt_bar_ori(expt_bar_ori > 1e4) = bar_ori;
    end
    
    cur_block_set = find(included_type & ~expt_binoc' & expt_Fr == 1 & expt_bar_ori == bar_ori);
    cur_block_set(ismember(cur_block_set,ignore_blocks)) = [];
else
    cur_block_set = bar_expts;
    cur_block_set(ismember(cur_block_set,ignore_blocks)) = [];
end
n_blocks = length(cur_block_set);

%%
if ismember(Expt_num,single_bar_expts)
    sp_dx = 0.125/4;
else
base_sp_dx = Expts{cur_block_set(1)}.Stimvals.dw;
sp_dx = base_sp_dx/spatial_usfac/scale_fac;
end
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
all_trial_wi = [];
all_trial_blocknums = [];
all_trial_start_times = [];
all_trial_end_times = [];
all_bin_edge_pts = [];

trial_toffset = zeros(length(cur_block_set),1);
cur_spkind_offset = 0;
cur_toffset = 0;
for ee = 1:n_blocks;
    fprintf('Expt %d Block %d of %d\n',Expt_num,ee,n_blocks);
    cur_block = cur_block_set(ee);
    
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
    all_trial_start_times = cat(1,all_trial_start_times,trial_start_times(use_trials)' + cur_toffset);
    all_trial_end_times = cat(1,all_trial_end_times,trial_end_times(use_trials)' + cur_toffset);
    all_trial_blocknums = cat(1,all_trial_blocknums,ones(length(use_trials),1)*ee);
    
    if strcmp(Expt_name,'G093')
        trial_wi = [Expts{cur_block}.Trials(:).wi];
        trial_wi = trial_wi(id_inds);
        all_trial_wi = cat(1,all_trial_wi,trial_wi(use_trials)');
    end
    
    if ~ismember(Expt_num,single_bar_expts)
        fname = sprintf('%s/stims/Expt%d_stim',data_dir,cur_block);
        load(fname);
        buffer_pix = floor((expt_npix(cur_block) - full_nPix)/2);
        cur_use_pix = (1:full_nPix) + buffer_pix;
    end
    
    n_trials = length(use_trials);
    for tt = 1:n_trials
        cur_stim_times = Expts{cur_block}.Trials(use_trials(tt)).Start(:)/1e4;
        if ~ismember(Expt_num,single_bar_expts)
            n_frames = size(left_stim_mats{use_trials(tt)},1);
        else
            cur_bar_Op = [Expts{cur_block}.Trials(use_trials(tt)).Op];
            n_frames = length(cur_bar_Op);
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
        
        if n_frames > min_trial_dur/dt
            use_frames = min(length(cur_stim_times),n_frames);
            all_stim_times = [all_stim_times; cur_stim_times + cur_toffset];
            all_t_axis = [all_t_axis; cur_t_axis + cur_toffset];
            all_t_bin_edges = [all_t_bin_edges; cur_t_edges + cur_toffset];
            all_tsince_start = [all_tsince_start; cur_tsince_start];
            all_blockvec = [all_blockvec; ones(size(cur_t_axis))*ee];
            all_trialvec = [all_trialvec; ones(size(cur_t_axis))*(tt + trial_cnt)];
            all_bin_edge_pts = [all_bin_edge_pts; length(all_t_bin_edges)];
        end
    end
    if Expt_name(1) == 'M'
        trial_toffset(ee) = all_t_bin_edges(end);
        cur_toffset = trial_toffset(ee);
    end
    trial_cnt = trial_cnt + n_trials;
end

%% DEFINE DATA USED FOR ANALYSIS
used_inds = find(all_tsince_start >= beg_buffer & (trial_dur-all_tsince_start) >= end_buffer);
if strcmp(Expt_name,'G093')
    un_wi_vals = unique(all_trial_wi);
    use_wi_trials = find(all_trial_wi == un_wi_vals(2));
    used_inds = used_inds(ismember(all_trialvec(used_inds),use_wi_trials));
end
NT = length(used_inds);

%% PROCESS EYE TRACKING DATA
if Expt_name(1) == 'G'
    trial_toffset = zeros(length(cur_block_set),1);
end
[all_eye_vals,all_eye_ts,all_eye_speed,et_params] = process_ET_data(all_t_axis,all_blockvec,cur_block_set,Expt_name,trial_toffset);
interp_eye_speed = interp1(all_eye_ts,all_eye_speed,all_t_axis);

%compute corrected eye data in bar-oriented frame
[corrected_eye_vals,corrected_eye_vals_interp]  = get_corrected_ET_data(Expts(cur_block_set),all_eye_vals,all_eye_ts,...
    all_t_axis,all_blockvec,bar_ori*ones(length(cur_block_set),1),used_inds);

[saccades,et_params] = detect_saccades(corrected_eye_vals,all_eye_speed,all_eye_ts,et_params);

if Expt_name(1) == 'M' && ~ismember(Expt_num,single_bar_expts)
    par_thresh = 5;
    orth_thresh = 1.25;
else
    par_thresh = 4;
    orth_thresh = 1;
end
[out_of_range] = detect_bad_fixation(corrected_eye_vals_interp,all_trialvec,used_inds,par_thresh,orth_thresh);
fract_out = length(out_of_range)/length(used_inds);
fprintf('Eliminating %.4f of data out of window\n',fract_out);
used_inds(ismember(used_inds,out_of_range)) = [];
NT = length(used_inds);

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
interp_sac_start_inds(bad_sacs) = [];
interp_sac_stop_inds(bad_sacs) = [];

%% CREATE SACCADE PREDICTOR MATS
saccade_start_inds = find(ismember(used_inds,interp_sac_start_inds));
used_saccade_set = find(ismember(interp_sac_start_inds,used_inds));
%nearest index in the used data set of the saccade stop time
saccade_stop_inds = round(interp1(used_inds,1:length(used_inds),interp_sac_stop_inds(used_saccade_set)))';

saccade_stop_inds(isnan(saccade_stop_inds)) = length(used_inds);

sac_amps = [saccades(:).amplitude];
is_micro = sac_amps(used_saccade_set) < 1;
big_sacs = find(~is_micro);
micro_sacs = find(is_micro);

saccade_trial_inds = all_trialvec(used_inds(saccade_start_inds));

%% DEFINE FIXATION POINTS
trial_start_inds = [1; find(diff(all_trialvec(used_inds)) ~= 0) + 1];
trial_end_inds = [find(diff(all_trialvec(used_inds)) ~= 0); NT];

fix_start_inds = sort([trial_start_inds; saccade_stop_inds]);
fix_stop_inds = sort([trial_end_inds; saccade_start_inds]);
fix_durs = (fix_stop_inds-fix_start_inds)*dt;
fix_start_inds(fix_durs<=0) = [];
fix_stop_inds(fix_durs <= 0) = [];
fix_durs(fix_durs <= 0) = [];
n_fixs = length(fix_start_inds);

%push the effects of saccades forward in time
sac_shift = round(0.05/dt);
pfix_start_inds = fix_start_inds;
pfix_stop_inds = fix_stop_inds;
for i = 1:length(fix_start_inds)
    next_trial = trial_start_inds(find(trial_start_inds >= fix_start_inds(i),1,'first'));
    if next_trial > fix_start_inds(i) + sac_shift
        pfix_start_inds(i) = fix_start_inds(i) + sac_shift;
    end
    next_trial = trial_start_inds(find(trial_start_inds >= fix_stop_inds(i),1,'first'));
    if next_trial > fix_stop_inds(i) + sac_shift
        pfix_stop_inds(i) = fix_stop_inds(i) + sac_shift;
    end
end

%for all times within the (forward-projected) saccade, use 'prior'
%state-transition model
use_prior = zeros(NT,1);
for i = 1:n_fixs-1
    use_prior((pfix_stop_inds(i)+1):pfix_start_inds(i+1)) = 1;
end
fix_ids = nan(NT,1);
pfix_ids = nan(NT,1);
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):(fix_stop_inds(ii));
    fix_ids(cur_inds) = ii;
    cur_inds = pfix_start_inds(ii):(pfix_stop_inds(ii));
    pfix_ids(cur_inds) = ii;
end

trial_ids = nan(NT,1);
for ii = 1:n_trials
    cur_inds = trial_start_inds(ii):trial_end_inds(ii);
    trial_ids(cur_inds) = ii;
end

%%
% smooth_eyepos = corrected_eye_vals_interp;
% smooth_eyepos(isnan(smooth_eyepos)) = 0;
% 
% %smooth out fast transients in eye signal
% eye_smooth_sig = round(0.025/dt);
% interp_inds = [];
% for ii = 1:n_fixs
%     cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
%     if length(cur_inds) > eye_smooth_sig*5;
%         smooth_eyepos(used_inds(cur_inds),1) = jmm_smooth_1d_cor(smooth_eyepos(used_inds(cur_inds),1),eye_smooth_sig,2);
%         smooth_eyepos(used_inds(cur_inds),2) = jmm_smooth_1d_cor(smooth_eyepos(used_inds(cur_inds),2),eye_smooth_sig,2);
%         smooth_eyepos(used_inds(cur_inds),3) = jmm_smooth_1d_cor(smooth_eyepos(used_inds(cur_inds),3),eye_smooth_sig,2);
%         smooth_eyepos(used_inds(cur_inds),4) = jmm_smooth_1d_cor(smooth_eyepos(used_inds(cur_inds),4),eye_smooth_sig,2);
%     end
%     interp_inds = [interp_inds; cur_inds'];
% end
% 
% interp_inds = unique(interp_inds);
% smooth_eyepos(used_inds,:) = interp1(used_inds(interp_inds),smooth_eyepos(used_inds(interp_inds),:),used_inds);
% 
% %subtract off within-trial median
% for ii = 1:length(trial_start_inds)
%     cur_inds = trial_start_inds(ii):trial_end_inds(ii);
%     smooth_eyepos(used_inds(cur_inds),:) = bsxfun(@minus,smooth_eyepos(used_inds(cur_inds),:),median(smooth_eyepos(used_inds(cur_inds),:)));
% end
% 
% smooth_eyepos = smooth_eyepos(used_inds,:);

%% remove block-wise median from coil measured position
for bb = 1:length(cur_block_set)
%    binds = find(all_blockvec == bb);
   binds = used_inds(all_blockvec(used_inds) == bb);
   corrected_eye_vals_interp(binds,:) = bsxfun(@minus,corrected_eye_vals_interp(binds,:),nanmedian(corrected_eye_vals_interp(binds,:)));
end
%% Construct estimated eye position
clear it_fix* drift*
cd(anal_dir)

if Expt_name(1) == 'G'
    if bar_ori == 0
        old_data_name = 'monoc_eyecorr_hbar2';
        data_name = 'monoc_eyecorr_hbar_highres3';
    else
        old_data_name = 'monoc_eyecorr_vbar2';
        data_name = 'monoc_eyecorr_vbar_highres3';
    end
else
    data_name = 'monoc_eyecorr_highres3';
end
% if bar_ori == 0
%     old_data_name = 'monoc_eyecorr_hbar';
%     data_name = 'monoc_eyecorr_hbar_highres';
% else
%     old_data_name = 'monoc_eyecorr_vbar';
%     data_name = 'monoc_eyecorr_vbar_highres';
% end

fprintf('Loading %s\n',data_name);
load(data_name,'drift_*','et_tr_set','best_*');
% load(old_data_name,'it_fix_post_mean','it_fix_post_std');

% best_fix_cor = it_fix_post_mean(end,:)*spatial_usfac/old_usfac;
% best_fix_std = it_fix_post_std(end,:)*spatial_usfac/old_usfac;

fin_fix_corr = nan(NT,1);
fin_fix_std = nan(NT,1);
fin_fix_corr(~isnan(fix_ids)) = best_fix_cor(fix_ids(~isnan(fix_ids)));
fin_fix_corr = interp1(find(~isnan(fix_ids)),fin_fix_corr(~isnan(fix_ids)),1:NT);
fin_fix_std(~isnan(fix_ids)) = best_fix_std(fix_ids(~isnan(fix_ids)));
fin_fix_std = interp1(find(~isnan(fix_ids)),fin_fix_std(~isnan(fix_ids)),1:NT);

fin_fix_corr = fin_fix_corr*sp_dx;
fin_fix_std = fin_fix_std*sp_dx;

fin_drift_corr = drift_post_mean(end,:)*sp_dx;
fin_drift_std = drift_post_std(end,:)*sp_dx;
fix_inds = [];
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
    fix_inds = [fix_inds cur_inds];
%     if length(cur_inds) > sac_shift
%         fin_drift_corr(cur_inds(1:end-sac_shift+1)) = fin_drift_corr(cur_inds(sac_shift:end));
%         fin_drift_std(cur_inds(1:end-sac_shift+1)) = fin_drift_std(cur_inds(sac_shift:end));
%     end
end
for ii = 1:length(trial_start_inds)
    cur_inds = trial_start_inds(ii):trial_end_inds(ii);
    fin_drift_corr(cur_inds(1:end-sac_shift)) = fin_drift_corr(cur_inds(sac_shift+1:end));
    fin_drift_std(cur_inds(1:end-sac_shift)) = fin_drift_std(cur_inds(sac_shift+1:end));
end
fin_drift_corr = interp1(find(~isnan(fix_ids)),fin_drift_corr(~isnan(fix_ids)),1:NT);
fin_drift_std = interp1(find(~isnan(fix_ids)),fin_drift_std(~isnan(fix_ids)),1:NT);

fin_tot_corr = fin_fix_corr(:) + fin_drift_corr(:);
fin_tot_std = sqrt(fin_fix_std(:).^2 + fin_drift_std(:).^2);


%% construct estimated position with eye-coil info (if using)
% data_name ='./monoc_eyecorr_Cprior';
if spatial_usfac == 4
data_name ='./monoc_eyecorr_highres5_Cprior';
end
fprintf('Loading %s\n',data_name);
load(data_name,'drift_*','best_*');

Pfin_fix_corr = nan(NT,1);
Pfin_fix_std = nan(NT,1);
Pfin_fix_corr(~isnan(fix_ids)) = best_fix_cor(fix_ids(~isnan(fix_ids)));
Pfin_fix_corr = interp1(find(~isnan(fix_ids)),Pfin_fix_corr(~isnan(fix_ids)),1:NT);
Pfin_fix_std(~isnan(fix_ids)) = best_fix_std(fix_ids(~isnan(fix_ids)));
Pfin_fix_std = interp1(find(~isnan(fix_ids)),Pfin_fix_std(~isnan(fix_ids)),1:NT);

Pfin_fix_corr = Pfin_fix_corr*sp_dx;
Pfin_fix_std = Pfin_fix_std*sp_dx;

Pfin_drift_corr = drift_post_mean(end,:)*sp_dx;
Pfin_drift_std = drift_post_std(end,:)*sp_dx;
for ii = 1:length(trial_start_inds)
    cur_inds = trial_start_inds(ii):trial_end_inds(ii);
    Pfin_drift_corr(cur_inds(1:end-sac_shift)) = Pfin_drift_corr(cur_inds(sac_shift+1:end));
    Pfin_drift_std(cur_inds(1:end-sac_shift)) = Pfin_drift_std(cur_inds(sac_shift+1:end));
end
Pfin_drift_corr = interp1(find(~isnan(fix_ids)),Pfin_drift_corr(~isnan(fix_ids)),1:NT);
Pfin_drift_std = interp1(find(~isnan(fix_ids)),Pfin_drift_std(~isnan(fix_ids)),1:NT);

Pfin_tot_corr = Pfin_fix_corr(:) + Pfin_drift_corr(:);
Pfin_tot_std = sqrt(Pfin_fix_std(:).^2 + Pfin_drift_std(:).^2);

%%
n_trials = length(unique(all_trialvec));
sim_sac_expts = find(~expt_has_ds(cur_block_set));
sim_sac_trials = find(ismember(all_trial_blocknums,sim_sac_expts))';

%% CREATE EXAMPLE TACE
% G086: USING TRIAL 179 for example Fig.

close all

print_on = false;
fig_width = 4.86; %3.27 4.86 6.83
rel_height = 0.6;
use_both_coils = true;
use_ep = true;
yl = [-0.4 0.4];

H = figure();
% for tt = 1:n_trials
% for tt = 121
% for tt = sim_sac_trials
% for tt = [30 31 74 113 121 179 223 239 436 522 556 563 630 658 688 736 746 757 759 913 930 949] %G086
% for tt = [31 113 121 179 239 556 757 949] %G086
% for tt = [52 427 524 756 780 856 893] %M270
% for tt = [42 90 93 150 171 172 193 224 327 331 520 543 727 753 913 959 995 1153] %M275
for tt = [172 331 913 995 1153] %M275
% tt = 574
    fprintf('Trial %d of %d\n',tt,n_trials);
    uu = find(all_trialvec(used_inds) == tt);
    if ~isempty(uu)
        bt = all_t_axis(used_inds(uu(1)));
        et = all_t_axis(used_inds(uu(end)));
        dur = et-bt;
        if dur > 3
            cur_sac_inds = find(ismember(saccade_start_inds,uu));
            rel_sac_start_times = all_t_axis(used_inds(saccade_start_inds(cur_sac_inds))) - bt;
            rel_sac_end_times = all_t_axis(used_inds(saccade_stop_inds(cur_sac_inds))) - bt;            
            
            hold on
            h1=shadedErrorBar(all_t_axis(used_inds(uu))-bt,fin_tot_corr(uu),fin_tot_std(uu),{'color','k'},0);
            if use_ep
                hh1=shadedErrorBar(all_t_axis(used_inds(uu))-bt,Pfin_tot_corr(uu),Pfin_tot_std(uu),{'color','m'},0);
            end
            h2=plot(all_t_axis(used_inds(uu))-bt,corrected_eye_vals_interp(used_inds(uu),2),'r');
            if use_both_coils
                h3=plot(all_t_axis(used_inds(uu))-bt,corrected_eye_vals_interp(used_inds(uu),4),'b');
            end
            %             h2=plot(all_t_axis(used_inds(uu))-bt,smooth_eyepos(uu,2),'k--');
%             h3=plot(all_t_axis(used_inds(uu))-bt,smooth_eyepos(uu,4),'r--');
            line([0 dur],[0 0],'color','k','linestyle','--');
            xlim([0 dur]);
            ylim(yl);
            xlabel('Time (s)');
            ylabel('Orthoganol position (deg)');
            title(sprintf('Trial %d',tt));
            for ii = 1:length(rel_sac_start_times)
                line(rel_sac_start_times([ii ii]),yl,'color',[0.2 0.8 0.2],'linestyle','--');
                %                 line(rel_sac_end_times([ii ii]),yl,'color',[0.2 0.8 0.2],'linestyle','--');
            end
            
            if print_on
                if use_both_coils
                    if use_ep
                        legend([h1.mainLine hh1.mainLine h2 h3],{'No prior','With prior','Left coil','Right coil'})
                    else
                        legend([h1.mainLine h2 h3],{'Inferred','Left coil','Right coil'})
                    end
                else
                    legend([h1.mainLine h2],{'Inferred','Coil'})
                end
                delete(h1.edge);
                if use_ep
                   delete(hh1.edge); 
                end
                figufy(H);
                fname = [fig_dir sprintf('Example_trace_%s_T%d',Expt_name,tt)];
                exportfig(H,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
                close(H);
            else
                pause
                clf(H);
            end
        end
    end
end


%% MICROSACCADE AMP COMPARISON
% measured_seq = smooth_eyepos(:,[2 4]);
measured_seq = corrected_eye_vals_interp(used_inds,[2 4]);

sac_buff_inds = 2;
min_fix_dur = 0.1;

saccade_prefix = fix_ids(saccade_start_inds);
saccade_postfix = fix_ids(saccade_stop_inds);
saccade_prefix_dur = nan(size(saccade_start_inds));
saccade_postfix_dur = nan(size(saccade_start_inds));
saccade_prefix_dur(~isnan(saccade_prefix)) = fix_durs(saccade_prefix(~isnan(saccade_prefix)));
saccade_postfix_dur(~isnan(saccade_postfix)) = fix_durs(saccade_postfix(~isnan(saccade_postfix)));

too_short = find(saccade_prefix_dur  < min_fix_dur | saccade_postfix_dur < min_fix_dur);
long_enough = setdiff(1:length(saccade_start_inds),too_short);

start_pts = saccade_start_inds(long_enough) - sac_buff_inds;
end_pts = saccade_stop_inds(long_enough) + sac_buff_inds;
start_pts(start_pts < 1) = 1; end_pts(end_pts > length(used_inds)) = length(used_inds);
m_pre_pos = measured_seq(start_pts,:);
m_post_pos = measured_seq(end_pts,:);
m_delta_pos = (m_post_pos - m_pre_pos);


inferred_pre_pos = fin_tot_corr(start_pts,:);
inferred_post_pos = fin_tot_corr(end_pts,:);
inferred_delta_pos = (inferred_post_pos - inferred_pre_pos);

saccade_blocks = all_blockvec(used_inds(saccade_start_inds));

use_micros = ismember(long_enough',micro_sacs) & ~isnan(m_delta_pos(:,1)) & ~isnan(inferred_delta_pos);
use_nonmicros = ~ismember(long_enough',micro_sacs) & ~isnan(m_delta_pos(:,1)) & ~isnan(inferred_delta_pos);
% use_micros = ismember(saccade_blocks(long_enough),sim_sac_expts) & ismember(long_enough',micro_sacs) & ~isnan(m_delta_pos(:,1)) & ~isnan(inferred_delta_pos);

%%
% measured_seq = smooth_eyepos(:,[2 4]);
measured_seq = corrected_eye_vals_interp(used_inds,[2 4]);
min_fix_dur = 0.1;
long_fix_dur = 0.15;

inferred_drift = nan(size(fin_tot_corr));
measured_drift = nan(length(fin_tot_corr),2);
inferred_fix_avg = nan(n_fixs,1);
measured_fix_avg = nan(n_fixs,2);
fix_ind_vec = nan(length(fin_tot_corr),1);
inf_drift_slope = nan(n_fixs,1);
meas_drift_slope = nan(n_fixs,1);
fix_inds = [];
pfix_inds = [];
long_fix_inds = [];
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
    cur_inf = fin_tot_corr(cur_inds);
    inferred_fix_avg(ii) = nanmedian(fin_tot_corr(cur_inds));
    inferred_drift(cur_inds) = cur_inf-inferred_fix_avg(ii);
    
    measured_fix_avg(ii,:) = nanmedian(measured_seq(cur_inds,:));
    measured_drift(cur_inds,:) = bsxfun(@minus,measured_seq(cur_inds,:),measured_fix_avg(ii,:));
    
    %     if fix_durs(ii) > min_fix_dur
    %         temp = robustfit(cur_inds*dt,inferred_drift(cur_inds));
    %         inf_drift_slope(ii) = temp(2);
    %
    %         temp = robustfit(cur_inds*dt,measured_drift(cur_inds,1));
    %         meas_drift_slope(ii) = temp(2);
    %     end
    
    fix_inds = [fix_inds cur_inds];
    if length(cur_inds)*dt >= long_fix_dur
        long_fix_inds = [long_fix_inds cur_inds];
    end
    fix_ind_vec(cur_inds) = ii;
    if length(cur_inds) > sac_shift
    cur_inds(1:sac_shift) = [];
    end
    pfix_inds = [pfix_inds cur_inds];
end
fix_dur_vec = nan(length(fin_tot_corr),1);
fix_dur_vec(~isnan(fix_ind_vec)) = fix_durs(fix_ind_vec(~isnan(fix_ind_vec)));

used_fixs = find(fix_durs > long_fix_dur);
used_fix_inds = find(fix_dur_vec >= long_fix_dur);


close all

%% Plot figure comparing drift
coil_drift = measured_drift(used_fix_inds,1);
% coil_drift = mean(measured_drift(used_fix_inds,:),2);

close all
sqrt_scale = true;
n_bins = 50;
xr = [-0.1 0.1];
h1=figure;
if sqrt_scale
    DensityPlot_jmm(inferred_drift(used_fix_inds),coil_drift,'ynormal','xrange',xr,'yrange',xr,'sd',[3 3],'sqrtsc');
else
    [h,det] = DensityPlot_jmm(inferred_drift(used_fix_inds),coil_drift,'ynormal','xrange',xr,'yrange',xr,'sd',[3 3]);
    eps = -0.5;
    lZ = log10(det.z);
    lZ(lZ < eps) = eps;
    imagesc(det.x(1,:),det.y(:,1),lZ);
    set(gca,'ydir','normal');
    caxis([eps 4]);
end

line([-0.5 0.5],[-0.5 0.5],'color','w');
xlabel('Inferred amplitude (deg)');
ylabel('Measured amplitude (deg)');
xlim(xr); ylim(xr);
set(gca,'xtick',-0.2:0.05:0.2,'ytick',-0.2:0.05:0.2);
axis square

b_ax = linspace(xr(1),xr(2),n_bins);

INF_nx = histc(inferred_drift(used_fix_inds),b_ax);
INF_nx = INF_nx/sum(INF_nx);

MEAS_nx = histc(coil_drift,b_ax);
MEAS_nx = MEAS_nx/sum(MEAS_nx);

DIFF_nx = histc(inferred_drift(used_fix_inds) - coil_drift,b_ax);
DIFF_nx = DIFF_nx/sum(DIFF_nx);

h5=figure; hold on
stairs(b_ax,INF_nx,'b');
stairs(b_ax,MEAS_nx,'r');
stairs(b_ax,DIFF_nx,'k');
xlim(xr); 
set(gca,'xtick',-0.2:0.05:0.2,'ytick',[]);
% set(gca,'yscale','log');
% ylim([0.0001 2]);

%%
fig_width = 3.27; 
rel_height = 0.85;

figufy(h1);
fname = [fig_dir 'drift_density_LP.pdf'];
exportfig(h1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(h1);

figufy(h5);
fname = [fig_dir 'drift_dists_LP.pdf'];
exportfig(h5,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(h5);

%% plot figure comparing microsacs
% coil_delta_pos = m_delta_pos(use_micros,1);
coil_delta_pos = mean(m_delta_pos(use_micros,:),2);

sqrt_scale = false;
n_bins = 50;
xr = [-0.4 0.4];

close all

h1=figure;
if sqrt_scale
    [hs,det] = DensityPlot_jmm(inferred_delta_pos(use_micros),coil_delta_pos,'ynormal','xrange',xr,'yrange',xr,'sd',[3 3],'sqrtsc');
else
    DensityPlot_jmm(inferred_delta_pos(use_micros),coil_delta_pos,'ynormal','xrange',xr,'yrange',xr,'sd',[3 3]);
end
% plot(inferred_delta_pos(use_micros),m_delta_pos(use_micros,1),'.','markersize',4);
% xlim([-0.5 0.5]); ylim([-0.5 0.5]);
line([-0.5 0.5],[-0.5 0.5],'color','w');
xlabel('Inferred amplitude (deg)');
ylabel('Measured amplitude (deg)');
xlim(xr); ylim(xr);
set(gca,'xtick',-0.4:0.1:0.4,'ytick',-0.4:0.1:0.4);
axis square


b_ax = linspace(-0.4,0.4,n_bins+1);

INF_nx = histc(inferred_delta_pos(use_micros),b_ax);
INF_nx = INF_nx/sum(INF_nx);
% h2=figure; hold on
% stairs(b_ax,INF_nx,'b');
% yl = ylim();
% line([0 0],yl,'color','k')
% xlim(b_ax([1 end])); 
% title('Inferred position');

MEAS_nx = histc(coil_delta_pos,b_ax);
MEAS_nx = MEAS_nx/sum(MEAS_nx);
% h3=figure; hold on
% stairs(b_ax,MEAS_nx,'r');
% yl = ylim();
% line([0 0],yl,'color','k')
% xlim(b_ax([1 end])); 
% title('Measured position');

DIFF_nx = histc(inferred_delta_pos(use_micros) - coil_delta_pos,b_ax);
DIFF_nx = DIFF_nx/sum(DIFF_nx);
% h4=figure; hold on
% stairs(b_ax,DIFF_nx,'k');
% yl = ylim();
% line([0 0],yl,'color','k')
% xlim(b_ax([1 end])); 
% title('Difference');

h5=figure; hold on
stairs(b_ax,INF_nx,'b');
stairs(b_ax,MEAS_nx,'r');
stairs(b_ax,DIFF_nx,'k');
xlim(b_ax([1 end])); 
xlim(xr); 
set(gca,'xtick',-0.4:0.1:0.4,'ytick',[]);

%% PRINT FILTER BEFORE/AFTER COMPARISONS
% fig_width = 3.27; 
fig_width = 3.27; 
rel_height = 0.85;

figufy(h1);
fname = [fig_dir 'microsac_density3.pdf'];
exportfig(h1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(h1);

% figufy(h2);
% fname = [fig_dir 'microsac_inferred.pdf'];
% exportfig(h2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(h2);
% 
% figufy(h3);
% fname = [fig_dir 'microsac_measured.pdf'];
% exportfig(h3,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(h3);

% figufy(h4);
% fname = [fig_dir 'microsac_diff.pdf'];
% exportfig(h4,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(h4);

figufy(h5);
fname = [fig_dir 'microsac_alldiff2.pdf'];
exportfig(h5,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(h5);

%% plot figure comparing absolute positions
usable_inds = fix_inds;
coil_meas = measured_seq(usable_inds,1);
% coil_meas = mean(measured_seq(usable_inds,:),2);

close all
n_bins = 50;
xr = [-0.4 0.4];
sqrt_scale = false;

h1=figure;
if sqrt_scale
    DensityPlot_jmm(fin_tot_corr(usable_inds),coil_meas,'ynormal','xrange',xr,'yrange',xr,'sd',[3 3],'sqrtsc');
else
    DensityPlot_jmm(fin_tot_corr(usable_inds),coil_meas,'ynormal','xrange',xr,'yrange',xr,'sd',[3 3]);    
end
% plot(fin_tot_corr(usable_inds),measured_seq(usable_inds,1),'.','markersize',1);
line([-0.5 0.5],[-0.5 0.5],'color','w');
xlabel('Inferred position (deg)','fontsize',12);
ylabel('Measured position (deg)','fontsize',12);
set(gca,'fontsize',10,'fontname','arial');
set(gca,'xtick',-0.4:0.1:0.4,'ytick',-0.4:0.1:0.4);
axis square
xlim(xr); ylim(xr);
fillPage(gcf,'papersize',[5 5]);
fname = sprintf('Full_dens_ori%d',bar_ori);

b_ax = linspace(xr(1),xr(2),n_bins+1);
INF_nx = histc(fin_tot_corr(usable_inds)',b_ax);
INF_nx = INF_nx/sum(INF_nx);
% figure; hold on
% stairs(b_ax,INF_nx,'b');
% yl = ylim();
% line([0 0],yl,'color','k')
% xlim(b_ax([1 end])); set(gca,'xtick',[],'ytick',[]);
% box off;
% title('Inferred position','fontsize',10);

h3=figure; hold on
MEAS_nx = histc(coil_meas,b_ax);
MEAS_nx = MEAS_nx/sum(MEAS_nx);
stairs(b_ax,MEAS_nx,'r');
xlim([-0.4 0.4]);
yl = ylim();
line([0 0],yl,'color','k')
xlim(b_ax([1 end])); 
set(gca,'xtick',[],'ytick',[]);
title('Measured position','fontsize',10);

DIFF_nx = histc(fin_tot_corr(usable_inds) - coil_meas,b_ax);
DIFF_nx = DIFF_nx/sum(DIFF_nx);
% % figure; hold on
% stairs(b_ax,DIFF_nx,'k');
% yl = ylim();
% % line([0 0],yl,'color','k')
% % xlim(b_ax([1 end])); set(gca,'xtick',[],'ytick',[]);
% % fillPage(gcf,'papersize',[4 4]);
% % box off;
% % title('Difference','fontsize',10);

h5=figure; hold on
stairs(b_ax,INF_nx,'b');
stairs(b_ax,MEAS_nx,'r');
stairs(b_ax,DIFF_nx,'k');
set(gca,'xtick',-0.4:0.1:0.4,'ytick',[]);
xlim(xr);
%%
figufy(h5);
fname = [fig_dir 'alleyepos_alldiff2.pdf'];
exportfig(h5,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(h5);

figufy(h1);
fname = [fig_dir 'alleyepos_dens2.pdf'];
exportfig(h1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(h1);

%%
figufy(h3);
fname = [fig_dir 'measured_ovdiff.pdf'];
exportfig(h3,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(h3);

%%
n_bins = 100;
unc_ax = linspace(0,0.04,n_bins);
unc_cdist = histc(fin_tot_std(fix_inds),unc_ax);
unc_cdist = unc_cdist/max(unc_cdist);

h1 = figure();
plot(unc_ax,unc_cdist,'k','linewidth',1);

fig_width = 3.27; 
rel_height = 0.85;
figufy(h1);
fname = [fig_dir 'inf_unc_dist.pdf'];
exportfig(h1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(h1);

%% ROBUST STATS ESTIMATOR FIG
xx = linspace(-1.25,1.25,200);
y_inf = ksdensity(fin_tot_corr,xx);
y_meas = ksdensity(measured_seq(:,1),xx);

r_reg = std(fin_tot_corr);
r_rob = robust_std_dev(fin_tot_corr);

figure
% plot(xx,y_inf,xx,y_meas,'r');
plot(xx,y_inf);
hold on
plot(xx,normpdf(xx,0,r_reg),'k');
plot(xx,normpdf(xx,0,r_rob),'r');
% set(gca,'yscale','log')
% ylim([1e-4 10]); 
xlim(xx([1 end]));



%% CREATE EXAMPLE OF 2-STAGE PROCEDUR
% G086: SST: 436 913   AT: (LOTS)

close all

print_on = true;
fig_width = 4.86; %3.27 4.86 6.83
rel_height = 0.6;
use_both_coils = false;
yl = [-0.4 0.4];

H = figure();
% for tt = 1:n_trials
% for tt = sim_sac_trials
for tt = [30 31 74 113 179] %G086
% for tt = [179] %G086
% for tt = [52 427 524 756 780 856 893] %M270
% for tt = [89 228 530 574 591 812 813 962 992 1132] %M275
% tt = 574
    fprintf('Trial %d of %d\n',tt,n_trials);
    uu = find(all_trialvec(used_inds) == tt);
    if ~isempty(uu)
        bt = all_t_axis(used_inds(uu(1)));
        et = all_t_axis(used_inds(uu(end)));
        dur = et-bt;
        if dur > 3
            cur_sac_inds = find(ismember(saccade_start_inds,uu));
            rel_sac_start_times = all_t_axis(used_inds(saccade_start_inds(cur_sac_inds))) - bt;
            rel_sac_end_times = all_t_axis(used_inds(saccade_stop_inds(cur_sac_inds))) - bt;            
            
            hold on
            hh1=shadedErrorBar(all_t_axis(used_inds(uu))-bt,fin_fix_corr(uu),fin_fix_std(uu),{'color','m'},0);
            hh2=shadedErrorBar(all_t_axis(used_inds(uu))-bt,fin_drift_corr(uu),fin_drift_std(uu),{'color',[0.2 0.8 0.2]},0);
            h1=shadedErrorBar(all_t_axis(used_inds(uu))-bt,fin_tot_corr(uu),fin_tot_std(uu),{'color','k'},0);
     
            h2=plot(all_t_axis(used_inds(uu))-bt,corrected_eye_vals_interp(used_inds(uu),2),'r');
            if use_both_coils
                h3=plot(all_t_axis(used_inds(uu))-bt,corrected_eye_vals_interp(used_inds(uu),4),'b');
            end
            
            line([0 dur],[0 0],'color','k','linestyle','--');
            xlim([0 dur]);
            ylim(yl);
            xlabel('Time (s)');
            ylabel('Orthoganol position (deg)');
            title(sprintf('Trial %d',tt));
            for ii = 1:length(rel_sac_start_times)
                line(rel_sac_start_times([ii ii]),yl,'color',[0.2 0.8 0.2],'linestyle','--');
                %                 line(rel_sac_end_times([ii ii]),yl,'color',[0.2 0.8 0.2],'linestyle','--');
            end
            
            if print_on
                delete(h1.edge);
                delete(hh1.edge);
                delete(hh2.edge);
                figufy(H);
                fname = [fig_dir sprintf('Example_2stage_%s_T%d',Expt_name,tt)];
                exportfig(H,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
                close(H);
            else
                pause
                clf(H);
            end
        end
    end
end


%%
clear it_fix* drift*
cd(anal_dir)

data_name = 'monoc_eyecorr_hbar_highres3_tight';

fprintf('Loading %s\n',data_name);
load(data_name,'drift_*','et_tr_set','best_*');

fin_fix_corr = nan(NT,1);
fin_fix_std = nan(NT,1);
fin_fix_corr(~isnan(fix_ids)) = best_fix_cor(fix_ids(~isnan(fix_ids)));
fin_fix_corr = interp1(find(~isnan(fix_ids)),fin_fix_corr(~isnan(fix_ids)),1:NT);
fin_fix_std(~isnan(fix_ids)) = best_fix_std(fix_ids(~isnan(fix_ids)));
fin_fix_std = interp1(find(~isnan(fix_ids)),fin_fix_std(~isnan(fix_ids)),1:NT);

fin_fix_corr = fin_fix_corr*sp_dx;
fin_fix_std = fin_fix_std*sp_dx;

fin_drift_corr = drift_post_mean(end,:)*sp_dx;
fin_drift_std = drift_post_std(end,:)*sp_dx;
fix_inds = [];
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
    fix_inds = [fix_inds cur_inds];
%     if length(cur_inds) > sac_shift
%         fin_drift_corr(cur_inds(1:end-sac_shift+1)) = fin_drift_corr(cur_inds(sac_shift:end));
%         fin_drift_std(cur_inds(1:end-sac_shift+1)) = fin_drift_std(cur_inds(sac_shift:end));
%     end
end
for ii = 1:length(trial_start_inds)
    cur_inds = trial_start_inds(ii):trial_end_inds(ii);
    fin_drift_corr(cur_inds(1:end-sac_shift)) = fin_drift_corr(cur_inds(sac_shift+1:end));
    fin_drift_std(cur_inds(1:end-sac_shift)) = fin_drift_std(cur_inds(sac_shift+1:end));
end
fin_drift_corr = interp1(find(~isnan(fix_ids)),fin_drift_corr(~isnan(fix_ids)),1:NT);
fin_drift_std = interp1(find(~isnan(fix_ids)),fin_drift_std(~isnan(fix_ids)),1:NT);

tight_fin_tot_corr = fin_fix_corr(:) + fin_drift_corr(:);
tight_fin_tot_std = sqrt(fin_fix_std(:).^2 + fin_drift_std(:).^2);



data_name = 'monoc_eyecorr_hbar_highres3_loose';

fprintf('Loading %s\n',data_name);
load(data_name,'drift_*','et_tr_set','best_*');

fin_fix_corr = nan(NT,1);
fin_fix_std = nan(NT,1);
fin_fix_corr(~isnan(fix_ids)) = best_fix_cor(fix_ids(~isnan(fix_ids)));
fin_fix_corr = interp1(find(~isnan(fix_ids)),fin_fix_corr(~isnan(fix_ids)),1:NT);
fin_fix_std(~isnan(fix_ids)) = best_fix_std(fix_ids(~isnan(fix_ids)));
fin_fix_std = interp1(find(~isnan(fix_ids)),fin_fix_std(~isnan(fix_ids)),1:NT);

fin_fix_corr = fin_fix_corr*sp_dx;
fin_fix_std = fin_fix_std*sp_dx;

fin_drift_corr = drift_post_mean(end,:)*sp_dx;
fin_drift_std = drift_post_std(end,:)*sp_dx;
fix_inds = [];
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
    fix_inds = [fix_inds cur_inds];
%     if length(cur_inds) > sac_shift
%         fin_drift_corr(cur_inds(1:end-sac_shift+1)) = fin_drift_corr(cur_inds(sac_shift:end));
%         fin_drift_std(cur_inds(1:end-sac_shift+1)) = fin_drift_std(cur_inds(sac_shift:end));
%     end
end
for ii = 1:length(trial_start_inds)
    cur_inds = trial_start_inds(ii):trial_end_inds(ii);
    fin_drift_corr(cur_inds(1:end-sac_shift)) = fin_drift_corr(cur_inds(sac_shift+1:end));
    fin_drift_std(cur_inds(1:end-sac_shift)) = fin_drift_std(cur_inds(sac_shift+1:end));
end
fin_drift_corr = interp1(find(~isnan(fix_ids)),fin_drift_corr(~isnan(fix_ids)),1:NT);
fin_drift_std = interp1(find(~isnan(fix_ids)),fin_drift_std(~isnan(fix_ids)),1:NT);

loose_fin_tot_corr = fin_fix_corr(:) + fin_drift_corr(:);
loose_fin_tot_std = sqrt(fin_fix_std(:).^2 + fin_drift_std(:).^2);

%% CREATE EXAMPLE TACE
% G086: SST: 436 913   AT: (LOTS)

close all

print_on = true;
fig_width = 4.86; %3.27 4.86 6.83
rel_height = 0.6;
yl = [-0.5 0.5];
offset = 0.1;

H = figure();
% for tt = [30 31 74 113 121 179 223 239 436 522 556 563 630 658 688 736 746 757 759 913 930 949] %G086
for tt = [179 688] %G086
    fprintf('Trial %d of %d\n',tt,n_trials);
    uu = find(all_trialvec(used_inds) == tt);
    if ~isempty(uu)
        bt = all_t_axis(used_inds(uu(1)));
        et = all_t_axis(used_inds(uu(end)));
        dur = et-bt;
        if dur > 3
            cur_sac_inds = find(ismember(saccade_start_inds,uu));
            rel_sac_start_times = all_t_axis(used_inds(saccade_start_inds(cur_sac_inds))) - bt;
            rel_sac_end_times = all_t_axis(used_inds(saccade_stop_inds(cur_sac_inds))) - bt;            
            
            hold on
            h1=shadedErrorBar(all_t_axis(used_inds(uu))-bt,fin_tot_corr(uu),fin_tot_std(uu),{'color','k'},0);
            h2=shadedErrorBar(all_t_axis(used_inds(uu))-bt,loose_fin_tot_corr(uu)-offset,loose_fin_tot_std(uu),{'color','m'},0);
            h3=shadedErrorBar(all_t_axis(used_inds(uu))-bt,tight_fin_tot_corr(uu)+offset,tight_fin_tot_std(uu),{'color','g'},0);
            plot(all_t_axis(used_inds(uu))-bt,corrected_eye_vals_interp(used_inds(uu),2),'r');
            line([0 dur],[0 0],'color','k','linestyle','--');
            xlim([0 dur]);
            ylim(yl);
            xlabel('Time (s)');
            ylabel('Orthoganol position (deg)');
            title(sprintf('Trial %d',tt));
            for ii = 1:length(rel_sac_start_times)
                line(rel_sac_start_times([ii ii]),yl,'color',[0.2 0.8 0.2],'linestyle','--');
                %                 line(rel_sac_end_times([ii ii]),yl,'color',[0.2 0.8 0.2],'linestyle','--');
            end
            
            if print_on
                delete(h1.edge);
                delete(h2.edge);
                 delete(h3.edge);
               figufy(H);
                fname = [fig_dir sprintf('LooseTight_Ex_trace_T%d',Expt_name,tt)];
                exportfig(H,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
                close(H);
            else
                pause
                clf(H);
            end
        end
    end
end

%%
ep = 0.04;
unc_ax = linspace(0,ep,100);
unc_dist = hist(fin_tot_std(fin_tot_std < ep),unc_ax);
unc_dist = unc_dist/sum(unc_dist);


unc_dist2 = hist(tight_fin_tot_std(tight_fin_tot_std < ep),unc_ax);
unc_dist2 = unc_dist2/sum(unc_dist2);

unc_dist3 = hist(loose_fin_tot_std(loose_fin_tot_std < ep),unc_ax);
unc_dist3 = unc_dist3/sum(unc_dist3);

figure
plot(unc_ax,unc_dist,'k',unc_ax,unc_dist2,'r',unc_ax,unc_dist3,'b')


