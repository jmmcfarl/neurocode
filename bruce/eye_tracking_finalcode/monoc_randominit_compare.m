clear all
close all
addpath('~/James_scripts/bruce/eye_tracking/');
dir_prefix = '/home/james';
% dir_prefix = '/Volumes/james';

fig_dir = [dir_prefix '/Analysis/bruce/ET_final/'];
single_bar_expts = [232 235 239];

%%
Expt_name ='G086';
Expt_num = str2num(Expt_name(2:end));

bar_ori = 0;

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

if ~ismember(Expt_num,single_bar_expts)
    sp_dx = 0.0565/spatial_usfac/scale_fac;
else
    sp_dx = 0.125/4;
end

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
for bb = 1:length(cur_block_set)
%    binds = find(all_blockvec == bb);
   binds = used_inds(all_blockvec(used_inds) == bb);
   corrected_eye_vals_interp(binds,:) = bsxfun(@minus,corrected_eye_vals_interp(binds,:),nanmedian(corrected_eye_vals_interp(binds,:)));
end

%%
clear it_fix* drift*
cd(anal_dir)

if bar_ori == 0
    old_data_name = 'monoc_eyecorr_hbar2';
    data_name = 'monoc_eyecorr_hbar_highres3';
else
    old_data_name = 'monoc_eyecorr_vbar2';
    data_name = 'monoc_eyecorr_vbar_highres3';
end

fprintf('Loading %s\n',data_name);
load(data_name,'drift_*','et_tr_set','best_*');
load(old_data_name,'it_fix_post_mean','it_fix_post_std');

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
for ii = 1:length(trial_start_inds)
    cur_inds = trial_start_inds(ii):trial_end_inds(ii);
    fin_drift_corr(cur_inds(1:end-sac_shift)) = fin_drift_corr(cur_inds(sac_shift+1:end));
    fin_drift_std(cur_inds(1:end-sac_shift)) = fin_drift_std(cur_inds(sac_shift+1:end));
end
fin_drift_corr = interp1(find(~isnan(fix_ids)),fin_drift_corr(~isnan(fix_ids)),1:NT);
fin_drift_std = interp1(find(~isnan(fix_ids)),fin_drift_std(~isnan(fix_ids)),1:NT);

fin_tot_corr = fin_fix_corr(:) + fin_drift_corr(:);
fin_tot_std = sqrt(fin_fix_std(:).^2 + fin_drift_std(:).^2);


%%
clear it_fix* drift*
cd(anal_dir)

if bar_ori == 0
    old_data_name = 'monoc_eyecorr_hbar_randinit2';
    data_name = 'monoc_eyecorr_hbar_randinit_highres3';
else
    old_data_name = 'monoc_eyecorr_vbar_randinit2';
    data_name = 'monoc_eyecorr_vbar_randinit_highres3';
end

fprintf('Loading %s\n',data_name);
load(data_name,'drift_*','et_tr_set','best_*');
load(old_data_name,'it_fix_post_mean','it_fix_post_std');

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
for ii = 1:length(trial_start_inds)
    cur_inds = trial_start_inds(ii):trial_end_inds(ii);
    fin_drift_corr(cur_inds(1:end-sac_shift)) = fin_drift_corr(cur_inds(sac_shift+1:end));
    fin_drift_std(cur_inds(1:end-sac_shift)) = fin_drift_std(cur_inds(sac_shift+1:end));
end
fin_drift_corr = interp1(find(~isnan(fix_ids)),fin_drift_corr(~isnan(fix_ids)),1:NT);
fin_drift_std = interp1(find(~isnan(fix_ids)),fin_drift_std(~isnan(fix_ids)),1:NT);

rfin_tot_corr = fin_fix_corr(:) + fin_drift_corr(:);
rfin_tot_std = sqrt(fin_fix_std(:).^2 + fin_drift_std(:).^2);

%%
n_trials = length(unique(all_trialvec));
sim_sac_expts = find(~expt_has_ds(cur_block_set));
sim_sac_trials = find(ismember(all_trial_blocknums,sim_sac_expts))';
%%

%%
print_on = true;
yl = [-0.4 0.4];
fig_width = 4.86; %3.27 4.86 6.83
rel_height = 0.6;

close all
n_trials = length(unique(all_trialvec));
% for tt = 1:n_trials
H = figure; 
% for tt = [30 31 74 113 121 179 223 239 436 522 556 563 630 658 688 736 746 757 759 913 930 949] %G086
% tt = 74;
tt = 179;
for tt = [30] %G086
% for tt = [96 137 154 179 376 409]
    uu = find(all_trialvec(used_inds) == tt);
    if ~isempty(uu)
        bt = all_t_axis(used_inds(uu(1)));
        et = all_t_axis(used_inds(uu(end)));
        dur = et-bt;
        if dur > 3.5
            cur_sac_inds = find(ismember(saccade_start_inds,uu));
            rel_sac_start_times = all_t_axis(used_inds(saccade_start_inds(cur_sac_inds))) - bt;
            rel_sac_end_times = all_t_axis(used_inds(saccade_stop_inds(cur_sac_inds))) - bt;            

            hold on;
            h1=shadedErrorBar(all_t_axis(used_inds(uu))-bt,fin_tot_corr(uu),fin_tot_std(uu),{'color','m'});
            h2=shadedErrorBar(all_t_axis(used_inds(uu))-bt,rfin_tot_corr(uu),rfin_tot_std(uu),{'color','k'});
            h3=plot(all_t_axis(used_inds(uu))-bt,corrected_eye_vals_interp(used_inds(uu),2),'r','linewidth',1);
%             h4=plot(all_t_axis(used_inds(uu))-bt,init_eyepos(uu),'g','linewidth',1)
            xlim([0 dur]);
            ylim(yl);
            xlabel('Time (s)','fontsize',10);
            ylabel('Orthoganol position (deg)','fontsize',10);
            title(sprintf('Trial %d',tt));
            for ii = 1:length(rel_sac_start_times)
                line(rel_sac_start_times([ii ii]),yl,'color',[0.2 0.8 0.2],'linestyle','--');
            end
            
            if print_on
                delete(h1.edge);
                delete(h2.edge);
                figufy(H);
                fname = [fig_dir sprintf('Example_randinit_%s_T%d',Expt_name,tt)];
                exportfig(H,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
                close(H);
            else
                pause
                clf(H);
            end
        end
    end
end

