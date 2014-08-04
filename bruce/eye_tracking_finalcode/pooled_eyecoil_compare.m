clear all
close all

use_simsac_only = false;

addpath(genpath('~/James_scripts/chronux/spectral_analysis/'));
% addpath('~/James_scripts/bruce/eye_tracking/');
fin_dir = '~/Analysis/bruce/ET_final/';
% fin_dir = '/Volumes/james/Analysis/bruce/ET_final/';

Expt_list = {'G085','G086','G087','G088','G089','G091','G093','G095'};
ori_list = [zeros(size(Expt_list))];

%add vertical oris
ori_list = [ori_list 90*ones(size(Expt_list))];
Expt_list = repmat(Expt_list,1,2);

%there's no vertical bar for G095
Expt_list(end) = [];
ori_list(end) = [];

Expt_list = cat(2,Expt_list,{'M266','M270','M275','M277'});
ori_list = cat(2,ori_list,[80 60 135 70]);

% Expt_list = cat(2,Expt_list,{'M232','M235','M239'});
% ori_list = cat(2,ori_list,[50 30 130]);

single_bar_expts = [232 235 239];
%%
use_sac_kerns = 1;
use_coils = [0 0]; %[L R]

flen = 12;
use_nPix = 16;
min_trial_dur = 0.75;
% spatial_usfac = 2;
spatial_usfac = 4;

for EE = 1:length(Expt_list)
    
    clear expt_* included_type
    
    Expt_num = str2num(Expt_list{EE}(2:end));
    Expt_name = Expt_list{EE};
    bar_ori = ori_list(EE);
    
    fprintf('Analyzing expt %s, ori %d, %d/%d\n',Expt_name,ori_list(EE),EE,length(Expt_list));
    
    data_dir = ['~/Data/bruce/' Expt_name];
%     data_dir = ['/Volumes/james/Data/bruce/' Expt_name];
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
    
    anal_dir = ['~/Analysis/bruce/' Expt_name '/ET_final/'];
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
    
    if ~any(strcmpi(Expt_name,{'M232','M235','M239'}))
        sim_sac_expts = find(~expt_has_ds(cur_block_set));
        imback_gs_expts = find(expt_has_ds(cur_block_set) & expt_imback(cur_block_set)');
        grayback_gs_expts = find(expt_has_ds(cur_block_set) & ~expt_imback(cur_block_set)');
    else
        sim_sac_expts = [];
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
    
    %%
    if use_simsac_only
        ss_set = find(ismember(all_blockvec(used_inds),sim_sac_expts));
    else
        ss_set = 1:length(used_inds);
    end
    new_NT = length(ss_set);
    if ~isempty(ss_set)
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
        
        %% SUBTRACT BLOCK-WISE MEDIAN EYE POSITION
        for bb = 1:length(cur_block_set)
            binds = used_inds(all_blockvec(used_inds) == bb);
            corrected_eye_vals_interp(binds,:) = bsxfun(@minus,corrected_eye_vals_interp(binds,:),nanmedian(corrected_eye_vals_interp(binds,:)));
        end
        
        %%
        proc_eyepos = corrected_eye_vals_interp;
                
        %subtract off within-trial median
        for ii = 1:length(trial_start_inds)
            cur_inds = trial_start_inds(ii):trial_end_inds(ii);
            proc_eyepos(used_inds(cur_inds),:) = bsxfun(@minus,proc_eyepos(used_inds(cur_inds),:),median(proc_eyepos(used_inds(cur_inds),:)));
        end
        
        proc_eyepos = proc_eyepos(used_inds(ss_set),:);
        raw_eyepos = corrected_eye_vals_interp(used_inds(ss_set),:);
        all_proc_eyepos{EE} = proc_eyepos(:,[2 4]);
        all_raw_eyepos{EE} = raw_eyepos(:,[2 4]);
        
        %%
        clear it_fix* drift*
        cd(anal_dir)
        
        if Expt_name(1) == 'G'
            if bar_ori == 0
                old_data_name ='./monoc_eyecorr_hbar2.mat';
                data_name ='./monoc_eyecorr_hbar_highres3.mat';
                mod_data_name = './monoc_eyecorr_hbar_mods.mat';
            elseif bar_ori == 90
                old_data_name ='./monoc_eyecorr_vbar2.mat';
                data_name ='./monoc_eyecorr_vbar_highres3.mat';
                mod_data_name = './monoc_eyecorr_vbar_mods.mat';
            end
        else
            old_data_name ='./monoc_eyecorr2_Cprior.mat';
            data_name ='./monoc_eyecorr_highres4.mat';
            mod_data_name = './monoc_eyecorr_mods.mat';
        end
        fprintf('Loading %s\n',data_name);
        load(old_data_name,'it_R2*','dit_R2*');
        load(data_name,'drift_*','best_fix*','et_tr_set');
        fprintf('Loading %s\n',mod_data_name);
        load(mod_data_name,'all_mod_SU*');
        
        fin_fix_corr = nan(NT,1);
        fin_fix_std = nan(NT,1);
%         fin_fix_corr(~isnan(fix_ids)) = it_fix_post_mean(end,fix_ids(~isnan(fix_ids)));
        fin_fix_corr(~isnan(fix_ids)) = best_fix_cor(fix_ids(~isnan(fix_ids)));
        fin_fix_corr = interp1(find(~isnan(fix_ids)),fin_fix_corr(~isnan(fix_ids)),1:NT);
%         fin_fix_std(~isnan(fix_ids)) = it_fix_post_std(end,fix_ids(~isnan(fix_ids)));
        fin_fix_std(~isnan(fix_ids)) = best_fix_std(fix_ids(~isnan(fix_ids)));
        fin_fix_std = interp1(find(~isnan(fix_ids)),fin_fix_std(~isnan(fix_ids)),1:NT);
        
        fin_fix_corr = fin_fix_corr*sp_dx;
        fin_fix_std = fin_fix_std*sp_dx;
        
        fin_drift_corr = drift_post_mean(end,:)*sp_dx;
        fin_drift_std = drift_post_std(end,:)*sp_dx;
        
%         for ii = 1:n_fixs
%             cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
%             if length(cur_inds) > sac_shift
%                 fin_drift_corr(cur_inds(1:end-sac_shift+1)) = fin_drift_corr(cur_inds(sac_shift:end));
%                 fin_drift_std(cur_inds(1:end-sac_shift+1)) = fin_drift_std(cur_inds(sac_shift:end));
%             end
%         end
for ii = 1:length(trial_start_inds)
    cur_inds = trial_start_inds(ii):trial_end_inds(ii);
    fin_drift_corr(cur_inds(1:end-sac_shift)) = fin_drift_corr(cur_inds(sac_shift+1:end));
    fin_drift_std(cur_inds(1:end-sac_shift)) = fin_drift_std(cur_inds(sac_shift+1:end));
end
fin_drift_corr = interp1(find(~isnan(fix_ids)),fin_drift_corr(~isnan(fix_ids)),1:NT);
 fin_drift_std = interp1(find(~isnan(fix_ids)),fin_drift_std(~isnan(fix_ids)),1:NT);
       
        fin_tot_corr = fin_fix_corr(:) + fin_drift_corr(:);
        fin_tot_std = sqrt(fin_fix_std(:).^2 + fin_drift_std(:).^2);
        
        fin_tot_corr = fin_tot_corr(ss_set);
        fin_tot_std = fin_tot_std(ss_set);
        
        all_inf_eyepos{EE} = fin_tot_corr;
        all_inf_eyestd{EE} = fin_tot_std;
        
        %%
        
        eye_data(EE).init_R2 = it_R2(1,:);
        eye_data(EE).fin_R2 = dit_R2(end,:);
        eye_data(EE).tr_set = et_tr_set;
        
        su_inds = find(all_mod_SU(et_tr_set) > 0);
        clear drift_R2_LOO fix_R2_L00
        [~,n_drift_its,~] = size(dit_R2_LOO);
        n_sus = length(su_inds);
        for ii = 1:n_sus
            if Expt_name(1) == 'G'
                cur_ind = ii;
            else
                cur_ind = su_inds(ii);
            end
            drift_R2_LOO(ii) = squeeze(dit_R2_LOO(cur_ind,end,et_tr_set(su_inds(ii))));
            fix_R2_L00(ii) = squeeze(it_R2(1,et_tr_set(su_inds(ii))));
        end
        eye_data(EE).init_R2_SU = fix_R2_L00;
        eye_data(EE).fin_R2_SU_LOO = drift_R2_LOO;
        eye_data(EE).fin_R2_SU = dit_R2(end,et_tr_set(su_inds));
        
        %%
%         measured_seq = smooth_eyepos(:,[2 4]);
        measured_seq = raw_eyepos(:,[2 4]);
%         measured_seq = mean(raw_eyepos(:,[2 4]),2);
        
        min_fix_dur = 0.1;
        long_fix_dur = 0.15;
        
        inferred_drift = nan(size(fin_tot_corr));
        measured_drift = nan(length(fin_tot_corr),2);
        inferred_fix_avg = nan(n_fixs,1);
        measured_fix_avg = nan(n_fixs,2);
        fix_ind_vec = nan(length(fin_tot_corr),1);
        fix_inds = [];
        for ii = 1:n_fixs
            cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
            if ismember(cur_inds(1),ss_set)
                cur_inds = find(ismember(ss_set,cur_inds));
                cur_inf = fin_tot_corr(cur_inds);
                inferred_fix_avg(ii) = nanmedian(fin_tot_corr(cur_inds));
                inferred_drift(cur_inds) = cur_inf-inferred_fix_avg(ii);
                
                measured_fix_avg(ii,:) = nanmedian(measured_seq(cur_inds,:));
                measured_drift(cur_inds,:) = bsxfun(@minus,measured_seq(cur_inds,:),measured_fix_avg(ii,:));
                
                fix_ind_vec(cur_inds) = ii;
                fix_inds = [fix_inds; cur_inds(:)];
            end
        end
        fix_dur_vec = nan(length(fin_tot_corr),1);
        fix_dur_vec(~isnan(fix_ind_vec)) = fix_durs(fix_ind_vec(~isnan(fix_ind_vec)));
        
        usable_data = ~isnan(measured_seq(:,1)) & ~isnan(fin_tot_corr);
        
        udata = ~isnan(fix_ind_vec) & usable_data;
        [eye_data(EE).drift_corrs,eye_data(EE).drif_pvals] = corr(measured_drift(udata,:),inferred_drift(udata),'type','spearman');
        [eye_data(EE).drift_corrs_LR,eye_data(EE).drif_pvals_LR] = corr(measured_drift(udata,1),measured_drift(udata,2),'type','spearman');
        [eye_data(EE).tot_corrs,eye_data(EE).tot_pvals] = corr(measured_seq(udata,:),fin_tot_corr(udata),'type','spearman');
        [eye_data(EE).tot_corrs_LR,eye_data(EE).tot_pvals_LR] = corr(measured_seq(udata,1),measured_seq(udata,1),'type','spearman');
        
        udata = fix_dur_vec >= min_fix_dur & usable_data;
        [eye_data(EE).drift_corrs_mfd,eye_data(EE).drif_pvals_mfd] = corr(measured_drift(udata,:),inferred_drift(udata),'type','spearman');
        [eye_data(EE).drift_corrs_mfd_LR,eye_data(EE).drif_pvals_mfd_LR] = corr(measured_drift(udata,1),measured_drift(udata,2),'type','spearman');
        [eye_data(EE).tot_corrs_mfd,eye_data(EE).tot_pvals_mfd] = corr(measured_seq(udata,:),fin_tot_corr(udata),'type','spearman');
        [eye_data(EE).tot_corrs_mfd_LR,eye_data(EE).tot_pvals_mfd_LR] = corr(measured_seq(udata,1),measured_seq(udata,2),'type','spearman');
        
        long_udata = fix_dur_vec >= long_fix_dur & usable_data;
        [eye_data(EE).drift_corrs_mfdl,eye_data(EE).drif_pvals_mfdl] = corr(measured_drift(long_udata,:),inferred_drift(long_udata),'type','spearman');
        [eye_data(EE).tot_corrs_mfdl,eye_data(EE).tot_pvals_mfdl] = corr(measured_seq(long_udata,:),fin_tot_corr(long_udata),'type','spearman');
        
        udata = ~isnan(measured_fix_avg(:,1)) & ~isnan(inferred_fix_avg);
        [eye_data(EE).fix_corrs,eye_data(EE).fix_pvals] = corr(measured_fix_avg(udata,:),inferred_fix_avg(udata),'type','spearman');
        udata = fix_durs >= min_fix_dur & ~isnan(measured_fix_avg(:,1)) & ~isnan(inferred_fix_avg);
        [eye_data(EE).fix_corrs_mfd,eye_data(EE).fix_pvals_mfd] = corr(measured_fix_avg(udata,:),inferred_fix_avg(udata),'type','spearman');
        udata = fix_durs >= long_fix_dur & ~isnan(measured_fix_avg(:,1)) & ~isnan(inferred_fix_avg);
        [eye_data(EE).fix_corrs_mfdl,eye_data(EE).fix_pvals_mfdl] = corr(measured_fix_avg(udata,:),inferred_fix_avg(udata),'type','spearman');
        
        eye_data(EE).inf_tot_std = robust_std_dev(fin_tot_corr(fix_inds));
        eye_data(EE).meas_tot_std = robust_std_dev(measured_seq(fix_inds,:));
        eye_data(EE).meas_tot_std_proc = robust_std_dev(proc_eyepos(fix_inds,:));
        eye_data(EE).tot_err = robust_std_dev(bsxfun(@minus,measured_seq(fix_inds,:),fin_tot_corr(fix_inds)));
        eye_data(EE).tot_err_proc = robust_std_dev(bsxfun(@minus,proc_eyepos(fix_inds,:),fin_tot_corr(fix_inds)));
        eye_data(EE).inf_drift_std = robust_std_dev(inferred_drift(long_udata));
        eye_data(EE).inf_drift_var = nanvar(inferred_drift(fix_inds));
        eye_data(EE).inf_tot_var = nanvar(fin_tot_corr(fix_inds));
        eye_data(EE).meas_drift_std = robust_std_dev(measured_drift(long_udata,:));
        eye_data(EE).drift_err = robust_std_dev(bsxfun(@minus,measured_drift(long_udata,:),inferred_drift(long_udata)));
        
        eye_data(EE).inf_median = nanmedian(fin_tot_corr(fix_inds));
        eye_data(EE).meas_median = nanmedian(measured_seq(fix_inds,:));
        eye_data(EE).med_unc = nanmedian(fin_tot_std(fix_inds));
        
        eye_data(EE).avg_fix_dur = nanmean(fix_durs);
        %%
        sac_buff_inds = 2;
        
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
       
        bad = find(diff(start_pts) == 0 | diff(end_pts) == 0);
        end_pts(bad) = []; start_pts(bad) = [];
        
        uset = find(ismember(start_pts,ss_set) & ismember(end_pts,ss_set));
        long_enough = long_enough(uset);
        
       
        start_pts = find(ismember(ss_set,start_pts(uset)));
        end_pts = find(ismember(ss_set,end_pts(uset)));
        m_pre_pos = measured_seq(start_pts,:);
        m_post_pos = measured_seq(end_pts,:);
        m_delta_pos = (m_post_pos - m_pre_pos);
        
        inferred_pre_pos = fin_tot_corr(start_pts,:);
        inferred_post_pos = fin_tot_corr(end_pts,:);
        inferred_delta_pos = (inferred_post_pos - inferred_pre_pos);
        
        use_micros = ismember(long_enough',micro_sacs) & ~isnan(m_delta_pos(:,1)) & ~isnan(inferred_delta_pos);
        use_nonmicros = ~ismember(long_enough',micro_sacs) & ~isnan(m_delta_pos(:,1)) & ~isnan(inferred_delta_pos);
        
        [eye_data(EE).sac_corrs,eye_data(EE).sac_pvals] = corr(m_delta_pos(use_nonmicros,:),inferred_delta_pos(use_nonmicros),'type','spearman');
        [eye_data(EE).msac_corrs,eye_data(EE).msac_pvals] = corr([m_delta_pos(use_micros,:)],inferred_delta_pos(use_micros),'type','spearman');
        [eye_data(EE).sac_corrs_LR,eye_data(EE).sac_pvals_LR] = corr(m_delta_pos(use_nonmicros,1),m_delta_pos(use_nonmicros,2),'type','spearman');
        [eye_data(EE).msac_corrs_LR,eye_data(EE).msac_pvals_LR] = corr(m_delta_pos(use_micros,1),m_delta_pos(use_micros,2),'type','spearman');
        
        eye_data(EE).inf_sac_median = nanmedian(abs(inferred_delta_pos(use_nonmicros)));
        eye_data(EE).meas_sac_median = nanmedian(abs(m_delta_pos(use_nonmicros,:)));
        eye_data(EE).inf_msac_median = nanmedian(abs(inferred_delta_pos(use_micros)));
        eye_data(EE).meas_msac_median = nanmedian(abs(m_delta_pos(use_micros,:)));
        
        eye_data(EE).inferred_delta_pos = inferred_delta_pos;
        eye_data(EE).measured_delta_pos = m_delta_pos;
        eye_data(EE).use_micros = use_micros;
        eye_data(EE).use_nonmicros = use_nonmicros;
        
        %%
        eye_data(EE).data_dur = length(used_inds);
        eye_data(EE).n_tr_units = length(et_tr_set);
        if Expt_name(1) == 'G'
            eye_data(EE).n_tr_SUs = sum(et_tr_set > 96);
        else
            eye_data(EE).n_tr_SUs = sum(et_tr_set > 24);
        end
        
        %%
        n_bins = 50;
        xr = [-0.1 0.1];
        if Expt_name(1) == 'G'
            coil_drift = measured_drift(:,1); %left coil drift
            coil_msac = m_delta_pos(:,1); %left coil sacs
            coil_abs = measured_seq(:,1); %left coil meas
        else
            coil_drift = mean(measured_drift,2); %avg of left and right coils
            coil_msac = mean(m_delta_pos,2); %avg of left and right coil msacs
            coil_abs = mean(measured_seq,2); %avg of left and right coil meas
        end
        [h,eye_data(EE).drift_jDist] = DensityPlot_jmm(inferred_drift(long_udata),coil_drift(long_udata),'ynormal','xrange',xr,'yrange',xr,'sd',[3 3],'noplot');
        dax = linspace(xr(1),xr(2),n_bins+1);
        inf_drift_marg = histc(inferred_drift(long_udata),dax);
        inf_drift_marg = inf_drift_marg/sum(inf_drift_marg);
        meas_drift_marg = histc(coil_drift(long_udata),dax);
        meas_drift_marg = meas_drift_marg/sum(meas_drift_marg);
        eye_data(EE).inf_drift_marg = inf_drift_marg;
        eye_data(EE).meas_drift_marg = meas_drift_marg;
        eye_data(EE).drift_ax = dax;
        
        n_bins = 50;
        xr = [-0.4 0.4];
        [hs,eye_data(EE).msac_jDist] = DensityPlot_jmm(inferred_delta_pos(use_micros),coil_msac(use_micros),'ynormal','xrange',xr,'yrange',xr,'sd',[3 3],'noplot');
        max = linspace(xr(1),xr(2),n_bins + 1);
        inf_msac_marg = histc(inferred_delta_pos(use_micros),max);
        inf_msac_marg = inf_msac_marg/sum(inf_msac_marg);
        meas_msac_marg = histc(coil_msac(use_micros),max);
        meas_msac_marg = meas_msac_marg/sum(meas_msac_marg);
        eye_data(EE).inf_msac_marg = inf_msac_marg;
        eye_data(EE).meas_drift_marg = meas_drift_marg;
        eye_data(EE).micro_ax = max;
        
        n_bins = 50;
        xr = [-0.4 0.4];
        [hs,eye_data(EE).tot_jDist] = DensityPlot_jmm(fin_tot_corr(fix_inds),coil_abs(fix_inds),'ynormal','xrange',xr,'yrange',xr,'sd',[3 3],'noplot');
        tax = linspace(xr(1),xr(2),n_bins + 1);
        inf_tot_marg = histc(fin_tot_corr(fix_inds),tax);
        inf_tot_marg = inf_tot_marg/sum(inf_tot_marg);
        meas_tot_marg = histc(coil_abs(fix_inds),tax);
        meas_tot_marg = meas_tot_marg/sum(meas_tot_marg);
        eye_data(EE).inf_tot_marg = inf_tot_marg;
        eye_data(EE).meas_tot_marg = meas_tot_marg;
        eye_data(EE).tot_ax = max;
        
        %%
%         params.Fs = 1/dt;
%         params.tapers = [5 9];
%         params.trialave = 1;
%         win = 100;
%         
%         [C,phi,S12,S1,S2,f]=coherencysegc(measured_seq(:,1),fin_tot_corr(:),win,params);        
%         
%         eye_data(EE).f = f;
%         eye_data(EE).meas_spec(:,1) = log10(S1);
%         eye_data(EE).inf_spec = log10(S2);
%         eye_data(EE).coh(:,1) = C;
%                 
%         [C,phi,S12,S1,S2,f]=coherencysegc(measured_seq(:,2),fin_tot_corr(:),win,params);        
%         eye_data(EE).meas_spec(:,2) = log10(S1);
%         eye_data(EE).coh(:,2) = C;
    end
end

%%
cd(fin_dir)
dname = 'all_pooled_eye_data_nosmooth_bsub_hres3';
% dname = 'all_pooled_eye_data_nosmooth_bsub_hres4';
if use_simsac_only
    dname = [dname '_simsac'];
end
save(dname,'eye_data','Expt_list','ori_list','all_inf_eyepos','all_inf_eyestd','all_proc_eyepos','all_raw_eyepos');

%%
fig_dir = '/home/james/Analysis/bruce/ET_final/';

udata = find(cellfun(@(x)length(x),all_inf_eyepos) > 0);

tot_err = reshape([eye_data(udata).tot_err],[2 length(udata)]);
tot_err_proc = reshape([eye_data(udata).tot_err_proc],[4 length(udata)]);
tot_corrs_mfd = reshape([eye_data(udata).tot_corrs_mfd],[2 length(udata)]);
drift_corrs_mfdl = reshape([eye_data(udata).drift_corrs_mfdl],[2 length(udata)]);
inf_tot_std = [eye_data(udata).inf_tot_std];
inf_tot_var = [eye_data(udata).inf_tot_var];
inf_drift_std = [eye_data(udata).inf_drift_std];
inf_drift_var = [eye_data(udata).inf_drift_var];
meas_tot_std = reshape([eye_data(udata).meas_tot_std],[2 length(udata)]);
meas_tot_std_proc = reshape([eye_data(udata).meas_tot_std_proc],[4 length(udata)]);
meas_drift_std = reshape([eye_data(udata).meas_drift_std],[2 length(udata)]);
msac_corrs = reshape([eye_data(udata).msac_corrs],[2 length(udata)]);
med_unc = [eye_data(udata).med_unc];
jbeExpts = find(strncmpi(Expt_list(udata),'G',1));
lemExpts = find(strncmpi(Expt_list(udata),'M',1));

% init_R2 = {eye_data(udata).init_R2};
% all_tr_set = {eye_data(udata).tr_set};
% fin_R2 = {eye_data(udata).fin_R2};
% init_R2_SU = {eye_data(udata).init_R2_SU};
% fin_R2_SU = {eye_data(udata).fin_R2_SU};
% fin_R2_SU_LOO = {eye_data(udata).fin_R2_SU_LOO};
% med_init_R2 = cellfun(@(X)nanmedian(X),init_R2);
% med_fin_R2 = cellfun(@(X)nanmedian(X),fin_R2);

jbe_hori = jbeExpts(ori_list(udata(jbeExpts)) == 0);
jbe_vert = jbeExpts(ori_list(udata(jbeExpts)) == 90);

%%
for ii = 1:length(eye_data)
    use_micros = eye_data(ii).use_micros;
    inf_msac_SD(ii) = robust_std_dev(eye_data(ii).inferred_delta_pos(use_micros));
    meas_msac_SD(ii,:) = robust_std_dev(eye_data(ii).measured_delta_pos(use_micros,:));
    msac_diff = bsxfun(@minus,eye_data(ii).measured_delta_pos(use_micros,:),eye_data(ii).inferred_delta_pos(use_micros));
    err_msac_SD(ii,:) = robust_std_dev(msac_diff);
end
    
%% BOXPLOT OF INF-MEAS ERRORS
X = [tot_err(1,jbe_hori)'; tot_err_proc(2,jbe_hori)'; tot_err(1,jbe_vert)'; tot_err_proc(2,jbe_vert)'];
G = [ones(length(jbe_hori),1); 2*ones(length(jbe_vert),1); 3*ones(length(jbe_hori),1); 4*ones(length(jbe_vert),1)];
boxplot(X,G,'labels',{'Vertical','Vertical-proc','Horitontal','Horizontal-proc'});
ylabel('Error SD (deg)');

%% BOXPLOT OF INFERRED AND MEASURED SDs
X = [inf_tot_std(1,jbe_hori)'; meas_tot_std(1,jbe_hori)'; meas_tot_std_proc(2,jbe_hori)'; ...
    inf_tot_std(1,jbe_vert)'; meas_tot_std(1,jbe_vert)'; meas_tot_std_proc(2,jbe_vert)'];
G = [ones(length(jbe_hori),1); 2*ones(length(jbe_hori),1); 3*ones(length(jbe_hori),1); ...
    4*ones(length(jbe_vert),1); 5*ones(length(jbe_vert),1); 6*ones(length(jbe_vert),1)];
boxplot(X,G,'labels',{'Vert-inf','Vert-meas','Vert-Pmeas','Hor-inf','Hor-meas','Hor-Pmeas'});
ylabel('Eye position SD (deg)');

%% BOXPLOT OF INFERRED AND MEASURED SDs with LEM
X = [inf_tot_std(1,jbe_hori)'; meas_tot_std(1,jbe_hori)'; meas_tot_std_proc(2,jbe_hori)'; ...
    inf_tot_std(1,jbe_vert)'; meas_tot_std(1,jbe_vert)'; meas_tot_std_proc(2,jbe_vert)'; ...
    inf_tot_std(1,lemExpts)'; meas_tot_std(1,lemExpts)'; meas_tot_std_proc(2,lemExpts)'; meas_tot_std(2,lemExpts)'; meas_tot_std_proc(4,lemExpts)'];
G = [ones(length(jbe_hori),1); 2*ones(length(jbe_hori),1); 3*ones(length(jbe_hori),1); ...
    4*ones(length(jbe_vert),1); 5*ones(length(jbe_vert),1); 6*ones(length(jbe_vert),1); ...
    7*ones(length(lemExpts),1); 8*ones(length(lemExpts),1); 9*ones(length(lemExpts),1); 10*ones(length(lemExpts),1); 11*ones(length(lemExpts),1);];
boxplot(X,G,'labels',{'Vert-inf','Vert-meas','Vert-Pmeas','Hor-inf','Hor-meas','Hor-Pmeas','lem-inf','lem-meas1','lem-Pmeas1','lem-meas2','lem-Pmeas2'});
ylabel('Eye position SD (deg)');

%% BOXPLOT Microsaccade correlation 
X = [msac_corrs(1,jbe_hori)'; msac_corrs(1,jbe_vert)'; msac_corrs(1,lemExpts)'; msac_corrs(2,lemExpts)'];
G = [ones(length(jbe_hori),1); 2*ones(length(jbe_vert),1); 3*ones(length(lemExpts),1); 4*ones(length(lemExpts),1)];
boxplot(X,G,'labels',{'Vertical','Horitontal','lem-1','lem-2'});
ylabel('Microsaccade correlation');

%% BOXPLOT OF DRIFT CORRELATION
X = [drift_corrs_mfdl(1,jbe_hori)'; drift_corrs_mfdl(1,jbe_vert)'; drift_corrs_mfdl(1,lemExpts)'; drift_corrs_mfdl(2,lemExpts)'];
G = [ones(length(jbe_hori),1); 2*ones(length(jbe_vert),1); 3*ones(length(lemExpts),1); 4*ones(length(lemExpts),1)];
boxplot(X,G,'labels',{'Vertical','Horitontal','lem-1','lem-2'});
ylabel('Microsaccade correlation');


%%
dname = 'all_pooled_eye_data_nosmooth_bsub_hres3';
load(dname)

udata = find(cellfun(@(x)length(x),all_inf_eyepos) > 0);

tot_err = reshape([eye_data(udata).tot_err],[2 length(udata)]);
drift_corrs_mfdl = reshape([eye_data(udata).drift_corrs_mfdl],[2 length(udata)]);
inf_tot_std = [eye_data(udata).inf_tot_std];
inf_tot_var = [eye_data(udata).inf_tot_var];
inf_drift_std = [eye_data(udata).inf_drift_std];
inf_drift_var = [eye_data(udata).inf_drift_var];
meas_tot_std = reshape([eye_data(udata).meas_tot_std],[2 length(udata)]);
msac_corrs = reshape([eye_data(udata).msac_corrs],[2 length(udata)]);
jbeExpts = find(strncmpi(Expt_list(udata),'G',1));
lemExpts = find(strncmpi(Expt_list(udata),'M',1));
jbe_hori = jbeExpts(ori_list(udata(jbeExpts)) == 0);
jbe_vert = jbeExpts(ori_list(udata(jbeExpts)) == 90);

dname = 'all_pooled_eye_data_nosmooth_bsub_hres3_simsac';
load(dname)
udata = find(cellfun(@(x)length(x),all_inf_eyepos) > 0);

ss_tot_err = reshape([eye_data(udata).tot_err],[2 length(udata)]);
ss_drift_corrs_mfdl = reshape([eye_data(udata).drift_corrs_mfdl],[2 length(udata)]);
ss_inf_tot_std = [eye_data(udata).inf_tot_std];
ss_inf_tot_var = [eye_data(udata).inf_tot_var];
ss_inf_drift_std = [eye_data(udata).inf_drift_std];
ss_inf_drift_var = [eye_data(udata).inf_drift_var];
ss_meas_tot_std = reshape([eye_data(udata).meas_tot_std],[2 length(udata)]);
ss_msac_corrs = reshape([eye_data(udata).msac_corrs],[2 length(udata)]);
ss_jbeExpts = find(strncmpi(Expt_list(udata),'G',1));
ss_lemExpts = find(strncmpi(Expt_list(udata),'M',1));
ss_jbe_hori = ss_jbeExpts(ori_list(udata(ss_jbeExpts)) == 0);
ss_jbe_vert = ss_jbeExpts(ori_list(udata(ss_jbeExpts)) == 90);


%% BOXPLOTS COMPARING ALL DATA VS SS-ONLY DATA
close all
X = [inf_tot_std(1,jbe_hori)'; ss_inf_tot_std(1,ss_jbe_hori)'; ...
    meas_tot_std(1,jbe_hori)'; ss_meas_tot_std(1,ss_jbe_hori)';];
G = [ones(length(jbe_hori),1); 2*ones(length(ss_jbe_hori),1); ...
    3*ones(length(jbe_hori),1); 4*ones(length(ss_jbe_hori),1)];
h1 = figure();
boxplot(X,G,'labels',{'Inf','Inf-SS','Meas','Meas-SS'});
ylabel('Eye position SD (deg)');
ylim([0 0.22])

X = [msac_corrs(1,jbe_hori)'; ss_msac_corrs(1,ss_jbe_hori)'; ...
    drift_corrs_mfdl(1,jbe_hori)'; ss_drift_corrs_mfdl(1,ss_jbe_hori)';];
G = [ones(length(jbe_hori),1); 2*ones(length(jbe_hori),1); ...
    3*ones(length(ss_jbe_hori),1); 4*ones(length(ss_jbe_hori),1)];
h2 = figure();
boxplot(X,G,'labels',{'Msac','Msac-SS','Drift','Drift-SS'});
ylabel('Eye position SD (deg)');
ylim([0 1]);

fig_width = 3.27; rel_height = 1.5;
figufy(h1);
fname = [fig_dir 'SD_Simsac_compare.pdf'];
exportfig(h1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(h1);

figufy(h2);
fname = [fig_dir 'Corr_Simsac_compare.pdf'];
exportfig(h2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(h2);

%%
all_R2_SU_imp = [];
all_R2_SU_imp_LOO = [];
all_tot_imp = [];
for ii = 1:length(init_R2)
    avg_init_R2(ii) = median(init_R2{ii}(all_tr_set{ii}));
    avg_fin_R2(ii) = median(fin_R2{ii}(all_tr_set{ii}));
    avg_tot_imp(ii) = median(fin_R2{ii}(all_tr_set{ii})./init_R2{ii}(all_tr_set{ii}));
    
    all_tot_imp = cat(2,all_tot_imp,fin_R2{ii}(all_tr_set{ii})./init_R2{ii}(all_tr_set{ii}));
    
    avg_R2_SU_imp(ii) = median(fin_R2_SU{ii}./init_R2_SU{ii});
    avg_R2_SU_imp_LOO(ii) = median(fin_R2_SU{ii}./init_R2_SU{ii});
    
    all_R2_SU_imp = cat(2,all_R2_SU_imp,fin_R2_SU{ii}./init_R2_SU{ii});
    all_R2_SU_imp_LOO = cat(2,all_R2_SU_imp_LOO,fin_R2_SU_LOO{ii}./init_R2_SU{ii});
end

all_percent_diff_LOO = (all_R2_SU_imp - all_R2_SU_imp_LOO)./all_R2_SU_imp_LOO;
%%
X = [avg_init_R2(jbe_hori)'; avg_fin_R2(jbe_hori)'; avg_init_R2(jbe_vert)'; avg_fin_R2(jbe_vert)'];
G = [ones(length(jbe_hori),1); 2*ones(length(jbe_hori),1);3*ones(length(jbe_vert),1);4*ones(length(jbe_vert),1) ];
boxplot(X,G,'labels',{'Vertical pre','Vertical post','Horitontal pre','Horitontal post'});
ylabel('Microsaccade correlation');

%%
figure
X = [avg_tot_imp(jbe_hori)'; avg_tot_imp(jbe_vert)'];
G = [ones(length(jbe_hori),1); 2*ones(length(jbe_vert),1)];
boxplot(X,G,'labels',{'Vertical','Horitontal'});
ylabel('Microsaccade correlation');

%% MODEL IMP VS MEAS/INF EYE SD
mS = 6;

h = figure;
subplot(2,1,1)
hold on
plot(inf_tot_std(jbe_hori), avg_tot_imp(jbe_hori),'o','markersize',mS);
plot(inf_tot_std(jbe_vert), avg_tot_imp(jbe_vert),'ro','markersize',mS);
xlabel('Inferred SD (deg)');
ylabel('R2 improvement');
legend('Vertical','Horizontal');

subplot(2,1,2)
hold on
plot(meas_tot_std(1,jbe_hori), avg_tot_imp(jbe_hori),'o','markersize',mS);
plot(meas_tot_std(1,jbe_vert), avg_tot_imp(jbe_vert),'ro','markersize',mS);
xlabel('Measured SD (deg)');
ylabel('R2 improvement');
legend('Vertical','Horizontal');

% fig_width = 3.27; rel_height = 1.5;
% figufy(h);
% fname = [fig_dir 'Modimp_vs_eyeSD.pdf'];
% exportfig(h,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(h);

%%
meas_spec(:,2,jbe_hori([3 4])) = nan;
inf_coh(:,2,jbe_hori([3 4])) = nan;

f = eye_data(1).f;
meas_spec = reshape([eye_data(udata).meas_spec],[length(f) 2 length(udata)]);
inf_spec = reshape([eye_data(udata).inf_spec],[length(f) 1 length(udata)]);
inf_coh = reshape([eye_data(udata).coh],[length(f) 2 length(udata)]);
logF = logspace(log10(0.05),log10(30),500);

avg_spec_hori = squeeze(mean(meas_spec(:,:,jbe_hori),3));
std_spec_hori = squeeze(std(meas_spec(:,:,jbe_hori),[],3));
avg_spec_hori = interp1(f,avg_spec_hori,logF);
std_spec_hori = interp1(f,std_spec_hori,logF);
avg_spec_vert = squeeze(mean(meas_spec(:,:,jbe_vert),3));
std_spec_vert = squeeze(std(meas_spec(:,:,jbe_vert),[],3));
avg_spec_vert = interp1(f,avg_spec_vert,logF);
std_spec_vert = interp1(f,std_spec_vert,logF);

avg_inf_hori = squeeze(nanmean(inf_spec(:,:,jbe_hori),3));
std_inf_hori = squeeze(nanstd(inf_spec(:,:,jbe_hori),[],3));
avg_inf_hori = interp1(f,avg_inf_hori,logF);
std_inf_hori = interp1(f,std_inf_hori,logF);
avg_inf_vert = squeeze(nanmean(inf_spec(:,:,jbe_vert),3));
std_inf_vert = squeeze(nanstd(inf_spec(:,:,jbe_vert),[],3));
avg_inf_vert = interp1(f,avg_inf_vert,logF);
std_inf_vert = interp1(f,std_inf_vert,logF);

h1= figure; hold on
shadedErrorBar(logF,avg_spec_hori(:,1),std_spec_hori(:,1),{'color','k'});
shadedErrorBar(logF,avg_spec_vert(:,1),std_spec_vert(:,1),{'color','r'});
shadedErrorBar(logF,avg_inf_hori,std_inf_hori,{'color','b'});
shadedErrorBar(logF,avg_inf_vert,std_inf_vert,{'color','g'});
set(gca,'xscale','log')
xlim([0.05 30])
xlabel('Frequency (Hz)');
ylabel('Log Power');


avg_coh_hori = squeeze(mean(inf_coh(:,1,jbe_hori),3));
std_coh_hori = squeeze(std(inf_coh(:,1,jbe_hori),[],3));
avg_coh_hori = interp1(f,avg_coh_hori,logF);
std_coh_hori = interp1(f,std_coh_hori,logF);
avg_coh_vert = squeeze(nanmean(inf_coh(:,1,jbe_vert),3));
std_coh_vert = squeeze(nanstd(inf_coh(:,1,jbe_vert),[],3));
avg_coh_vert = interp1(f,avg_coh_vert,logF);
std_coh_vert = interp1(f,std_coh_vert,logF);

h2= figure; hold on
shadedErrorBar(logF,avg_coh_hori,std_coh_hori);
shadedErrorBar(logF,avg_coh_vert,std_coh_vert,{'color','r'});
set(gca,'xscale','log')
xlim([0.05 30])
xlabel('Frequency (Hz)');
ylabel('Coherence');

%%
% fig_width = 3.27; rel_height = 0.8;
% figufy(h1);
% fname = [fig_dir 'Eye_spec.pdf'];
% exportfig(h1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(h1);
% 
% figufy(h2);
% fname = [fig_dir 'Eye_coh.pdf'];
% exportfig(h2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(h2);

%%
meas_spec(:,2,jbe_hori([3 4])) = nan;
inf_coh(:,2,jbe_hori([3 4])) = nan;

figure; hold on
errorbar(f,squeeze(mean(meas_spec(:,2,jbe_hori),3)),squeeze(std(meas_spec(:,2,jbe_hori),[],3)),'k');
errorbar(f,squeeze(mean(meas_spec(:,1,jbe_hori),3)),squeeze(std(meas_spec(:,1,jbe_hori),[],3)));
set(gca,'xscale','log');
xlim([0.05 30])

figure; hold on
errorbar(f,squeeze(mean(inf_coh(:,2,jbe_hori),3)),squeeze(std(inf_coh(:,2,jbe_hori),[],3)),'k');
errorbar(f,squeeze(mean(inf_coh(:,1,jbe_hori),3)),squeeze(std(inf_coh(:,1,jbe_hori),[],3)));
set(gca,'xscale','log');
xlim([0.05 30])
