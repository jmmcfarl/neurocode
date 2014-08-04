clear all
% close all
addpath('~/James_scripts/bruce/eye_tracking/');

% Expt_list = {'G085','G086','G087','G088','G089','G091','G093','G095'};
Expt_list = {'M266','M270','M275','M277','M281','M287','M289','M294'};
% ori_list = ones(size(Expt_list))*0;

ori_list = [80 60 135 70 140 90 160 40];
single_bar_expts = [232 235 239];

for EE = 1:length(Expt_list)
    clear expt_* included_type
    
    Expt_num = str2num(Expt_list{EE}(2:end));
    Expt_name = Expt_list{EE};
    bar_ori = ori_list(EE);
    
    fprintf('Analyzing expt %s, ori %d, %d/%d\n',Expt_name,ori_list(EE),EE,length(Expt_list));
    
    if Expt_num < 280
    data_dir = ['~/Data/bruce/' Expt_name];
    elseif Expt_num <289
    data_dir = ['/media/NTlab_data2/Data/bruce/' Expt_name];
    else
    data_dir = ['/media/NTlab_data3/Data/bruce/' Expt_name];
    end
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
    min_trial_dur = 0.75;

    beg_buffer = 0.2;
    end_buffer = 0.05;
    if ~ismember(Expt_num,single_bar_expts)
        trial_dur = 4;
    else
        trial_dur = 2;
    end
    
    use_right_eye = false;
    
    n_use_blocks = Inf;
        
    if Expt_name(1) == 'G'
        n_probes = 96;
    else 
        n_probes = 24;
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
        expt_names = cell(1,length(Expts));
        expt_dds = nan(1,length(Expts));
        expt_bar_ori = nan(1,length(Expts));
        expt_sac_dir = nan(1,length(Expts));
        expt_Fr = nan(1,length(Expts));
        expt_imback = nan(1,length(Expts));
        included_type = false(1,length(Expts));
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
    
    %% load overall su data
    % LOAD REFCLUSTERS
    cluster_dir = ['~/Analysis/bruce/' Expt_name '/clustering'];
    fname = [cluster_dir '/final_cluster.mat'];
    if exist(fname,'file')
        load(fname);
        SU_numbers = unique(SU_ID_mat(~isnan(SU_ID_mat)));
        for ii = 1:length(SU_numbers)
            SU_tot_nblocks = sum(SU_ID_mat(:) == SU_numbers(ii));
        end
        fprintf('%d SUs Clustered\n',length(SU_numbers));
        
    else
        disp('No Cluster data found.');
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
    all_spk_times = cell(n_probes,1);
    all_clust_ids = cell(n_probes,1);
    all_spk_inds = cell(n_probes,1);
    
    trial_toffset = zeros(length(cur_block_set),1);
    cur_spkind_offset = 0;
    cur_toffset = 0;
    for ee = 1:n_blocks;
        fprintf('Expt %d Block %d of %d\n',Expt_num,ee,n_blocks);
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
            cur_spkind_offset = cur_spkind_offset+ round(32e3*(max(trial_end_times)-min(trial_start_times) + 50));
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
        
%% BIN SPIKES FOR MU AND SU
% LOAD REFCLUSTERS
fname = [cluster_dir '/final_cluster.mat'];
if exist(fname,'file')
    load(fname);
    SU_numbers = unique(SU_ID_mat(~isnan(SU_ID_mat)));
    for ii = 1:length(SU_numbers)
        SU_tot_nblocks = sum(SU_ID_mat(:) == SU_numbers(ii));
    end
    fprintf('%d SUs Clustered\n',length(SU_numbers));
    
else
    disp('No Cluster data found.');
end

%for SU probes
all_binned_sua = nan(length(all_t_axis),length(SU_numbers));
cur_used_blocks = find(ismember(SU_target_blocks,cur_block_set));
SU_ID_mat = SU_ID_mat(cur_used_blocks,:);
[CC,BB] = meshgrid(1:length(SU_clust_data),1:length(cur_used_blocks));

su_probes = nan(1,length(SU_numbers));
all_su_spk_times = cell(length(SU_numbers),1);
all_su_spk_inds = cell(length(SU_numbers),1);
clear SU_block_probes
for ss = 1:length(SU_numbers)
    used_clust_set = unique(CC(SU_ID_mat==SU_numbers(ss))); %set of clusters used to capture this SU
    SU_block_probes(ss,:) = nan(1,length(cur_block_set));
    cur_su_spk_times = [];
    cur_su_spk_inds = [];
    cur_blocks = [];
    for cc = 1:length(used_clust_set)
        cur_clust = used_clust_set(cc);
        cur_probe = SU_clust_data(cur_clust).probe_num;
        cur_clust_label = SU_clust_data(cur_clust).cluster_label;
        cur_blocks = [cur_blocks; find(SU_ID_mat(:,cur_clust) == SU_numbers(ss))];
        SU_block_probes(ss,cur_blocks) = cur_probe;
        
        all_su_inds = all_clust_ids{cur_probe} == cur_clust_label;
        cur_su_spk_times = all_spk_times{cur_probe}(all_su_inds);
        cur_su_spk_inds = all_spk_inds{cur_probe}(all_su_inds);
        spk_block_inds = round(interp1(all_t_axis,all_blockvec,cur_su_spk_times));
        cur_su_spk_times = cur_su_spk_times(ismember(spk_block_inds,cur_blocks));   
        cur_su_spk_inds = cur_su_spk_inds(ismember(spk_block_inds,cur_blocks));
        
        all_su_spk_times{ss} = cat(1,all_su_spk_times{ss},cur_su_spk_times(:));
        all_su_spk_inds{ss} = cat(1,all_su_spk_inds{ss},cur_su_spk_inds(:));
    end
    if ~isempty(all_su_spk_times{ss})
    cur_suahist = histc(all_su_spk_times{ss},all_t_bin_edges);
    cur_suahist(all_bin_edge_pts) = [];
    cur_id_set = ismember(all_blockvec,cur_blocks);
    all_binned_sua(cur_id_set,ss) = cur_suahist(cur_id_set);
    su_probes(ss) = mode(SU_block_probes(ss,~isnan(SU_block_probes(ss,:))));
    end
end

double_spike_buffer = 3; %number of samples (in either direction) to exclude double spikes from adjacent-probe SUs
all_binned_mua = nan(length(all_t_axis),n_probes);
for cc = 1:n_probes
    %this is the set of blocks where this probe had an SU, and the
    %correspodning SU numbers
    cur_set = find(SU_block_probes == cc);
    if ~isempty(cur_set)
        [cur_SS,cur_BB] = ind2sub([length(SU_numbers) length(cur_used_blocks)],cur_set);
    else
        cur_SS = [];
    end
    unique_su_nums = unique(cur_SS);
    cur_mua_inds = find(all_clust_ids{cc} >= 1);
    
    %remove spikes from isolated SUs on the same probe from the MU
    for ss = 1:length(unique_su_nums)
        cur_mua_inds(ismember(all_spk_inds{cc}(cur_mua_inds),all_su_spk_inds{unique_su_nums(ss)})) = [];
    end

    nearby_probes = [cc-1 cc+1]; nearby_probes(nearby_probes < 1 | nearby_probes > n_probes) = [];
    cur_set = find(ismember(SU_block_probes,nearby_probes));
    if ~isempty(cur_set)
        [cur_SS,cur_BB] = ind2sub([length(SU_numbers) length(cur_used_blocks)],cur_set);
    else
        cur_SS = [];
    end
    unique_su_nums = unique(cur_SS); %set of SU numbers picked up on adjacent probes
    if ~isempty(unique_su_nums)
        double_spikes = [];
        for ss = 1:length(unique_su_nums)
            cur_blocked_inds = bsxfun(@plus,all_su_spk_inds{unique_su_nums(ss)},-double_spike_buffer:double_spike_buffer);
            double_spikes = [double_spikes; find(ismember(all_spk_inds{cc}(cur_mua_inds),cur_blocked_inds))];
        end
        fprintf('Eliminating %d of %d double spikes in MUA\n',length(double_spikes),length(cur_mua_inds));
        cur_mua_inds(double_spikes) = [];
    end
    [cur_spkhist,cur_bin_inds] = histc(all_spk_times{cc}(cur_mua_inds),all_t_bin_edges);
    cur_spkhist(all_bin_edge_pts) = [];
    all_binned_mua(:,cc) = cur_spkhist;
end
    
       
    %% CREATE EVENT PREDICTORS FOR REAL AND SIM SACCADES (POOLING MICRO AND MACRO SACS)
    Xblock = zeros(length(all_stim_times),n_blocks);
    for i = 1:n_blocks
        cur_set = find(all_blockvec==i);
        Xblock(cur_set,i) = 1;
    end
    
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
    clear cur_sua_data
    for ii = 1:length(su_probes)
        cur_sua_data(ii).probe_num = su_probes(ii);
        cur_sua_data(ii).su_num = SU_numbers(ii);
        cur_sua_data(ii).num_blocks = sum(~isnan(SU_block_probes(ii,:)),2);
        cur_sua_data(ii).n_spikes = nansum(all_binned_sua(used_inds,ii));
        cur_sua_data(ii).mean_rate = nanmean(all_binned_sua(used_inds,ii));
    end
    clear cur_mua_data
    for ii = 1:n_probes
        cur_mua_data(ii).probe_num = ii;
        cur_mua_data(ii).n_spikes = nansum(all_binned_mua(used_inds,ii));
        cur_mua_data(ii).mean_rate = nanmean(all_binned_mua(used_inds,ii));
    end
    Expt_sua_data{EE}.Expt_num = Expt_num;
    Expt_sua_data{EE}.bar_ori = bar_ori;
    Expt_sua_data{EE}.sua_data = cur_sua_data;
    Expt_sua_data{EE}.mua_data = cur_mua_data;
end

%%
cd ~/Analysis/bruce/ET_final/
if bar_ori == 0
    sname = 'SUA_sum_data_hor';
else
    sname = 'SUA_sum_data_ver';
end
% sname = 'SUA_sum_data_lem';
save(sname,'Expt_sua_data');
