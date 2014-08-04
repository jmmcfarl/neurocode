function [clusterDetails,spike_features,sum_fig] = detect_and_cluster_init(sfile,params,use_chs)

if nargin < 3 || isempty(use_chs)
    use_chs = nan;
end

%% DEFAULT PARAMETERS
if nargin < 2 || isempty(params)
    params = struct();
end
if ~isfield(params,'filt_cutoff') %high-pass for filtering raw V for spike detection
    params.filt_cutoff = [100 nan]; 
end
if ~isfield(params,'add_Vmean') %whether or not to add in across-electrode average FullV
    params.add_Vmean = 0;
end
if ~isfield(params,'thresh_sign') %whether to detect on peaks or valleys
    params.thresh_sign = -1;
end
if ~isfield(params,'target_rate') %target spike detection rate. 
    params.target_rate = 50;
end
if ~isfield(params,'spk_pts')
    params.spk_pts = [-12:27]; %set of time points relative to trigger to use for classification
end
if ~isfield(params,'outlier_thresh')
    params.outlier_thresh = 7; %threshold on Mahal distance (sqrt) to count a point as an outlier
end
if ~isfield(params,'verbose')
    params.verbose = 2; %display text
end
if ~isfield(params,'use_best_only')
    params.use_best_only = 0; %use only the best spike waveform for calculating template features.
end
if ~isfield(params,'cluster_bias')
    params.cluster_bias = 0.85; %bias to control type 2 errors for classifying SUs
end
if ~isfield(params,'summary_plot')
    params.summary_plot = 1; 
end
if ~isfield(params,'reg_lambda')
    params.reg_lambda = 0; %regularization on Cov Mats for EM (Matlab based, doesnt seem to help)
end
if ~isfield(params,'n_pcs')
    params.n_pcs = 4; %number of PCs
end
if ~isfield(params,'n_tdims')
    params.n_timds = 4; %number of voltage dimensions
end
if ~isfield(params,'try_features')
    %[PCs voltage energy templates]
    params.try_features = [1 2 4]; %which features to try clustering with.
end
if ~isfield(params,'min_Pcomp')
    min_Pcomp = 0.005; %minimum cluster probability (basically minimum firing rate)
end
if ~isfield(params,'max_n_retriggers')
    params.max_n_retriggers = 3; %max number of times to try retriggering
end
if ~isfield(params,'noise_thresh') 
    params.noise_thresh = 3.5; %number of sigma (robust estimate of bkgnd noise) to use as max trigger value
end
if ~isfield(params,'max_back_comps')
    params.max_back_comps = 3; %maximum number of Gaussians to try modeling background noise with
end

%% LOAD VOLTAGE SIGNAL

%loads in (high-pass filtered) voltage signal
if isstr(sfile)
    [V,Vtime,Fs] = Load_FullV(sfile, params.add_Vmean, params.filt_cutoff,use_chs);
else
    V = sfile.V(:,use_chs);
    Vtime = sfile.Vtime;
    Fs = sfile.Fs;
end
%% DETECT SPIKES
if strcmp(params.target_rate,'median')
    target_Nspks = 'median';
else
    target_Nspks = params.target_rate*length(V)/Fs;
end

if length(use_chs) == 1
    trig_ch = 1;
elseif length(use_chs) == 2
    if use_chs(1) == 1
        trig_ch = 1 ;
    else
        trig_ch = 2;
    end
else
    trig_ch = 2;
end

[spk_id, trig_thresh,noise_sigma] = triggerSpikes(V(:,trig_ch),params.thresh_sign,target_Nspks);

%check if identified trigger threshold is too high above the estimated
%noise level. If so, lower threshold and retrigger
if trig_thresh/noise_sigma >= params.noise_thresh
    fprintf('Lowering trigger threshold to %.2f sigma\n',params.noise_thresh);
    new_trig = params.noise_thresh*noise_sigma;
    [spk_id, trig_thresh] = triggerSpikes(V(:,trig_ch),params.thresh_sign,target_Nspks,new_trig);
end
spk_id(spk_id <= abs(params.spk_pts(1)) | spk_id >= length(V)-params.spk_pts(end)) = []; %get rid of spikes at the edges

%extract spike snippets
Spikes = getSpikeSnippets(V,Vtime,spk_id,params.spk_pts,trig_ch);

% artifact detection
artifact_ids = find_spike_artifacts(Spikes,params);
Spikes.V(artifact_ids,:,:) = [];
Spikes.times(artifact_ids) = [];
Spikes.trig_vals(artifact_ids) = [];
spk_id(artifact_ids) = []; %bug-fix 12/17/2013
if params.verbose > 0
    fprintf('Removed %d potential artifacts\n',length(artifact_ids));
end

[N_spks, N_samps, N_chs] = size(Spikes.V);

%% GMM AUTOCLUSTERING USING MULTIPLE FEATURES AND MULTIPLE INITIALIZATIONS. 
params.trig_ch = trig_ch;
[clusterDetails, spike_xy, spike_features] = autocluster_init(Spikes,params);
if ischar(sfile)
    clusterDetails.rawV_file = sfile;
end
clusterDetails.trig_thresh = trig_thresh;
clusterDetails.trig_ch = trig_ch;
clusterDetails.use_chs = use_chs;
clusterDetails.recDur = length(V)/Fs;
clusterDetails.params = params;
clusterDetails.times = Spikes.times;
clusterDetails.spk_inds = spk_id(:);

clusterDetails = compute_cluster_stats(clusterDetails,Spikes,spike_features);

%% Check if trigger threshold is too high
if clusterDetails.iso_dists(1) > 2 %if there is a reasonable SU
    needToRetrigger = true;
    n_retriggers = 0;
    while needToRetrigger && n_retriggers < params.max_n_retriggers
        su_inds = find(clusterDetails.spike_clusts == 2);
        [nn,xx] = hist(Spikes.trig_vals(su_inds),500);
        nn = smooth(nn,50)/sum(nn); %smoothed spike amp density
        prev_rate = sum(clusterDetails.n_spks)/clusterDetails.recDur;
        new_rate = prev_rate*2; %try doubling trigger rate
        if params.thresh_sign > 0
            nn = flipud(nn);
        end
        if nn(end) < max(nn)/5 || new_rate > 250
            needToRetrigger = false;
            break;
        end
        fprintf('Possible SU cutoff detected, trying retriggering\n');
        new_params = params;
        new_params.target_rate = new_rate;
        fixed = 0; %retrigger and fit new model params
        [clusterDetails,spike_features,spike_xy,Spikes] = apply_clustering(sfile,clusterDetails,new_params,fixed);
        clusterDetails = compute_cluster_stats(clusterDetails,Spikes,spike_features);
        
        %may need additional components to model background with new
        %triggering
        cur_n_back_comps = clusterDetails.Ncomps - 1;
        while cur_n_back_comps < params.max_back_comps
            cur_n_back_comps = cur_n_back_comps + 1;
            fprintf('Trying background split %d of %d\n',cur_n_back_comps,params.max_back_comps);
            
            [spike_xy,xyproj_mat] = Project_GMMfeatures(spike_features, clusterDetails.gmm_fit,clusterDetails.cluster_labels);
            %implement a sign convention that positive values of the first dimension
            %correspond to higher spike amplitudes
            bb = corr(spike_xy(:,1),abs(Spikes.trig_vals(:)));
            if bb < 0
                spike_xy(:,1) = -spike_xy(:,1);
                xyproj_mat(:,1) = -xyproj_mat(:,1);
            end
            [new_GMM_obj, new_distance,new_comp_idx, new_clust_labels] = ...
                try_backgnd_splitting(clusterDetails.gmm_fit,Spikes.V,spike_features,spike_xy,xyproj_mat,clusterDetails.comp_idx,clusterDetails.cluster_labels,params);
            
            fprintf('Orig: %.3f New: %.3f \n',clusterDetails.dprime,new_distance);
            if new_distance > clusterDetails.dprime
                clusterDetails.gmm_fit = new_GMM_obj;
                clusterDetails.comp_idx = new_comp_idx;
                clusterDetails.cluster_labels = new_clust_labels;
                clusterDetails = compute_cluster_stats(clusterDetails,Spikes,spike_features);                
            end
        end

        n_retriggers = n_retriggers + 1;
    end
end


%% CHECK FOR ADDITIONAL SUs BY FITTING MODELS TO THE BACKGROUND SPIKES RECURSIVELY
max_n_sus = 5;
if clusterDetails.iso_dists(1) > 2
    cur_n_SUs = 1;
    new_cluster = clusterDetails;
    while cur_n_SUs <= max_n_sus
        fprintf('Trying to fit %d SUs\n',cur_n_SUs+1);
        
        cur_back_spikes = find(new_cluster.spike_clusts == 1); %current background spikes
        [back_GMM, back_dist,back_comp_idx,back_clust_labels,back_cluster_stats] = ...
            GMM_fit(Spikes.V(cur_back_spikes,:),spike_features(cur_back_spikes,:),2,params);
        if isobject(back_GMM)
            cur_n_back_comps = 1;
            while cur_n_back_comps < params.max_back_comps
                cur_n_back_comps = cur_n_back_comps + 1;
                fprintf('Trying background split %d of %d\n',cur_n_back_comps,params.max_back_comps);
                
                [temp_spike_xy,temp_xyproj_mat] = Project_GMMfeatures(spike_features(cur_back_spikes,:), back_GMM,back_clust_labels);
                %implement a sign convention that positive values of the first dimension
                %correspond to higher spike amplitudes
                bb = corr(temp_spike_xy(:,1),abs(Spikes.trig_vals(cur_back_spikes)));
                if bb < 0
                    temp_spike_xy(:,1) = -temp_spike_xy(:,1);
                    temp_xyproj_mat(:,1) = -temp_xyproj_mat(:,1);
                end
                [newback_GMM_obj, newback_distance,newback_comp_idx, newback_clust_labels] = ...
                    try_backgnd_splitting(back_GMM,Spikes.V(cur_back_spikes,:),spike_features(cur_back_spikes,:),...
                    temp_spike_xy,temp_xyproj_mat,back_comp_idx,back_clust_labels,params);
                
                fprintf('Orig: %.3f New: %.3f \n',back_dist,newback_distance);
                if newback_distance > back_dist
                    back_GMM = newback_GMM_obj;
                    back_dist = newback_distance;
                    back_comp_idx = newback_comp_idx;
                    back_clust_labels = newback_clust_labels;
                end
            end
            
            %             init_comp_idx = clusterDetails.comp_idx;
            %             init_comp_idx = init_comp_idx + length(back_clust_labels);
            %             init_comp_idx(cur_back_spikes) = back_comp_idx;
            
            init_comp_idx = clusterDetails.comp_idx;
            uids = init_comp_idx <= 0;
            buids = back_comp_idx <= 0;
            cuids = setdiff(1:N_spks,[uids; cur_back_spikes(buids)]);
            init_comp_idx = init_comp_idx + length(back_clust_labels);
            init_comp_idx(cur_back_spikes) = back_comp_idx;
            [~,~,init_comp_idx(cuids)] = unique(init_comp_idx(cuids)); %bug fixed 11-15-13
            
            init_cluster_labels = clusterDetails.cluster_labels;
            init_cluster_labels(init_cluster_labels == 1) = [];
            init_cluster_labels = [back_clust_labels init_cluster_labels + 1];
            
            n_comps = length(init_cluster_labels);
            [back_GMM, back_dist,back_comp_idx,back_clust_labels,back_cluster_stats] = ...
                GMM_fit(Spikes.V,spike_features,n_comps,params,[],init_comp_idx,init_cluster_labels);
            
            [new_Lratio,new_iso_dist] = compute_cluster_Lratio(spike_features,back_GMM,back_comp_idx,back_clust_labels);
            n_good_SUs = sum(new_iso_dist > 2 & new_Lratio < 1e3);
            if n_good_SUs > cur_n_SUs
                clusterDetails.gmm_fit = back_GMM;
                clusterDetails.comp_idx = back_comp_idx;
                clusterDetails.cluster_labels = back_clust_labels;
                clusterDetails = compute_cluster_stats(clusterDetails,Spikes,spike_features);  
                new_cluster = clusterDetails;
            else
                break;
            end
        end
        cur_n_SUs = cur_n_SUs + 1;
    end
end

%% CREATE SUMMARY FIGURE
if params.summary_plot > 0
    
    sum_fig = create_summary_cluster_fig(clusterDetails,Spikes,spike_xy,clusterDetails.params);
    
end
