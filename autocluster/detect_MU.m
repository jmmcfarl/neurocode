function [clusterDetails] = detect_MU(sfile,params,use_chs,trig_thresh)

if nargin < 3 || isempty(use_chs)
    use_chs = nan;
end
if nargin < 4
    trig_thresh = [];
    use_fix_thresh = false;
else
    use_fix_thresh = true;
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
if ~isfield(params,'noise_thresh') 
    params.noise_thresh = 3.5; %number of sigma (robust estimate of bkgnd noise) to use as max trigger value
end
if ~isfield(params,'try_features')
    %[PCs voltage energy templates]
    params.try_features = [1]; %which features to try clustering with.
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

if use_fix_thresh
[spk_id] = triggerSpikes_fixthresh(V(:,trig_ch),params.thresh_sign,trig_thresh);    
else
[spk_id, trig_thresh,noise_sigma] = triggerSpikes(V(:,trig_ch),params.thresh_sign,target_Nspks);
end

%check if identified trigger threshold is too high above the estimated
%noise level. If so, lower threshold and retrigger
if ~use_fix_thresh
if trig_thresh/noise_sigma >= params.noise_thresh
    fprintf('Lowering trigger threshold to %.2f sigma\n',params.noise_thresh);
    new_trig = params.noise_thresh*noise_sigma;
    [spk_id, trig_thresh] = triggerSpikes(V(:,trig_ch),params.thresh_sign,target_Nspks,new_trig);
end
end
spk_id(spk_id <= abs(params.spk_pts(1)) | spk_id >= length(V)-params.spk_pts(end)) = []; %get rid of spikes at the edges

%extract spike snippets
Spikes = getSpikeSnippets(V,Vtime,spk_id,params.spk_pts,trig_ch);

% artifact detection
artifact_ids = find_spike_artifacts(Spikes,params);
Spikes.V(artifact_ids,:,:) = [];
Spikes.times(artifact_ids) = [];
Spikes.trig_vals(artifact_ids) = [];
if params.verbose > 0
    fprintf('Removed %d potential artifacts\n',length(artifact_ids));
end

[N_spks, N_samps, N_chs] = size(Spikes.V);


%%
comp_idx = ones(N_spks,1);
cluster_labels = [1];
[cluster_labels, cluster_stats] = relabel_clusters(Spikes.V,comp_idx,cluster_labels);
clusterDetails.mean_spike = cluster_stats.mean_spike;
clusterDetails.std_spike = cluster_stats.std_spike;
clusterDetails.comp_idx = int16(comp_idx);
clusterDetails.cluster_labels = cluster_labels;
clusterDetails.failed = 0;
clusterDetails.Ncomps = 1;

clusterDetails.trig_thresh = trig_thresh;
clusterDetails.trig_ch = trig_ch;
clusterDetails.use_chs = use_chs;
clusterDetails.recDur = length(V)/Fs;
clusterDetails.params = params;
clusterDetails.times = Spikes.times;
clusterDetails.spk_inds = spk_id(:);

uids = clusterDetails.comp_idx > 0;
spike_clusts = int16(nan(size(clusterDetails.comp_idx)));
spike_clusts(uids) = (clusterDetails.cluster_labels(clusterDetails.comp_idx(uids)));
clusterDetails.spike_clusts = spike_clusts(:);


