function [state_seq,fhmm] = hsmm_uds_get_fixthresh_uds_state_seq_mix(sig_features,Fs,UDS_segs)

% use fixed threshold crossing method to infer state sequence
%
% Input:
%       hmm: structure with the initialized parameters of an hmm
%       emiss: observation sequence
%
% Output:
%       hmm: structure of parameters with the following fields
%           p: dimension of emissions
%           T: number of time samples
%           Fs: sample frequency of emissions

fhmm.Nsegs = size(UDS_segs,1);
UDS_samples = [];
for i = 1:size(UDS_segs,1)
   UDS_samples = [UDS_samples UDS_segs(i,1):UDS_segs(i,2)]; 
end

%% initialize hsmm
options = statset('MaxIter',200,'TolFun',1e-8);
obj = gmdistribution.fit(sig_features(UDS_samples,:),2,'Options',options);
[~,state_order] = sort(obj.mu(:,1));
state_means = obj.mu(state_order,:);
state_vars = squeeze(obj.Sigma(:,:,state_order));
priors = obj.PComponents(state_order);

%set the threshold as the minimum of the density estimate in the range
%between the two state means
thresh_range = linspace(state_means(1),state_means(2),1e4);
np_dens = ksdensity(sig_features(UDS_samples,:),thresh_range);
[pk_amps,pk_locs] = findpeaks(-np_dens);
if ~isempty(pk_amps)
    [~,biggest] = max(pk_amps);
    thresh_loc = pk_locs(biggest);
else
    [~,thresh_loc] = min(np_dens);
end
threshold1 = thresh_range(thresh_loc);

%set the threshold as the first value of the observation where the
%probability of the up state exceeds that of the down state
state_lliks(:,1) = log(priors(1)) - 1/2*log(2*pi*state_vars(1)) - ...
    (thresh_range-state_means(1)).^2/(2*state_vars(1));
state_lliks(:,2) = log(priors(2)) - 1/2*log(2*pi*state_vars(2)) - ...
    (thresh_range-state_means(2)).^2/(2*state_vars(2));
threshold2 = thresh_range(find(state_lliks(:,2) > state_lliks(:,1),1,'first'));
if isempty(threshold2)
    threshold2 = mean(state_means);
end

% create hmm structure
fhmm.state_vars = state_vars;
fhmm.state_means =state_means;
fhmm.threshold = 0.5*(threshold1+threshold2);
fhmm.priors = priors;

%set all instances where the observation is greater than threshold to the
%up state
for i = 1:fhmm.Nsegs   
    cur_data = sig_features(UDS_segs(i,1):UDS_segs(i,2));
    state_seq{i} = ones(size(cur_data));
    state_seq{i}(cur_data > fhmm.threshold) = 2;  
end
