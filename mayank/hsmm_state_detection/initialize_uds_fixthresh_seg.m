function [hmm] = initialize_uds_fixthresh_seg(emiss,Fs,UDS_segs)

%%initialize hidden semi-markov model.  meanparams determines whether any
%%of the observations is nonstationary.  currently the function only
%%supports a maximum of one nonstationary emission parameter.  this must be
%%the first column of emiss.  Currently only allows 2 hidden states
hmm.drive_letter = 'G';
addpath(strcat(hmm.drive_letter,':\WC_Germany\hmm_state_detect'))

hmm.p=size(emiss,2);
hmm.T=size(emiss,1);
hmm.Fs = Fs; %sampling frequency of the observations
hmm.K = 2;
hmm.Nsegs = size(UDS_segs,1);
hmmFs = 2016/40; %default sampling frequency for HMM analysis

%these are the indices for the segments of data with UDS
UDS_samples = [];
for i = 1:hmm.Nsegs
    UDS_segind_new_start = round((UDS_segs(i,1)-1)*Fs/hmmFs)+1;
    enddiff = length(emiss)-UDS_segind_new_start+1;
    cT = UDS_segs(i,2)-UDS_segs(i,1)+2;
    curT(i) = min(round(Fs/hmmFs*cT),enddiff);
    UDS_seg_inds{i} = UDS_segind_new_start:...
        UDS_segind_new_start+curT(i)-1;
    UDS_seg_inds{i}(UDS_seg_inds{i} > length(emiss)) = [];
    UDS_segs_new(i,:) = UDS_seg_inds{i}([1 end]); 
    UDS_samples = [UDS_samples UDS_seg_inds{i}];
end

hmm.UDS_segs = UDS_segs_new;
hmm.Nsegs = size(UDS_segs_new,1);
hmm.uds_dur = length(UDS_samples)/Fs;

%% GMM initializations
options = statset('MaxIter',200,'TolFun',1e-8);
obj = gmdistribution.fit(emiss(UDS_samples,:),hmm.K,'Options',options);
[dummy,state_order] = sort(obj.mu(:,1));
state_means = obj.mu(state_order,:);
state_vars = squeeze(obj.Sigma(:,:,state_order));
priors = obj.PComponents(state_order);

%set the threshold as the first value of the observation where the
%probability of the up state exceeds that of the down state
thresh_range = linspace(state_means(1),state_means(2),1e4);
state_lliks(:,1) = log(priors(1)) - 1/2*log(2*pi*state_vars(1)) - ...
    (thresh_range-state_means(1)).^2/(2*state_vars(1));
state_lliks(:,2) = log(priors(2)) - 1/2*log(2*pi*state_vars(2)) - ...
    (thresh_range-state_means(2)).^2/(2*state_vars(2));
threshold = thresh_range(find(state_lliks(:,2) > state_lliks(:,1),1,'first'));

if isempty(threshold)
    threshold = mean(state_means);
end

%% create hmm structure

hmm.state_vars = state_vars;
hmm.state_means =state_means;
hmm.threshold = threshold;
hmm.priors = priors;

