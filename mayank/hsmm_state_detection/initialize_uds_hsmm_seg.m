function [hmm] = initialize_uds_hsmm_seg(emiss,Fs,meanparams,UDS_segs,drive_letter)

%%initialize hidden semi-markov model.  meanparams determines whether any
%%of the observations is nonstationary.  currently the function only
%%supports a maximum of one nonstationary emission parameter.  this must be
%%the first column of emiss.  Currently only allows 2 hidden states
hmm.drive_letter = drive_letter;
addpath(strcat(hmm.drive_letter,':\WC_Germany\hmm_state_detect'))

hmm.p=size(emiss,2); %dimensionality of emissions vector
hmm.T=size(emiss,1); %number of time samples
hmm.Fs = Fs; %sample frequency (Hz)
hmm.K = 2; %number of hidden states
hmm.UDS_segs = UDS_segs; %segments containing UDS
hmm.Nsegs = size(UDS_segs,1); %number of discrete UDS segments

%these are the indices for the segments of data with UDS
UDS_samples = [];
for i = 1:hmm.Nsegs
    UDS_seg_inds{i} = [UDS_segs(i,1):UDS_segs(i,2)];
    UDS_samples = [UDS_samples UDS_segs(i,1):UDS_segs(i,2)];
    time{i} = (UDS_segs(i,1):UDS_segs(i,2))/Fs;
end
hmm.uds_dur = length(UDS_samples)/hmm.Fs; %total duration of UDS (seconds)

meantype = meanparams.meantype;
%if meantype is set to variable, initialize state mean function parameters
if strcmp(meantype,'variable')
    hmm.windowSize = meanparams.windowSize;
    hmm.windowSlide = meanparams.windowSlide;
    %set up a low-pass filter for the meanfunctions to help denoise the
    %initialization
    slid_avg_hcf = 0.02; 
    niqf = Fs/2;
    [slid_b,slid_a] = butter(2,slid_avg_hcf/niqf,'low');
end

%% choose random starting values for A and Pi from dir priors
%use a uniform prior on pi parameter (this is inconsequential, but
%is kept for generality)
pi_alpha = 1; 
Pi_alpha = ones(1,hmm.K)*pi_alpha; 
for i = 1:hmm.Nsegs
    Pi(i,:) = dirichletrnd(Pi_alpha);
end

%intialize the transition matrix to have expected state durations of 1s
%under the geometric distribution
A_other = 1/Fs/(hmm.K-1);
A_self = 1-A_other;
for i = 1:hmm.K
    for j = 1:hmm.K
        if i==j
            A(i,j)=A_self;
        else
            A(i,j)=A_other;
        end
    end
end
%% if we are using a nonstationary emission parameter, estimate its time
%% varying state means
if strcmp(meantype,'variable')
    ov_meanfuns = [];
    for i = 1:hmm.Nsegs
        [state_means,state_t] = ...
            locate_state_means(emiss(UDS_seg_inds{i},1),hmm.windowSize,hmm.windowSlide,Fs);
        state_t = state_t + UDS_segs(i,1)/Fs;
        state_t  = [time{i}(1) state_t (time{i}(end)+1)];
        state_means = [state_means(:,1) state_means state_means(:,end)];
        meanfuns{i} = interp1(state_t,state_means',time{i}); %interpolate state mean functions onto the time grid
        %low-pass filter the state meanfunctions
        for k = 1:hmm.K
            meanfuns{i}(:,k) = filtfilt(slid_b,slid_a,meanfuns{i}(:,k));
        end
        ov_meanfuns = [ov_meanfuns; meanfuns{i}];
    end
end

%% GMM initializations
options = statset('MaxIter',200,'TolFun',1e-8);

if strcmp(meantype,'fixed')
    %if using fixed state means fit a K state mixture model and set the
    %state with the larger average as the up state
    obj = gmdistribution.fit(emiss(UDS_samples,:),hmm.K,'Options',options);
    [~,state_order] = sort(obj.mu(:,1));
    dstate_fixmean = obj.mu(state_order(1),:);
    ustate_fixmean = obj.mu(state_order(2),:);
    dstate_var = squeeze(obj.Sigma(:,:,state_order(1)));
    ustate_var = squeeze(obj.Sigma(:,:,state_order(2)));
else %for time-varying state means we need to estimate state-conditional covariance matrices
    
    %fit a mixture model to the difference from the down-state avg from the
    %sliding window density estimate
    temp_emiss = emiss(UDS_samples,:);
    temp_emiss(:,1) = emiss(UDS_samples,1)-ov_meanfuns(:,1);
    temp_emiss(any(isnan(temp_emiss),2),:) = [];
    obj = gmdistribution.fit(temp_emiss,hmm.K,'Options',options);
    dstate = find(abs(obj.mu(:,1)) == min(abs(obj.mu(:,1)))); %the down state is the state with mean closest to 0
    dstate_var = squeeze(obj.Sigma(:,:,dstate));
    if hmm.p > 1
        dstate_fixmean = obj.mu(dstate,2:end);
    end
    
    %repeat for the up state
    temp_emiss = emiss(UDS_samples,:);
    temp_emiss(:,1) = emiss(UDS_samples,1)-ov_meanfuns(:,2);
    temp_emiss(any(isnan(temp_emiss),2),:) = [];
    obj = gmdistribution.fit(temp_emiss,hmm.K,'Options',options);
    ustate = find(abs(obj.mu(:,1)) == min(abs(obj.mu(:,1))));
    ustate_var = squeeze(obj.Sigma(:,:,ustate));
    if hmm.p > 1
        ustate_fixmean = obj.mu(ustate,2:end);
    end
end
        
%% create hmm structure
hmm.A = A; %transition probability matrix
hmm.Pi = Pi; %initial state distribution

if hmm.p > 1 || strcmp(meantype,'fixed') %for fixed state means
    hmm.state(1).fixedmean = dstate_fixmean;
    hmm.state(2).fixedmean = ustate_fixmean;
end
if strcmp(meantype,'variable') %for variable state mean functions
    for i = 1:hmm.Nsegs
        hmm.state(1).meanfun{i} = meanfuns{i}(:,1);
        hmm.state(2).meanfun{i} = meanfuns{i}(:,2);
    end
end
hmm.state(1).var = dstate_var;
hmm.state(2).var = ustate_var;

hmm.covtype = 'full'; %default, full covariance matrix

%state duration distribution models
hmm.state(1).dur_type = 'inv_gauss';
hmm.state(2).dur_type = 'inv_gauss';
   
hmm.meantype = meantype;

