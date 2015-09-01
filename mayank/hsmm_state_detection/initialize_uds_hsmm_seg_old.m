function [hmm] = initialize_uds_hsmm_seg_old(emiss,Fs,meanparams,UDS_segs,drive_letter)

%%initialize hidden semi-markov model.  meanparams determines whether any
%%of the observations is nonstationary.  currently the function only
%%supports a maximum of one nonstationary emission parameter.  this must be
%%the first column of emiss.  Currently only allows 2 hidden states
hmm.drive_letter = drive_letter;

addpath(strcat(hmm.drive_letter,':\Code\fullBNT-1.0.4\kPMstats\'))
addpath(strcat(hmm.drive_letter,':\Code\fullBNT-1.0.4\netlab3.3\'))
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
mix=gmm(hmm.p,hmm.K,'full');
gmm_options(3) = 1e-10; %tolerance
gmm_options(5) = 1;%reset cov if singular values
gmm_options(14) = 200; %max iterations

options = statset('MaxIter',200);

if strcmp(meantype,'fixed')
    mix=gmm(hmm.p,hmm.K,'full');
    gmm_options(3) = 1e-10; %tolerance
    gmm_options(5) = 1;%reset cov if singular values
    gmm_options(14) = 200; %max iterations
    mix = gmmem(mix, emiss(UDS_samples,:), gmm_options);
    [dummy,state_order] = sort(mix.centres(:,1));
    dstate_fixmean = mix.centres(state_order(1),:);
    ustate_fixmean = mix.centres(state_order(2),:);
    dstate_var = squeeze(mix.covars(:,:,state_order(1)));
    ustate_var = squeeze(mix.covars(:,:,state_order(2)));
else
    %for difference from down state meanfunction
    temp_emiss = emiss(UDS_samples,:);
    temp_emiss(:,1) = emiss(UDS_samples,1)-ov_meanfuns(:,1);
    temp_emiss(any(isnan(temp_emiss),2),:) = [];
    mix_d = gmmem(mix, temp_emiss, gmm_options);
    
    obj = gmdistribution.fit(temp_emiss,hmm.K,'options',options);
    
    dstate = find(abs(mix_d.centres(:,1))==min(abs(mix_d.centres(:,1)))); %find mixture component with 0 mean
    dstate_var = squeeze(mix_d.covars(:,:,dstate));
    if hmm.p > 1
        dstate_fixmean = mix_d.centres(dstate,2:end);
    end
    
    temp_emiss = emiss(UDS_samples,:);
    temp_emiss(:,1) = emiss(UDS_samples,1)-ov_meanfuns(:,2);
     temp_emiss(any(isnan(temp_emiss),2),:) = [];
   mix_u = gmmem(mix, temp_emiss, gmm_options);
    ustate = find(abs(mix_u.centres(:,1))==min(abs(mix_u.centres(:,1)))); %find mixture component with 0 mean
    ustate_var = squeeze(mix_u.covars(:,:,ustate));
    if hmm.p > 1
        ustate_fixmean = mix_u.centres(ustate,2:end);
    end
end

%% create hmm structure
hmm.A = A;
hmm.Pi = Pi;

if hmm.p > 1 || strcmp(meantype,'fixed')
    hmm.state(1).fixedmean = dstate_fixmean;
    hmm.state(2).fixedmean = ustate_fixmean;
end
if strcmp(meantype,'variable')
    for i = 1:hmm.Nsegs
        hmm.state(1).meanfun{i} = meanfuns{i}(:,1);
        hmm.state(2).meanfun{i} = meanfuns{i}(:,2);
    end
end
hmm.state(1).var = dstate_var;
hmm.state(2).var = ustate_var;

hmm.covtype = 'full'; %default

hmm.state(1).dur_type = 'inv_gauss';
hmm.state(2).dur_type = 'gamma';
   
hmm.meantype = meantype;

