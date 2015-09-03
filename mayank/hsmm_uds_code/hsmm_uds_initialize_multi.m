function [hmm] = hsmm_uds_initialize_multi(emiss,Fs,params)

% initialize hidden semi-markov model for inferring UP-DOWN states. Only
% for bivariate signal features
%
% Input:
%       emiss:      column vector of signal feature ('emissions' or 'observation sequence'
%       to be used for inferring UDS
%       Fs:    sample frequency of emiss vector
%       params:     parameter structure containing the following fields
%            UDS_segs: desired down-sampling factor for computing spectrogram
%               meantype: type of state-conditional Gaussian mean parameters.
%               Can be 'variable' or 'fixed'. Default is 'variable'
%            movingwin: if using variable state-conditional means, this
%               determines how the data is windowed for sliding-window density
%               estimation. [windowSize windowSlide]. In seconds. Default:
%               [50 1].

% Output:
%        hmm: structure containing information about the initialized HMM

if size(emiss,2) > size(emiss,1)
    error('Dimensions of emiss are incorrect');
end
p=size(emiss,2); %dimensionality of emissions vector (only works for d=1)
hmm.T=size(emiss,1); %number of time samples
hmm.Fs = Fs; %sample frequency (Hz)
hmm.K = 2; %number of hidden states
hmm.p = p;

if isfield(params,'UDS_segs')
    UDS_segs = params.UDS_segs;
else
    UDS_segs = [1 hmm.T];
end
hmm.UDS_segs = UDS_segs; %segments containing UDS
hmm.Nsegs = size(UDS_segs,1); %number of discrete UDS segments

%these are the indices for the segments of data with UDS
UDS_samples = [];
for i = 1:hmm.Nsegs
    UDS_seg_inds{i} = [UDS_segs(i,1):UDS_segs(i,2)]; %cell array containing the indices of UDS data separately within each UDS segment
    UDS_samples = [UDS_samples UDS_segs(i,1):UDS_segs(i,2)]; %vector of all indices of UDS data
    time{i} = (UDS_segs(i,1):UDS_segs(i,2))/Fs; %cell array of time axes for each UDS segment
end
hmm.uds_dur = length(UDS_samples)/hmm.Fs; %total duration of UDS (seconds)

if isfield(params,'meantype')
    meantype = params.meantype;
else
    meantype = 'variable';
end
%if meantype is set to variable, initialize state mean function parameters
if strcmp(meantype,'variable')
    if isfield(params,'movingwin')
        hmm.windowSize = params.movingwin(1);
        hmm.windowSlide = params.movingwin(2);
    else
        hmm.windowSize = 50;
        hmm.windowSlide = 1;
    end
    %set up a low-pass filter for the meanfunctions to help denoise the
    %initialization
    slid_avg_hcf = 0.02; %in Hz
    niqf = Fs/2;
    [slid_b,slid_a] = butter(2,slid_avg_hcf/niqf,'low');
end

%set the distribution on the initial state probability to be uniform (for
%each UDS segment
for i = 1:hmm.Nsegs
    Pi(i,:) = [0.5 0.5];
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

    %loop over dimensions
    for d = 1:p
        ov_meanfuns = [];
        all_meandiffs = [];
        for i = 1:hmm.Nsegs
            [cur_peaks{i},state_t{i},emissdist,emissrange] = hsmm_uds_init_state_means(emiss(UDS_seg_inds{i},d),[hmm.windowSize hmm.windowSlide],Fs);
            all_meandiffs = [all_meandiffs; diff(cur_peaks{i},1,2)];
        end
        avg_meandiff = nanmean(all_meandiffs);
        for i = 1:hmm.Nsegs
            no_uppeak = isnan(cur_peaks{i}(:,2));
            no_downpeak = isnan(cur_peaks{i}(:,1));
            cur_peaks{i}(no_uppeak,2) = cur_peaks{i}(no_uppeak,1) + avg_meandiff;
            cur_peaks{i}(no_downpeak,1) = cur_peaks{i}(no_downpeak,2) - avg_meandiff;
            
            state_t{i} = state_t{i} + UDS_segs(i,1)/Fs;
            state_t{i} = [time{i}(1) state_t{i} (time{i}(end)+1)];
            cur_peaks{i} = [cur_peaks{i}(1,:); cur_peaks{i}; cur_peaks{i}(end,:)];
            meanfuns{i}(d,:,:) = interp1(state_t{i},cur_peaks{i},time{i});
            %low-pass filter the state meanfunctions
            for k = 1:hmm.K
                meanfuns{i}(d,:,k) = filtfilt(slid_b,slid_a,meanfuns{i}(d,:,k));
            end
            ov_meanfuns = cat(2,ov_meanfuns,meanfuns{i});
        end
    end
end

%% OPTIONAL DENSITY PLOT
% [~,temp_t,emissdist_d1,emissrange_d1] = hsmm_uds_init_state_means(emiss(:,1),[hmm.windowSize hmm.windowSlide],Fs);
% [~,temp_t,emissdist_d2,emissrange_d2] = hsmm_uds_init_state_means(emiss(:,2),[hmm.windowSize hmm.windowSlide],Fs);
% 
% figure
% subplot(2,1,1)
% pcolor(temp_t,emissrange_d1,log(emissdist_d1'));shading flat
% hold on
% caxis([-4 -0.5])
% for i = 1:hmm.Nsegs
%     plot(time{i},meanfuns_d1{i}(:,1),'w','linewidth',2)
%     plot(time{i},meanfuns_d1{i}(:,2),'w','linewidth',2)
% end
% xlabel('Time (s)','fontsize',16)
% ylabel('Amplitude (z)','fontsize',16)
% title('Sliding window density estimation','fontsize',16)
% subplot(2,1,2)
% pcolor(temp_t,emissrange_d2,log(emissdist_d2'));shading flat
% hold on
% caxis([-4 -0.5])
% for i = 1:hmm.Nsegs
%     plot(time{i},meanfuns_d2{i}(:,1),'w','linewidth',2)
%     plot(time{i},meanfuns_d2{i}(:,2),'w','linewidth',2)
% end
% xlabel('Time (s)','fontsize',16)
% ylabel('Amplitude (z)','fontsize',16)
% title('Sliding window density estimation','fontsize',16)

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
    for d = 1:p
        temp_emiss(:,d) = emiss(UDS_samples,d)-squeeze(ov_meanfuns(d,:,1))';
    end
    temp_emiss(any(isnan(temp_emiss),2),:) = [];
    obj = gmdistribution.fit(temp_emiss,hmm.K,'Options',options);
    dstate = find(abs(obj.mu(:,1)) == min(abs(obj.mu(:,1)))); %the down state is the state with mean closest to 0
    dstate_var = squeeze(obj.Sigma(:,:,dstate));
    
    %repeat for the up state
    for d = 1:p
        temp_emiss(:,d) = emiss(UDS_samples,d)-squeeze(ov_meanfuns(d,:,2))';
    end
    temp_emiss(any(isnan(temp_emiss),2),:) = [];
    obj = gmdistribution.fit(temp_emiss,hmm.K,'Options',options);
    ustate = find(abs(obj.mu(:,1)) == min(abs(obj.mu(:,1))));
    ustate_var = squeeze(obj.Sigma(:,:,ustate));
end

%% create hmm structure
hmm.A = A; %transition probability matrix
hmm.Pi = Pi; %initial state distribution

if strcmp(meantype,'fixed') %for fixed state means
    hmm.state(1).fixedmean = dstate_fixmean;
    hmm.state(2).fixedmean = ustate_fixmean;
end
if strcmp(meantype,'variable') %for variable state mean functions
    for i = 1:hmm.Nsegs
        for d = 1:p
           hmm.state(1).meanfun{i} = meanfuns{i}(:,:,1);
           hmm.state(2).meanfun{i} = meanfuns{i}(:,:,2);
        end
    end
end
hmm.state(1).var = dstate_var;
hmm.state(2).var = ustate_var;
hmm.meantype = meantype;

