function [hmm] = initialize_uds_hsmm_seg_dualvar(emiss,Fs,meanparams,UDS_segs,drive_letter)

%%initialize hidden semi-markov model.  meanparams determines whether any
%%of the observations is nonstationary.  currently the function only
%%supports a maximum of one nonstationary emission parameter.  this must be
%%the first column of emiss.  Currently only allows 2 hidden states
hmm.drive_letter = drive_letter;
addpath(strcat(hmm.drive_letter,':\WC_Germany\hmm_state_detect'))

hmm.p=size(emiss,2);
hmm.T=size(emiss,1);
hmm.Fs = Fs;
hmm.K = 2;
hmm.UDS_segs = UDS_segs;
hmm.Nsegs = size(UDS_segs,1);

%these are the indices for the segments of data with UDS
UDS_samples = [];
for i = 1:hmm.Nsegs
    UDS_seg_inds{i} = [UDS_segs(i,1):UDS_segs(i,2)];
    UDS_samples = [UDS_samples UDS_segs(i,1):UDS_segs(i,2)];
    time{i} = (UDS_segs(i,1):UDS_segs(i,2))/Fs;
end
hmm.uds_dur = length(UDS_samples)/hmm.Fs;

meantype = meanparams.meantype;
if strcmp(meantype,'variable')
    hmm.windowSize = meanparams.windowSize;
    hmm.windowSlide = meanparams.windowSlide;
    slid_avg_hcf = 0.02;
    niqf = Fs/2;
    [slid_b,slid_a] = butter(2,slid_avg_hcf/niqf,'low');
end

%% choose random starting values for A and Pi from dir priors
pi_alpha = 1;
Pi_alpha = ones(1,hmm.K)*pi_alpha;
for i = 1:hmm.Nsegs
    Pi(i,:) = dirichletrnd(Pi_alpha);
end
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
    ov_meanfuns2 = [];
    for i = 1:hmm.Nsegs
        [state_means,state_t] = ...
            locate_state_means(emiss(UDS_seg_inds{i},1),hmm.windowSize,hmm.windowSlide,Fs);
        state_t = state_t + UDS_segs(i,1)/Fs;
        state_t  = [time{i}(1) state_t (time{i}(end)+1)];
        state_means = [state_means(:,1) state_means state_means(:,end)];
        meanfuns{i} = interp1(state_t,state_means',time{i});
        for k = 1:hmm.K
            meanfuns{i}(:,k) = filtfilt(slid_b,slid_a,meanfuns{i}(:,k));
        end
        ov_meanfuns = [ov_meanfuns; meanfuns{i}];

        [state_means,state_t] = ...
            locate_state_means(emiss(UDS_seg_inds{i},2),hmm.windowSize,hmm.windowSlide,Fs);
        state_t = state_t + UDS_segs(i,1)/Fs;
        state_t  = [time{i}(1) state_t (time{i}(end)+1)];
        state_means = [state_means(:,1) state_means state_means(:,end)];
        meanfuns2{i} = interp1(state_t,state_means',time{i});
        for k = 1:hmm.K
            meanfuns2{i}(:,k) = filtfilt(slid_b,slid_a,meanfuns2{i}(:,k));
        end
        ov_meanfuns2 = [ov_meanfuns2; meanfuns2{i}];
    end
end

%% GMM initializations
options = statset('MaxIter',200,'TolFun',1e-8);

%for difference from down state meanfunction
temp_emiss = emiss(UDS_samples,:);
temp_emiss(:,1) = emiss(UDS_samples,1)-ov_meanfuns(:,1);
temp_emiss(:,2) = emiss(UDS_samples,2)-ov_meanfuns2(:,1);
obj_d = gmdistribution.fit(temp_emiss,hmm.K,'Options',options);
dstate = find(abs(obj_d.mu(:,1))==min(abs(obj_d.mu(:,1)))); %find mixture component with 0 mean
dstate_var = squeeze(obj_d.Sigma(:,:,dstate));

temp_emiss = emiss(UDS_samples,:);
temp_emiss(:,1) = emiss(UDS_samples,1)-ov_meanfuns(:,2);
temp_emiss(:,2) = emiss(UDS_samples,2)-ov_meanfuns2(:,2);
obj_u = gmdistribution.fit(temp_emiss,hmm.K,'Options',options);
ustate = find(abs(obj_u.mu(:,1))==min(abs(obj_u.mu(:,1)))); %find mixture component with 0 mean
ustate_var = squeeze(obj_u.Sigma(:,:,ustate));

%% create hmm structure
hmm.A = A;
hmm.Pi = Pi;

if strcmp(meantype,'variable')
    for i = 1:hmm.Nsegs
        hmm.state(1).meanfun{i} = meanfuns{i}(:,1);
        hmm.state(2).meanfun{i} = meanfuns{i}(:,2);
        hmm.state(1).meanfun2{i} = meanfuns2{i}(:,1);
        hmm.state(2).meanfun2{i} = meanfuns2{i}(:,2);
    end
end
hmm.state(1).var = dstate_var;
hmm.state(2).var = ustate_var;

hmm.covtype = 'full'; %default

hmm.state(1).dur_type = 'inv_gauss';
hmm.state(2).dur_type = 'gamma';
   
hmm.meantype = meantype;

