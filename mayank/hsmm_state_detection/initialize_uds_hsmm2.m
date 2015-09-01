function [hmm] = initialize_uds_hsmm(emiss,Fs,meanparams,UDS_segs)

%%initialize hidden semi-markov model.  meanparams determines whether any
%%of the observations is nonstationary.  currently the function only
%%supports a maximum of one nonstationary emission parameter.  this must be
%%the first column of emiss.  Currently only allows 2 hidden states
hmm.drive_letter = 'G';

addpath(strcat(hmm.drive_letter,':\Code\fullBNT-1.0.4\kPMstats\'))
addpath(strcat(hmm.drive_letter,':\Code\fullBNT-1.0.4\netlab3.3\'))
addpath(strcat(hmm.drive_letter,':\WC_Germany\hmm_state_detect'))

hmm.p=size(emiss,2);
hmm.T=size(emiss,1);
hmm.Fs = Fs;
hmm.K = 2;
hmm.UDS_segs = UDS_segs;
hmm.Nsegs = size(UDS_segs,1);

UDS_samples = [];
for i = 1:hmm.Nsegs
    UDS_seg_inds{i} = [UDS_samples UDS_segs(i,1):UDS_segs(i,2)];
    UDS_samples = [UDS_samples UDS_segs(i,1):UDS_segs(i,2)];
    time{i} = (UDS_segs(i,1):UDS_segs(i,2))/Fs;
end

meantype = meanparams.meantype;
if strcmp(meantype,'variable')
    hmm.windowSize = meanparams.windowSize;
    hmm.windowSlide = meanparams.windowSlide;
    hmm.hcf = meanparams.state_mean_hcf;
end

%% choose random starting values for A and Pi from dir priors
pi_alpha = 1;
Pi_alpha = ones(1,hmm.K)*pi_alpha;
Pi = dirichletrnd(Pi_alpha);
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
        state_t  = [time{i}(1) state_t time{i}(end)];
        state_means = [state_means(:,1) state_means state_means(:,end)];
        meanfuns{i} = interp1(state_t,state_means',time{i});
        ov_meanfuns = [ov_meanfuns meanfuns{i}];
    end
end

%% GMM initializations
mix=gmm(hmm.p,hmm.K,'full');
gmm_options(3) = 1e-10; %tolerance
gmm_options(5) = 1;%reset cov if singular values
gmm_options(14) = 100; %max iterations

if strcmp(meantype,'fixed')
    mix=gmm(hmm.p,hmm.K,'full');
    gmm_options(3) = 1e-10; %tolerance
    gmm_options(5) = 1;%reset cov if singular values
    gmm_options(14) = 100; %max iterations
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
    mix_d = gmmem(mix, temp_emiss, gmm_options);
    dstate = find(abs(mix_d.centres(:,1))==min(abs(mix_d.centres(:,1)))); %find mixture component with 0 mean
    dstate_var = squeeze(mix_d.covars(:,:,dstate));
    if hmm.p > 1
        dstate_fixmean = mix_d.centres(dstate,2:end);
    end
    
    temp_emiss = emiss(UDS_samples,:);
    temp_emiss(:,1) = emiss(UDS_samples,1)-ov_meanfuns(:,2);
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
    hmm.state(1).fixmean = dstate_fixmean;
    hmm.state(2).fixmean = ustate_fixmean;
end
if strcmp(meantype,'variable')
    for i = 1:Nsegs
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

