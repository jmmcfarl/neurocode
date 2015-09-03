clear all
close all

cd G:\WC_Germany\persistent_revised\
load pers_revised_dir_2010

%%

raw_Fs = 2016;
Fs_desired = 40; %(in Hz)
dsf = round(raw_Fs/Fs_desired);
Fs = raw_Fs/dsf;
niqf = raw_Fs/2;
lf_lcf = 0.05; %low cut-off for the low-freq filter
lf_hcf = 2; %high cut-off for the low-freq filter
hf_lcf = 20; %low cut-off for high freq filter
hf_hcf = 100; %high cut-off for high freq filter
hf_smooth = 0.05; %std dev of gaussian smoothing for hf RMS
windowSize = 60;
windowSlide = 10;

%% initializations
%compute filter coefficients
[b_lf,a_lf] = butter(2,[lf_lcf/niqf lf_hcf/niqf]);
[b_hf,a_hf] = butter(2,[hf_lcf/niqf hf_hcf/niqf]);
cll = [-0.5 2.4];
cl = [0 200];
for d = 1:length(dir_array)
    cd(dir_array{d})
    load used_data wcv_minus_spike lf8
    
    %compute low-freq observations
%     obs_lf = filtfilt(b_lf,a_lf,wcv_minus_spike);
%     obs_lf = zscore(downsample(obs_lf,dsf));
%     obs_lf = obs_lf(:);
%     %compute high-freq RMS observations
%     obs_hf = filtfilt(b_hf,a_hf,wcv_minus_spike);
%     obs_hf = sqrt(jmm_smooth_1d_cor(obs_hf.^2,round(hf_smooth*raw_Fs))); %smoothed RMS power
%     obs_hf = log(obs_hf-min(obs_hf)+1); %log transform the rms power to normalize
% %     obs_hf = log(obs_hf); %log transform the rms power to normalize
%     obs_hf = zscore(downsample(obs_hf,dsf));
%     obs_hf = obs_hf(:);
%     time = (1:length(obs_lf))/Fs;
    
%     [state_means,state_t,obs_dist,obsrange,uni_times] = ...
%         locate_state_means(obs_lf,windowSize,windowSlide,Fs);
%     state_t  = [0 state_t max(time)];
%     state_means = [state_means(:,1) state_means state_means(:,end)];
%     lf_meanfuns = interp1(state_t,state_means',time);
%     
%     obs_lf = obs_lf - lf_meanfuns(:,2);
    
%     mix=gmm(2,2,'full');
%     gmm_options(3) = 1e-10; %tolerance
%     gmm_options(5) = 1;%reset cov if singular values
%     gmm_options(14) = 200; %max iterations
%     mix = gmmem(mix,[obs_lf obs_hf],gmm_options);
%     hist3_plot(obs_lf,obs_hf,[40 40],[],[],cll,1), hold on
%     gaussplot(mix.centres(1,:),sqrt(mix.covars(:,:,1)),2,'k')
%     gaussplot(mix.centres(2,:),sqrt(mix.covars(:,:,2)),2,'w')
%     ylabel('HF Power (z)','FontSize',14)
%     xlabel('LF Amp (z)','FontSize',14)
%     t_name = ['G:\WC_Germany\persistent_2010\joint_dist_lf_hf\mp_log_' f_names{d}];
%     print('-dpng',t_name), close
%       hist3_plot(obs_lf,obs_hf,[40 40],[],[],cl,0), hold on
%     gaussplot(mix.centres(1,:),sqrt(mix.covars(:,:,1)),2,'k')
%     gaussplot(mix.centres(2,:),sqrt(mix.covars(:,:,2)),2,'w')
%     ylabel('HF Power (z)','FontSize',14)
%     xlabel('LF Amp (z)','FontSize',14)
%     t_name = ['G:\WC_Germany\persistent_2010\joint_dist_lf_hf\mp_' f_names{d}];
%     print('-dpng',t_name), close  
    
    
    obs_lf = filtfilt(b_lf,a_lf,lf8);
    obs_lf = zscore(downsample(obs_lf,dsf));
    obs_lf = obs_lf(:);
    %compute high-freq RMS observations
    obs_hf = filtfilt(b_hf,a_hf,lf8);
    obs_hf = sqrt(jmm_smooth_1d_cor(obs_hf.^2,round(hf_smooth*raw_Fs))); %smoothed RMS power
    obs_hf = log(obs_hf-min(obs_hf)+1); %log transform the rms power to normalize
%     obs_hf = log(obs_hf); %log transform the rms power to normalize
    obs_hf = zscore(downsample(obs_hf,dsf));
    obs_hf = obs_hf(:);
    
%         [state_means,state_t,obs_dist,obsrange,uni_times] = ...
%         locate_state_means(obs_lf,windowSize,windowSlide,Fs);
%     state_t  = [0 state_t max(time)];
%     state_means = [state_means(:,1) state_means state_means(:,end)];
%     lf_meanfuns = interp1(state_t,state_means',time);
%     
%     obs_lf = obs_lf - lf_meanfuns(:,2);

    
%     mix=gmm(2,2,'full');
%     gmm_options(3) = 1e-10; %tolerance
%     gmm_options(5) = 1;%reset cov if singular values
%     gmm_options(14) = 200; %max iterations
%     mix = gmmem(mix,[obs_lf obs_hf],gmm_options);
%     hist3_plot(obs_lf,obs_hf,[40 40],[],[],cl,0), hold on
%     gaussplot(mix.centres(1,:),sqrt(mix.covars(:,:,1)),2,'k')
%     gaussplot(mix.centres(2,:),sqrt(mix.covars(:,:,2)),2,'w')
%     ylabel('HF Power (z)','FontSize',14)
%     xlabel('LF Amp (z)','FontSize',14)
%     t_name = ['G:\WC_Germany\persistent_2010\joint_dist_lf_hf\lf8_' f_names{d}];
%     print('-dpng',t_name), close
%      hist3_plot(obs_lf,obs_hf,[40 40],[],[],cll,1), hold on
%     gaussplot(mix.centres(1,:),sqrt(mix.covars(:,:,1)),2,'k')
%     gaussplot(mix.centres(2,:),sqrt(mix.covars(:,:,2)),2,'w')
%     ylabel('HF Power (z)','FontSize',14)
%     xlabel('LF Amp (z)','FontSize',14)
%     t_name = ['G:\WC_Germany\persistent_2010\joint_dist_lf_hf\lf8_log_' f_names{d}];
%     print('-dpng',t_name), close
   
end