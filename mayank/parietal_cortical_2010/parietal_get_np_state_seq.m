function [smooth_state_seq,raw_state_seq,Fs] = parietal_get_np_state_seq(raw_data,raw_Fs,f_names)

addpath('G:\Code\smoothing\software')
addpath('G:\Code\FullBNT-1.0.4\KPMstats\')
addpath('G:\Code\FullBNT-1.0.4\netlab3.3')
addpath('G:\WC_Germany\new_stellate_analysis\')
addpath('G:\WC_Germany\hsmm_state_detection')

Fs_desired = 50; %(in Hz)
dsf = round(raw_Fs/Fs_desired);
Fs = raw_Fs/dsf;
niqf = raw_Fs/2;
hf_lcf = 20; %low cut-off for high freq filter
hf_hcf = 100; %high cut-off for high freq filter
hf_smooth = 0.075; %default 0.05 std dev of gaussian smoothing for hf RMS
num_states = 2; %number of hidden states (must be 2 in current implementation)

min_state_dur = round(Fs*0.05);

%% initializations
%compute filter coefficients
[b_hf,a_hf] = butter(2,[hf_lcf/niqf hf_hcf/niqf]);

obs_hf = filtfilt(b_hf,a_hf,raw_data);
obs_hf = sqrt(jmm_smooth_1d_cor(obs_hf.^2,round(hf_smooth*raw_Fs))); %smoothed RMS power
obs_hf = zscore(obs_hf);
obs_hf = log(obs_hf-min(obs_hf)+1); %log transform the rms power to normalize
% obs_hf = log(obs_hf); %log transform the rms power to normalize
obs_hf = zscore(downsample(obs_hf,dsf));
obs_hf = obs_hf(:);

T = length(obs_hf);
t_axis = (1:T)/Fs;

%%
mix=gmm(1,num_states,'full');
gmm_options(3) = 1e-20; %tolerance
gmm_options(5) = 1;%reset cov if singular values
gmm_options(14) = 500; %max iterations

mix = gmmem(mix, obs_hf, gmm_options);
centres = mix.centres;
[dummy,state_order] = sort(centres);
dmean = centres(state_order(1));
dcovar = squeeze(mix.covars(:,:,state_order(1)));
umean = centres(state_order(2));
ucovar = squeeze(mix.covars(:,:,state_order(2)));
priors = mix.priors(state_order);
x = linspace(-4,4,10e3);
g1 = priors(1)/sqrt(2*pi*dcovar)*exp(-(x-dmean).^2/(2*dcovar));
g2 = priors(2)/sqrt(2*pi*ucovar)*exp(-(x-umean).^2/(2*ucovar));
[dummy,dmean_loc] = max(g1);
threshold = dmean_loc + find(g1(dmean_loc:end) < g2(dmean_loc:end),1,'first');
threshold = x(threshold);

ksdensity(obs_hf), hold on
plot(x,g1,'r',x,g2,'k')
yl = ylim;
line([threshold threshold],yl,'Color','k')
t_names = ['G:\WC_Germany\parietal_cortical_2010\test_meanfuns\np_dist_' f_names];
print('-dpng',t_names), close

%% 
up_trans = find(obs_hf(1:end-1) < threshold & obs_hf(2:end) > threshold);
down_trans = find(obs_hf(1:end-1) > threshold & obs_hf(2:end) < threshold);
up_trans(up_trans > down_trans(end)) = [];
down_trans(down_trans < up_trans(1)) = [];
state_seq = ones(size(obs_hf));
for i = 1:length(up_trans)
    state_seq(up_trans(i):down_trans(i)) = 2;
end

[smooth_state_seq,orig_dur,reject_dur,sojourn_times] = thresh_state_smooth(state_seq,Fs,50,50);
raw_state_seq = state_seq;
