
clear all
close all
cd G:\WC_Germany\parietal_cortical_2010\
load parietal_cortical_2010
addpath('G:\WC_Germany\parietal_cortical_2010\')
addpath('G:\WC_Germany\persistent_2010\')

raw_Fs = 2016;
niqf = raw_Fs/2;
dsf = 8;
Fsd = raw_Fs/dsf;

lcf = 0.05;
hcf = 40;

lcf2 = 2;
hcf2 = 20;

lcf3 = 2;
hcf3 = 5;
% niqf2 = Fsd/2;
% [b,a] = butter(2,[lcf3 hcf3]/niqf2);

winSize = 10;

%%
d = find_struct_field_vals(sess_data,'name','2007-04-17_A')

cdir = sess_data(d).directory;
cdir(1) = 'G';
disp(sprintf('session %d',d))
cd(cdir);

load used_data lf8 wcv_minus_spike 

% load hsmm_state_seq_seg_lf_indivdesynch_dp
% mp_state_seq = hsmm_bbstate_seq;
% load hsmm_state_seq8_seg_lf_indivdesynch_dp
% lfp_state_seq = hsmm_state_seq8;
% 
% used_state_seq = mp_state_seq;

lf8_f = get_lf_features(lf8,raw_Fs,Fsd,[lcf 20]);
lf8_f2 = get_lf_features_caus(lf8,raw_Fs,Fsd,[lcf2 hcf2]);
lf8_f3 = get_lf_features(lf8,raw_Fs,Fsd,[lcf3 hcf3]);
% lf8_f3a = get_lf_features_acaus(lf8,raw_Fs,Fsd,[lcf3 hcf3]);
% lf8_f3c = get_lf_features_caus(lf8,raw_Fs,Fsd,[lcf3 hcf3]);
t = (1:length(lf8_f))/Fsd;
wcv_f = get_lf_features(wcv_minus_spike,raw_Fs,Fsd,[lcf hcf]);
wcv_f2 = get_lf_features_caus(wcv_minus_spike,raw_Fs,Fsd,[lcf2 40]);



dt = round(winSize*Fsd);
ct = 1;
et = ct + dt;
T = round(t(end)*Fsd);
while ct < T
    plot(t(ct:et),lf8_f(ct:et),'r','linewidth',1), hold on
    plot(t(ct:et),lf8_f2(ct:et)/3-2.5,'r','linewidth',1)
    plot(t(ct:et),wcv_f(ct:et),'b','linewidth',1)
    plot(t(ct:et),wcv_f2(ct:et)/3-4,'b','linewidth',1)
    xlim([t(ct) t(et)])
    ylim([-5 3])
    ct = ct + dt;
    et = ct + dt;
    pause
    clf
end

%%
% min_freq = 1.5; max_freq = 5; delta_j = 0.1;
% k0 = 6; %wavenumber for morlet wavelet
% fourier_factor = (4*pi)/(k0 + sqrt(2 + k0^2));
% min_scale = 1/max_freq/fourier_factor;
% max_scale = 1/min_freq/fourier_factor;
% n_scales = round(log2(max_scale/min_scale)/delta_j);
% [wav_trans,periods,scales,COI] = wavelet(lf8_f,1/Fsd,0,delta_j,min_scale,n_scales);
% wfreqs = 1./periods;
% raw_scalogram = abs(wav_trans).^2;
% raw_phase = angle(wav_trans);
% [peak_pow,peak_loc] = max(raw_scalogram);
% peak_freq = wfreqs(peak_loc);
% peak_phase = raw_phase(peak_loc,1:length(t));
