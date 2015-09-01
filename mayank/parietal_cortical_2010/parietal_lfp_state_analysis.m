clear all
close all
cd G:\WC_Germany\parietal_cortical_2010
load parietal_cortical_2010

raw_fs = 2016;
niqf = raw_fs/2;
dsf = 8;
fsd = raw_fs/dsf;
lcf = 0.1;
hcf = 120;
[b,a] = butter(2,[lcf/niqf hcf/niqf]);
[b2,a2] = butter(2,[lcf/niqf 2/niqf]);
[b3,a3] = butter(2,[20/niqf 100/niqf]);
for d = 1:length(sess_data)
    
    cd(sess_data(d).directory)
    pwd
    
    load used_data lf8
    load hsmm_state_seq
   
    lf8_f = filtfilt(b,a,lf8);
    lf8_f = zscore(downsample(lf8_f,dsf));
    lf8_lf = filtfilt(b2,a2,lf8);
    lf8_lf = zscore(downsample(lf8_lf,dsf));
    lf8_hf = filtfilt(b3,a3,lf8);
    lf8_hf = jmm_smooth_1d_cor(lf8_hf.^2,round(0.05*raw_fs));
    lf8_hf = sqrt(lf8_hf);
    lf8_hf = log(lf8_hf - min(lf8_hf)+1);
    lf8_hf = zscore(downsample(lf8_hf,dsf));
    state_t = (1:length(hsmm_state_seqn))/Fs;
    t = (1:length(lf8_f))/fsd;
    
    state_t = [0 state_t];
    state_seq = [hsmm_state_seqn(1);hsmm_state_seqn];
    
    state_seq = interp1(state_t,state_seq,t);
    state_seq = logical(state_seq-1);
    
    
    min_freq = 3; max_freq = 100; delta_j = 0.1;
    k0 = 6; %wavenumber for morlet wavelet
    fourier_factor = (4*pi)/(k0 + sqrt(2 + k0^2));
    min_scale = 1/max_freq/fourier_factor;
    max_scale = 1/min_freq/fourier_factor;
    n_scales = round(log2(max_scale/min_scale)/delta_j);
    [wav_trans,periods,scales,COI] = wavelet(lf8_f,1/fsd,0,delta_j,min_scale,n_scales);
    wfreqs = 1./periods;
    inv_scales = 1./scales';
    raw_scalogram = abs(wav_trans).^2;
    n = size(raw_scalogram,2);

    sm_scalogram = smoothwavelet(inv_scales(:,ones(1,n)).*raw_scalogram,1/fsd,delta_j,scales,0.6);
    sm_scalogram_d = downsample(log(sm_scalogram),4)';
    sm_scalogram_d = zscore(sm_scalogram_d);
    
    fd = downsample(wfreqs,4);
    [class,err,POSTERIOR,logp,coeff] = classify([],sm_scalogram_d,state_seq);
    L = coeff(1,2).linear
    temp = sm_scalogram_d*L;
%     figure, ksdensity(temp)
%     figure
    plot(fd,L,'b')
    
    band1 = find(wfreqs < 10);
    band2 = find(wfreqs > 10 & wfreqs < 40);
    band3 = find(wfreqs > 40);
    pow_band1 = zscore(mean(log(raw_scalogram(band1,:))));
    pow_band2 = zscore(mean(log(raw_scalogram(band2,:))));
    pow_band3 = zscore(mean(log(raw_scalogram(band3,:))));
    
    dsf2 = 4;
    td = downsample(t,dsf2);
    state_seq_d = downsample(state_seq,dsf2);
    sm_scalogram_d = downsample(sm_scalogram,dsf2);
    
    raw_scalogram = zscore(raw_scalogram');
    
    raw_scalogram = zscore(log(raw_scalogram - repmat(min(raw_scalogram),length(lf8_f),1)+0.1));
    [class,err,POSTERIOR,logp,coeff] = classify([],sm_scalogram_d,state_seq_d);
    
end