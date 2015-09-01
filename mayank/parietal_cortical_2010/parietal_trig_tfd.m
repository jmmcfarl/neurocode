clear all
close all

addpath('F:\Code\smoothing\software')
addpath('F:\Code\FullBNT-1.0.4\KPMstats\')
addpath('F:\Code\FullBNT-1.0.4\netlab3.3')
addpath('F:\WC_Germany\new_stellate_analysis\')
addpath('F:\WC_Germany\hsmm_state_detection')
addpath('F:\WC_Germany\parietal_cortical_2010\')
addpath('F:\Code\WC_anal\general\')
addpath('F:\Code\wavelet_tools\')

cd F:\WC_Germany\parietal_cortical_2010\
load parietal_cortical_2010
load F:\WC_Germany\parietal_cortical_2010\desynch_times_mp_lf8_2_24_09

% parietal = find_struct_field_vals(sess_data,'region','parietal');
% sess_data = sess_data(parietal);
% desynch_start_times = desynch_start_times(parietal);
% desynch_stop_times = desynch_stop_times(parietal);

%get rid of interneurons
interneurons = find_struct_field_vals(sess_data,'cell_type','interneuron');
sess_data(interneurons) = [];
desynch_times(interneurons) = [];

% %get rid of sessions without lfp
% bad_sessions = [];
% for i = 1:length(sess_data)
%     if ismember(0,sess_data(i).gains)
%         bad_sessions = [bad_sessions i];
%     end
% end
% sess_data(bad_sessions) = [];
% desynch_start_times(bad_sessions) = [];
% desynch_stop_times(bad_sessions) = [];

thomas_el = find_struct_field_vals(sess_data,'thom_elec',1);
sess_data = sess_data(thomas_el);
desynch_times = desynch_times(thomas_el);

parietal = find_struct_field_vals(sess_data,'region','parietal');
frontal = find_struct_field_vals(sess_data,'region','frontal');
prefrontal = find_struct_field_vals(sess_data,'region','prefrontal');

n = length(sess_data);

raw_Fs = 2016;
dsf = 16;
niqf = raw_Fs/2;
Fsd = raw_Fs/dsf;
Fs_bb = raw_Fs/8;
bb_dsf = Fs_bb/Fsd;

[b,a] = butter(2,[0.05/niqf 45/niqf]);
% [b2,a2] = butter(2,[0.05/niqf 40/niqf]);

%look 0.5s after the transition you want and 1.5s before (reversed for down
%transitions)
forwardlag = round(1.5*Fsd);
backlag = round(0.5*Fsd);
coi_extra = round(0.25*Fsd);
lags = -backlag:forwardlag;

min_freq = 4; max_freq = 40; delta_j = 0.1;
k0 = 6; %wavenumber for morlet wavelet
fourier_factor = (4*pi)/(k0 + sqrt(2 + k0^2));
min_scale = 1/max_freq/fourier_factor;
max_scale = 1/min_freq/fourier_factor;
n_scales = round(log2(max_scale/min_scale)/delta_j);

avg_trig_tfd = zeros(n,n_scales+1,length(lags));
trig_lf8 = zeros(n,length(lags));
trig_mp = zeros(n,length(lags));

up_lag = zeros(n,n_scales+1);
modulation = zeros(n,n_scales+1);

n = length(sess_data);

for d = 1:n
    cdir = sess_data(d).directory;
    cdir(1) = 'F';
    cd(cdir)
    pwd
    
    load used_data lf8 lf4 lf5 wcv_minus_spike
    
    lf8_f = filtfilt(b,a,lf8);
    lf8_f = downsample(lf8_f,dsf)/sess_data(d).gains(8);
    lf5_f = filtfilt(b,a,lf5);
    lf5_f = downsample(lf5_f,dsf)/sess_data(d).gains(5);
    if sess_data(d).thom_elec
        lf4_f = filtfilt(b,a,lf4);
        lf4_f = downsample(lf4_f,dsf)/sess_data(d).gains(4);
    end
%     wcv_f = filtfilt(b,a,wcv_minus_spike);
%     wcv_f = zscore(downsample(wcv_f,dsf));

    
    load hsmm_state_seq_seg_lf_pert15   
    load hsmm_state_seq8_seg_lf_pert15
%     load hsmm_state_seq4_lf_pert15
    mp_state_seq =  hsmm_bbstate_seq;
    lfp_state_seq = hsmm_bbstate_seq8;
            
    %compute overall scalogram
    [wav_trans,periods,scales,COI] = wavelet(lf8_f,1/Fsd,0,delta_j,min_scale,n_scales);
    wfreqs = 1./periods;  
    raw_scalogram = abs(wav_trans).^2;
    inv_scales = 1./scales';
    tp = size(raw_scalogram,2);
    raw_scalogram = smoothwavelet(inv_scales(:,ones(1,tp)).*raw_scalogram,1/Fsd,delta_j,scales);
    raw_scalogram = log(raw_scalogram + 1e-5);
    mean_power8 = mean(raw_scalogram,2);  std_power8 = std(raw_scalogram,[],2);

%     [wav_trans,periods,scales,COI] = wavelet(lf5_f,1/Fsd,0,delta_j,min_scale,n_scales);
%     wfreqs = 1./periods;
%     raw_scalogram = abs(wav_trans).^2;
%     raw_scalogram = smoothwavelet(inv_scales(:,ones(1,tp)).*raw_scalogram,1/Fsd,delta_j,scales);
%     raw_scalogram = log(raw_scalogram + 1e-5);
%     mean_power5 = mean(raw_scalogram,2);
%     std_power5 = std(raw_scalogram,[],2);

    if sess_data(d).thom_elec
        [wav_trans,periods,scales,COI] = wavelet(lf4_f,1/Fsd,0,delta_j,min_scale,n_scales);
        wfreqs = 1./periods;
        raw_scalogram = abs(wav_trans).^2;
        raw_scalogram = smoothwavelet(inv_scales(:,ones(1,tp)).*raw_scalogram,1/Fsd,delta_j,scales);
        raw_scalogram = log(raw_scalogram + 1e-5);
        mean_power4 = mean(raw_scalogram,2);
        std_power4 = std(raw_scalogram,[],2);
    end
    
%     [wav_trans,periods,scales,COI] = wavelet(wcv_f,1/Fsd,0,delta_j,min_scale,n_scales);
%     wfreqs = 1./periods;
%     raw_scalogram = abs(wav_trans).^2;        
%     raw_scalogram = log(raw_scalogram + 1e-5);
%     mean_power = mean(raw_scalogram,2);
%     std_power = std(raw_scalogram,[],2);

    %compute used state transition indices
    [up_state_ind,down_state_ind] = get_used_state_trans_seg(lfp_state_seq,Fsd,hmm8.Fs,hmm8.UDS_segs,2*Fsd);
    up_state_ind = round(up_state_ind/bb_dsf);
    down_state_ind = round(down_state_ind/bb_dsf);
    up_state_ind(up_state_ind <1) = 1;
    down_state_ind(down_state_ind < 1) = 1;
    
    %% cycle through all used state transitions
    used_trans_ind = up_state_ind;
    
    n_trans = size(used_trans_ind,1);
%     trig_tfd = nan(n_trans,n_scales+1,length(lags));
    trig_tfd8 = nan(n_trans,n_scales+1,length(lags));
%     trig_tfd5 = nan(n_trans,n_scales+1,length(lags));
%     trig_coh8 = nan(n_trans,n_scales+1,length(lags));
%     trig_coh5 = nan(n_trans,n_scales+1,length(lags));   
    if sess_data(d).thom_elec
        trig_tfd4 = nan(n_trans,n_scales+1,length(lags));
%         trig_coh4 = nan(n_trans,n_scales+1,length(lags));
        trig_coh84 = nan(n_trans,n_scales+1,length(lags));
    end
    for i = 1:n_trans
        %take a segment around the transition and add a buffer to
        %accomodate the COI effects
        cur_seg = used_trans_ind(i,1)-backlag-coi_extra:used_trans_ind(i,1)+forwardlag+coi_extra;
        [wav_trans8,periods,scales,COI] = wavelet(lf8_f(cur_seg),1/Fsd,0,delta_j,min_scale,n_scales);
        wfreqs = 1./periods;
        raw_scalogram8 = abs(wav_trans8).^2;
%         [wav_trans5,periods,scales,COI] = wavelet(lf5_f(cur_seg),1/Fsd,0,delta_j,min_scale,n_scales);
%         raw_scalogram5 = abs(wav_trans5).^2;
        if sess_data(d).thom_elec
              %%% if you want to check against random state transitions
                rand_trans = ceil(rand*n_trans);
                cur_seg = used_trans_ind(rand_trans,1)-backlag-coi_extra:used_trans_ind(rand_trans,1)+forwardlag+coi_extra;
                %%% 
                
            [wav_trans4,periods,scales,COI] = wavelet(lf4_f(cur_seg),1/Fsd,0,delta_j,min_scale,n_scales);
            raw_scalogram4 = abs(wav_trans4).^2;
        end
%         if you want to check against random state transitions
%                 rand_trans = ceil(rand*n_trans);
%                 cur_seg = used_trans_ind(rand_trans,1)-backlag-coi_extra:used_trans_ind(rand_trans,1)+forwardlag+coi_extra;
%         [wav_trans,periods,scales,COI] = wavelet(wcv_f(cur_seg),1/Fsd,0,delta_j,min_scale,n_scales);
%         wfreqs = 1./periods;
%         raw_scalogram = abs(wav_trans).^2;
        inv_scales = 1./scales';
        tp = size(raw_scalogram8,2);
        sm_scalogram8 = smoothwavelet(inv_scales(:,ones(1,tp)).*raw_scalogram8,1/Fsd,delta_j,scales);
%         sm_scalogram5 = smoothwavelet(inv_scales(:,ones(1,tp)).*raw_scalogram5,1/Fsd,delta_j,scales);
        if sess_data(d).thom_elec
            sm_scalogram4 = smoothwavelet(inv_scales(:,ones(1,tp)).*raw_scalogram4,1/Fsd,delta_j,scales);
        end
%         sm_scalogram = smoothwavelet(inv_scales(:,ones(1,tp)).*raw_scalogram,1/Fsd,delta_j,scales);
        
%         %compute coheragrams
%         cross_scalogram = wav_trans.*conj(wav_trans8);
%         sm_cross_scalogram = smoothwavelet(inv_scales(:,ones(1,tp)).*cross_scalogram,1/Fsd,delta_j,scales);
%         wav_cohere8 = abs(sm_cross_scalogram).^2./(sm_scalogram.*sm_scalogram8);
%         cross_scalogram = wav_trans.*conj(wav_trans5);
%         sm_cross_scalogram = smoothwavelet(inv_scales(:,ones(1,tp)).*cross_scalogram,1/Fsd,delta_j,scales);
%         wav_cohere5 = abs(sm_cross_scalogram).^2./(sm_scalogram.*sm_scalogram5);
        if sess_data(d).thom_elec
            cross_scalogram = wav_trans8.*conj(wav_trans4);
            sm_cross_scalogram = smoothwavelet(inv_scales(:,ones(1,tp)).*cross_scalogram,1/Fsd,delta_j,scales);
            wav_cohere84 = abs(sm_cross_scalogram).^2./(sm_scalogram8.*sm_scalogram4);
%             cross_scalogram = wav_trans.*conj(wav_trans4);
%             sm_cross_scalogram = smoothwavelet(inv_scales(:,ones(1,tp)).*cross_scalogram,1/Fsd,delta_j,scales);
%             wav_cohere4 = abs(sm_cross_scalogram).^2./(sm_scalogram.*sm_scalogram4);
        end
        
        %get rid of COI buffer
%         wav_cohere8 = wav_cohere8(:,coi_extra+1:end-coi_extra);
%         wav_cohere5 = wav_cohere5(:,coi_extra+1:end-coi_extra);
%         used_scalogram = raw_scalogram(:,coi_extra+1:end-coi_extra);
        used_scalogram8 = sm_scalogram8(:,coi_extra+1:end-coi_extra);
%         used_scalogram5 = sm_scalogram5(:,coi_extra+1:end-coi_extra);
        if sess_data(d).thom_elec
            used_scalogram4 = sm_scalogram4(:,coi_extra+1:end-coi_extra);
%             wav_cohere4 = wav_cohere4(:,coi_extra+1:end-coi_extra);
            wav_cohere84 = wav_cohere84(:,coi_extra+1:end-coi_extra);
        end
        
        %for up transitions
        next_trans = used_trans_ind(i,2)-used_trans_ind(i,1);
        
        %%% for random transitions
        rand_next_trans = used_trans_ind(rand_trans,2)-used_trans_ind(rand_trans,1);
        next_trans = min(next_trans,rand_next_trans);
        %%%
        
        if next_trans < forwardlag
%             used_scalogram(:,backlag+next_trans+1:end) = nan;
            used_scalogram8(:,backlag+next_trans+1:end) = nan;
%             used_scalogram5(:,backlag+next_trans+1:end) = nan;
%             wav_cohere8(:,backlag+next_trans+1:end) = nan;
%             wav_cohere5(:,backlag+next_trans+1:end) = nan;
            if sess_data(d).thom_elec
                used_scalogram4(:,backlag+next_trans+1:end) = nan;
%                 wav_cohere4(:,backlag+next_trans+1:end) = nan;              
                wav_cohere84(:,backlag+next_trans+1:end) = nan;              
            end
            %             used_lf8seg(backlag+next_trans+1:end) = nan;
        end
        
%         %% for down transitions
%         if i > 1
%             prev_trans = used_trans_ind(i,1)-used_trans_ind(i-1,2);
%             if prev_trans < backlag
%                used_scalogram(:,1:(backlag - prev_trans)) = nan; 
%                used_lf8seg(1:(backlag-prev_trans)) = nan;
%             end
%         end
        
%         trig_tfd(i,:,:) = used_scalogram;
        trig_tfd8(i,:,:) = used_scalogram8;
%         trig_tfd5(i,:,:) = used_scalogram8;
%         trig_coh8(i,:,:) = wav_cohere8;
%         trig_coh5(i,:,:) = wav_cohere5;
        if sess_data(d).thom_elec
            trig_tfd4(i,:,:) = used_scalogram4;
%             trig_coh4(i,:,:) = wav_cohere4;           
            trig_coh84(i,:,:) = wav_cohere84;           
        end
%         trig_lf8(i,:) = used_lf8seg;
    end
%     trig_tfd = log(trig_tfd + 1e-5);
%     trig_tfd = shiftdim(trig_tfd,1);
%     trig_tfd = trig_tfd - repmat(mean_power,[1,length(lags),n_trans]);
%     trig_tfd = trig_tfd./repmat(std_power,[1,length(lags),n_trans]);
%     trig_tfd = shiftdim(trig_tfd,2);
    trig_tfd8 = log(trig_tfd8 + 1e-5);
    trig_tfd8 = shiftdim(trig_tfd8,1);
    trig_tfd8 = trig_tfd8 - repmat(mean_power8,[1,length(lags),n_trans]);
    trig_tfd8 = trig_tfd8./repmat(std_power8,[1,length(lags),n_trans]);
    trig_tfd8 = shiftdim(trig_tfd8,2);
%     trig_tfd5 = log(trig_tfd5 + 1e-5);
%     trig_tfd5 = shiftdim(trig_tfd5,1);
%     trig_tfd5 = trig_tfd5 - repmat(mean_power5,[1,length(lags),n_trans]);
%     trig_tfd5 = trig_tfd5./repmat(std_power5,[1,length(lags),n_trans]);
%     trig_tfd5 = shiftdim(trig_tfd5,2);
%     trig_coh8 = atanh(trig_coh8);
%     trig_coh5 = atanh(trig_coh5);
    if sess_data(d).thom_elec
        trig_tfd4 = log(trig_tfd4 + 1e-5);
        trig_tfd4 = shiftdim(trig_tfd4,1);
        trig_tfd4 = trig_tfd4 - repmat(mean_power4,[1,length(lags),n_trans]);
        trig_tfd4 = trig_tfd4./repmat(std_power4,[1,length(lags),n_trans]);
        trig_tfd4 = shiftdim(trig_tfd4,2);
%         trig_coh4 = atanh(trig_coh4);
        trig_coh84 = atanh(trig_coh84);
    end
    
%     mean_trig_tfd = squeeze(nanmean(trig_tfd));
%     mean_trig_tfd8 = squeeze(nanmean(trig_tfd8));
%     mean_trig_tfd5 = squeeze(nanmean(trig_tfd5));
%     mean_trig_coh8 = squeeze(nanmean(trig_coh8));
%     mean_trig_coh5 = squeeze(nanmean(trig_coh8));
%     if sess_data(d).thom_elec
%         mean_trig_tfd4 = squeeze(nanmean(trig_tfd4));
%         mean_trig_coh4 = squeeze(nanmean(trig_coh4));
%         mean_trig_coh84 = squeeze(nanmean(trig_coh84));
%     end
%     betafit = zeros(4,length(wfreqs));
%     fit_lags = find(lags/Fsd > -0.4 & lags/Fsd < 0.2);
%     for j = 1:length(wfreqs)
%         cur_data = mean_trig_tfd(j,fit_lags);
%         init_guess = [2*std(cur_data) -50 20 0];
%         [betafit(:,j), resid] = nlinfit(lags(fit_lags),cur_data,@my_sigmoid,init_guess);
%         r_sq(j) = 1-sum(resid.^2)/sum((cur_data-mean(cur_data)).^2);
%     end
%     rel_lags = betafit(2,:)/Fsd-0.1;
%     sig_mod = betafit(1,:);
%     sig_mod(r_sq < 0.7) = nan;
%     rel_lags(r_sq < 0.7) = nan;
%     freq_up_lag(d,:) = rel_lags;
%     freq_modulation(d,:) = sig_mod;   
%     clear betafit
    
%     avg_trig_tfd(d,:,:) = squeeze(nanmean(trig_tfd));
    avg_trig_tfd8(d,:,:) = squeeze(nanmean(trig_tfd8));
%     avg_trig_tfd5(d,:,:) = squeeze(nanmean(trig_tfd5));
%     avg_trig_coh8(d,:,:) = squeeze(nanmean(trig_coh8));
%     avg_trig_coh5(d,:,:) = squeeze(nanmean(trig_coh5));
    if sess_data(d).thom_elec
        avg_trig_tfd4(d,:,:) = squeeze(nanmean(trig_tfd4));
        avg_trig_coh84(d,:,:) = squeeze(nanmean(trig_coh84));
%         avg_trig_coh4(d,:,:) = squeeze(nanmean(trig_coh4));
    end
end
%%

cd F:\WC_Germany\parietal_cortical_2010
save up8_trig_wvlet_specgram_thom_rand avg_trig* lags Fsd wfreqs 

%%
n = length(sess_data);
parietal = find_struct_field_vals(sess_data,'region','parietal');
frontal = find_struct_field_vals(sess_data,'region','frontal');
prefrontal = find_struct_field_vals(sess_data,'region','prefrontal');

sup = find_struct_field_vals(sess_data,'layer','23');
deep = setdiff(1:n,sup);
pfc = setdiff(1:n,parietal);

sup_par = intersect(sup,parietal);
deep_par = intersect(deep,parietal);
sup_pfc = intersect(sup,pfc);
deep_pfc = intersect(deep,pfc);
figure
pcolor(lags/Fsd,wfreqs,squeeze(nanmean(avg_trig_coh84)));shading flat

cur_set = sup_par;


% figure
% errorbar(wfreqs,nanmean(freq_up_lag),nanstd(freq_up_lag)/sqrt(n))

% figure
% errorbar(wfreqs,nanmean(freq_modulation),nanstd(freq_modulation)/sqrt(n))