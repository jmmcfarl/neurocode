clear all
close all
%%
cd C:\WC_Germany\sven_thomas_combined\
load ./combined_dir_nd.mat	
addpath('C:\Code\WC_anal\general\')
addpath('C:\WC_Germany\Overall_EC\')
addpath('C:\Code\Chronux\spectral_analysis\continuous\')
addpath('C:\WC_Germany\hsmm_state_detection\')
uset = sort([l3mec l3lec l3mec_np l3lec_np]);
all_cells = 1:length(combined_dir);
l3mec = find(ismember(all_cells(uset),l3mec));
l3lec = find(ismember(all_cells(uset),l3lec));
l3mec_np = find(ismember(all_cells(uset),l3mec_np));
l3lec_np = find(ismember(all_cells(uset),l3lec_np));
combined_dir = combined_dir(uset);
hpc_mua = hpc_mua(uset);
hpc_lfp = hpc_lfp(uset);
ctx_lfp = ctx_lfp(uset);

%% DEFINE PARAMETERS
raw_Fs = 2016;
dsf = 8;
Fsd = raw_Fs/dsf;
niqf = raw_Fs/2;
hcf = 40;
lcf = 0.05;
lcf_hf = 15;
hcf_hf = 80;
hcf_sm = 0.025;

mua_rate_sm = 0.05;

params.Fs = Fsd;
params.err = [2 0.05];
params.tapers = [7 13];
W = 0.1;
params.fpass = [0 20];
winlength = 50;
winslide = 50;
movingwin = [winlength winslide];

f_i = linspace(0,20,500);
%%
for d = 1:length(combined_dir)
    d
    cur_dir = combined_dir{d};
    cd(cur_dir)
    load ./used_data lf8 lf7 wcv_minus_spike 
    if ctx_lfp(d) == 7
        lf8 = lf7;
    end
    [desynch_times,desynch_inds,P_lf8,f,t] = locate_desynch_times_individual_v2(lf8);

    [lf8_lf,t_axis] = get_lf_features(lf8,raw_Fs,Fsd,[lcf hcf]);
    wcv_lf = get_lf_features(wcv_minus_spike,raw_Fs,Fsd,[lcf hcf]);
    
    if hpc_lfp(d) == 3
        load ./used_data lf3 lf5
        if ismember(d,old_data_inds)
            lf3 = lf3 + lf5; %redefine LF3 wrt gnd
        end
        hpc_lf = get_lf_features(lf3,raw_Fs,Fsd,[lcf hcf]);
        hpc_hf = get_hf_features(lf3,raw_Fs,Fsd,[lcf_hf hcf_hf],hcf_sm);
    else
        load ./used_data lf2
        hpc_lf = get_lf_features(lf2,raw_Fs,Fsd,[lcf hcf]);
        hpc_hf = get_hf_features(lf2,raw_Fs,Fsd,[lcf_hf hcf_hf],hcf_sm);
    end
    if ~isnan(hpc_mua(d))
        load ./mua_data2
        load ./sync_times.mat
        synctt = synct;
        hpc_mua_times = mua_times{hpc_mua(d)};
        hpc_mua_times(hpc_mua_times < synctt(1) | hpc_mua_times > synctt(end)) = [];
        temp_mua_rate = hist(hpc_mua_times,synct)*raw_Fs;
        temp_mua_rate = jmm_smooth_1d_cor(temp_mua_rate,mua_rate_sm);
        mua_rate = downsample(temp_mua_rate,dsf);
        mua_rate = zscore(mua_rate);
%         [log_mua,offset(d)] = log_transform_sig(mua_rate);
%         mua_rate = zscore(log_mua);

        if length(mua_rate) > length(t_axis)
            mua_rate = mua_rate(1:length(t_axis));
        end
    end
    
    if ~isempty(desynch_times)
        desynch_start = round(interp1(t_axis,1:length(t_axis),desynch_times(:,1)));
        desynch_stop = round(interp1(t_axis,1:length(t_axis),desynch_times(:,2)));
        synch_start = desynch_stop;
        synch_stop = desynch_start;
        if desynch_start(1) ~= 1
           synch_start = [1; synch_start]; 
        end
        if desynch_stop(end) ~= length(t_axis)
            synch_stop = [synch_stop; length(t_axis)];
        end
        seg_inds = [synch_start(:) synch_stop(:)];
    else
        desynch_start = [];
        desynch_stop = [];
        seg_inds = [1 length(t_axis)];
    end
    seg_durs = diff(seg_inds,[],2)/Fsd;
    too_short = find(seg_durs < winlength);
    seg_inds(too_short,:) = [];
    
    [Cw8,Phiw8,Sw8,Pww_88,f,ConfC,dof] = coherencyc_unequal_length_trials_jmm([wcv_lf lf8_lf],movingwin,params,seg_inds);
    [Cw2,Phiw2,Sw2,Pww_22,f,ConfC,dof] = coherencyc_unequal_length_trials_jmm([wcv_lf hpc_lf],movingwin,params,seg_inds);
    [C82,Phi82,S82,P88_22,f,ConfC,dof] = coherencyc_unequal_length_trials_jmm([lf8_lf hpc_lf],movingwin,params,seg_inds);
    [Cw2h,Phiw2h,Sw2h,Pww_22h,f,ConfC,dof] = coherencyc_unequal_length_trials_jmm([wcv_lf hpc_hf],movingwin,params,seg_inds);
    [C82h,Phi82h,S82h,P88_22h,f,ConfC,dof] = coherencyc_unequal_length_trials_jmm([lf8_lf hpc_hf],movingwin,params,seg_inds);
    if ~isnan(hpc_mua(d))
        [Cm8,Phim8,Sm8,Pmm_88,f,ConfC,dof] = coherencyc_unequal_length_trials_jmm([mua_rate' lf8_lf],movingwin,params,seg_inds);
        [Cwm,Phiwm,Swm,Pww_mm,f,ConfC,dof] = coherencyc_unequal_length_trials_jmm([wcv_lf mua_rate'],movingwin,params,seg_inds);
        [Cm2h,Phim2h,Sm2h,Pmm_22h,f,ConfC,dof] = coherencyc_unequal_length_trials_jmm([mua_rate' hpc_hf],movingwin,params,seg_inds);
    end
    
    S8w = conj(Sw8);
    S28 = conj(S82);
    S2w = conj(Sw2);
    S2h8 = conj(S82h);
    S2hw = conj(Sw2h);
    
%     aCw8_cor_temp = atanh(Cw8) - 1/(dof-2);
    aCw8_cor_temp = atanh(Cw8);
    aCw8_cor(d,:) = interp1(f,aCw8_cor_temp,f_i);
%     aCw2_cor_temp = atanh(Cw2) - 1/(dof-2);
     aCw2_cor_temp = atanh(Cw2);
   aCw2_cor(d,:) = interp1(f,aCw2_cor_temp,f_i);
%     aCw2h_cor_temp = atanh(Cw2h) - 1/(dof-2);
     aCw2h_cor_temp = atanh(Cw2h);
   aCw2h_cor(d,:) = interp1(f,aCw2h_cor_temp,f_i);
%     aC82_cor_temp = atanh(C82) - 1/(dof-2);
    aC82_cor_temp = atanh(C82);
    aC82_cor(d,:) = interp1(f,aC82_cor_temp,f_i);
%     aC82h_cor_temp = atanh(C82h) - 1/(dof-2);
    aC82h_cor_temp = atanh(C82h);
    aC82h_cor(d,:) = interp1(f,aC82h_cor_temp,f_i);
    
    lPww(d,:) = interp1(f,log10(Pww_88(:,1)),f_i);
    lP88(d,:) = interp1(f,log10(Pww_88(:,2)),f_i);
    lP22h(d,:) = interp1(f,log10(P88_22h(:,2)),f_i);
    lP22(d,:) = interp1(f,log10(P88_22(:,2)),f_i);
    
    if ~isnan(hpc_mua(d))
        S2hm = conj(Sm2h);
        S8m = conj(Sm8);
        Smw = conj(Swm);
        aCwm_cor_temp = atanh(Cwm) - 1/(dof-2);
        aCwm_cor(d,:) = interp1(f,aCwm_cor_temp,f_i);
        aCm8_cor_temp = atanh(Cm8) - 1/(dof-2);
        aCm8_cor(d,:) = interp1(f,aCm8_cor_temp,f_i);
        aCm2h_cor_temp = atanh(Cm2h) - 1/(dof-2);
        aCm2h_cor(d,:) = interp1(f,aCm2h_cor_temp,f_i);
        lPmm(d,:) = interp1(f,log10(Pww_mm(:,2)),f_i);
    else
       aCwm_cor(d,:) = nan(1,length(f_i));
       aCm8_cor(d,:) = nan(1,length(f_i));
       aCm2h_cor(d,:) = nan(1,length(f_i));
       lPmm(d,:) = nan(1,length(f_i));
    end
    
    num = abs(S2w.*Pww_88(:,2) - S28.*S8w).^2;
    denom = (Pww_22(:,2).*Pww_88(:,2) - abs(S28).^2).*(Pww_22(:,1).*Pww_88(:,2) - abs(Sw8).^2);
    partial_w2(d,:) = interp1(f,sqrt(num./denom),f_i);
    
    num = abs(S2hw.*Pww_88(:,2) - S2h8.*S8w).^2;
    denom = (Pww_22h(:,2).*Pww_88(:,2) - abs(S2h8).^2).*(Pww_22h(:,1).*Pww_88(:,2) - abs(Sw8).^2);
    partial_w2h(d,:) = interp1(f,sqrt(num./denom),f_i);
    
    num = abs(S28.*Pww_88(:,1) - S2w.*Sw8).^2;
    denom = (P88_22(:,2).*Pww_88(:,1) - abs(S2w).^2).*(P88_22(:,1).*Pww_88(:,1) - abs(S8w).^2);
    partial_82(d,:) = interp1(f,sqrt(num./denom),f_i);
    
    num = abs(S2h8.*Pww_88(:,1) - S2hw.*Sw8).^2;
    denom = (Pww_22h(:,2).*Pww_88(:,1) - abs(S2hw).^2).*(Pww_88(:,2).*Pww_88(:,1) - abs(S8w).^2);
    partial_82h(d,:) = interp1(f,sqrt(num./denom),f_i);
    
    if ~isnan(hpc_mua(d))
        num = abs(Smw.*Pww_88(:,2) - Sm8.*S8w).^2;
        denom = (Pww_mm(:,2).*Pww_88(:,2) - abs(Sm8).^2).*(Pww_mm(:,1).*Pww_88(:,2) - abs(Sw8).^2);
        partial_wm(d,:) = interp1(f,sqrt(num./denom),f_i);
        
        num = abs(Sm8.*Pww_88(:,1) - Smw.*Sw8).^2;
        denom = (Pww_mm(:,2).*Pww_88(:,1) - abs(Smw).^2).*(Pww_88(:,2).*Pww_88(:,1) - abs(S8w).^2);
        partial_8m(d,:) = interp1(f,sqrt(num./denom),f_i);
    else
       partial_wm(d,:) = nan(1,length(f_i));
       partial_8m(d,:) = nan(1,length(f_i));
    end
    correction(d) = 1/(dof-2);   
end

%%
cd C:\WC_Germany\sven_thomas_combined\
save combined_coh_mua_np aC* f_i lP* partial* correction
%%
n_recs = size(lPww,1);
n_mec = length(l3mec);
n_lec = length(l3lec);
used_hpc_lfp = find(~isnan(hpc_lfp));
mec_hpc = l3mec(ismember(l3mec,used_hpc_lfp));
lec_hpc = l3lec(ismember(l3lec,used_hpc_lfp));
lred = [0.6 0.3 0.3];
lblue = [0.5 0.5 0.9];
Cw8_cor = tanh(aCw8_cor);

figure;hold on
plot(f_i,mean(Cw8_cor(mec_hpc,:)),'r','linewidth',2)
plot(f_i,mean(Cw8_cor(lec_hpc,:)),'b','linewidth',2)
plot(f_i,mean(partial_w2h(mec_hpc,:)),'color',lred,'linewidth',2)
plot(f_i,mean(partial_w2h(lec_hpc,:)),'color',lblue,'linewidth',2)
legend('MEC-Ctx LFP','LEC-Ctx LFP','MEC-Hpc LFP','LEC-Hpc LFP')
legend('MEC-Ctx LFP','LEC-Ctx LFP')
shadedErrorBar(f_i,mean(Cw8_cor(mec_hpc,:)),std(Cw8_cor(mec_hpc,:))./sqrt(length(mec_hpc)),{'r'});
shadedErrorBar(f_i,mean(Cw8_cor(lec_hpc,:)),std(Cw8_cor(lec_hpc,:))./sqrt(length(lec_hpc)),{'b'});
shadedErrorBar(f_i,mean(partial_w2h(mec_hpc,:)),std(partial_w2h(mec_hpc,:))./sqrt(length(mec_hpc)),{'color',lred});
shadedErrorBar(f_i,mean(partial_w2h(lec_hpc,:)),std(partial_w2h(lec_hpc,:))./sqrt(length(lec_hpc)),{'color',lblue});
xlim([0 1.5])
xlabel('Frequency (Hz)','fontsize',14)
ylabel('Coherence','fontsize',14)

%%
used_hpc_lfp = 1:n_recs;
mec_hpc = l3mec(ismember(l3mec,used_hpc_lfp));
lec_hpc = l3lec(ismember(l3lec,used_hpc_lfp));
lred = [0.6 0.3 0.3];
lblue = [0.5 0.5 0.9];

figure;hold on
plot(f_i,mean(partial_w2h(mec_hpc,:)),'r','linewidth',2)
plot(f_i,mean(partial_w2h(lec_hpc,:)),'b','linewidth',2)
plot(f_i,mean(partial_82h(mec_hpc,:)),'color',lred,'linewidth',2)
plot(f_i,mean(partial_82h(lec_hpc,:)),'color',lblue,'linewidth',2)
legend('MEC-Hpc LFP','LEC-Hpc LFP','Ctx-Hpc LFP (MEC)','Ctx-Hpc LFP (LEC)')
shadedErrorBar(f_i,mean(partial_w2h(mec_hpc,:)),std(partial_w2h(mec_hpc,:))./sqrt(length(mec_hpc)),{'r'});
shadedErrorBar(f_i,mean(partial_w2h(lec_hpc,:)),std(partial_w2h(lec_hpc,:))./sqrt(length(lec_hpc)),{'b'});
shadedErrorBar(f_i,mean(partial_82h(mec_hpc,:)),std(partial_82h(mec_hpc,:))./sqrt(length(mec_hpc)),{'color',lred});
shadedErrorBar(f_i,mean(partial_82h(lec_hpc,:)),std(partial_82h(lec_hpc,:))./sqrt(length(lec_hpc)),{'color',lblue});
xlim([0 2])
xlabel('Frequency (Hz)','fontsize',14)
ylabel('Partial coherence','fontsize',14)


%%
used_mua = find(~isnan(hpc_mua));
gray = [0.2 0.2 0.2];
Cwm_cor = tanh(aCwm_cor);
C8m_cor = tanh(aCm8_cor);

figure;hold on
shadedErrorBar(f_i,mean(Cwm_cor(used_mua,:)),std(Cwm_cor(used_mua,:))./sqrt(length(used_mua)),{'g'});
shadedErrorBar(f_i,mean(C8m_cor(used_mua,:)),std(C8m_cor(used_mua,:))./sqrt(length(used_mua)),{'color',gray});
xlim([0 2])
xlabel('Frequency (Hz)','fontsize',14)
ylabel('Coherence','fontsize',14)

%%
used_mua = find(~isnan(hpc_mua));
figure;hold on; set(gca,'fontname','arial','fontsize',14,'fontname','arial')
plot(f_i,mean(partial_wm(used_mua,:)),'r','linewidth',2)
plot(f_i,mean(partial_8m(used_mua,:)),'color',[0.2 0.2 0.2],'linewidth',2);
legend('MEC MP','Ctx LFP')
shadedErrorBar(f_i,mean(partial_wm(used_mua,:)),std(partial_wm(used_mua,:))./sqrt(length(used_mua)),{'r'});
shadedErrorBar(f_i,mean(partial_8m(used_mua,:)),std(partial_8m(used_mua,:))./sqrt(length(used_mua)),{'color',[0.2 0.2 0.2]});
xlim([0 1.])
xlabel('Frequency (Hz)','fontsize',16,'fontname','arial')
ylabel('Coherence','fontsize',16,'fontname','arial')
%%
Cwm_cor = tanh(aCwm_cor);
C8m_cor = tanh(aCm8_cor);
uds_freqs = find(f_i > 0.1 & f_i < 1);
[Cwm_peak,Cwm_peakf] = max(Cwm_cor(:,uds_freqs),[],2);
[C8m_peak,C8m_peakf] = max(C8m_cor(:,uds_freqs),[],2);

[Cwm_par_peak,Cwm_par_peakf] = max(partial_wm(:,uds_freqs),[],2);
[C8m_par_peak,C8m_par_peakf] = max(partial_8m(:,uds_freqs),[],2);

