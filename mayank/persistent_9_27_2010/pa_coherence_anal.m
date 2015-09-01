clear all
close all
%%
load G:\WC_Germany\overall_EC\overall_EC_dir.mat
addpath('G:\Code\WC_anal\general\')
addpath('G:\WC_Germany\Overall_EC\')
addpath('G:\Code\Chronux\spectral_analysis\continuous\')
addpath('G:\WC_Germany\hsmm_state_detection\\')

drive_letter = 'G';

used_data = [l3mec_p l3lec_p];
sess_data = sess_data(used_data);

%%
dsf = 16;
Fsd = 2016/dsf;
niqf = 2016/2;
% [b,a] = butter(2,[0.05/niqf 40/niqf]);
[b,a] = butter(2,40/niqf,'low');

params.Fs = Fsd;
params.err = [2 0.05];
params.tapers = [7 13];
W = 0.1;
params.fpass = [0 8];
winlength = 50;
winslide = 50;
movingwin = [winlength winslide];

f_i = linspace(0,8,500);

for d = 1:length(sess_data)
    
    cdir = sess_data(d).directory;
    cdir(1) = 'G';
    disp(sprintf('session %d',d))
    cd(cdir);
    s_name = strcat(sess_data(d).region,'_l',sess_data(d).layer,'_',sess_data(d).name);
    load used_data lf8 lf3 wcv_minus_spike lf5 lf2
% lf3 = lf2;
lf2r = lf2-lf5;
lf3g = lf3+lf5;

    wcv_d = filtfilt(b,a,wcv_minus_spike);
    wcv_d = downsample(wcv_d,dsf);
    lf8_d = filtfilt(b,a,lf8);
    lf8_d = downsample(lf8_d,dsf);
    lf5_d = filtfilt(b,a,lf5);
    lf5_d = downsample(lf5_d,dsf);
    lf3_d = filtfilt(b,a,lf3);
    lf3_d = downsample(lf3_d,dsf);
    lf3_gd = filtfilt(b,a,lf3g);
    lf3_gd = downsample(lf3_gd,dsf);
    lf2_d = filtfilt(b,a,lf2);
    lf2_d = downsample(lf2_d,dsf);
    lf2_rd = filtfilt(b,a,lf2r);
    lf2_rd = downsample(lf2_rd,dsf);
    t_axis = (1:length(wcv_d))/Fsd;
    
    lf3_hf = get_hf_features(lf3,2016,Fsd,[20 80],0.05);
    lf5_hf = get_hf_features(lf5,2016,Fsd,[20 80],0.05);
    
%     wcv_d = zscore(wcv_d);
%     lf8_d = zscore(lf8_d);
%     lf3_d = zscore(lf3_d);
%     lf2_d = zscore(lf2_d);
%     lf5_d = zscore(lf5_d);

    t_axis = (1:length(lf3_d))/Fsd;
    load pa_hsmm_state_seq8
    [new_seg_inds] = resample_uds_seg_inds(hmm8.UDS_segs,hmm8.Fs,Fsd,length(lf8_d));
    
%     [Cw8,Phiw8,Sw8,Pww_88,f,ConfC,dof] = coherencyc_unequal_length_trials_fixedW([wcv_d lf8_d],W,params,new_seg_inds);
%     [Cw3,Phiw3,Sw3,Pww_33,~,ConfC,dof] = coherencyc_unequal_length_trials_fixedW([wcv_d lf3_d],W,params,new_seg_inds);
%     [C83,Phi83,S83,P88_33,f,ConfC,dof] = coherencyc_unequal_length_trials_fixedW([lf8_d lf3_d],W,params,new_seg_inds);
    [Cw8,Phiw8,Sw8,Pww_88,f,ConfC,dof] = coherencyc_unequal_length_trials_jmm([wcv_d lf8_d],movingwin,params,new_seg_inds);
    [Cw3,Phiw3,Sw3,Pww_33,~,ConfC,dof] = coherencyc_unequal_length_trials_jmm([wcv_d lf3_d],movingwin,params,new_seg_inds);
    [Cw5,Phiw5,Sw5,Pww_55,~,ConfC,dof] = coherencyc_unequal_length_trials_jmm([wcv_d lf5_d],movingwin,params,new_seg_inds);
    [C83,Phi83,S83,P88_33,f,ConfC,dof] = coherencyc_unequal_length_trials_jmm([lf8_d lf3_d],movingwin,params,new_seg_inds);
    [C83g,Phi83g,S83g,P88_33g,f,ConfC,dof] = coherencyc_unequal_length_trials_jmm([lf8_d lf3_gd],movingwin,params,new_seg_inds);
    [C82,Phi82,S82,P88_22,f,ConfC,dof] = coherencyc_unequal_length_trials_jmm([lf8_d lf2_d],movingwin,params,new_seg_inds);
    [C82r,Phi82r,S82r,P88_22r,f,ConfC,dof] = coherencyc_unequal_length_trials_jmm([lf8_d lf2_rd],movingwin,params,new_seg_inds);
    [C85,Phi85,S85,P88_55,f,ConfC,dof] = coherencyc_unequal_length_trials_jmm([lf8_d lf5_d],movingwin,params,new_seg_inds);
    [C53,Phi53,S53,P55_33,f,ConfC,dof] = coherencyc_unequal_length_trials_jmm([lf5_d lf3_d],movingwin,params,new_seg_inds);
    [C32,Phi32,S32,P33_22,f,ConfC,dof] = coherencyc_unequal_length_trials_jmm([lf3_d lf2_d],movingwin,params,new_seg_inds);
%     [Cw3h,Phiw3h,Sw3h,Pww_33h,~,ConfC,dof] = coherencyc_unequal_length_trials_jmm([wcv_d lf3_hf],movingwin,params,new_seg_inds);
%     [C83h,Phi83h,S83h,P88_33h,f,ConfC,dof] = coherencyc_unequal_length_trials_jmm([lf8_d lf3_hf],movingwin,params,new_seg_inds);
%     [C85h,Phi85h,S85h,P88_55h,f,ConfC,dof] = coherencyc_unequal_length_trials_jmm([lf8_d lf5_hf],movingwin,params,new_seg_inds);
%     [Cw5h,Phiw5h,Sw5h,Pww_55h,~,ConfC,dof] = coherencyc_unequal_length_trials_jmm([wcv_d lf5_hf],movingwin,params,new_seg_inds);
 
    S38 = conj(S83);
    S35 = conj(S53);
    S58 = conj(S85);
%     S3h8 = conj(S83h);
%     S5h8 = conj(S85h);
    S8w = conj(Sw8);
    S3w = conj(Sw3);
    S5w = conj(Sw5);
%     S3hw = conj(Sw3h);
%     S5hw = conj(Sw5h);
    
    aCw8_cor_temp = atanh(Cw8) - 1/(dof-2);
    aCw8_cor(d,:) = interp1(f,aCw8_cor_temp,f_i);
    Pw8(d,:) = interp1(f,Phiw8,f_i);
    
    aCw3_cor_temp = atanh(Cw3) - 1/(dof-2);
    aCw3_cor(d,:) = interp1(f,aCw3_cor_temp,f_i);
    Pw3(d,:) = interp1(f,Phiw3,f_i);

    aCw5_cor_temp = atanh(Cw5) - 1/(dof-2);
    aCw5_cor(d,:) = interp1(f,aCw5_cor_temp,f_i);
    Pw5(d,:) = interp1(f,Phiw5,f_i);

%     aCw3h_cor_temp = atanh(Cw3h) - 1/(dof-2);
%     aCw3h_cor(d,:) = interp1(f,aCw3h_cor_temp,f_i);
%     Pw3h(d,:) = interp1(f,Phiw3h,f_i);

%     aCw5h_cor_temp = atanh(Cw5h) - 1/(dof-2);
%     aCw5h_cor(d,:) = interp1(f,aCw5h_cor_temp,f_i);
%     Pw5h(d,:) = interp1(f,Phiw5h,f_i);

    aC83_cor_temp = atanh(C83) - 1/(dof-2);
    aC83_cor(d,:) = interp1(f,aC83_cor_temp,f_i);
    P83(d,:) = interp1(f,Phi83,f_i);

    aC83g_cor_temp = atanh(C83g) - 1/(dof-2);
    aC83g_cor(d,:) = interp1(f,aC83g_cor_temp,f_i);
    P83g(d,:) = interp1(f,Phi83g,f_i);

    aC82_cor_temp = atanh(C82) - 1/(dof-2);
    aC82_cor(d,:) = interp1(f,aC82_cor_temp,f_i);
    P82(d,:) = interp1(f,Phi82,f_i);

    aC82r_cor_temp = atanh(C82r) - 1/(dof-2);
    aC82r_cor(d,:) = interp1(f,aC82r_cor_temp,f_i);
    P82r(d,:) = interp1(f,Phi82r,f_i);

    aC53_cor_temp = atanh(C53) - 1/(dof-2);
    aC53_cor(d,:) = interp1(f,aC53_cor_temp,f_i);
    P53(d,:) = interp1(f,Phi53,f_i);
    
    aC32_cor_temp = atanh(C32) - 1/(dof-2);
    aC32_cor(d,:) = interp1(f,aC32_cor_temp,f_i);
    P32(d,:) = interp1(f,Phi32,f_i);

    %     aC83h_cor_temp = atanh(C83h) - 1/(dof-2);
%     aC83h_cor(d,:) = interp1(f,aC83h_cor_temp,f_i);
%     P83h(d,:) = interp1(f,Phi83h,f_i);
% 
%     aC85h_cor_temp = atanh(C85h) - 1/(dof-2);
%     aC85h_cor(d,:) = interp1(f,aC85h_cor_temp,f_i);
%     P85h(d,:) = interp1(f,Phi85h,f_i);

    lPww(d,:) = interp1(f,log10(Pww_88(:,1)),f_i);
    lP88(d,:) = interp1(f,log10(Pww_88(:,2)),f_i);
    lP55(d,:) = interp1(f,log10(P55_33(:,1)),f_i);
    lP33(d,:) = interp1(f,log10(Pww_33(:,2)),f_i);
    lP33g(d,:) = interp1(f,log10(P88_33g(:,2)),f_i);
    lP22(d,:) = interp1(f,log10(P88_22(:,2)),f_i);
    lP22r(d,:) = interp1(f,log10(P88_22r(:,2)),f_i);
%     lP33h(d,:) = interp1(f,log10(Pww_33h(:,2)),f_i);
%     lP55h(d,:) = interp1(f,log10(Pww_55h(:,2)),f_i);
    
    num = abs(S3w.*Pww_88(:,2) - S38.*S8w).^2;
    denom = (Pww_33(:,2).*Pww_88(:,2) - abs(S38).^2).*(Pww_33(:,1).*Pww_88(:,2) - abs(Sw8).^2);
    partial_w3(d,:) = interp1(f,num./denom,f_i);

    num = abs(S3w.*Pww_55(:,2) - S35.*S5w).^2;
    denom = (Pww_33(:,2).*Pww_55(:,2) - abs(S35).^2).*(Pww_33(:,1).*Pww_55(:,2) - abs(Sw5).^2);
    partial_w3_5(d,:) = interp1(f,num./denom,f_i);

    num = abs(S5w.*Pww_33(:,2) - S53.*S3w).^2;
    denom = (Pww_55(:,2).*Pww_33(:,2) - abs(S53).^2).*(Pww_55(:,1).*Pww_33(:,2) - abs(Sw3).^2);
    partial_w5_3(d,:) = interp1(f,num./denom,f_i);

    num = abs(S5w.*Pww_88(:,2) - S58.*S8w).^2;
    denom = (Pww_55(:,2).*Pww_88(:,2) - abs(S58).^2).*(Pww_55(:,1).*Pww_88(:,2) - abs(Sw8).^2);
    partial_w5(d,:) = interp1(f,num./denom,f_i);

    num = abs(S8w.*Pww_55(:,2) - S85.*S5w).^2;
    denom = (Pww_88(:,2).*Pww_55(:,2) - abs(S85).^2).*(Pww_88(:,1).*Pww_55(:,2) - abs(Sw5).^2);
    partial_w8(d,:) = interp1(f,num./denom,f_i);

%     num = abs(S3hw.*Pww_88(:,2) - S3h8.*S8w).^2;
%     denom = (Pww_33h(:,2).*Pww_88(:,2) - abs(S3h8).^2).*(Pww_33h(:,1).*Pww_88(:,2) - abs(Sw8).^2);
%     partial_w3h(d,:) = interp1(f,num./denom,f_i);
% 
%     num = abs(S5hw.*Pww_88(:,2) - S5h8.*S8w).^2;
%     denom = (Pww_55h(:,2).*Pww_88(:,2) - abs(S5h8).^2).*(Pww_55h(:,1).*Pww_88(:,2) - abs(Sw8).^2);
%     partial_w5h(d,:) = interp1(f,num./denom,f_i);

    numer = abs(S8w.*Pww_33(:,2) - S83.*S3w).^2;
    denomer = (Pww_88(:,2).*Pww_33(:,2) - abs(S83).^2).*(Pww_88(:,1).*Pww_33(:,2) - abs(Sw3).^2);
    partial_w8(d,:) = interp1(f,numer./denomer,f_i);
    
    correction(d) = 1/(dof-2);
    
    used_freqs = find(f_i < 1);
%     smooth_spectrum = jmm_smooth_1d_cor(lP88(d,used_freqs),3);
%     [temp,temploc] = findpeaks(smooth_spectrum,'npeaks',1);
%     if ~isempty(temp)
%         max_P88(d)= temp; max_P88_loc(d) = temploc;
%     end
    smooth_spectrum = jmm_smooth_1d_cor(lP55(d,used_freqs),3);
    [temp,temploc] = findpeaks(smooth_spectrum,'npeaks',1);
    if ~isempty(temp)
        max_P55(d)= temp; max_P55_loc(d) = temploc;
    end
%     peak_P83(d) = P83(d,max_P88_loc(d));
    peak_P53(d) = P53(d,max_P55_loc(d));

end

%%
cd G:\WC_Germany\persistent_9_27_2010\
save pa_coherence_data_6 aC* f_i Pw8 Pw3* P83* P82* lP* partial* peak* max_*
%%
mec = 1:22;
lec = 23:36;
bad_lf3 = [7 9 11 12 13 16 17 25 31 32 35];
good_lf3 = setdiff(1:36,bad_lf3);
good_lf3_mec = good_lf3(ismember(good_lf3,mec));
good_lf3_lec = good_lf3(ismember(good_lf3,lec));

%%
for i = 1:36
    subplot(2,1,1)
    plot(f_i,P83(i,:),'b')
    hold on
    plot(f_i,P83g(i,:),'c')
    plot(f_i,P82(i,:),'k')
    plot(f_i,P82r(i,:),'g')
    xlim([0 1])
    ylim([-pi pi])
    i
    subplot(2,1,2)
    plot(f_i,10.^lP88(i,:),'r')
    hold on
    plot(f_i,10.^lP33(i,:),'b')
    plot(f_i,10.^lP33g(i,:),'c')
    plot(f_i,10.^lP22(i,:),'k')
    plot(f_i,10.^lP22r(i,:),'g')
    xlim([0 1])
    pause
    clf
end

%%
Y = pdist(peak_P53',@(Xi,Xj) abs(circ_dist2(Xi,Xj)));
Z = linkage(Y,'average');
n_clusters = 2;
T = cluster(Z,n_clusters);
figure
plot(peak_P53,temp,'.','markersize',14)
hold on
plot(peak_P53(T==1),temp(T==1),'r.','markersize',14)
plot(peak_P53(lec),temp(lec),'ko','linewidth',2)
xlabel('Relative Phase (radians','fontsize',14)
ylabel('Coherence','fontsize',14)

%%
Cw8 = tanh(aCw8_cor);

% figure
% hold on
% h = errorbar(f_i,nanmean(Cw8(mec,:)),nanstd(Cw8(mec,:))/sqrt(length(mec)));
% errorbar_tick(h,.001,'units');
% h = errorbar(f_i,nanmean(Cw8(lec,:)),nanstd(Cw8(lec,:))/sqrt(length(lec)),'color',[0.2 0.6 0.2]);
% errorbar_tick(h,.001,'units');
% h = errorbar(f_i,nanmean(partial_w3(mec,:)),nanstd(partial_w3(mec,:))/sqrt(length(mec)),'c');
% errorbar_tick(h,.001,'units');
% h = errorbar(f_i,nanmean(partial_w3(lec,:)),nanstd(partial_w3(lec,:))/sqrt(length(lec)),'g');
% errorbar_tick(h,.001,'units');
% xlim([0 1])

figure
set(gca,'fontname','arial')
hold on
h = errorbar(f_i,nanmean(Cw8(mec,:)),nanstd(Cw8(mec,:))/sqrt(length(mec)));
errorbar_tick(h,.001,'units');
h = errorbar(f_i,nanmean(Cw8(lec,:)),nanstd(Cw8(lec,:))/sqrt(length(lec)),'color',[0.2 0.6 0.2]);
errorbar_tick(h,.001,'units');
h = errorbar(f_i,nanmean(partial_w3(good_lf3_mec,:)),nanstd(partial_w3(good_lf3_mec,:))/sqrt(length(good_lf3_mec)),'c');
errorbar_tick(h,.001,'units');
h = errorbar(f_i,nanmean(partial_w3(good_lf3_lec,:)),nanstd(partial_w3(good_lf3_lec,:))/sqrt(length(good_lf3_lec)),'g');
errorbar_tick(h,.001,'units');
xlim([0 1])
legend('MEC-LF8','LEC-LF8','MEC-LF3','LEC-LF3')
xlabel('Frequency (Hz)','fontsize',14,'fontname','arial')
ylabel('Coherence','fontsize',14,'fontname','arial')

figure
hold on
h = errorbar(f_i,nanmean(Cw8(mec,:)),nanstd(Cw8(mec,:))/sqrt(length(mec)));
errorbar_tick(h,.001,'units');
h = errorbar(f_i,nanmean(Cw8(lec,:)),nanstd(Cw8(lec,:))/sqrt(length(lec)),'color',[0.2 0.6 0.2]);
errorbar_tick(h,.001,'units');
h = errorbar(f_i,nanmean(partial_w3_5(good_lf3_mec,:)),nanstd(partial_w3_5(good_lf3_mec,:))/sqrt(length(good_lf3_mec)),'c');
errorbar_tick(h,.001,'units');
h = errorbar(f_i,nanmean(partial_w3_5(good_lf3_lec,:)),nanstd(partial_w3_5(good_lf3_lec,:))/sqrt(length(good_lf3_lec)),'g');
errorbar_tick(h,.001,'units');
xlim([0 1])
legend('MEC-LF8','LEC-LF8','MEC-LF3','LEC-LF3')
xlabel('Frequency (Hz)','fontsize',14)
ylabel('Coherence','fontsize',14)

figure
hold on
h = errorbar(f_i,nanmean(Cw8(mec,:)),nanstd(Cw8(mec,:))/sqrt(length(mec)));
errorbar_tick(h,.001,'units');
h = errorbar(f_i,nanmean(Cw8(lec,:)),nanstd(Cw8(lec,:))/sqrt(length(lec)),'color',[0.2 0.6 0.2]);
errorbar_tick(h,.001,'units');
h = errorbar(f_i,nanmean(partial_w5_3(good_lf3_mec,:)),nanstd(partial_w5_3(good_lf3_mec,:))/sqrt(length(good_lf3_mec)),'c');
errorbar_tick(h,.001,'units');
h = errorbar(f_i,nanmean(partial_w5_3(good_lf3_lec,:)),nanstd(partial_w5_3(good_lf3_lec,:))/sqrt(length(good_lf3_lec)),'g');
errorbar_tick(h,.001,'units');
xlim([0 1])
legend('MEC-LF8','LEC-LF8','MEC-LF3','LEC-LF3')
xlabel('Frequency (Hz)','fontsize',14)
ylabel('Coherence','fontsize',14)

figure
hold on
h = errorbar(f_i,nanmean(Cw8(mec,:)),nanstd(Cw8(mec,:))/sqrt(length(mec)));
errorbar_tick(h,.001,'units');
h = errorbar(f_i,nanmean(Cw8(lec,:)),nanstd(Cw8(lec,:))/sqrt(length(lec)),'color',[0.2 0.6 0.2]);
errorbar_tick(h,.001,'units');
h = errorbar(f_i,nanmean(partial_w5(good_lf3_mec,:)),nanstd(partial_w5(good_lf3_mec,:))/sqrt(length(good_lf3_mec)),'c');
errorbar_tick(h,.001,'units');
h = errorbar(f_i,nanmean(partial_w5(good_lf3_lec,:)),nanstd(partial_w5(good_lf3_lec,:))/sqrt(length(good_lf3_lec)),'g');
errorbar_tick(h,.001,'units');
xlim([0 1])
legend('MEC-LF8','LEC-LF8','MEC-LF3','LEC-LF3')
xlabel('Frequency (Hz)','fontsize',14)
ylabel('Coherence','fontsize',14)

figure
hold on
h = errorbar(f_i,nanmean(Cw8(mec,:)),nanstd(Cw8(mec,:))/sqrt(length(mec)));
errorbar_tick(h,.001,'units');
h = errorbar(f_i,nanmean(Cw8(lec,:)),nanstd(Cw8(lec,:))/sqrt(length(lec)),'color',[0.2 0.6 0.2]);
errorbar_tick(h,.001,'units');
h = errorbar(f_i,nanmean(partial_w3h(good_lf3_mec,:)),nanstd(partial_w3h(good_lf3_mec,:))/sqrt(length(good_lf3_mec)),'c');
errorbar_tick(h,.001,'units');
h = errorbar(f_i,nanmean(partial_w3h(good_lf3_lec,:)),nanstd(partial_w3h(good_lf3_lec,:))/sqrt(length(good_lf3_lec)),'g');
errorbar_tick(h,.001,'units');
xlim([0 1])
legend('MEC-LF8','LEC-LF8','MEC-LF3HF','LEC-LF3HF')
xlabel('Frequency (Hz)','fontsize',14)
ylabel('Coherence','fontsize',14)

figure
hold on
h = errorbar(f_i,nanmean(Cw8(mec,:)),nanstd(Cw8(mec,:))/sqrt(length(mec)));
errorbar_tick(h,.001,'units');
h = errorbar(f_i,nanmean(Cw8(lec,:)),nanstd(Cw8(lec,:))/sqrt(length(lec)),'color',[0.2 0.6 0.2]);
errorbar_tick(h,.001,'units');
h = errorbar(f_i,nanmean(partial_w5h(good_lf3_mec,:)),nanstd(partial_w5h(good_lf3_mec,:))/sqrt(length(good_lf3_mec)),'c');
errorbar_tick(h,.001,'units');
h = errorbar(f_i,nanmean(partial_w5h(good_lf3_lec,:)),nanstd(partial_w5h(good_lf3_lec,:))/sqrt(length(good_lf3_lec)),'g');
errorbar_tick(h,.001,'units');
xlim([0 1])
legend('MEC-LF8','LEC-LF8','MEC-LF3HF','LEC-LF3HF')
xlabel('Frequency (Hz)','fontsize',14)
ylabel('Coherence','fontsize',14)

used_freqs = find(f_i < 1);
[peak_Cw8,peak_Cw8_loc] = max(Cw8(:,used_freqs),[],2);
[peak_partw3,peak_partw3_loc] = max(partial_w3(:,used_freqs),[],2);