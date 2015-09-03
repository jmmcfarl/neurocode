clear all
close all
%%
load F:\WC_Germany\overall_EC\overall_EC_dir.mat
addpath('F:\Code\WC_anal\general\')
addpath('F:\WC_Germany\Overall_EC\')
addpath('F:\Code\Chronux\spectral_analysis\continuous\')
addpath('F:\WC_Germany\hsmm_state_detection\\')

drive_letter = 'F';

%%
dsf = 32;
Fsd = 2016/dsf;
niqf = 2016/2;
[b,a] = butter(2,[0.05/niqf 10/niqf]);

params.Fs = Fsd;
params.err = [2 0.05];
params.tapers = [10 19];
W = 0.1;
params.fpass = [0 8];
% winlength = 25;
% winslide = 25;
% movingwin = [winlength winslide];

f_i = linspace(W,10,250);

for d = 1:length(sess_data)
    
    cdir = sess_data(d).directory;
    cdir(1) = 'F';
    disp(sprintf('session %d',d))
    cd(cdir);
    s_name = strcat(sess_data(d).region,'_l',sess_data(d).layer,'_',sess_data(d).name);
    load used_data lf8 lf2 lf5 lf3 wcv_minus_spike

    wcv_d = filtfilt(b,a,wcv_minus_spike);
    wcv_d = downsample(wcv_d,dsf)/sess_data(d).gains(1);
    lf8_d = filtfilt(b,a,lf8);
    lf8_d = downsample(lf8_d,dsf)/sess_data(d).gains(8);
    lf3_d = filtfilt(b,a,lf3);
    lf3_d = downsample(lf3_d,dsf)/sess_data(d).gains(3);
    lf2_d = filtfilt(b,a,lf2);
    lf2_d = downsample(lf2_d,dsf)/sess_data(d).gains(2);
    lf5_d = filtfilt(b,a,lf5);
    lf5_d = downsample(lf5_d,dsf)/sess_data(d).gains(5);
%     lf2_r = lf2_d - lf5_d;
    lf2_r = lf2/sess_data(d).gains(2)-lf5/sess_data(d).gains(5);
    lf2_r = downsample(filtfilt(b,a,lf2_r),dsf);    
    t_axis = (1:length(wcv_d))/Fsd;
    
    t_axis = (1:length(lf3_d))/Fsd;
    load ec_hmm_state_seq8
    [new_seg_inds] = resample_uds_seg_inds(hmm8.UDS_segs,hmm8.Fs,Fsd,length(lf8_d));
    
    [Cw8,Phiw8,Sw8,Pww_88,f,ConfC,dof] = coherencyc_unequal_length_trials_fixedW([wcv_d lf8_d],W,params,new_seg_inds);
    [Cw3,Phiw3,Sw3,Pww_33,~,ConfC,dof] = coherencyc_unequal_length_trials_fixedW([wcv_d lf3_d],W,params,new_seg_inds);
    [C83,Phi83,S83,P88_33,f,ConfC,dof] = coherencyc_unequal_length_trials_fixedW([lf8_d lf3_d],W,params,new_seg_inds);
    [Cw2r,Phiw2r,Sw2r,Pww_22r,~,ConfC,dof] = coherencyc_unequal_length_trials_fixedW([wcv_d lf2_r],W,params,new_seg_inds);
    [C82r,Phi82r,S82r,P88_22r,f,ConfC,dof] = coherencyc_unequal_length_trials_fixedW([lf8_d lf2_r],W,params,new_seg_inds);
%     [Cw2,Phiw2,Sw2,Pww_22,~,ConfC,dof] = coherencyc_unequal_length_trials_fixedW([wcv_d lf2_d],W,params,new_seg_inds);
%     [C82,Phi82,S82,P88_22,f,ConfC,dof] = coherencyc_unequal_length_trials_fixedW([lf8_d lf2_d],W,params,new_seg_inds);
%     [C52r,Phi52r,S52r,P55_22r,f,ConfC,dof] = coherencyc_unequal_length_trials_fixedW([lf5_d lf2_r],W,params,new_seg_inds);
%     [C53,Phi53,S53,P55_33,f,ConfC,dof] = coherencyc_unequal_length_trials_fixedW([lf5_d lf3_d],W,params,new_seg_inds);
 
    S38 = conj(S83);
    S2r8 = conj(S82r);
%     S28 = conj(S82);
    S8w = conj(Sw8);
    S3w = conj(Sw3);
    S2rw = conj(Sw2r);
    
    aCw8_cor_temp = atanh(Cw8) - 1/(dof-2);
    aCw8_cor(d,:) = interp1(f,aCw8_cor_temp,f_i);
    Pw8(d,:) = interp1(f,Phiw8,f_i);
    
    aCw3_cor_temp = atanh(Cw3) - 1/(dof-2);
    aCw3_cor(d,:) = interp1(f,aCw3_cor_temp,f_i);
    Pw3(d,:) = interp1(f,Phiw3,f_i);
    
    aCw2r_cor_temp = atanh(Cw2r) - 1/(dof-2);
    aCw2r_cor(d,:) = interp1(f,aCw2r_cor_temp,f_i);
    Pw2r(d,:) = interp1(f,Phiw2r,f_i);
    
    aC83_cor_temp = atanh(C83) - 1/(dof-2);
    aC83_cor(d,:) = interp1(f,aC83_cor_temp,f_i);
    P83(d,:) = interp1(f,Phi83,f_i);
    
    aC82r_cor_temp = atanh(C82r) - 1/(dof-2);
    aC82r_cor(d,:) = interp1(f,aC82r_cor_temp,f_i);
    P82r(d,:) = interp1(f,Phi82r,f_i);
    
%     aC52r_cor_temp = atanh(C52r) - 1/(dof-2);
%     aC52r_cor(d,:) = interp1(f,aC52r_cor_temp,f_i);
%     P52r(d,:) = interp1(f,Phi52r,f_i);
%     
%     aC53_cor_temp = atanh(C53) - 1/(dof-2);
%     aC53_cor(d,:) = interp1(f,aC53_cor_temp,f_i);
%     P53(d,:) = interp1(f,Phi53,f_i);
    
    lPww(d,:) = interp1(f,log10(Pww_88(:,1)),f_i);
    lP88(d,:) = interp1(f,log10(Pww_88(:,2)),f_i);
%     lP55(d,:) = interp1(f,log10(P55_33(:,1)),f_i);
    lP33(d,:) = interp1(f,log10(Pww_33(:,2)),f_i);
%     lP22(d,:) = interp1(f,log10(Pww_22(:,2)),f_i);
    lP22r(d,:) = interp1(f,log10(Pww_22r(:,2)),f_i);
    
    num = abs(S3w.*Pww_88(:,2) - S38.*S8w).^2;
    denom = (Pww_33(:,2).*Pww_88(:,2) - abs(S38).^2).*(Pww_33(:,1).*Pww_88(:,2) - abs(Sw8).^2);
    partial_w3(d,:) = interp1(f,num./denom,f_i);
%     
%     num = abs(S2w.*Pww_88(:,2) - S28.*S8w).^2;
%     denom = (Pww_22(:,2).*Pww_88(:,2) - abs(S28).^2).*(Pww_22(:,1).*Pww_88(:,2) - abs(Sw8).^2);
%     partial_w2(d,:) = interp1(f,num./denom,f_i);
%     
    num = abs(S2rw.*Pww_88(:,2) - S2r8.*S8w).^2;
    denom = (Pww_22r(:,2).*Pww_88(:,2) - abs(S2r8).^2).*(Pww_22r(:,1).*Pww_88(:,2) - abs(Sw8).^2);
    partial_w2r(d,:) = interp1(f,num./denom,f_i);
%     
    numer = abs(S8w.*Pww_33(:,2) - S83.*S3w).^2;
    denomer = (Pww_88(:,2).*Pww_33(:,2) - abs(S83).^2).*(Pww_88(:,1).*Pww_33(:,2) - abs(Sw3).^2);
    partial_w8(d,:) = interp1(f,numer./denomer,f_i);
    
end

%%
mec = find_struct_field_vals(sess_data,'region','MEC');
layer3 = find_struct_field_vals(sess_data,'layer','3');
layer2 = find_struct_field_vals(sess_data,'layer','2');
layer23 = find_struct_field_vals(sess_data,'layer','23');
lec = find_struct_field_vals(sess_data,'region','LEC');
l3mec = intersect(mec,layer3);
l2mec = intersect(mec,layer2);
l3lec = intersect(lec,layer3);
l2lec = intersect(lec,layer2);
l3mec(24:end) = [];
l23mec = intersect(mec,layer23);
l23mec = unique([l2mec l3mec l23mec]);

frange = find(f_i > 0.25 & f_i < 0.6);
for ii = 1:size(P83,1)
    avg_P83(ii) = circ_mean(P83(ii,frange));
end
correct_lf3phase = find(avg_P83 > -0.4);
l3mec_c = l3mec(avg_P83(l3mec) > -0.4);
l3mec_w = l3mec(avg_P83(l3mec) < -0.4);
l3lec_c = l3lec(avg_P83(l3lec) > -0.4);
l3lec_w = l3lec(avg_P83(l3lec) < -0.4);

max_partial_w3 = max(partial_w3,[],2);
max_partial_w2r = max(partial_w2r,[],2);
max_partial_w8 = max(partial_w8,[],2);
max_C83 = max(aC83_cor,[],2);

for ii = 1:size(lPww,1)
    [max_Pww(ii),max_Pww_loc(ii)] = nanmax(lPww(ii,:));
    [max_P88(ii),max_P88_loc(ii)] = nanmax(lP88(ii,:));
    [max_P33(ii),max_P33_loc(ii)] = nanmax(lP33(ii,:));
end
max_Pww_loc = f_i(max_Pww_loc);
max_P88_loc = f_i(max_P88_loc);
max_P33_loc = f_i(max_P33_loc);

partial_w3_f2 = partial_w3(:,4);
partial_w3_f4 = partial_w3(:,10);
partial_w8_f2 = partial_w8(:,4);
partial_w8_f4 = partial_w8(:,10);

cd F:\WC_Germany\overall_EC\
save overall_EC_coherence_data_gains aC* f_i Pw8 Pw3 P83 P82r Pw2r lP* partial* max_*
%%
figure
h = errorbar(f_i,nanmean(partial_w3(l3mec_c,:)),nanstd(partial_w3(l3mec_c,:))/sqrt(length(l3mec_c)));
errorbar_tick(h,.01,'units');
hold on
h=errorbar(f_i,nanmean(partial_w3(l3lec_c,:)),nanstd(partial_w3(l3lec_c,:))/sqrt(length(l3lec_c)),'r');
errorbar_tick(h,.01,'units')
xlim([0 2])

figure
h = errorbar(f_i,nanmean(partial_w8(l3mec,:)),nanstd(partial_w8(l3mec,:))/sqrt(length(l3mec)));
errorbar_tick(h,.01,'units');
hold on
h=errorbar(f_i,nanmean(partial_w8(l3lec,:)),nanstd(partial_w8(l3lec,:))/sqrt(length(l3lec)),'r');
errorbar_tick(h,.01,'units')
xlim([0 2])

figure
h = errorbar(f_i,nanmean(aCw8_cor(l3mec,:)),nanstd(aCw8_cor(l3mec,:))/sqrt(length(l3mec)));
errorbar_tick(h,.01,'units');
hold on
h=errorbar(f_i,nanmean(aCw8_cor(l3lec,:)),nanstd(aCw8_cor(l3lec,:))/sqrt(length(l3lec)),'r');
errorbar_tick(h,.01,'units')
xlim([0 2])

figure
h = errorbar(f_i,nanmean(aCw3_cor(l3mec,:)),nanstd(aCw3_cor(l3mec,:))/sqrt(length(l3mec)));
errorbar_tick(h,.01,'units');
hold on
h=errorbar(f_i,nanmean(aCw3_cor(l3lec,:)),nanstd(aCw3_cor(l3lec,:))/sqrt(length(l3lec)),'r');
errorbar_tick(h,.01,'units')
xlim([0 2])

figure
h = errorbar(f_i,nanmean(aC83_cor(l3mec,:)),nanstd(aC83_cor(l3mec,:))/sqrt(length(l3mec)));
errorbar_tick(h,.01,'units');
hold on
h=errorbar(f_i,nanmean(aC83_cor(l3lec,:)),nanstd(aC83_cor(l3lec,:))/sqrt(length(l3lec)),'r');
errorbar_tick(h,.01,'units')
xlim([0 2])

figure
h = errorbar(f_i,nanmean(aC83_cor(l3mec,:)),nanstd(aC83_cor(l3mec,:))/sqrt(length(l3mec)));
errorbar_tick(h,.01,'units');
hold on
h=errorbar(f_i,nanmean(aCw3_cor(l3mec,:)),nanstd(aCw3_cor(l3mec,:))/sqrt(length(l3mec)),'r');
errorbar_tick(h,.01,'units')
h=errorbar(f_i,nanmean(aCw8_cor(l3mec,:)),nanstd(aCw8_cor(l3mec,:))/sqrt(length(l3mec)),'c');
errorbar_tick(h,.01,'units')
xlim([0 2])

figure
h = errorbar(f_i,nanmean(aC83_cor(l3lec,:)),nanstd(aC83_cor(l3lec,:))/sqrt(length(l3lec)));
errorbar_tick(h,.01,'units');
hold on
h=errorbar(f_i,nanmean(aCw3_cor(l3lec,:)),nanstd(aCw3_cor(l3lec,:))/sqrt(length(l3lec)),'r');
errorbar_tick(h,.01,'units')
h=errorbar(f_i,nanmean(aCw8_cor(l3lec,:)),nanstd(aCw8_cor(l3lec,:))/sqrt(length(l3lec)),'c');
errorbar_tick(h,.01,'units')
xlim([0 2])

figure
h=errorbar(f_i,nanmean(lPww(l3mec,:)),nanstd(lPww(l3mec,:))/sqrt(length(l3mec)));
errorbar_tick(h,.01,'units');
hold on
h=errorbar(f_i,nanmean(lPww(l3lec,:)),nanstd(lPww(l3lec,:))/sqrt(length(l3lec)),'g');
errorbar_tick(h,.01,'units');
xlim([0 2])
legend('L3MEC','L3LEC')
xlabel('Frequency (Hz)')
ylabel('Power')

figure
h=errorbar(f_i,nanmean(lP88),nanstd(lP88)/sqrt(size(lP88,1)),'r');
errorbar_tick(h,.01,'units');
hold on
h=errorbar(f_i,nanmean(lP33(correct_lf3phase,:)),nanstd(lP33(correct_lf3phase,:))/sqrt(length(correct_lf3phase)),'k');
errorbar_tick(h,.01,'units');
xlim([0 2])
legend('LF8','LF3')
xlabel('Frequency (Hz)')
ylabel('Power')


