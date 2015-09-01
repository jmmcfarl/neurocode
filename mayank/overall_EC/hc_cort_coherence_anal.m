clear all
close all
%%
load G:\WC_Germany\wc_database.mat
addpath('G:\Code\WC_anal\general\')
addpath('G:\WC_Germany\Overall_EC\')
addpath('G:\Code\Chronux\spectral_analysis\continuous\')
addpath('G:\WC_Germany\hsmm_state_detection\\')

drive_letter = 'G';

%%
dsf = 32;
Fsd = 2016/dsf;
niqf = 2016/2;
[b,a] = butter(2,[0.05/niqf 10/niqf]);

params.Fs = Fsd;
params.err = [2 0.05];
params.tapers = [5 9];
W = 0.05;
params.fpass = [0 8];
winlength = 50;
winslide = 50;
win = [winlength winslide];

f_i = linspace(W,10,250);

for d = 95:length(wc_db)
    
    cdir = wc_db(d).directory;
    cdir(1) = 'G';
    disp(sprintf('session %d',d))
    cd(cdir);
    s_name = strcat(wc_db(d).region,'_',wc_db(d).celltype,'_',wc_db(d).date);
    load used_data lf8 lf3 wcv_minus_spike

    wcv_d = filtfilt(b,a,wcv_minus_spike);
    wcv_d = downsample(wcv_d,dsf);
    lf8_d = filtfilt(b,a,lf8);
    lf8_d = downsample(lf8_d,dsf);
    lf3_d = filtfilt(b,a,lf3);
    lf3_d = downsample(lf3_d,dsf);
    
    t_axis = (1:length(wcv_d))/Fsd;
%         load ec_hmm_state_seq
%     [new_seg_inds] = resample_uds_seg_inds(hmm.UDS_segs,hmm.Fs,Fsd,length(wcv_d));

%     [Cw8,Phiw8,Sw8,Pww_88,f,ConfC,dof] = coherencyc_unequal_length_trials_fixedW([wcv_d lf8_d],W,params,new_seg_inds);
%     [Cw3,Phiw3,Sw3,Pww_33,f,ConfC,dof] = coherencyc_unequal_length_trials_fixedW([wcv_d lf3_d],W,params,new_seg_inds);
%     [C83,Phi83,S83,P88_33,f,ConfC,dof] = coherencyc_unequal_length_trials_fixedW([lf8_d lf3_d],W,params,new_seg_inds);
%      [Cw8,Phiw8,Sw8,Pww_88,f,ConfC,dof] = coherencyc_unequal_length_trials([wcv_d lf8_d],movingwin,params,new_seg_inds);
%     [Cw3,Phiw3,Sw3,Pww_33,f,ConfC,dof] = coherencyc_unequal_length_trials([wcv_d lf3_d],movingwin,params,new_seg_inds);
%     [C83,Phi83,S83,P88_33,f,ConfC,dof] = coherencyc_unequal_length_trials([lf8_d lf3_d],movingwin,params,new_seg_inds);
 %     [Cw8,Phiw8,Sw8,Pww_88,f,ConfC,dof] = coherencyc_unequal_length_trials_fixedW([wcv_d lf8_d],W,params,new_seg_inds);
%     [Cw3,Phiw3,Sw3,Pww_33,f,ConfC,dof] = coherencyc_unequal_length_trials_fixedW([wcv_d lf3_d],W,params,new_seg_inds);
%     [C83,Phi83,S83,P88_33,f,ConfC,dof] = coherencyc_unequal_length_trials_fixedW([lf8_d lf3_d],W,params,new_seg_inds);

    [Cw8,Phiw8,Sw8,Sww,S88,f,ConfC] = coherencysegc(wcv_d, lf8_d,win,params);
    [Cw3,Phiw3,Sw3,Sww,S33,f,ConfC] = coherencysegc(wcv_d, lf3_d,win,params);
    [C83,Phi83,S83,S88,S33,f,ConfC] = coherencysegc(lf8_d, lf3_d,win,params);
    S38 = conj(S83);
    S8w = conj(Sw8);
    S3w = conj(Sw3);

    Pw8(d,:) = interp1(f,Phiw8,f_i);
    Pw3(d,:) = interp1(f,Phiw3,f_i);
    P83(d,:) = interp1(f,Phi83,f_i);
    lSww(d,:) = interp1(f,log10(Sww),f_i);
    lS88(d,:) = interp1(f,log10(S88),f_i);
    lS33(d,:) = interp1(f,log10(S33),f_i);
    
    num = abs(S3w.*S88 - S38.*S8w).^2;
    denom = (S33.*S88 - abs(S38).^2).*(Sww.*S88 - abs(Sw8).^2);
    partial_w3(d,:) = interp1(f,num./denom,f_i); 
 
    numer = abs(S8w.*S33 - S83.*S3w).^2;
    denomer = (S88.*S33 - abs(S83).^2).*(Sww.*S33 - abs(Sw3).^2);
    partial_w8(d,:) = interp1(f,numer./denomer,f_i); 
    
end

%%
cd G:\WC_Germany\overall_EC\
save hc_cort_coherence_data f_i Pw8 Pw3 P83 lS* partial*

%%
ca1 = find_struct_field_vals(wc_db,'region','ca1');
ca3 = find_struct_field_vals(wc_db,'region','ca3');
pyr = find_struct_field_vals(wc_db,'celltype','pyr');
int = find_struct_field_vals(wc_db,'celltype','int');
dg = find_struct_field_vals(wc_db,'region','dg');

ca1pyr = intersect(ca1,pyr);
ca1int = intersect(ca1,int);

frange = find(f_i > 0.2 & f_i < 0.6);
avg_P83 = mean(P83(:,frange),2);
ca1pyr_c = ca1pyr(avg_P83(ca1pyr) > -0.5);
ca1int_c = ca1int(avg_P83(ca1int) > -0.5);
dg_c = dg(avg_P83(dg) > -0.5);
ca3_c = ca3(avg_P83(ca3) > -0.5);


figure
h = errorbar(f_i,nanmean(partial_w3(ca1pyr_c,:)),nanstd(partial_w3(ca1pyr_c,:))/sqrt(length(ca1pyr_c)));
errorbar_tick(h,.01,'units');
hold on
h=errorbar(f_i,nanmean(partial_w3(ca1int_c,:)),nanstd(partial_w3(ca1int_c,:))/sqrt(length(ca1int_c)),'r');
errorbar_tick(h,.01,'units')
h=errorbar(f_i,nanmean(partial_w3(ca3_c,:)),nanstd(partial_w3(ca3_c,:))/sqrt(length(ca3_c)),'k');
errorbar_tick(h,.01,'units')
h=errorbar(f_i,nanmean(partial_w3(dg_c,:)),nanstd(partial_w3(dg_c,:))/sqrt(length(dg_c)),'g');
errorbar_tick(h,.01,'units')
xlim([0 5])

figure
h = errorbar(f_i,nanmean(partial_w8(ca1pyr,:)),nanstd(partial_w8(ca1pyr,:))/sqrt(length(ca1pyr)));
errorbar_tick(h,.01,'units');
hold on
h=errorbar(f_i,nanmean(partial_w8(ca1int,:)),nanstd(partial_w8(ca1int,:))/sqrt(length(ca1int)),'r');
errorbar_tick(h,.01,'units')
h=errorbar(f_i,nanmean(partial_w8(ca3,:)),nanstd(partial_w8(ca3,:))/sqrt(length(ca3)),'k');
errorbar_tick(h,.01,'units')
h=errorbar(f_i,nanmean(partial_w8(dg,:)),nanstd(partial_w8(dg,:))/sqrt(length(dg)),'g');
errorbar_tick(h,.01,'units')


figure
h = errorbar(f_i,nanmean(partial_w3(ca1pyr,:)),nanstd(partial_w3(ca1pyr,:))/sqrt(length(ca1pyr)));
errorbar_tick(h,.01,'units');
hold on
h=errorbar(f_i,nanmean(partial_w3(ca1int,:)),nanstd(partial_w3(ca1int,:))/sqrt(length(ca1int)),'r');
errorbar_tick(h,.01,'units')
h=errorbar(f_i,nanmean(partial_w3(ca3,:)),nanstd(partial_w3(ca3,:))/sqrt(length(ca3)),'k');
errorbar_tick(h,.01,'units')
h=errorbar(f_i,nanmean(partial_w3(dg,:)),nanstd(partial_w3(dg,:))/sqrt(length(dg)),'g');
errorbar_tick(h,.01,'units')
xlim([0 5])

figure
h = errorbar(f_i,nanmean(partial_w3(ca1pyr_c,:)),nanstd(partial_w3(ca1pyr_c,:))/sqrt(length(ca1pyr_c)));
errorbar_tick(h,.01,'units');
hold on
h=errorbar(f_i,nanmean(partial_w3(ca1int_c,:)),nanstd(partial_w3(ca1int_c,:))/sqrt(length(ca1int_c)),'r');
errorbar_tick(h,.01,'units')
h=errorbar(f_i,nanmean(partial_w3(ca3_c,:)),nanstd(partial_w3(ca3_c,:))/sqrt(length(ca3_c)),'k');
errorbar_tick(h,.01,'units')
h=errorbar(f_i,nanmean(partial_w3(dg_c,:)),nanstd(partial_w3(dg_c,:))/sqrt(length(dg_c)),'g');
errorbar_tick(h,.01,'units')
xlim([0 2])

figure
h=errorbar(f_i,nanmean(lSww(ca1pyr,:)),nanstd(lSww(ca1pyr,:))/sqrt(length(ca1pyr)));
errorbar_tick(h,.01,'units');
hold on
h=errorbar(f_i,nanmean(lSww(ca1int,:)),nanstd(lSww(ca1int,:))/sqrt(length(ca1int)),'r');
errorbar_tick(h,.01,'units');
h=errorbar(f_i,nanmean(lSww(ca3,:)),nanstd(lSww(ca3,:))/sqrt(length(ca3)),'k');
errorbar_tick(h,.01,'units');
h=errorbar(f_i,nanmean(lSww(dg,:)),nanstd(lSww(dg,:))/sqrt(length(dg)),'g');
errorbar_tick(h,.01,'units');
xlim([0 2])
xlabel('Frequency (Hz)')
ylabel('Power')
