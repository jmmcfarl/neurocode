clear all
close all
%%
load E:\WC_Germany\overall_EC\overall_allcells_dir
addpath('E:\Code\WC_anal\general\')
addpath('E:\WC_Germany\Overall_EC\')
addpath('E:\Code\Chronux\spectral_analysis\continuous\')
addpath('E:\WC_Germany\hsmm_state_detection\\')

drive_letter = 'E';

load overall_allcells_coherence_data_gains
load overall_allcells_coherence_clustered

%% restrict analysis to good hippocampal LFPs
sess_data_hipp = sess_data(good_hipp_lf2r);

partial_w8 = partial_w8(good_hipp_lf2r,:);
% partial_w2 = partial_w2(used_hipp_lfp,:);
partial_wh = partial_w2r(good_hipp_lf2r,:);
% partial_w3 = partial_w3(used_hipp_lfp,:);
% P83 = P83(good_hipp_lf2r,:);
% P82 = P82(used_hipp_lfp,:);
P8h = P82r(good_hipp_lf2r,:);
% Pw8 = Pw8(good_hipp_lf2r,:);
% Pw3 = Pw3(used_hipp_lfp,:);
% Pw2 = Pw2(used_hipp_lfp,:);
Pwh = Pw2r(good_hipp_lf2r,:);
% lPww = lPww(good_hipp_lf2r,:);
% lP88 = lP88(good_hipp_lf2r,:);
% lP33 = lP33(used_hipp_lfp,:);
% lP22 = lP22(used_hipp_lfp,:);
lPhh = lP22r(good_hipp_lf2r,:);

aC8h_cor = aC82r_cor(good_hipp_lf2r,:);
aCw8_cor = aCw8_cor(good_hipp_lf2r,:);
aCwh_cor = aCw2r_cor(good_hipp_lf2r,:);

used_f = find(f_i < 1);
[max_aC8h,max_aC8h_loc] = max(aC8h_cor(:,used_f),[],2);
[max_aCw8,max_aCw8_loc] = max(aCw8_cor(:,used_f),[],2);
[max_aCwh,max_aCwh_loc] = max(aCwh_cor(:,used_f),[],2);

%%
mec = find_struct_field_vals(sess_data_hipp,'region','MEC');
layer3 = find_struct_field_vals(sess_data_hipp,'layer','3');
layer2 = find_struct_field_vals(sess_data_hipp,'layer','2');
layer23 = find_struct_field_vals(sess_data_hipp,'layer','23');
layer5 = find_struct_field_vals(sess_data_hipp,'layer','56');
lec = find_struct_field_vals(sess_data_hipp,'region','LEC');
l3mec = intersect(mec,layer3);
l2mec = intersect(mec,layer2);
l3lec = intersect(lec,layer3);
l3mec(24:end) = [];
l5mec = intersect(mec,layer5);
l23mec = intersect(mec,layer23);
l23mec = unique([l2mec l3mec l23mec]);
ca1 = find_struct_field_vals(sess_data_hipp,'region','ca1');
ca3 = find_struct_field_vals(sess_data_hipp,'region','ca3');
dg = find_struct_field_vals(sess_data_hipp,'region','dg');
pyr = find_struct_field_vals(sess_data_hipp,'cell_type','pyr');
ca1p = intersect(ca1,pyr);
ca1i = setdiff(ca1,pyr);

%%
figure
h = errorbar(f_i,nanmean(partial_wh(l3mec,:)),nanstd(partial_wh(l3mec,:))/sqrt(length(l3mec)));
errorbar_tick(h,.01,'units');
hold on
h=errorbar(f_i,nanmean(partial_wh(l3lec,:)),nanstd(partial_wh(l3lec,:))/sqrt(length(l3lec)),'r');
errorbar_tick(h,.01,'units')
xlim([0 1.5])
title('Partial MP-hipp LFP Coh')
xlabel('Frequency (Hz)','fontsize',14)
ylabel('Partial coherence','fontsize',14)
legend('L3MEC','L3LEC')

figure
h = errorbar(f_i,nanmean(partial_wh(l3mec,:)),nanstd(partial_wh(l3mec,:))/sqrt(length(l3mec)));
errorbar_tick(h,.01,'units');
hold on
h=errorbar(f_i,nanmean(partial_wh(l2mec,:)),nanstd(partial_wh(l2mec,:))/sqrt(length(l2mec)),'color',[0.7 0.1 0.2]);
errorbar_tick(h,.01,'units')
h=errorbar(f_i,nanmean(partial_wh(l5mec,:)),nanstd(partial_wh(l5mec,:))/sqrt(length(l5mec)),'color',[0.2 0.5 0.4]);
errorbar_tick(h,.01,'units')
h=errorbar(f_i,nanmean(partial_wh(l3lec,:)),nanstd(partial_wh(l3lec,:))/sqrt(length(l3lec)),'r');
errorbar_tick(h,.01,'units')
h=errorbar(f_i,nanmean(partial_wh(ca1p,:)),nanstd(partial_wh(ca1p,:))/sqrt(length(ca1p)),'k');
errorbar_tick(h,.01,'units')
h=errorbar(f_i,nanmean(partial_wh(ca1i,:)),nanstd(partial_wh(ca1i,:))/sqrt(length(ca1i)),'g');
errorbar_tick(h,.01,'units')
h=errorbar(f_i,nanmean(partial_wh(dg,:)),nanstd(partial_wh(dg,:))/sqrt(length(dg)),'c');
errorbar_tick(h,.01,'units')
h=errorbar(f_i,nanmean(partial_wh(ca3,:)),nanstd(partial_wh(ca3,:))/sqrt(length(ca3)),'y');
errorbar_tick(h,.01,'units')
xlim([0 1.5])
title('Partial MP-hipp LFP Coh')
legend('L3MEC','L2MEC','L5MEC','LEC','CA1p','CA1i','DG','Ca3')
xlabel('Frequency (Hz)','fontsize',14)
ylabel('Partial coherence','fontsize',14)

figure
h = errorbar(f_i,nanmean(lPww(l3mec,:)),nanstd(lPww(l3mec,:))/sqrt(length(l3mec))); 
errorbar_tick(h,.01,'units');
hold on
h = errorbar(f_i,nanmean(lPww(l3lec,:)),nanstd(lPww(l3mec,:))/sqrt(length(l3lec)),'g'); 
errorbar_tick(h,.01,'units');
xlim([0 1.5])
title('MP Power spectra')
xlabel('Frequency (Hz)','fontsize',14)
ylabel('Power (AU)','fontsize',14)
legend('L3MEC','L3LEC')

figure
h = errorbar(f_i,nanmean(lP88),nanstd(lP88)/sqrt(size(lP88,1)),'r'); 
errorbar_tick(h,.01,'units');
hold on
h = errorbar(f_i,nanmean(lPhh),nanstd(lPhh)/sqrt(size(lPhh,1)),'k'); 
errorbar_tick(h,.01,'units');
xlim([0 1.5])
title('LFP Power spectra')
xlabel('Frequency (Hz)','fontsize',14)
ylabel('Power (AU)','fontsize',14)
legend('Cortical LFP','Hippocampal LFP')


for i = 1:length(f_i)
    mavg_Pwh(i) = circ_mean(Pwh(l3mec,i));
    mstd_Pwh(i) = circ_std(Pwh(l3mec,i));
    lavg_Pwh(i) = circ_mean(Pwh(l3lec,i));
    lstd_Pwh(i) = circ_std(Pwh(l3lec,i));
end

%%
mec = find_struct_field_vals(sess_data,'region','MEC');
layer3 = find_struct_field_vals(sess_data,'layer','3');
layer2 = find_struct_field_vals(sess_data,'layer','2');
layer23 = find_struct_field_vals(sess_data,'layer','23');
layer5 = find_struct_field_vals(sess_data,'layer','56');
lec = find_struct_field_vals(sess_data,'region','LEC');
l3mec = intersect(mec,layer3);
l2mec = intersect(mec,layer2);
l3lec = intersect(lec,layer3);
l3mec(24:end) = [];
l5mec = intersect(mec,layer5);
l23mec = intersect(mec,layer23);
l23mec = unique([l2mec l3mec l23mec]);
ca1 = find_struct_field_vals(sess_data,'region','ca1');
ca3 = find_struct_field_vals(sess_data,'region','ca3');
dg = find_struct_field_vals(sess_data,'region','dg');
pyr = find_struct_field_vals(sess_data,'cell_type','pyr');
ca1p = intersect(ca1,pyr);
ca1i = setdiff(ca1,pyr);

%%

figure
h = errorbar(f_i,nanmean(partial_w8(l3mec,:)),nanstd(partial_w8(l3mec,:))/sqrt(length(l3mec)));
errorbar_tick(h,.01,'units');
hold on
h=errorbar(f_i,nanmean(partial_w8(l2mec,:)),nanstd(partial_w8(l2mec,:))/sqrt(length(l2mec)),'color',[0.7 0.1 0.2]);
errorbar_tick(h,.01,'units')
h=errorbar(f_i,nanmean(partial_w8(l5mec,:)),nanstd(partial_w8(l5mec,:))/sqrt(length(l5mec)),'color',[0.2 0.5 0.4]);
errorbar_tick(h,.01,'units')
h=errorbar(f_i,nanmean(partial_w8(l3lec,:)),nanstd(partial_w8(l3lec,:))/sqrt(length(l3lec)),'r');
errorbar_tick(h,.01,'units')
h=errorbar(f_i,nanmean(partial_w8(ca1p,:)),nanstd(partial_w8(ca1p,:))/sqrt(length(ca1p)),'k');
errorbar_tick(h,.01,'units')
h=errorbar(f_i,nanmean(partial_w8(ca1i,:)),nanstd(partial_w8(ca1i,:))/sqrt(length(ca1i)),'g');
errorbar_tick(h,.01,'units')
h=errorbar(f_i,nanmean(partial_w8(dg,:)),nanstd(partial_w8(dg,:))/sqrt(length(dg)),'c');
errorbar_tick(h,.01,'units')
h=errorbar(f_i,nanmean(partial_w8(ca3,:)),nanstd(partial_w8(ca3,:))/sqrt(length(ca3)),'y');
errorbar_tick(h,.01,'units')
xlim([0 1.5])
title('Partial MP-Cort LFP Coh')
legend('L3MEC','L2MEC','L5MEC','LEC','CA1p','CA1i','DG','Ca3')
xlabel('Frequency (Hz)','fontsize',14)
ylabel('Partial coherence','fontsize',14)


figure
h = errorbar(f_i,nanmean(partial_w8(l3mec,:)),nanstd(partial_w8(l3mec,:))/sqrt(length(l3mec)));
errorbar_tick(h,.01,'units');
hold on
h=errorbar(f_i,nanmean(partial_w8(l3lec,:)),nanstd(partial_w8(l3lec,:))/sqrt(length(l3lec)),'r');
errorbar_tick(h,.01,'units')
xlim([0 1.5])
title('Partial MP-cortical LFP Coh')
xlabel('Frequency (Hz)','fontsize',14)
ylabel('Partial coherence','fontsize',14)
legend('L3MEC','L3LEC')


for i = 1:length(f_i)
    mavg_Pw8(i) = circ_mean(Pw8(l3mec,i));
    mstd_Pw8(i) = circ_std(Pw8(l3mec,i));
    lavg_Pw8(i) = circ_mean(Pw8(l3lec,i));
    lstd_Pw8(i) = circ_std(Pw8(l3lec,i));
end

%%
figure
h = errorbar(f_i,mavg_Pw8,mstd_Pw8/sqrt(length(l3mec)));
errorbar_tick(h,.01,'units');
hold on
h=errorbar(f_i,lavg_Pw8,lstd_Pw8/sqrt(length(l3lec)),'r');
errorbar_tick(h,.01,'units')
xlim([0 1.5])
title('MP-CORT LFP Relative Phase')

figure
h = errorbar(f_i,mavg_Pwh,mstd_Pwh/sqrt(length(l3mec)));
errorbar_tick(h,.01,'units');
hold on
h=errorbar(f_i,lavg_Pwh,lstd_Pwh/sqrt(length(l3lec)),'r');
errorbar_tick(h,.01,'units')
xlim([0 1.5])
title('MP-HIPP LFP Relative Phase')
