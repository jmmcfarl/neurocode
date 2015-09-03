clear all
close all

addpath('G:\WC_Germany\parietal_cortical_2010\')
addpath('G:\WC_Germany\hsmm_state_detection\')
addpath('G:\WC_Germany\hsmm_uds_code\')

cd G:\WC_Germany\parietal_cortical_2010\
load parietal_cortical_2010
load G:\WC_Germany\parietal_cortical_2010\desynch_times_individual

%get rid of interneurons
interneurons = find_struct_field_vals(sess_data,'cell_type','interneuron');
sess_data(interneurons) = [];

frontal = find_struct_field_vals(sess_data,'region','frontal');
prefrontal = find_struct_field_vals(sess_data,'region','prefrontal');
parietal = find_struct_field_vals(sess_data,'region','parietal');
thom_el = find_struct_field_vals(sess_data,'thom_elec',1);
thom_par = thom_el(find(ismember(thom_el,parietal)));
thom_pfc = setdiff(thom_el,thom_par);
parietal = find_struct_field_vals(sess_data,'region','parietal');
frontal = setdiff(1:length(sess_data),parietal);
superficial = find_struct_field_vals(sess_data,'layer','23');
deep = setdiff(1:length(sess_data),superficial);

sess_data = sess_data(thom_el);

n = length(sess_data);
for d = 1:n
%     d=2; %d=3
    cdir = sess_data(d).directory;
    cdir(1) = 'G';
    cd(cdir)
    pwd
    load ./hsmm_state_seq4_seg_lf_3_31_2011_v3
    load ./fm_state_seq4_lf_3_31_2011_v3
    
    %%
    Fs_bb = 252;
    
    min_state_dur = 0;
    max_state_dur = 30;
    dur_range = min_state_dur:0.025:max_state_dur;
    
    hmm_state_durations4 = hsmm_uds_compute_state_durations(hmm_bbstate_seq4,Fs_bb);
    hsmm_state_durations4 = hsmm_uds_compute_state_durations(hsmm_bbstate_seq4,Fs_bb);
    fm_state_durations4 = hsmm_uds_compute_state_durations(fm_state_seq4,Fs_bb);
    fmz_state_durations4 = hsmm_uds_compute_state_durations(fm_state_seqz4,Fs_bb);
    for i = 1:hmm4.K
        %estimate the empirical state duration pmf
        hmm_emp_pmf4(d,i,:) = hist(hmm_state_durations4{i},dur_range);
        hmm_emp_pmf4(d,i,:) = hmm_emp_pmf4(d,i,:)/sum(hmm_emp_pmf4(d,i,:));       
        hsmm_emp_pmf4(d,i,:) = hist(hsmm_state_durations4{i},dur_range);
        hsmm_emp_pmf4(d,i,:) = hsmm_emp_pmf4(d,i,:)/sum(hsmm_emp_pmf4(d,i,:));       
        fm_emp_pmf4(d,i,:) = hist(fm_state_durations4{i},dur_range);
        fm_emp_pmf4(d,i,:) = fm_emp_pmf4(d,i,:)/sum(fm_emp_pmf4(d,i,:));        
        fmz_emp_pmf4(d,i,:) = hist(fmz_state_durations4{i},dur_range);
        fmz_emp_pmf4(d,i,:) = fmz_emp_pmf4(d,i,:)/sum(fmz_emp_pmf4(d,i,:));
    end
    
%     
%      load ./hsmm_state_seq_seg_lf_3_31_2011_v2
%     load ./fm_state_seq_lf_3_31_2011_v2
%       hmm_state_durations = hsmm_uds_compute_state_durations(hmm_bbstate_seq,Fs_bb);
%     hsmm_state_durations = hsmm_uds_compute_state_durations(hsmm_bbstate_seq,Fs_bb);
%     fm_state_durations = hsmm_uds_compute_state_durations(fm_state_seq,Fs_bb);
%     fmz_state_durations = hsmm_uds_compute_state_durations(fm_state_seqz,Fs_bb);
%     for i = 1:hmm.K
%         %estimate the empirical state duration pmf
%         hmm_emp_pmf(d,i,:) = hist(hmm_state_durations{i},dur_range);
%         hmm_emp_pmf(d,i,:) = hmm_emp_pmf(d,i,:)/sum(hmm_emp_pmf(d,i,:));       
%         hsmm_emp_pmf(d,i,:) = hist(hsmm_state_durations{i},dur_range);
%         hsmm_emp_pmf(d,i,:) = hsmm_emp_pmf(d,i,:)/sum(hsmm_emp_pmf(d,i,:));       
%         fm_emp_pmf(d,i,:) = hist(fm_state_durations{i},dur_range);
%         fm_emp_pmf(d,i,:) = fm_emp_pmf(d,i,:)/sum(fm_emp_pmf(d,i,:));        
%         fmz_emp_pmf(d,i,:) = hist(fmz_state_durations{i},dur_range);
%         fmz_emp_pmf(d,i,:) = fmz_emp_pmf(d,i,:)/sum(fmz_emp_pmf(d,i,:));
%     end
 
end

%%
short_thresh = 0.15;
for d = 1:n
    temp = find(dur_range > short_thresh,1,'first');
    n_short_hmm4(d,1) = sum(hmm_emp_pmf4(d,1,1:temp),3);
    n_short_hmm4(d,2) = sum(hmm_emp_pmf4(d,2,1:temp),3);  
    n_short_hsmm4(d,1) = sum(hsmm_emp_pmf4(d,1,1:temp),3);
    n_short_hsmm4(d,2) = sum(hsmm_emp_pmf4(d,2,1:temp),3);
    n_short_fm4(d,1) = sum(fm_emp_pmf4(d,1,1:temp),3);
    n_short_fm4(d,2) = sum(fm_emp_pmf4(d,2,1:temp),3);
    n_short_fmz4(d,1) = sum(fmz_emp_pmf4(d,1,1:temp),3);
    n_short_fmz4(d,2) = sum(fmz_emp_pmf4(d,2,1:temp),3);
end
    
%%
figure
h = errorbar(dur_range,squeeze(nanmean(hmm_emp_pmf4(:,1,:))),squeeze(nanstd(hmm_emp_pmf4(:,1,:)))/sqrt(21));
errorbar_tick(h,.001,'units')
hold on
h = errorbar(dur_range,squeeze(nanmean(hsmm_emp_pmf4(:,1,:))),squeeze(nanstd(hsmm_emp_pmf4(:,1,:)))/sqrt(21),'r');
errorbar_tick(h,.001,'units')
h = errorbar(dur_range,squeeze(nanmean(fm_emp_pmf4(:,1,:))),squeeze(nanstd(fm_emp_pmf4(:,1,:)))/sqrt(21),'k');
errorbar_tick(h,.001,'units')
h = errorbar(dur_range,squeeze(nanmean(fmz_emp_pmf4(:,1,:))),squeeze(nanstd(fmz_emp_pmf4(:,1,:)))/sqrt(21),'g');
errorbar_tick(h,.001,'units')
xlim([0.02 10])
ylim([5e-5 5e-2])
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('DOWN state duration (s)','fontsize',14)
ylabel('Relative frequency','fontsize',14)
legend('HMM','HSMM','SMM-TC','Np-TC')

figure
h = errorbar(dur_range,squeeze(nanmean(hmm_emp_pmf4(:,2,:))),squeeze(nanstd(hmm_emp_pmf4(:,2,:)))/sqrt(21));
errorbar_tick(h,.001,'units')
hold on
h = errorbar(dur_range,squeeze(nanmean(hsmm_emp_pmf4(:,2,:))),squeeze(nanstd(hsmm_emp_pmf4(:,2,:)))/sqrt(21),'r');
errorbar_tick(h,.001,'units')
h = errorbar(dur_range,squeeze(nanmean(fm_emp_pmf4(:,2,:))),squeeze(nanstd(fm_emp_pmf4(:,2,:)))/sqrt(21),'k');
errorbar_tick(h,.001,'units')
h = errorbar(dur_range,squeeze(nanmean(fmz_emp_pmf4(:,2,:))),squeeze(nanstd(fmz_emp_pmf4(:,2,:)))/sqrt(21),'g');
errorbar_tick(h,.001,'units')
xlim([0.02 4])
ylim([5e-5 7e-2])
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('UP state duration (s)','fontsize',14)
ylabel('Relative frequency','fontsize',14)
legend('HMM','HSMM','SMM-TC','Np-TC')

%%
% figure
% h = errorbar(dur_range,squeeze(nanmean(hmm_emp_pmf(:,1,:))),squeeze(nanstd(hmm_emp_pmf(:,1,:)))/3);
% errorbar_tick(h,.001,'units')
% hold on
% h = errorbar(dur_range,squeeze(nanmean(hsmm_emp_pmf(:,1,:))),squeeze(nanstd(hsmm_emp_pmf(:,1,:)))/3,'r');
% errorbar_tick(h,.001,'units')
% h = errorbar(dur_range,squeeze(nanmean(fm_emp_pmf(:,1,:))),squeeze(nanstd(fm_emp_pmf(:,1,:)))/3,'k');
% errorbar_tick(h,.001,'units')
% h = errorbar(dur_range,squeeze(nanmean(fmz_emp_pmf(:,1,:))),squeeze(nanstd(fmz_emp_pmf(:,1,:)))/3,'g');
% errorbar_tick(h,.001,'units')
% xlim([0.02 10])
% ylim([5e-5 5e-2])
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% xlabel('DOWN state duration (s)','fontsize',14)
% ylabel('Relative frequency','fontsize',14)
% legend('HMM','HSMM','SMM-TC','Np-TC')
% 
% figure
% h = errorbar(dur_range,squeeze(nanmean(hmm_emp_pmf(:,2,:))),squeeze(nanstd(hmm_emp_pmf(:,2,:)))/3);
% errorbar_tick(h,.001,'units')
% hold on
% h = errorbar(dur_range,squeeze(nanmean(hsmm_emp_pmf(:,2,:))),squeeze(nanstd(hsmm_emp_pmf(:,2,:)))/3,'r');
% errorbar_tick(h,.001,'units')
% h = errorbar(dur_range,squeeze(nanmean(fm_emp_pmf(:,2,:))),squeeze(nanstd(fm_emp_pmf(:,2,:)))/3,'k');
% errorbar_tick(h,.001,'units')
% h = errorbar(dur_range,squeeze(nanmean(fmz_emp_pmf(:,2,:))),squeeze(nanstd(fmz_emp_pmf(:,2,:)))/3,'g');
% errorbar_tick(h,.001,'units')
% xlim([0.02 4])
% ylim([5e-5 7e-2])
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% xlabel('UP state duration (s)','fontsize',14)
% ylabel('Relative frequency','fontsize',14)
% legend('HMM','HSMM','SMM-TC','Np-TC')
% 
