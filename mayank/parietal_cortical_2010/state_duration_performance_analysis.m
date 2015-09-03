clear all
close all

addpath('G:\WC_Germany\parietal_cortical_2010\')
addpath('G:\WC_Germany\hsmm_state_detection\')
addpath('G:\WC_Germany\hsmm_uds_code\')

cd G:\WC_Germany\parietal_cortical_2010\
load parietal_cortical_2010
load G:\WC_Germany\parietal_cortical_2010\desynch_times_mp_lf8


frontal = find_struct_field_vals(sess_data,'region','frontal');
prefrontal = find_struct_field_vals(sess_data,'region','prefrontal');
parietal = find_struct_field_vals(sess_data,'region','parietal');

% sess_data = sess_data(parietal);
% desynch_times = desynch_times(parietal);

%get rid of interneurons
interneurons = find_struct_field_vals(sess_data,'cell_type','interneuron');
sess_data(interneurons) = [];

raw_Fs = 2016;
dsf = 8;
Fsd = raw_Fs/dsf;

frontal = find_struct_field_vals(sess_data,'region','frontal');
prefrontal = find_struct_field_vals(sess_data,'region','prefrontal');
parietal = find_struct_field_vals(sess_data,'region','parietal');
thom_el = find_struct_field_vals(sess_data,'thom_elec',1);
thom_par = thom_el(find(ismember(thom_el,parietal)));
thom_pfc = setdiff(thom_el,thom_par);
sess_data = sess_data(thom_el);

n = length(sess_data);

min_state_dur = 0;
max_state_dur = 30;
dur_range = min_state_dur:0.05:max_state_dur;

for d = 1:n
    
        cdir = sess_data(d).directory;
    cdir(1) = 'G';
    cd(cdir)
    pwd
    
    load ./hsmm_state_seq4_seg_lf_4_5_2011
    load ./fm_state_seq4_lf_4_11_2011

    hmm_state_durations4 = hsmm_uds_compute_state_durations(hmm_bbstate_seq4,Fsd);
    hsmm_state_durations4 = hsmm_uds_compute_state_durations(hsmm_bbstate_seq4,Fsd);
    fm_state_durations4 = hsmm_uds_compute_state_durations(fm_state_seq4,Fsd);
    fmz_state_durations4 = hsmm_uds_compute_state_durations(fm_state_seqz4,Fsd);
    for i = 1:hsmm4.K
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
    
    shortest_hmm4(d,1) = min(hmm_state_durations4{1});
    shortest_hmm4(d,2) = min(hmm_state_durations4{2});
    
end

%%
short_thresh = 0.2;
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
g = [0 0.8 0];
figure
set(gca,'fontname','arial','fontsize',14)
h = errorbar(dur_range,squeeze(nanmean(hmm_emp_pmf4(:,1,:))),squeeze(nanstd(hmm_emp_pmf4(:,1,:)))/3);
errorbar_tick(h,.001,'units')
hold on
h = errorbar(dur_range,squeeze(nanmean(hsmm_emp_pmf4(:,1,:))),squeeze(nanstd(hsmm_emp_pmf4(:,1,:)))/3,'r');
errorbar_tick(h,.001,'units')
h = errorbar(dur_range,squeeze(nanmean(fm_emp_pmf4(:,1,:))),squeeze(nanstd(fm_emp_pmf4(:,1,:)))/3,'k');
errorbar_tick(h,.001,'units')
h = errorbar(dur_range,squeeze(nanmean(fmz_emp_pmf4(:,1,:))),squeeze(nanstd(fmz_emp_pmf4(:,1,:)))/3,'color',g);
errorbar_tick(h,.001,'units')
xlim([0.04 20])
ylim([1e-4 1e-1])
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('Duration (s)','fontname','arial','fontsize',16)
ylabel('Probability','fontname','arial','fontsize',16)
legend('HMM','HSMM','SMM-TC','Np-TC')

figure
h = errorbar(dur_range,squeeze(nanmean(hmm_emp_pmf4(:,2,:))),squeeze(nanstd(hmm_emp_pmf4(:,2,:)))/3);
errorbar_tick(h,.001,'units')
hold on
h = errorbar(dur_range,squeeze(nanmean(hsmm_emp_pmf4(:,2,:))),squeeze(nanstd(hsmm_emp_pmf4(:,2,:)))/3,'r');
errorbar_tick(h,.001,'units')
h = errorbar(dur_range,squeeze(nanmean(fm_emp_pmf4(:,2,:))),squeeze(nanstd(fm_emp_pmf4(:,2,:)))/3,'k');
errorbar_tick(h,.001,'units')
h = errorbar(dur_range,squeeze(nanmean(fmz_emp_pmf4(:,2,:))),squeeze(nanstd(fmz_emp_pmf4(:,2,:)))/3,'color',g);
errorbar_tick(h,.001,'units')
xlim([0.04 20])
ylim([1e-4 1e-1])
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('Duration (s)','fontname','arial','fontsize',16)
ylabel('Probability','fontname','arial','fontsize',16)
legend('HMM','HSMM','SMM-TC','Np-TC')

figure
plot(n_short_hsmm4(:,1),n_short_hsmm4(:,2),'o')
hold on
plot(n_short_hmm4(:,1),n_short_hmm4(:,2),'ro')
plot(n_short_fm4(:,1),n_short_fm4(:,2),'ko')
plot(n_short_fmz4(:,1),n_short_fmz4(:,2),'o','color',g)
xlabel('Fraction short DOWNs','fontname','arial','fontsize',16)
ylabel('Fraction short UPs','fontname','arial','fontsize',16)
legend('HSMM','HMM','SMM-TC','Np-TC')
X = [nanmean(n_short_hsmm4(:,1)); nanstd(n_short_hsmm4(:,1))/sqrt(21)];
Y = [nanmean(n_short_hsmm4(:,2)); nanstd(n_short_hsmm4(:,2))/sqrt(21)];
eplot(X',Y','b',2)
X = [nanmean(n_short_hmm4(:,1)); nanstd(n_short_hmm4(:,1))/sqrt(21)];
Y = [nanmean(n_short_hmm4(:,2)); nanstd(n_short_hmm4(:,2))/sqrt(21)];
eplot(X',Y','r',2)
X = [nanmean(n_short_fm4(:,1)); nanstd(n_short_fm4(:,1))/sqrt(21)];
Y = [nanmean(n_short_fm4(:,2)); nanstd(n_short_fm4(:,2))/sqrt(21)];
eplot(X',Y','k',2)
X = [nanmean(n_short_fmz4(:,1)); nanstd(n_short_fmz4(:,1))/sqrt(21)];
Y = [nanmean(n_short_fmz4(:,2)); nanstd(n_short_fmz4(:,2))/sqrt(21)];
eplot(X',Y',g,2)
