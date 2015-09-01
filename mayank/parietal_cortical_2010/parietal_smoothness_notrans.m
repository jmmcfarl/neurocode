%inputs:
% raw_data
% raw_Fs
clear all
close all

drive_letter = 'G';

addpath('G:\WC_Germany\hmm_state_detect\\')
addpath('G:\WC_Germany\hsmm_state_detection')
addpath('G:\WC_Germany\parietal_cortical_2010')
% addpath('G:\WC_Germany\new_stellate_analysis\')
addpath('G:\WC_Germany\persistent_2010\')

cd G:\WC_Germany\parietal_cortical_2010\
load parietal_cortical_2010
load desynch_times_individual

raw_Fs = 2016;
dsf = 8;
Fsd = raw_Fs/dsf;
Fs = 50.4;

% parietal = find_struct_field_vals(sess_data,'region','parietal');
% sess_data = sess_data(parietal);
interneurons = find_struct_field_vals(sess_data,'cell_type','interneuron');
sess_data(interneurons) = [];
desynch_times_mp(interneurons) = [];
desynch_times_lf8(interneurons) = [];
desynch_times_lf4(interneurons) = [];

frontal = find_struct_field_vals(sess_data,'region','frontal');
prefrontal = find_struct_field_vals(sess_data,'region','prefrontal');
parietal = find_struct_field_vals(sess_data,'region','parietal');
thom_el = find_struct_field_vals(sess_data,'thom_elec',1);
thom_par = thom_el(find(ismember(thom_el,parietal)));
thom_pfc = setdiff(thom_el,thom_par);

min_state_dur = 1;
max_state_dur = round(Fs*30);
dur_range = (1:max_state_dur)/Fs;

dur_hist_bincenters = (min_state_dur:max_state_dur)/Fs;

hcf_vals = [0.5 1 2 4 10];
smooth_vals = [0.025 0.05 0.1 0.2 0.3];

short_threshold = 0.2;

n = length(sess_data);

%% For LF4 
for d = 1:length(thom_el)
    d
    direct = sess_data(thom_el(d)).directory;
    direct(1) = 'G';
    cd(direct)
    load used_data lf4
    load smoothness_lf4
    for i = 1:length(hcf_vals)
        [hmm4_nt{i}] = parietal_get_hmm_lf_test_notrans(lf4,raw_Fs,hcf_vals(i),desynch_times_lf4{thom_el(d)},hmm_bbstate_seq4{i});
        up_dur_pmf4(d,i,:) = hist(state_durations4{d,i}{2},dur_hist_bincenters);
        up_dur_pmf4(d,i,:) = up_dur_pmf4(d,i,:)/sum(up_dur_pmf4(d,i,:));
        down_dur_pmf4(d,i,:) = hist(state_durations4{d,i}{1},dur_hist_bincenters);
        down_dur_pmf4(d,i,:) = down_dur_pmf4(d,i,:)/sum(down_dur_pmf4(d,i,:));
        min_up_dur4(d,i) = min(state_durations4{d,i}{2});
        min_down_dur4(d,i) = min(state_durations4{d,i}{1});
        mean_up_dur4(d,i) = mean(state_durations4{d,i}{2});
        mean_down_dur4(d,i) = mean(state_durations4{d,i}{1});
        meandiff = [];
        for s = 1:hmm4{i}.Nsegs
            meandiff = [meandiff; hmm4{i}.state(2).meanfun{s}-hmm4{i}.state(1).meanfun{s}];
        end
        %         meandiff = mean(meandiff);
        covar1 = hmm4{i}.state(1).var;
        covar2 = hmm4{i}.state(2).var;
        lf4_kl_lf(d,i) = gauss_kl_div(meandiff,covar1,covar2);
        
    end
    save smoothness_lf4_notrans hmm_bbstate_seq4 hmm4
    clear hmm*
    
    for i = 1:length(smooth_vals)
        [hmm_bbstate_seq4_hf{i},hmm4_hf{i},state_durations_hf4{d,i}] = parietal_get_hmm_hf_test(lf4,raw_Fs,smooth_vals(i),desynch_times_lf4{thom_el(d)});
        up_dur_hf_pmf4(d,i,:) = hist(state_durations_hf4{d,i}{2},dur_hist_bincenters);
        up_dur_hf_pmf4(d,i,:) = up_dur_hf_pmf4(d,i,:)/sum(up_dur_hf_pmf4(d,i,:));
        down_dur_hf_pmf4(d,i,:) = hist(state_durations_hf4{d,i}{1},dur_hist_bincenters);
        down_dur_hf_pmf4(d,i,:) = down_dur_hf_pmf4(d,i,:)/sum(down_dur_hf_pmf4(d,i,:));
        min_up_dur_hf4(d,i) = min(state_durations_hf4{d,i}{2});
        min_down_dur_hf4(d,i) = min(state_durations_hf4{d,i}{1});
        mean_up_dur_hf4(d,i) = mean(state_durations_hf4{d,i}{2});
        mean_down_dur_hf4(d,i) = mean(state_durations_hf4{d,i}{1});
        meandiff = [];
        for s = 1:hmm4_hf{i}.Nsegs
            meandiff = [meandiff; hmm4_hf{i}.state(2).meanfun{s}-hmm4_hf{i}.state(1).meanfun{s}];
        end
        %         meandiff = mean(meandiff);
        covar1 = hmm4_hf{i}.state(1).var;
        covar2 = hmm4_hf{i}.state(2).var;
        lf4_kl_hf(d,i) = gauss_kl_div(meandiff,covar1,covar2);
    end
    save smoothness_lf4_hf hmm_bbstate_seq4_hf hmm4_hf
    clear hmm*
    
    cd F:\WC_Germany\parietal_cortical_2010
    save smoothness_test_results_robust_dp_mp_lf8_lf4
end


%% For MP 
for d = thom_pfc 
    d
    direct = sess_data(d).directory;
    direct(1) = 'G';
    cd(direct)
    load used_data wcv_minus_spike
    for i = 1:length(hcf_vals)
        [hmm_bbstate_seq{i},hmm{i},state_durations{d,i}] = parietal_get_hmm_lf_test(wcv_minus_spike,raw_Fs,hcf_vals(i),desynch_times_mp{d});
        up_dur_pmf(d,i,:) = hist(state_durations{d,i}{2},dur_hist_bincenters);
        up_dur_pmf(d,i,:) = up_dur_pmf(d,i,:)/sum(up_dur_pmf(d,i,:));
        down_dur_pmf(d,i,:) = hist(state_durations{d,i}{1},dur_hist_bincenters);
        down_dur_pmf(d,i,:) = down_dur_pmf(d,i,:)/sum(down_dur_pmf(d,i,:));
        min_up_dur(d,i) = min(state_durations{d,i}{2});
        min_down_dur(d,i) = min(state_durations{d,i}{1});
        mean_up_dur(d,i) = mean(state_durations{d,i}{2});
        mean_down_dur(d,i) = mean(state_durations{d,i}{1});
        meandiff = [];
        for s = 1:hmm{i}.Nsegs
            meandiff = [meandiff; hmm{i}.state(2).meanfun{s}-hmm{i}.state(1).meanfun{s}];
        end
%         meandiff = mean(meandiff);
        covar1 = hmm{i}.state(1).var;
        covar2 = hmm{i}.state(2).var;
        mp_kl_lf(d,i) = gauss_kl_div(meandiff,covar1,covar2);
    end
    save smoothness_mp hmm_bbstate_seq hmm
    clear hmm*
    
    for i = 1:length(smooth_vals)
        [hmm_bbstate_seq_hf{i},hmm_hf{i},state_durations_hf{d,i}] = parietal_get_hmm_hf_test(wcv_minus_spike,raw_Fs,smooth_vals(i),desynch_times_mp{d});
        up_dur_hf_pmf(d,i,:) = hist(state_durations_hf{d,i}{2},dur_hist_bincenters);
        up_dur_hf_pmf(d,i,:) = up_dur_hf_pmf(d,i,:)/sum(up_dur_hf_pmf(d,i,:));
        down_dur_hf_pmf(d,i,:) = hist(state_durations_hf{d,i}{1},dur_hist_bincenters);
        down_dur_hf_pmf(d,i,:) = down_dur_hf_pmf(d,i,:)/sum(down_dur_hf_pmf(d,i,:));
        min_up_dur_hf(d,i) = min(state_durations_hf{d,i}{2});
        min_down_dur_hf(d,i) = min(state_durations_hf{d,i}{1});
        mean_up_dur_hf(d,i) = mean(state_durations_hf{d,i}{2});
        mean_down_dur_hf(d,i) = mean(state_durations_hf{d,i}{1});
        meandiff = [];
        for s = 1:hmm_hf{i}.Nsegs
            meandiff = [meandiff; hmm_hf{i}.state(2).meanfun{s}-hmm_hf{i}.state(1).meanfun{s}];
        end
%         meandiff = mean(meandiff);
        covar1 = hmm_hf{i}.state(1).var;
        covar2 = hmm_hf{i}.state(2).var;
        mp_kl_hf(d,i) = gauss_kl_div(meandiff,covar1,covar2);
    end  
    save smoothness_hf_mp hmm_bbstate_seq_hf hmm_hf
    clear hmm*
    
    cd F:\WC_Germany\parietal_cortical_2010
    save smoothness_test_results_robust_dp_mp_lf8_lf4
end

%%
cd F:\WC_Germany\parietal_cortical_2010
save smoothness_test_results_robust_dp_mp_lf8_lf4

%% FOR LF4
figure
cmap = colormap(jet(length(smooth_vals)));
subplot(2,1,1)
for i = 1:length(smooth_vals)
    ebar_mat(dur_hist_bincenters,up_dur_hf_pmf4(:,i,:),.001,cmap(i,:)), hold on
end
xlim([0 3])
set(gca,'yscale','log')
ylim([2e-4 8e-2])
subplot(2,1,2)
for i = 1:length(smooth_vals)
    ebar_mat(dur_hist_bincenters,down_dur_hf_pmf4(:,i,:),.001,cmap(i,:)), hold on
end
xlim([0 3])
set(gca,'yscale','log')
ylim([2e-4 8e-2])

figure
cmap = colormap(jet(length(hcf_vals)));
subplot(2,1,1)
for i = 1:length(hcf_vals)
    ebar_mat(dur_hist_bincenters,up_dur_pmf4(:,i,:),.001,cmap(i,:)), hold on
end
xlim([0 3])
set(gca,'yscale','log')
ylim([2e-4 5e-2])
subplot(2,1,2)
for i = 1:length(hcf_vals)
    ebar_mat(dur_hist_bincenters,down_dur_pmf4(:,i,:),.001,cmap(i,:)), hold on
end
xlim([0 3])
set(gca,'yscale','log')
ylim([2e-4 5e-2])

figure
errorbar(smooth_vals,mean(min_up_dur_hf4),std(min_up_dur_hf4)/sqrt(n)), hold on
errorbar(smooth_vals,mean(min_down_dur_hf4),std(min_down_dur_hf4)/sqrt(n),'r'), hold on

figure
errorbar(hcf_vals,mean(min_up_dur4),std(min_up_dur4)/sqrt(n)), hold on
errorbar(hcf_vals,mean(min_down_dur4),std(min_down_dur4)/sqrt(n),'r'), hold on

% figure
% errorbar(smooth_vals,mean(mean_up_dur_hf4),std(mean_up_dur_hf4)/sqrt(n)), hold on
% errorbar(smooth_vals,mean(mean_down_dur_hf4),std(mean_down_dur_hf4)/sqrt(n),'r'), hold on
% 
% figure
% errorbar(hcf_vals,mean(mean_up_dur4),std(mean_up_dur4)/sqrt(n)), hold on
% errorbar(hcf_vals,mean(mean_down_dur4),std(mean_down_dur4)/sqrt(n),'r'), hold on

figure
errorbar(smooth_vals,mean(lf4_kl_hf),std(lf4_kl_hf)/sqrt(n)), hold on
figure
errorbar(hcf_vals,mean(lf4_kl_lf),std(lf4_kl_lf)/sqrt(n),'r'), hold on

%% FOR LF8
figure
cmap = colormap(jet(length(smooth_vals)));
subplot(2,1,1)
for i = 1:length(smooth_vals)
    ebar_mat(dur_hist_bincenters,up_dur_hf_pmf8(:,i,:),.001,cmap(i,:)), hold on
end
xlim([0 3])
set(gca,'yscale','log')
ylim([2e-4 8e-2])
subplot(2,1,2)
for i = 1:length(smooth_vals)
    ebar_mat(dur_hist_bincenters,down_dur_hf_pmf8(:,i,:),.001,cmap(i,:)), hold on
end
xlim([0 3])
set(gca,'yscale','log')
ylim([2e-4 8e-2])

figure
cmap = colormap(jet(length(hcf_vals)));
subplot(2,1,1)
for i = 1:length(hcf_vals)
    ebar_mat(dur_hist_bincenters,up_dur_pmf8(:,i,:),.001,cmap(i,:)), hold on
end
xlim([0 3])
set(gca,'yscale','log')
ylim([2e-4 5e-2])
subplot(2,1,2)
for i = 1:length(hcf_vals)
    ebar_mat(dur_hist_bincenters,down_dur_pmf8(:,i,:),.001,cmap(i,:)), hold on
end
xlim([0 3])
set(gca,'yscale','log')
ylim([2e-4 5e-2])

figure
errorbar(smooth_vals,mean(min_up_dur_hf8),std(min_up_dur_hf8)/sqrt(n)), hold on
errorbar(smooth_vals,mean(min_down_dur_hf8),std(min_down_dur_hf8)/sqrt(n),'r'), hold on

figure
errorbar(hcf_vals,mean(min_up_dur8),std(min_up_dur8)/sqrt(n)), hold on
errorbar(hcf_vals,mean(min_down_dur8),std(min_down_dur8)/sqrt(n),'r'), hold on
% 
% figure
% errorbar(smooth_vals,mean(mean_up_dur_hf8),std(mean_up_dur_hf8)/sqrt(n)), hold on
% errorbar(smooth_vals,mean(mean_down_dur_hf8),std(mean_down_dur_hf8)/sqrt(n),'r'), hold on

% figure
% errorbar(hcf_vals,mean(mean_up_dur8),std(mean_up_dur8)/sqrt(n)), hold on
% errorbar(hcf_vals,mean(mean_down_dur8),std(mean_down_dur8)/sqrt(n),'r'), hold on

figure
errorbar(smooth_vals,mean(lf8_kl_hf),std(lf8_kl_hf)/sqrt(n)), hold on
figure
errorbar(hcf_vals,mean(lf8_kl_lf),std(lf8_kl_lf)/sqrt(n),'r'), hold on

%% FOR MP
figure
cmap = colormap(jet(length(smooth_vals)));
subplot(2,1,1)
for i = 1:length(smooth_vals)
    ebar_mat(dur_hist_bincenters,up_dur_hf_pmf(thom_pfc,i,:),.001,cmap(i,:)), hold on
end
xlim([0 3])
set(gca,'yscale','log')
ylim([2e-4 8e-2])
subplot(2,1,2)
for i = 1:length(smooth_vals)
    ebar_mat(dur_hist_bincenters,down_dur_hf_pmf(thom_pfc,i,:),.001,cmap(i,:)), hold on
end
xlim([0 3])
set(gca,'yscale','log')
ylim([2e-4 8e-2])

figure
cmap = colormap(jet(length(hcf_vals)));
subplot(2,1,1)
for i = 1:length(hcf_vals)
    ebar_mat(dur_hist_bincenters,up_dur_pmf(thom_pfc,i,:),.001,cmap(i,:)), hold on
end
xlim([0 3])
set(gca,'yscale','log')
ylim([2e-4 5e-2])
subplot(2,1,2)
for i = 1:length(hcf_vals)
    ebar_mat(dur_hist_bincenters,down_dur_pmf(thom_pfc,i,:),.001,cmap(i,:)), hold on
end
xlim([0 3])
set(gca,'yscale','log')
ylim([2e-4 5e-2])

figure
errorbar(smooth_vals,mean(min_up_dur_hf(thom_pfc,:)),std(min_up_dur_hf(thom_pfc,:))/sqrt(9)), hold on
errorbar(smooth_vals,mean(min_down_dur_hf(thom_pfc,:)),std(min_down_dur_hf(thom_pfc,:))/sqrt(9),'r'), hold on

figure
errorbar(hcf_vals,mean(min_up_dur(thom_pfc,:)),std(min_up_dur(thom_pfc,:))/sqrt(9)), hold on
errorbar(hcf_vals,mean(min_down_dur(thom_pfc,:)),std(min_down_dur(thom_pfc,:))/sqrt(9),'r'), hold on

figure
errorbar(smooth_vals,mean(mean_up_dur_hf),std(mean_up_dur_hf)/sqrt(n)), hold on
errorbar(smooth_vals,mean(mean_down_dur_hf),std(mean_down_dur_hf)/sqrt(n),'r'), hold on

figure
errorbar(hcf_vals,mean(mean_up_dur),std(mean_up_dur)/sqrt(n)), hold on
errorbar(hcf_vals,mean(mean_down_dur),std(mean_down_dur)/sqrt(n),'r'), hold on

figure
errorbar(smooth_vals,mean(mp_kl_hf),std(mp_kl_hf)/sqrt(n)), hold on
figure
errorbar(hcf_vals,mean(mp_kl_lf),std(mp_kl_lf)/sqrt(n),'r'), hold on

