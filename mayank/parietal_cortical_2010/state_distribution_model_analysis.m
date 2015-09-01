

clear all
% close all

addpath('G:\WC_Germany\parietal_cortical_2010\')
addpath('G:\WC_Germany\hsmm_state_detection\')

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

% dur_range = linspace(0,30,301);
% Fs_dur = 1/(dur_range(2)-dur_range(1));
% dur_range = dur_range + 1/Fs_dur/2;
% dur_range(end) = [];
%%
animal_id = [9 7 10 6 2 3 4 9 7 10 1 2 3 4 5 6 7 8 9 11 1];

for d = 1:n
    cdir = sess_data(d).directory;
    cdir(1) = 'G';
    cd(cdir)
    pwd
 
%         load fixmean_state_seqs_4_26_10_v2
%     [state_durationsf_mp] = compute_state_durations_seg(fixmean_state_seq,Fsd);
% 
if ismember(thom_el(d),thom_pfc)
%     load ./hsmm_state_seq_seg_lf_4_28_10_v3
    load ./hsmm_state_seq_seg_lf_4_5_2011
dur_range = hmm.dur_range;
dur_diff = dur_range(2) - dur_range(1);
cdur_range = [0 dur_range+dur_diff];
Fs_dur = 1/(dur_range(2)-dur_range(1));
    [state_durations_mp] = compute_state_durations_seg(hmm_bbstate_seq,Fsd);
    avg_up(d) = nanmean(state_durations_mp{2});
    avg_down(d) = nanmean(state_durations_mp{1});
    max_up(d) = nanmax(state_durations_mp{2});
    max_down(d) = nanmax(state_durations_mp{1});
    n_down_states = length(state_durations_mp{1});
    n_up_states = length(state_durations_mp{2});
    
    %% for down state
    [alpha,beta] = gamma_mlfit(state_durations_mp{1});
    gam_pmf(d,1,:) = gamma_pmf(dur_range,alpha,beta);
%     gam_cdf(d,1,:) = gammainc(beta*dur_range,alpha);
    
    [mu,lambda] = inverse_gauss_mlfit(state_durations_mp{1});
    ig_pmf(d,1,:) = inverse_gaussian_pmf(dur_range,mu,lambda);
    
%     ig_cdf(d,1,:) = normcdf(sqrt(lambda./hmm.dur_range).*(hmm.dur_range/mu-1),0,1)...
%         +exp(2*lambda/mu)*normcdf(-sqrt(lambda./hmm.dur_range).*(hmm.dur_range/mu+1),0,1);
    
%     lambda = 1/nanmean(state_durations_mp{1});
%     exp_pmf(d,1,:) = exponential_pmf(hmm.dur_range,lambda);
%     exp_cdf(d,1,:) = 1-exp(-lambda*hmm.dur_range);
    p_fit = 1/nanmean(state_durations_mp{1}*Fs_dur);
    geo_pmf(d,1,:) = geometric_pmf(1:length(dur_range),p_fit);
    
    emp_hist(d,1,:) = hist(state_durations_mp{1},dur_range);
    emp_hist(d,1,end) = 0;
    emp_pmf(d,1,:) = emp_hist(d,1,:)/sum(emp_hist(d,1,:));
    
    ks_gam(d,1) = ks_pmf(gam_pmf(d,1,:),emp_pmf(d,1,:));
    ks_ig(d,1) = ks_pmf(ig_pmf(d,1,:),emp_pmf(d,1,:));
    ks_exp(d,1) = ks_pmf(geo_pmf(d,1,:),emp_pmf(d,1,:));
    
    chi2_gam(d,1) = chi2pdof_pmf(dur_range,squeeze(emp_hist(d,1,:)),squeeze(gam_pmf(d,1,:))*n_down_states,2);
    chi2_ig(d,1) = chi2pdof_pmf(dur_range,squeeze(emp_hist(d,1,:)),squeeze(ig_pmf(d,1,:))*n_down_states,2);
    chi2_exp(d,1) = chi2pdof_pmf(dur_range,squeeze(emp_hist(d,1,:)),squeeze(geo_pmf(d,1,:))*n_down_states,1);
    
    kldiv_gam(d,1) = kldiv_emp_pmf(squeeze(emp_pmf(d,1,:)),squeeze(gam_pmf(d,1,:)));
    kldiv_ig(d,1) = kldiv_emp_pmf(squeeze(emp_pmf(d,1,:)),squeeze(ig_pmf(d,1,:)));
    kldiv_exp(d,1) = kldiv_emp_pmf(squeeze(emp_pmf(d,1,:)),squeeze(geo_pmf(d,1,:)));
    
    %for up state
    [alpha,beta] = gamma_mlfit(state_durations_mp{2});
    [gam_pmf(d,2,:)] = gamma_pmf(dur_range,alpha,beta);
    [mu,lambda] = inverse_gauss_mlfit(state_durations_mp{2});
    [ig_pmf(d,2,:)] = inverse_gaussian_pmf(dur_range,mu,lambda);
%     lambda = 1/mean(state_durations_mp{2});
%     exp_pmf(d,2,:) = exponential_pmf(hmm.dur_range,lambda);
    p_fit = 1/nanmean(state_durations_mp{2}*Fs_dur);
    geo_pmf(d,2,:) = geometric_pmf(1:length(dur_range),p_fit);
    
    emp_hist(d,2,:) = hist(state_durations_mp{2},dur_range);
    emp_hist(d,2,end) = 0;
    emp_pmf(d,2,:) = emp_hist(d,2,:)/sum(emp_hist(d,2,:));
    
    ks_gam(d,2) = ks_pmf(gam_pmf(d,2,:),emp_pmf(d,2,:));
    ks_ig(d,2) = ks_pmf(ig_pmf(d,2,:),emp_pmf(d,2,:));
    ks_exp(d,2) = ks_pmf(geo_pmf(d,2,:),emp_pmf(d,2,:));
    
    chi2_gam(d,2) = chi2pdof_pmf(dur_range,emp_hist(d,2,:),gam_pmf(d,2,:)*n_up_states,2);
    chi2_ig(d,2) = chi2pdof_pmf(dur_range,emp_hist(d,2,:),ig_pmf(d,2,:)*n_up_states,2);
    chi2_exp(d,2) = chi2pdof_pmf(dur_range,emp_hist(d,2,:),geo_pmf(d,2,:)*n_up_states,1);

    kldiv_gam(d,2) = kldiv_emp_pmf(squeeze(emp_pmf(d,2,:)),squeeze(gam_pmf(d,2,:)));
    kldiv_ig(d,2) = kldiv_emp_pmf(squeeze(emp_pmf(d,2,:)),squeeze(ig_pmf(d,2,:)));
    kldiv_exp(d,2) = kldiv_emp_pmf(squeeze(emp_pmf(d,2,:)),squeeze(geo_pmf(d,2,:)));
end    
    %     subplot(2,1,1)
    %     plot(hmm.dur_range,squeeze(emp_pmf(d,1,:))), hold on
    %     plot(hmm.dur_range,squeeze(gam_pmf(d,1,:)),'r')
    %      plot(hmm.dur_range,squeeze(ig_pmf(d,1,:)),'k')
    %     plot(hmm.dur_range,squeeze(exp_pmf(d,1,:)),'c')
    %     axis tight
    %    xlim([0 5])
    %        subplot(2,1,2)
    %     plot(hmm.dur_range,squeeze(emp_pmf(d,2,:))), hold on
    %     plot(hmm.dur_range,squeeze(gam_pmf(d,2,:)),'r')
    %      plot(hmm.dur_range,squeeze(ig_pmf(d,2,:)),'k')
    %     plot(hmm.dur_range,squeeze(exp_pmf(d,2,:)),'c')
    %     axis tight
    %    xlim([0 5])
    %     cname = strcat(sess_data(d).region,'_',sess_data(d).layer,'_',sess_data(d).name);
    %     t_names = ['F:\WC_Germany\parietal_cortical_2010\state_durs\mp_' cname];
    %     print('-dpng',t_names), close
    
    
    %
    %
    %     [state_durations] = compute_state_durations_seg(fixmean_state_seq,Fsd);
    %
    %     %for down state
    %     [alpha,beta] = gamma_mlfit(state_durations{1});
    %     [gamf_pmf(d,1,:)] = gamma_pmf(hmm.dur_range,alpha,beta);
    %     [mu,lambda] = inverse_gauss_mlfit(state_durations{1});
    %     [igf_pmf(d,1,:)] = inverse_gaussian_pmf(hmm.dur_range,mu,lambda);
    %     lambda = 1/mean(state_durations{1});
    %     expf_pmf(d,1,:) = exponential_pmf(hmm.dur_range,lambda);
    %
    %     empf_pmf(d,1,:) = hist(state_durations{1},hmm.dur_range);
    %     empf_pmf(d,1,end) = 0;
    %     empf_pmf(d,1,:) = empf_pmf(d,1,:)/sum(empf_pmf(d,1,:));
    %
    %     ks_gamf(d,1) = ks_pmf(gamf_pmf(d,1,:),empf_pmf(d,1,:));
    %     ks_igf(d,1) = ks_pmf(igf_pmf(d,1,:),empf_pmf(d,1,:));
    %     ks_expf(d,1) = ks_pmf(expf_pmf(d,1,:),empf_pmf(d,1,:));
    %
    %     %for up state
    %     [alpha,beta] = gamma_mlfit(state_durations{2});
    %     [gamf_pmf(d,2,:)] = gamma_pmf(hmm.dur_range,alpha,beta);
    %     [mu,lambda] = inverse_gauss_mlfit(state_durations{2});
    %     [igf_pmf(d,2,:)] = inverse_gaussian_pmf(hmm.dur_range,mu,lambda);
    %     lambda = 1/mean(state_durations{2});
    %     expf_pmf(d,2,:) = exponential_pmf(hmm.dur_range,lambda);
    %
    %     empf_pmf(d,2,:) = hist(state_durations{2},hmm.dur_range);
    %     empf_pmf(d,2,end) = 0;
    %     empf_pmf(d,2,:) = empf_pmf(d,2,:)/sum(empf_pmf(d,2,:));
    %
    %     ks_gamf(d,2) = ks_pmf(gamf_pmf(d,2,:),empf_pmf(d,2,:));
    %     ks_igf(d,2) = ks_pmf(igf_pmf(d,2,:),empf_pmf(d,2,:));
    %     ks_expf(d,2) = ks_pmf(expf_pmf(d,2,:),empf_pmf(d,2,:));
    
    
%     load ./hsmm_state_seq4_seg_lf_4_28_10_v3
load ./hsmm_state_seq4_seg_lf_4_5_2011
dur_range = hmm4.dur_range;
dur_diff = dur_range(2) - dur_range(1);
cdur_range = [0 dur_range+dur_diff];
Fs_dur = 1/(dur_range(2)-dur_range(1));
[state_durations] = compute_state_durations_seg(hmm_bbstate_seq4,Fsd);
    avg_up4(d) = nanmean(state_durations{2});
    avg_down4(d) = nanmean(state_durations{1});
    max_up4(d) = nanmax(state_durations{2});
    max_down4(d) = nanmax(state_durations{1});
    n_up_states4 = length(state_durations{2});
    n_down_states4 = length(state_durations{1});
  
    %for down state
    [alpha,beta] = gamma_mlfit(state_durations{1});
    [gam_pmf4(d,1,:)] = gamma_pmf(dur_range,alpha,beta);
    [mu,lambda] = inverse_gauss_mlfit(state_durations{1});
    [ig_pmf4(d,1,:)] = inverse_gaussian_pmf(dur_range,mu,lambda);
%     lambda = 1/mean(state_durations{1});
%     exp_pmf4(d,1,:) = exponential_pmf(hmm4.dur_range,lambda);
    p_fit = 1/nanmean(state_durations{1}*Fs_dur);
    geo_pmf4(d,1,:) = geometric_pmf(1:length(dur_range),p_fit);
    
    emp_hist4(d,1,:) = hist(state_durations{1},dur_range);
    emp_hist4(d,1,end) = 0;
    emp_pmf4(d,1,:) = emp_hist4(d,1,:)/sum(emp_hist4(d,1,:));
    
    ks_gam4(d,1) = ks_pmf(gam_pmf4(d,1,:),emp_pmf4(d,1,:));
    ks_ig4(d,1) = ks_pmf(ig_pmf4(d,1,:),emp_pmf4(d,1,:));
    ks_exp4(d,1) = ks_pmf(geo_pmf4(d,1,:),emp_pmf4(d,1,:));
    
    chi2_gam4(d,1) = chi2pdof_pmf(dur_range,emp_hist4(d,1,:),gam_pmf4(d,1,:)*n_down_states4,2);
    chi2_ig4(d,1) = chi2pdof_pmf(dur_range,emp_hist4(d,1,:),ig_pmf4(d,1,:)*n_down_states4,2);
    chi2_exp4(d,1) = chi2pdof_pmf(dur_range,emp_hist4(d,1,:),geo_pmf4(d,1,:)*n_down_states4,1);

    kldiv_gam4(d,1) = kldiv_emp_pmf(squeeze(emp_pmf4(d,1,:)),squeeze(gam_pmf4(d,1,:)));
    kldiv_ig4(d,1) = kldiv_emp_pmf(squeeze(emp_pmf4(d,1,:)),squeeze(ig_pmf4(d,1,:)));
    kldiv_exp4(d,1) = kldiv_emp_pmf(squeeze(emp_pmf4(d,1,:)),squeeze(geo_pmf4(d,1,:)));

    [~,down_bin_inds] = histc(state_durations{1},cdur_range);
    down_bin_inds(down_bin_inds==0) = [];
    ig_llik4(d,1) = sum(log(ig_pmf4(d,1,down_bin_inds)),3);
    ig_nllik4(d,1) = ig_llik4(d,1)/length(state_durations{1});
    gam_llik4(d,1) = sum(log(gam_pmf4(d,1,down_bin_inds)),3);
    gam_nllik4(d,1) = gam_llik4(d,1)/length(state_durations{1});
    geo_llik4(d,1) = sum(log(geo_pmf4(d,1,down_bin_inds)),3);
    geo_nllik4(d,1) = geo_llik4(d,1)/length(state_durations{1});
    
     ig_bic4(d,1) = -2*ig_llik4(d,1)+2*log(length(state_durations{1}));
    gam_bic4(d,1) = -2*gam_llik4(d,1)+2*log(length(state_durations{1}));
    geo_bic4(d,1) = -2*geo_llik4(d,1)+log(length(state_durations{1}));
   
    %for up state
    [alpha,beta] = gamma_mlfit(state_durations{2});
    [gam_pmf4(d,2,:)] = gamma_pmf(dur_range,alpha,beta);
    [mu,lambda] = inverse_gauss_mlfit(state_durations{2});
    [ig_pmf4(d,2,:)] = inverse_gaussian_pmf(dur_range,mu,lambda);
%     lambda = 1/mean(state_durations{2});
%     exp_pmf4(d,2,:) = exponential_pmf(hmm4.dur_range,lambda);
    p_fit = 1/nanmean(state_durations{2}*Fs_dur);
    geo_pmf4(d,2,:) = geometric_pmf(1:length(dur_range),p_fit);
    
    emp_hist4(d,2,:) = hist(state_durations{2},dur_range);
    emp_hist4(d,2,end) = 0;
    emp_pmf4(d,2,:) = emp_hist4(d,2,:)/sum(emp_hist4(d,2,:));
    
    ks_gam4(d,2) = ks_pmf(gam_pmf4(d,2,:),emp_pmf4(d,2,:));
    ks_ig4(d,2) = ks_pmf(ig_pmf4(d,2,:),emp_pmf4(d,2,:));
    ks_exp4(d,2) = ks_pmf(geo_pmf4(d,2,:),emp_pmf4(d,2,:));
  
    chi2_gam4(d,2) = chi2pdof_pmf(dur_range,emp_hist4(d,2,:),gam_pmf4(d,2,:)*n_up_states4,2);
    chi2_ig4(d,2) = chi2pdof_pmf(dur_range,emp_hist4(d,2,:),ig_pmf4(d,2,:)*n_up_states4,2);
    chi2_exp4(d,2) = chi2pdof_pmf(dur_range,emp_hist4(d,2,:),geo_pmf4(d,2,:)*n_up_states4,1);

    kldiv_gam4(d,2) = kldiv_emp_pmf(squeeze(emp_pmf4(d,2,:)),squeeze(gam_pmf4(d,2,:)));
    kldiv_ig4(d,2) = kldiv_emp_pmf(squeeze(emp_pmf4(d,2,:)),squeeze(ig_pmf4(d,2,:)));
    kldiv_exp4(d,2) = kldiv_emp_pmf(squeeze(emp_pmf4(d,2,:)),squeeze(geo_pmf4(d,2,:)));

    [~,up_bin_inds] = histc(state_durations{2},cdur_range);
    up_bin_inds(up_bin_inds==0) = [];
    ig_llik4(d,2) = sum(log(ig_pmf4(d,2,up_bin_inds)),3);
    ig_nllik4(d,2) = ig_llik4(d,2)/length(state_durations{2});
    gam_llik4(d,2) = sum(log(gam_pmf4(d,2,up_bin_inds)),3);
    gam_nllik4(d,2) = gam_llik4(d,2)/length(state_durations{2});
    geo_llik4(d,2) = sum(log(geo_pmf4(d,2,up_bin_inds)),3);
    geo_nllik4(d,2) = geo_llik4(d,2)/length(state_durations{2});

    ig_bic4(d,2) = -2*ig_llik4(d,2)+2*log(length(state_durations{2}));
    gam_bic4(d,2) = -2*gam_llik4(d,2)+2*log(length(state_durations{2}));
    geo_bic4(d,2) = -2*geo_llik4(d,2)+log(length(state_durations{2}));
    %         subplot(2,1,1)
    %     plot(hmm4.dur_range,squeeze(emp_pmf4(d,1,:))), hold on
    %     plot(hmm4.dur_range,squeeze(gam_pmf4(d,1,:)),'r')
    %      plot(hmm4.dur_range,squeeze(ig_pmf4(d,1,:)),'k')
    %     plot(hmm4.dur_range,squeeze(exp_pmf4(d,1,:)),'c')
    %     axis tight
    %    xlim([0 5])
    %        subplot(2,1,2)
    %     plot(hmm4.dur_range,squeeze(emp_pmf4(d,2,:))), hold on
    %     plot(hmm4.dur_range,squeeze(gam_pmf4(d,2,:)),'r')
    %      plot(hmm4.dur_range,squeeze(ig_pmf4(d,2,:)),'k')
    %     plot(hmm4.dur_range,squeeze(exp_pmf4(d,2,:)),'c')
    %     axis tight
    %    xlim([0 5])
    %     cname = strcat(sess_data(d).region,'_',sess_data(d).layer,'_',sess_data(d).name);
    %     t_names = ['F:\WC_Germany\parietal_cortical_2010\state_durs\lf4_' cname];
    %     print('-dpng',t_names), close
    
%     load hsmm_state_seq4_seg_hf_4_26_10_v2
%     [state_durationshf] = compute_state_durations_seg(hsmm_bbstate_seq4_hf,Fsd);
%     avg_up4hf(d) = nanmean(state_durationshf{2});
%     avg_down4hf(d) = nanmean(state_durationshf{1});
% 
%     %for down state
%     [alpha,beta] = gamma_mlfit(state_durationshf{1});
%     [gam_pmf4hf(d,1,:)] = gamma_pmf(hmm4.dur_range,alpha,beta);
%     [mu,lambda] = inverse_gauss_mlfit(state_durationshf{1});
%     [ig_pmf4hf(d,1,:)] = inverse_gaussian_pmf(hmm4.dur_range,mu,lambda);
% %     lambda = 1/mean(state_durationshf{1});
% %     exp_pmf4hf(d,1,:) = exponential_pmf(hmm4.dur_range,lambda);
%     p_fit = 1/mean(state_durationshf{1}*Fs);
%     geo_pmf4hf(d,1,:) = geometric_pmf(1:length(hmm4.dur_range),p_fit);
%     
%     emp_pmf4hf(d,1,:) = hist(state_durationshf{1},hmm4.dur_range);
%     emp_pmf4hf(d,1,end) = 0;
%     emp_pmf4hf(d,1,:) = emp_pmf4hf(d,1,:)/sum(emp_pmf4hf(d,1,:));
%     
%     ks_gam4hf(d,1) = ks_pmf(gam_pmf4hf(d,1,:),emp_pmf4hf(d,1,:));
%     ks_ig4hf(d,1) = ks_pmf(ig_pmf4hf(d,1,:),emp_pmf4hf(d,1,:));
%     ks_exp4hf(d,1) = ks_pmf(geo_pmf4hf(d,1,:),emp_pmf4hf(d,1,:));
%     
%     %for up state
%     [alpha,beta] = gamma_mlfit(state_durationshf{2});
%     [gam_pmf4hf(d,2,:)] = gamma_pmf(hmm4.dur_range,alpha,beta);
%     [mu,lambda] = inverse_gauss_mlfit(state_durationshf{2});
%     [ig_pmf4hf(d,2,:)] = inverse_gaussian_pmf(hmm4.dur_range,mu,lambda);
% %     lambda = 1/mean(state_durationshf{2});
% %     exp_pmf4hf(d,2,:) = exponential_pmf(hmm4.dur_range,lambda);
%       p_fit = 1/mean(state_durationshf{2}*Fs);
%     geo_pmf4hf(d,2,:) = geometric_pmf(1:length(hmm4.dur_range),p_fit);
%   
%     emp_pmf4hf(d,2,:) = hist(state_durationshf{2},hmm4.dur_range);
%     emp_pmf4hf(d,2,end) = 0;
%     emp_pmf4hf(d,2,:) = emp_pmf4hf(d,2,:)/sum(emp_pmf4hf(d,2,:));
%     
%     ks_gam4hf(d,2) = ks_pmf(gam_pmf4hf(d,2,:),emp_pmf4hf(d,2,:));
%     ks_ig4hf(d,2) = ks_pmf(ig_pmf4hf(d,2,:),emp_pmf4hf(d,2,:));
%     ks_exp4hf(d,2) = ks_pmf(geo_pmf4hf(d,2,:),emp_pmf4hf(d,2,:));
    
    %             subplot(2,1,1)
    %     plot(hmm4.dur_range,squeeze(emp_pmf4hf(d,1,:))), hold on
    %     plot(hmm4.dur_range,squeeze(gam_pmf4hf(d,1,:)),'r')
    %      plot(hmm4.dur_range,squeeze(ig_pmf4hf(d,1,:)),'k')
    %     plot(hmm4.dur_range,squeeze(exp_pmf4hf(d,1,:)),'c')
    %     axis tight
    %    xlim([0 5])
    %        subplot(2,1,2)
    %     plot(hmm4.dur_range,squeeze(emp_pmf4hf(d,2,:))), hold on
    %     plot(hmm4.dur_range,squeeze(gam_pmf4hf(d,2,:)),'r')
    %      plot(hmm4.dur_range,squeeze(ig_pmf4hf(d,2,:)),'k')
    %     plot(hmm4.dur_range,squeeze(exp_pmf4hf(d,2,:)),'c')
    %     axis tight
    %    xlim([0 5])
    %     cname = strcat(sess_data(d).region,'_',sess_data(d).layer,'_',sess_data(d).name);
    %     t_names = ['F:\WC_Germany\parietal_cortical_2010\state_durs\lf4hf_' cname];
    %     print('-dpng',t_names), close
    
%     load fixmean_state_seqs4_4_26_10_v2
%     [state_durationsf] = compute_state_durations_seg(fixmean_state_seq4,Fsd);
%     avg_up4f(d) = nanmean(state_durationsf{2});
%     avg_down4f(d) = nanmean(state_durationsf{1});
%     
%     %for down state
%     [alpha,beta] = gamma_mlfit(state_durationsf{1});
%     [gamf_pmf4(d,1,:)] = gamma_pmf(hmm4.dur_range,alpha,beta);
%     [mu,lambda] = inverse_gauss_mlfit(state_durationsf{1});
%     [igf_pmf4(d,1,:)] = inverse_gaussian_pmf(hmm4.dur_range,mu,lambda);
% %     lambda = 1/mean(state_durationsf{1});
% %     expf_pmf4(d,1,:) = exponential_pmf(hmm4.dur_range,lambda);
%       p_fit = 1/mean(state_durationsf{1}*Fs);
%     geof_pmf4(d,1,:) = geometric_pmf(1:length(hmm4.dur_range),p_fit);
%     
%     empf_pmf4(d,1,:) = hist(state_durationsf{1},hmm4.dur_range);
%     empf_pmf4(d,1,end) = 0;
%     empf_pmf4(d,1,:) = empf_pmf4(d,1,:)/sum(empf_pmf4(d,1,:));
%     
%     ks_gam4f(d,1) = ks_pmf(gamf_pmf4(d,1,:),empf_pmf4(d,1,:));
%     ks_ig4f(d,1) = ks_pmf(igf_pmf4(d,1,:),empf_pmf4(d,1,:));
%     ks_exp4f(d,1) = ks_pmf(geof_pmf4(d,1,:),empf_pmf4(d,1,:));
%     
%     %for up state
%     [alpha,beta] = gamma_mlfit(state_durationsf{2});
%     [gamf_pmf4(d,2,:)] = gamma_pmf(hmm4.dur_range,alpha,beta);
%     [mu,lambda] = inverse_gauss_mlfit(state_durationsf{2});
%     [igf_pmf4(d,2,:)] = inverse_gaussian_pmf(hmm4.dur_range,mu,lambda);
% %     lambda = 1/mean(state_durationsf{2});
% %     expf_pmf4(d,2,:) = exponential_pmf(hmm4.dur_range,lambda);
%        p_fit = 1/mean(state_durationsf{2}*Fs);
%     geof_pmf4(d,2,:) = geometric_pmf(1:length(hmm4.dur_range),p_fit);
%    
%     empf_pmf4(d,2,:) = hist(state_durationsf{2},hmm4.dur_range);
%     empf_pmf4(d,2,end) = 0;
%     empf_pmf4(d,2,:) = empf_pmf4(d,2,:)/sum(empf_pmf4(d,2,:));
%     
%     ks_gam4f(d,2) = ks_pmf(gamf_pmf4(d,2,:),empf_pmf4(d,2,:));
%     ks_ig4f(d,2) = ks_pmf(igf_pmf4(d,2,:),empf_pmf4(d,2,:));
%     ks_exp4f(d,2) = ks_pmf(geof_pmf4(d,2,:),empf_pmf4(d,2,:));
    
%     subplot(2,1,1)
%     plot(hmm4.dur_range,squeeze(empf_pmf4(d,1,:))), hold on
%     plot(hmm4.dur_range,squeeze(gamf_pmf4(d,1,:)),'r')
%     plot(hmm4.dur_range,squeeze(igf_pmf4(d,1,:)),'k')
%     plot(hmm4.dur_range,squeeze(expf_pmf4(d,1,:)),'c')
%     axis tight
%     xlim([0 5])
%     subplot(2,1,2)
%     plot(hmm4.dur_range,squeeze(empf_pmf4(d,2,:))), hold on
%     plot(hmm4.dur_range,squeeze(gamf_pmf4(d,2,:)),'r')
%     plot(hmm4.dur_range,squeeze(igf_pmf4(d,2,:)),'k')
%     plot(hmm4.dur_range,squeeze(expf_pmf4(d,2,:)),'c')
%     axis tight
%     xlim([0 5])
%     cname = strcat(sess_data(d).region,'_',sess_data(d).layer,'_',sess_data(d).name);
%     t_names = ['F:\WC_Germany\parietal_cortical_2010\state_durs\lf4fix_' cname];
%     print('-dpng',t_names), close
    
%     [dummy,ks_4f(d,1)] = kstest2(state_durations{1},state_durationsf{1});
%     [dummy,ks_4f(d,2)] = kstest2(state_durations{2},state_durationsf{2});
%     [dummy,ks_4mp(d,1),kss_4mp(d,1)] = kstest2(state_durations{1},state_durations_mp{1});
%     [dummy,ks_4mp(d,2),kss_4mp(d,2)] = kstest2(state_durations{2},state_durations_mp{2});
%     [dummy,ks_4fmp(d,2),kss_4fmp(d,1)] = kstest2(state_durationsf{1},state_durations_mp{1});
%     [dummy,ks_4fmp(d,2),kss_4fmp(d,2)] = kstest2(state_durationsf{2},state_durations_mp{2});
%     [dummy,ks_4fmpf(d,2),kss_4fmpf(d,1)] = kstest2(state_durationsf{1},state_durationsf_mp{1});
%     [dummy,ks_4fmpf(d,2),kss_4fmpf(d,2)] = kstest2(state_durationsf{2},state_durationsf_mp{2});
   
end


ig_aic4 = 2*2-2*ig_llik4;
gam_aic4 = 2*2-2*gam_llik4;
geo_aic4 = 2*1-2*geo_llik4;

%% 
thom_pfc = find(ismember(thom_el,thom_pfc));

%% FOR KS Stats
figure
plot(ks_exp4(:,1),ks_gam4(:,1),'.','markersize',16), hold on
plot(ks_exp4(:,2),ks_gam4(:,2),'r.','markersize',16)
% plot(ks_exp(thom_pfc,1),ks_gam(thom_pfc,1),'o','markersize',10), hold on
% plot(ks_exp(thom_pfc,2),ks_gam(thom_pfc,2),'ro','markersize',10)
line([0 0.5],[0 0.5],'color','k')
xlim([0 0.5]), ylim([0 0.5])

figure
plot(ks_exp4(:,1),ks_ig4(:,1),'.','markersize',16), hold on
plot(ks_exp4(:,2),ks_ig4(:,2),'r.','markersize',16)
% plot(ks_exp(thom_pfc,1),ks_ig(thom_pfc,1),'o','markersize',10), hold on
% plot(ks_exp(thom_pfc,2),ks_ig(thom_pfc,2),'ro','markersize',10)
line([0 0.5],[0 0.5],'color','k')
xlim([0 0.5]), ylim([0 0.5])

figure
plot(ks_gam4(:,1),ks_ig4(:,1),'.','markersize',16), hold on
plot(ks_gam4(:,2),ks_ig4(:,2),'r.','markersize',16)
% plot(ks_gam(thom_pfc,1),ks_ig(thom_pfc,1),'o','markersize',10), hold on
% plot(ks_gam(thom_pfc,2),ks_ig(thom_pfc,2),'ro','markersize',10)
line([0 0.5],[0 0.5],'color','k')
xlim([0 0.5]), ylim([0 0.5])

%% FOR KL DIV
figure
set(gca,'fontname','arial','fontsize',14)
plot(kldiv_exp4(:,1),kldiv_gam4(:,1),'.','markersize',16), hold on
plot(kldiv_exp4(:,2),kldiv_gam4(:,2),'r.','markersize',16)
line([0 2],[0 2],'color','k')
xlim([0.8 2]), ylim([0 2])
xlabel('KL GEO','fontsize',16,'fontname','arial')
ylabel('KL GAM','fontsize',16,'fontname','arial')

figure
set(gca,'fontname','arial','fontsize',14)
plot(kldiv_exp4(:,1),kldiv_ig4(:,1),'.','markersize',16), hold on
plot(kldiv_exp4(:,2),kldiv_ig4(:,2),'r.','markersize',16)
line([0 2],[0 2],'color','k')
xlim([0.8 2]), ylim([0 2])
xlabel('KL GEO','fontsize',16,'fontname','arial')
ylabel('KL IG','fontsize',16,'fontname','arial')

figure
set(gca,'fontname','arial','fontsize',14)
plot(kldiv_gam4(:,1),kldiv_ig4(:,1),'.','markersize',16), hold on
plot(kldiv_gam4(:,2),kldiv_ig4(:,2),'r.','markersize',16)
line([0 2],[0 2],'color','k')
xlim([0. 2]), ylim([0 2])
xlabel('KL GAM','fontsize',16,'fontname','arial')
ylabel('KL IG','fontsize',16,'fontname','arial')

%% FOR Likelihood
figure
set(gca,'fontname','arial','fontsize',14)
plot(geo_nllik4(:,1),gam_nllik4(:,1),'.','markersize',16), hold on
plot(geo_nllik4(:,2),gam_nllik4(:,2),'r.','markersize',16)
line([-6.5 -3.5],[-6.5 -3.5],'color','k')
xlim([-6.5 -4.5]), ylim([-6.5 -3.5])
xlabel('L GEO','fontsize',16,'fontname','arial')
ylabel('L GAM','fontsize',16,'fontname','arial')

figure
set(gca,'fontname','arial','fontsize',14)
plot(geo_nllik4(:,1),ig_nllik4(:,1),'.','markersize',16), hold on
plot(geo_nllik4(:,2),ig_nllik4(:,2),'r.','markersize',16)
line([-6.5 -3.5],[-6.5 -3.5],'color','k')
xlim([-6.5 -4.5]), ylim([-6.5 -3.5])
xlabel('L GEO','fontsize',16,'fontname','arial')
ylabel('L IG','fontsize',16,'fontname','arial')

figure
set(gca,'fontname','arial','fontsize',14)
plot(gam_nllik4(:,1),ig_nllik4(:,1),'.','markersize',16), hold on
plot(gam_nllik4(:,2),ig_nllik4(:,2),'r.','markersize',16)
line([-6.5 -3.5],[-6.5 -3.5],'color','k')
xlim([-6.5 -3.5]), ylim([-6.5 -3.5])
xlabel('L GAM','fontsize',16,'fontname','arial')
ylabel('L IG','fontsize',16,'fontname','arial')

%%
d = find_struct_field_vals(sess_data,'name','2007-04-17_C');
figure
set(gca,'fontname','arial','fontsize',14)
subplot(2,1,1)
bar(dur_range,squeeze(emp_pmf4(d,1,:))), hold on
plot(dur_range,squeeze(gam_pmf4(d,1,:)),'r')
plot(dur_range,squeeze(ig_pmf4(d,1,:)),'k')
plot(dur_range,squeeze(geo_pmf4(d,1,:)),'c')
axis tight
xlim([0 3])
xlabel('Duration (s)','fontsize',16,'fontname','arial')
ylabel('Relative Frequency','fontsize',16,'fontname','arial')

subplot(2,1,2)
bar(dur_range,squeeze(emp_pmf4(d,2,:))), hold on
plot(dur_range,squeeze(gam_pmf4(d,2,:)),'r')
plot(dur_range,squeeze(ig_pmf4(d,2,:)),'k')
plot(dur_range,squeeze(geo_pmf4(d,2,:)),'c')
axis tight
xlim([0 3])
xlabel('Duration (s)','fontsize',16,'fontname','arial')
ylabel('Relative Frequency','fontsize',16,'fontname','arial')


% figure
% subplot(2,1,1)
% bar(hmm4.dur_range,squeeze(emp_pmf4hf(d,1,:))), hold on
% plot(hmm4.dur_range,squeeze(gam_pmf4hf(d,1,:)),'r')
% plot(hmm4.dur_range,squeeze(ig_pmf4hf(d,1,:)),'k')
% plot(hmm4.dur_range,squeeze(exp_pmf4hf(d,1,:)),'c')
% axis tight
% xlim([0 3])
% subplot(2,1,2)
% bar(hmm4.dur_range,squeeze(emp_pmf4hf(d,2,:))), hold on
% plot(hmm4.dur_range,squeeze(gam_pmf4hf(d,2,:)),'r')
% plot(hmm4.dur_range,squeeze(ig_pmf4hf(d,2,:)),'k')
% plot(hmm4.dur_range,squeeze(exp_pmf4hf(d,2,:)),'c')
% axis tight
% xlim([0 3])


%%
% sm_fac = 3;
% for d = 1:n
%     semp_pmf4(d,1,:) = jmm_smooth_1d_cor(squeeze(emp_pmf4(d,1,:)),sm_fac);
%     semp_pmf4(d,2,:) = jmm_smooth_1d_cor(squeeze(emp_pmf4(d,2,:)),sm_fac);
%     sempf_pmf4(d,1,:) = jmm_smooth_1d_cor(squeeze(empf_pmf4(d,1,:)),sm_fac);
%     sempf_pmf4(d,2,:) = jmm_smooth_1d_cor(squeeze(empf_pmf4(d,2,:)),sm_fac);
%     semp_pmf4hf(d,1,:) = jmm_smooth_1d_cor(squeeze(emp_pmf4hf(d,1,:)),sm_fac);
%     semp_pmf4hf(d,2,:) = jmm_smooth_1d_cor(squeeze(emp_pmf4hf(d,2,:)),sm_fac);
% end
% figure
% ebar_mat(hmm4.dur_range,semp_pmf4(:,1,:)-semp_pmf4hf(:,1,:),.001,'b'), hold on
% ebar_mat(hmm4.dur_range,semp_pmf4(:,2,:)-semp_pmf4hf(:,2,:),.001,'r'), hold on
% xlim([0 2])

% figure
% subplot(2,1,1)
% ebar_mat(dur_range,emp_pmf4(:,1,:),.001,'b');
% xlim([0 3])
% hold on
% % ebar_mat(hmm4.dur_range,emp_pmf4hf(:,1,:),.001,'r');
% ebar_mat(dur_range,emp_pmf(thom_pfc,1,:),.001,'r');
% yl = ylim;
% ylim([0 yl(2)])
% subplot(2,1,2)
% ebar_mat(dur_range,emp_pmf4(:,2,:),.001,'b');
% xlim([0 3])
% hold on
% % ebar_mat(hmm4.dur_range,emp_pmf4hf(:,2,:),.001,'r');
% ebar_mat(dur_range,emp_pmf(thom_pfc,2,:),.001,'r');
% yl = ylim;
% ylim([0 yl(2)])

% figure
% subplot(2,1,1)
% ebar_mat(dur_range,emp_pmf4(:,1,:),.001,'b');
% xlim([0 3])
% hold on
% ebar_mat(dur_range,empf_pmf4(:,1,:),.001,'r');
% yl = ylim;
% ylim([0 yl(2)])
% subplot(2,1,2)
% ebar_mat(dur_range,emp_pmf4(:,2,:),.001,'b');
% xlim([0 3])
% hold on
% ebar_mat(dur_range,empf_pmf4(:,2,:),.001,'r');
% yl = ylim;
% ylim([0 yl(2)])

figure
set(gca,'fontname','arial','fontsize',14)
ebar_mat(dur_range,emp_pmf4(:,1,:),.001,'b');
xlim([0 3])
hold on
ebar_mat(dur_range,emp_pmf4(:,2,:),.001,'r');
yl = ylim;
ylim([0 yl(2)])
xlabel('Duration (s)','fontsize',16,'fontname','arial')
ylabel('Relative Frequency','fontsize',16,'fontname','arial')



%% within animal analysis
for i = 1:11
    cur_data = find(animal_id==i);
    ks_exp4_a(i,:) = mean(ks_exp4(cur_data,:));
    ks_gam4_a(i,:) = mean(ks_gam4(cur_data,:));
    ks_ig4_a(i,:) = mean(ks_ig4(cur_data,:));
end