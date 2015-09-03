clear all
close all

addpath('C:\WC_Germany\parietal_cortical_2010\')
addpath('C:\WC_Germany\hsmm_state_detection\')
cd C:\WC_Germany\parietal_cortical_2010\
load parietal_cortical_2010
load C:\WC_Germany\parietal_cortical_2010\desynch_times_individual


frontal = find_struct_field_vals(sess_data,'region','frontal');
prefrontal = find_struct_field_vals(sess_data,'region','prefrontal');
parietal = find_struct_field_vals(sess_data,'region','parietal');

% sess_data = sess_data(parietal);
% desynch_times = desynch_times(parietal);

%get rid of interneurons
interneurons = find_struct_field_vals(sess_data,'cell_type','interneuron');
sess_data(interneurons) = [];
desynch_times_lf4(interneurons) = [];

raw_Fs = 2016;
dsf = 12;
Fsd = raw_Fs/dsf;
% lf_lcf = 0.05;
% lf_hcf = 2;
% hf_lcf = 20;
% hf_hcf = 80;
% hf_smooth = 0.03;

% thomas_el = find_struct_field_vals(sess_data,'thom_elec',1);
% sess_data = sess_data(thomas_el);
% desynch_times = desynch_times_lf4(thomas_el);

frontal = find_struct_field_vals(sess_data,'region','frontal');
prefrontal = find_struct_field_vals(sess_data,'region','prefrontal');
parietal = find_struct_field_vals(sess_data,'region','parietal');
thom_el = find_struct_field_vals(sess_data,'thom_elec',1);
thom_par = thom_el(find(ismember(thom_el,parietal)));
thom_pfc = setdiff(thom_el,thom_par);
% sess_data = sess_data(thom_pfc);
% desynch_times_mp = desynch_times_mp(thom_pfc);
% desynch_times_lf4 = desynch_times_lf4(thom_pfc);
% desynch_times_lf8 = desynch_times_lf8(thom_pfc);
sess_data = sess_data(thom_el);


n = length(sess_data);
animal_id = [9 7 10 6 2 3 4 9 7 10 1 2 3 4 5 6 7 8 9 11 1];

for d = 1:n
    d
    cdir = sess_data(d).directory;
    cdir(1) = 'C';
    cd(cdir)
    %     load used_data wcv_minus_spike lf8 lf4
    %
    %     [wcv_hf,t_axis,Fs] = get_hf_features(wcv_minus_spike,raw_Fs,Fsd,[hf_lcf hf_hcf],hf_smooth);
    %     [lf8_hf,t_axis,Fs] = get_hf_features(lf8,raw_Fs,Fsd,[hf_lcf hf_hcf],hf_smooth);
    %     [wcv_lf,t_axis,Fs] = get_lf_features(wcv_minus_spike,raw_Fs,Fsd,[lf_lcf lf_hcf]);
    %     [lf8_lf,t_axis,Fs] = get_lf_features(lf8,raw_Fs,Fsd,[lf_lcf lf_hcf]);
    %
    %     if sess_data(d).thom_elec
    %         [lf4_hf,t_axis,Fs] = get_hf_features(lf4,raw_Fs,Fsd,[hf_lcf hf_hcf],hf_smooth);
    %         [lf4_lf,t_axis,Fs] = get_lf_features(lf4,raw_Fs,Fsd,[lf_lcf lf_hcf]);
    %     end
    
    %% a priori comparisons
    %     mix=gmm(1,2,'full');
    %     gmm_options(3) = 1e-10; %tolerance
    %     gmm_options(5) = 1;%reset cov if singular values
    %     gmm_options(14) = 200; %max iterations
    %     [mix, options, errlog] = gmmem(mix,wcv_lf,gmm_options);
    %     state_means = mix.centres;
    %     state_covars = squeeze(mix.covars);
    %     mp_pr_kl_lf(d) = gauss_kl_div(state_means(2)-state_means(1),state_covars(1),state_covars(2));
    %     mp_pr_max_lik_lf(d) = options(8)/length(wcv_lf);
    %
    %     mix=gmm(1,2,'full');
    %     gmm_options(3) = 1e-10; %tolerance
    %     gmm_options(5) = 1;%reset cov if singular values
    %     gmm_options(14) = 200; %max iterations
    %     [mix, options, errlog] = gmmem(mix,wcv_hf,gmm_options);
    %     state_means = mix.centres;
    %     state_covars = squeeze(mix.covars);
    %     mp_pr_kl_hf(d) = gauss_kl_div(state_means(2)-state_means(1),state_covars(1),state_covars(2));
    %     mp_pr_max_lik_hf(d) = options(8)/length(wcv_hf);
    %
    %     mix=gmm(1,2,'full');
    %     gmm_options(3) = 1e-10; %tolerance
    %     gmm_options(5) = 1;%reset cov if singular values
    %     gmm_options(14) = 200; %max iterations
    %     [mix, options, errlog] = gmmem(mix,lf8_lf,gmm_options);
    %     state_means = mix.centres;
    %     state_covars = squeeze(mix.covars);
    %     lf8_pr_kl_lf(d) = gauss_kl_div(state_means(2)-state_means(1),state_covars(1),state_covars(2));
    %     lf8_pr_max_lik_lf(d) = options(8)/length(lf8_lf);
    %
    %     mix=gmm(1,2,'full');
    %     gmm_options(3) = 1e-10; %tolerance
    %     gmm_options(5) = 1;%reset cov if singular values
    %     gmm_options(14) = 200; %max iterations
    %     [mix, options, errlog] = gmmem(mix,lf8_hf,gmm_options);
    %     state_means = mix.centres;
    %     state_covars = squeeze(mix.covars);
    %     lf8_pr_kl_hf(d) = gauss_kl_div(state_means(2)-state_means(1),state_covars(1),state_covars(2));
    %     lf8_pr_max_lik_hf(d) = options(8)/length(lf8_hf);
    %
    %     if sess_data(d).thom_elec
    %         mix=gmm(1,2,'full');
    %         gmm_options(3) = 1e-10; %tolerance
    %         gmm_options(5) = 1;%reset cov if singular values
    %         gmm_options(14) = 200; %max iterations
    %         [mix, options, errlog] = gmmem(mix,lf4_lf,gmm_options);
    %         state_means = mix.centres;
    %         state_covars = squeeze(mix.covars);
    %         lf4_pr_kl_lf(d) = gauss_kl_div(state_means(2)-state_means(1),state_covars(1),state_covars(2));
    %         lf4_pr_max_lik_lf(d) = options(8)/length(lf4_lf);
    %
    %         mix=gmm(1,2,'full');
    %         gmm_options(3) = 1e-10; %tolerance
    %         gmm_options(5) = 1;%reset cov if singular values
    %         gmm_options(14) = 200; %max iterations
    %         [mix, options, errlog] = gmmem(mix,lf4_hf,gmm_options);
    %         state_means = mix.centres;
    %         state_covars = squeeze(mix.covars);
    %         lf4_pr_kl_hf(d) = gauss_kl_div(state_means(2)-state_means(1),state_covars(1),state_covars(2));
    %         lf4_pr_max_lik_hf(d) = options(8)/length(lf4_hf);
    %     end
    
    %% post-hoc
    %     load hsmm_state_seq_seg_lf_4_10_10
    %     meandiff = [];
    %     for i = 1:hmm.Nsegs
    %         meandiff = [meandiff; hmm.state(2).meanfun{i}-hmm.state(1).meanfun{i}];
    %     end
    %     meandiff = mean(meandiff);
    %     covar1 = hmm.state(1).var;
    %     covar2 = hmm.state(2).var;
    %     mp_kl_lf(d) = gauss_kl_div(meandiff,covar1,covar2);
    %     mp_max_lik_lf(d) = hmm.max_lik/hmm.T;
    %
    %     load hsmm_state_seq_seg_hf_4_10_10
    %     meandiff = [];
    %     for i = 1:hmm.Nsegs
    %         meandiff = [meandiff; hmm_hf.state(2).meanfun{i}-hmm_hf.state(1).meanfun{i}];
    %     end
    %     meandiff = mean(meandiff);
    %     covar1 = hmm_hf.state(1).var;
    %     covar2 = hmm_hf.state(2).var;
    %     mp_kl_hf(d) = gauss_kl_div(meandiff,covar1,covar2);
    %     mp_max_lik_hf(d) = hmm_hf.max_lik/hmm_hf.T;
    %     clear hmm*
    
    %     load hsmm_state_seq8_seg_lf_4
    %     meandiff = [];
    %     for i = 1:hmm8.Nsegs
    %         meandiff = [meandiff; hmm8.state(2).meanfun{i}-hmm8.state(1).meanfun{i}];
    %     end
    %     meandiff = mean(meandiff);
    %     covar1 = hmm8.state(1).var;
    %     covar2 = hmm8.state(2).var;
    %     kl1 = gauss_kl_div(meandiff,covar1,covar2);
    %     kl2 = gauss_kl_div(-meandiff,covar2,covar1);
    %     lf8_kl_lf(d) = kl1+kl2;
    %     lf8_max_lik_lf(d) = hmm8.max_lik/hmm8.T;
    %
    %     load hsmm_state_seq8_seg_hf_4_10_10
    %     meandiff = [];
    %     for i = 1:hmm8.Nsegs
    %         meandiff = [meandiff; hmm8_hf.state(2).meanfun{i}-hmm8_hf.state(1).meanfun{i}];
    %     end
    %     meandiff = mean(meandiff);
    %     covar1 = hmm8_hf.state(1).var;
    %     covar2 = hmm8_hf.state(2).var;
    %       kl1 = gauss_kl_div(meandiff,covar1,covar2);
    %     kl2 = gauss_kl_div(-meandiff,covar2,covar1);
    %     lf8_kl_hf(d) = kl1+kl2;
    %     lf8_max_lik_hf(d) = hmm8_hf.max_lik/hmm8_hf.T;
    %     clear hmm*
    %
    %     load hsmm_state_seq8_seg_lfhf
    %     meandiff1 = [];
    %     meandiff2 = [];
    %     for i = 1:hmm8_lfhf.Nsegs
    %         meandiff1 = [meandiff1; hmm8_lfhf.state(2).meanfun{i}-hmm8_lfhf.state(1).meanfun{i}];
    %         meandiff2 = [meandiff2; hmm8_lfhf.state(2).meanfun2{i}-hmm8_lfhf.state(1).meanfun2{i}];
    %     end
    %     meandiff = [mean(meandiff1) mean(meandiff2)];
    %     covar1 = hmm8_lfhf.state(1).var;
    %     covar2 = hmm8_lfhf.state(2).var;
    %       kl1 = gauss_kl_div(meandiff,covar1,covar2);
    %     kl2 = gauss_kl_div(-meandiff,covar2,covar1);
    %     lf8_kl_lfhf(d) = kl1+kl2;
    %     lf8_max_lik_lfhf(d) = hmm8_lfhf.max_lik/hmm8_lfhf.T;
    %     clear hmm*
    %
    if sess_data(d).thom_elec
        load ./hsmm_state_seq4_seg_lf_4_5_2011
%         load hsmm_state_seq4_seg_lf_4_28_10_v3
        %         load hsmm_state_seq4_seg_lf_4_26_10_v2
        meandiff = [];
        for i = 1:hsmm4.Nsegs
            meandiff = [meandiff; hsmm4.state(2).meanfun{i}-hsmm4.state(1).meanfun{i}];
        end
        %         meandiff = mean(meandiff);
        covar1 = hsmm4.state(1).var;
        covar2 = hsmm4.state(2).var;
        kl1 = gauss_kl_div(meandiff',covar1,covar2);
        kl2 = gauss_kl_div(-meandiff',covar2,covar1);
        lf4_kl_lf(d) = kl1+kl2;
        lf4_bd_lf(d) = gauss_bd(meandiff',covar1,covar2);
%         lf4_he_lf(d) = mean(1 - sqrt(2*sqrt(covar1)*sqrt(covar2)/(covar1+covar2)).*exp(-0.25*meandiff/(covar1+covar2)));
%         lf4_max_lik_lf(d) = hmm4.max_lik/hmm4.T;

 
%                 load hsmm_state_seq_seg_lf_4_28_10_v3
%         %         load hsmm_state_seq4_seg_lf_4_26_10_v2
%         meandiff = [];
%         for i = 1:hmm.Nsegs
%             meandiff = [meandiff; hmm.state(2).meanfun{i}-hmm.state(1).meanfun{i}];
%         end
%         %         meandiff = mean(meandiff);
%         covar1 = hmm.state(1).var;
%         covar2 = hmm.state(2).var;
%         kl1 = gauss_kl_div(meandiff',covar1,covar2);
%         kl2 = gauss_kl_div(-meandiff',covar2,covar1);
%         mp_kl_lf(d) = kl1+kl2;
%         mp_max_lik_lf(d) = hmm4.max_lik/hmm4.T;

load hsmm_state_seq4_seg_hf_4_5_2011
%         load hsmm_state_seq4_seg_hf_4_28_10_v3
        %         load hsmm_state_seq4_seg_hf_4_26_10_v2
        meandiff = [];
        for i = 1:hsmm4_hf.Nsegs
            meandiff = [meandiff; hsmm4_hf.state(2).meanfun{i}-hsmm4_hf.state(1).meanfun{i}];
        end
        %         meandiff = mean(meandiff);
        covar1 = hsmm4_hf.state(1).var;
        covar2 = hsmm4_hf.state(2).var;
        kl1 = gauss_kl_div(meandiff',covar1,covar2);
        kl2 = gauss_kl_div(-meandiff',covar2,covar1);
        lf4_kl_hf(d) = kl1+kl2;
        lf4_bd_hf(d) = gauss_bd(meandiff',covar1,covar2);
%        lf4_he_hf(d) = mean(1 - sqrt(2*sqrt(covar1)*sqrt(covar2)/(covar1+covar2)).*exp(-0.25*meandiff/(covar1+covar2)));
        %         lf4_max_lik_hf(d) = hmm4_hf.max_lik/hmm4_hf.T;
        clear hmm*
        
load ./hsmm_state_seq4_seg_lfhf_4_5_2011
        %         load hsmm_state_seq4_seg_lfhf_4_28_10_v3_3
        meandiff = [];
        for i = 1:hsmm4_lfhf.Nsegs
            meandiff = [meandiff; hsmm4_lfhf.state(2).meanfun{i}-hsmm4_lfhf.state(1).meanfun{i}];
        end
        %         meandiff = [mean(meandiff1) mean(meandiff2)];
%         meandiff = [meandiff1 meandiff2];
        covar1 = hsmm4_lfhf.state(1).var;
        covar2 = hsmm4_lfhf.state(2).var;
        kl1 = gauss_kl_div(meandiff',covar1,covar2);
        kl2 = gauss_kl_div(-meandiff',covar2,covar1);
        lf4_kl_lfhf(d) = kl1+kl2;
        lf4_bd_lfhf(d) = gauss_bd(meandiff',covar1,covar2);
%         lf4_max_lik_lfhf(d) = hmm4_lfhf.max_lik/hmm4_lfhf.T;
        clear hmm*

    end
    
%     disp(sess_data(d).name)
%     disp(lf4_kl_lf(d))
%     disp(lf4_kl_hf(d))
%     disp(lf4_he_lf(d))
%     disp(lf4_he_hf(d))
%     pause
end

% thom_el = find_struct_field_vals(sess_data,'thom_elec',1);
% thom_par = thom_el(find(ismember(thom_el,parietal)));
% thom_pfc = setdiff(thom_el,thom_par);
% parietal = find_struct_field_vals(sess_data,'region','parietal');
% frontal = setdiff(1:length(sess_data),parietal);
% superficial = find_struct_field_vals(sess_data,'layer','23');
% deep = setdiff(1:length(sess_data),superficial);

% cd G:\WC_Germany\parietal_cortical_2010\

% save model_comparisons mp* lf8* lf4*
%%
%% within animal analysis
% for i = 1:11
%     cur_data = find(animal_id==i);
%     lf4_kl_lfhf_a(i) = mean(lf4_kl_lfhf(cur_data));
%     lf4_kl_lf_a(i) = mean(lf4_kl_lf(cur_data));
%     lf4_kl_hf_a(i) = mean(lf4_kl_hf(cur_data));
% end