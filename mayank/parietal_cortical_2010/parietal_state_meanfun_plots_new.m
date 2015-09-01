clear all
close all

addpath('G:\WC_Germany\parietal_cortical_2010\')

cd G:\WC_Germany\parietal_cortical_2010\
load parietal_cortical_2010
load G:\WC_Germany\parietal_cortical_2010\desynch_times_individual


frontal = find_struct_field_vals(sess_data,'region','frontal');
prefrontal = find_struct_field_vals(sess_data,'region','prefrontal');
parietal = find_struct_field_vals(sess_data,'region','parietal');

% sess_data = sess_data(parietal);
% desynch_times = desynch_times(parietal);

%get rid of interneurons
interneurons = find_struct_field_vals(sess_data,'cell_type','interneuron');
sess_data(interneurons) = [];
desynch_times_lf8(interneurons) = [];
desynch_times_lf4(interneurons) = [];
desynch_times_mp(interneurons) = [];


frontal = find_struct_field_vals(sess_data,'region','frontal');
prefrontal = find_struct_field_vals(sess_data,'region','prefrontal');
parietal = find_struct_field_vals(sess_data,'region','parietal');
thom_el = find_struct_field_vals(sess_data,'thom_elec',1);
thom_par = thom_el(find(ismember(thom_el,parietal)));
thom_pfc = setdiff(thom_el,thom_par);
sess_data = sess_data(thom_el);


raw_Fs = 2016;
dsf = 40;
Fsd = raw_Fs/dsf;
lf_lcf = 0.05;
lf_hcf = 2;
hf_lcf = 20; %low cut-off for the low-freq filter
hf_hcf = 80; %high cut-off for the low-freq filter
hf_smooth = 0.15;

n = length(sess_data);

for d = 1:n
    cdir = sess_data(d).directory;
    cdir(1) = 'G';
    cd(cdir)
    
    load used_data lf4 wcv_minus_spike
    
%     load hsmm_state_seq_seg_lf_4_28_10_v1
%     load hsmm_state_seq4_seg_lf_4_28_10_v3
load ./hsmm_state_seq4_seg_lf_4_5_2011    
    
%for LF4
[lf4_lf,t_axis,Fs] = get_lf_features(lf4,raw_Fs,Fsd,[lf_lcf lf_hcf]);
[obs_dist,obsrange,win_t] = compute_slide_dist(lf4_lf,hsmm4.windowSize,hsmm4.windowSlide,Fsd);
pcolor(win_t,obsrange,log(obs_dist'));shading flat
caxis([-4 0]),hold on
for i = 1:hmm4.Nsegs
    t_axis = (hsmm4.UDS_segs(i,1):hsmm4.UDS_segs(i,2))/Fsd;
    plot(t_axis,hsmm4.state(1).meanfun{i},'w','linewidth',2)
    plot(t_axis,hsmm4.state(2).meanfun{i},'w','linewidth',2)
end
cname = strcat(sess_data(d).region,'_',sess_data(d).layer,'_',sess_data(d).name);
t_names = ['G:\WC_Germany\parietal_cortical_2010\state_meanfuns_new\lf4_lf_' cname];
print('-dpng',t_names), close
    
    
    
    
    
    %% for HF
%     load hsmm_state_seq_seg_hf_4_28_10_v1
%     if sess_data(d).thom_elec
%         load hsmm_state_seq4_seg_hf_4_10_10
%     end
    
%     %for MP
%     [wcv_hf,t_axis,Fs] = get_hf_features(wcv_minus_spike,raw_Fs,Fsd,[hf_lcf hf_hcf],hf_smooth);
%     [obs_dist,obsrange,win_t] = compute_slide_dist(wcv_hf,hmm_hf.windowSize,hmm_hf.windowSlide,Fsd);
%     pcolor(win_t,obsrange,log(obs_dist'));shading flat
%     caxis([-4 0]),hold on
%     for i = 1:hmm_hf.Nsegs
%         t_axis = (hmm_hf.UDS_segs(i,1):hmm_hf.UDS_segs(i,2))/Fsd;
%         plot(t_axis,hmm_hf.state(1).meanfun{i},'w','linewidth',2)
%         plot(t_axis,hmm_hf.state(2).meanfun{i},'w','linewidth',2)
%     end
%     cname = strcat(sess_data(d).region,'_',sess_data(d).layer,'_',sess_data(d).name);
%     t_names = ['G:\WC_Germany\parietal_cortical_2010\state_meanfuns\mp_hf_' cname];
%     print('-dpng',t_names), close
    
%     %for LF8
%     [lf8_hf,t_axis,Fs] = get_hf_features(lf8,raw_Fs,Fsd,[hf_lcf hf_hcf],hf_smooth);
%     [obs_dist,obsrange,win_t] = compute_slide_dist(lf8_hf,hmm8_hf.windowSize,hmm8_hf.windowSlide,Fsd);
%     pcolor(win_t,obsrange,log(obs_dist'));shading flat
%     caxis([-4 0]),hold on
%     for i = 1:hmm8_hf.Nsegs
%         t_axis = (hmm8_hf.UDS_segs(i,1):hmm8_hf.UDS_segs(i,2))/Fsd;
%         plot(t_axis,hmm8_hf.state(1).meanfun{i},'w','linewidth',2)
%         plot(t_axis,hmm8_hf.state(2).meanfun{i},'w','linewidth',2)
%     end
%     cname = strcat(sess_data(d).region,'_',sess_data(d).layer,'_',sess_data(d).name);
%     t_names = ['F:\WC_Germany\parietal_cortical_2010\state_meanfuns\lf8_hf_' cname];
%     print('-dpng',t_names), close

%         for LF4
%             load hsmm_state_seq4_seg_hf_4_28_10_v2
load ./hsmm_state_seq4_seg_hf_4_5_2011.mat
        [lf4_hf,t_axis,Fs] = get_hf_features(lf4,raw_Fs,Fsd,[hf_lcf hf_hcf],hf_smooth);
        [obs_dist,obsrange,win_t] = compute_slide_dist(lf4_hf,hmm4_hf.windowSize,hmm4_hf.windowSlide,Fsd);
        pcolor(win_t,obsrange,log(obs_dist'));shading flat
        caxis([-4 0]),hold on
        for i = 1:hsmm4_hf.Nsegs
            t_axis = (hsmm4_hf.UDS_segs(i,1):hsmm4_hf.UDS_segs(i,2))/Fsd;
            plot(t_axis,hsmm4_hf.state(1).meanfun{i},'w','linewidth',2)
            plot(t_axis,hsmm4_hf.state(2).meanfun{i},'w','linewidth',2)
        end
        cname = strcat(sess_data(d).region,'_',sess_data(d).layer,'_',sess_data(d).name);
        t_names = ['G:\WC_Germany\parietal_cortical_2010\state_meanfuns_new\lf4_hf_' cname];
        print('-dpng',t_names), close
  
%     load hsmm_state_seq4_seg_lfhf_4_28_10_v2
load hsmm_state_seq4_seg_lfhf_4_5_2011
%     [lf4_hf,t_axis,Fs] = get_hf_features(lf4,raw_Fs,Fsd,[hf_lcf hf_hcf],hf_smooth);
%     [lf4_lf,t_axis,Fs] = get_lf_features(lf4,raw_Fs,Fsd,[lf_lcf lf_hcf]);
    [obs_dist,obsrange,win_t] = compute_slide_dist(lf4_lf,hsmm4_lfhf.windowSize,hsmm4_lfhf.windowSlide,Fsd);
    subplot(2,1,1)
    pcolor(win_t,obsrange,log(obs_dist'));shading flat
    caxis([-4 0]),hold on
    for i = 1:hsmm4_lfhf.Nsegs
        t_axis = (hsmm4_lfhf.UDS_segs(i,1):hsmm4_lfhf.UDS_segs(i,2))/Fsd;
        plot(t_axis,hsmm4_lfhf.state(1).meanfun{i}(:,1),'w','linewidth',2)
        plot(t_axis,hsmm4_lfhf.state(2).meanfun{i}(:,1),'w','linewidth',2)
    end
    [obs_dist,obsrange,win_t] = compute_slide_dist(lf4_hf,hsmm4_lfhf.windowSize,hsmm4_lfhf.windowSlide,Fsd);
    subplot(2,1,2)
    pcolor(win_t,obsrange,log(obs_dist'));shading flat
    caxis([-4 0]),hold on
    for i = 1:hsmm4_lfhf.Nsegs
        t_axis = (hsmm4_lfhf.UDS_segs(i,1):hsmm4_lfhf.UDS_segs(i,2))/Fsd;
        plot(t_axis,hsmm4_lfhf.state(1).meanfun{i}(:,2),'w','linewidth',2)
        plot(t_axis,hsmm4_lfhf.state(2).meanfun{i}(:,2),'w','linewidth',2)
    end
    cname = strcat(sess_data(d).region,'_',sess_data(d).layer,'_',sess_data(d).name);
    t_names = ['G:\WC_Germany\parietal_cortical_2010\state_meanfuns_new\lf4_lfhf_' cname];
    print('-dpng',t_names), close


end

cd G:\WC_Germany\parietal_cortical_2010\

