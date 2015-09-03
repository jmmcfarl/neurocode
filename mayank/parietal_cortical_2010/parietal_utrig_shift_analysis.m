clear all
close all

addpath('F:\WC_Germany\parietal_cortical_2010\')

cd F:\WC_Germany\parietal_cortical_2010\
load parietal_cortical_2010
load F:\WC_Germany\parietal_cortical_2010\desynch_times_mp_lf8

% %get rid of interneurons
% interneurons = find_struct_field_vals(sess_data,'cell_type','interneuron');
% sess_data(interneurons) = [];
% desynch_start_times(interneurons) = [];
% desynch_stop_times(interneurons) = [];

frontal = find_struct_field_vals(sess_data,'region','frontal');
prefrontal = find_struct_field_vals(sess_data,'region','prefrontal');
parietal = find_struct_field_vals(sess_data,'region','parietal');

% sess_data = sess_data(parietal);
% desynch_start_times = desynch_start_times(parietal);
% desynch_stop_times = desynch_stop_times(parietal);

%get rid of interneurons
interneurons = find_struct_field_vals(sess_data,'cell_type','interneuron');
sess_data(interneurons) = [];
desynch_times(interneurons) = [];

n = length(sess_data);

raw_Fs = 2016;
dsf = 8;
Fsd = raw_Fs/dsf;
niqf = raw_Fs/2;
wide_lcf = 0.05;
wide_hcf = 40;
high_lcf = 15;
high_hcf = 40;
hf_smooth = 0.015;

backlag = round(0.75*Fsd);
forwardlag = round(0.75*Fsd);
lags = -backlag:forwardlag;


% thomas_el = find_struct_field_vals(sess_data,'thom_elec',1);
% sess_data = sess_data(thomas_el);
% desynch_times = desynch_times(thomas_el);

for d = 1:n
    cdir = sess_data(d).directory;
    cdir(1) = 'F';
    cd(cdir)
    pwd

    %% load data
    load used_data lf5 lf8 lf4 wcv_minus_spike
    lf5_hf = get_hf_features(lf5,raw_Fs,Fsd,[high_lcf high_hcf],hf_smooth);
    lf8_hf = get_hf_features(lf8,raw_Fs,Fsd,[high_lcf high_hcf],hf_smooth);
    wcv_hf = get_hf_features(wcv_minus_spike,raw_Fs,Fsd,[high_lcf high_hcf],hf_smooth);
    lf5_f = get_lf_features(lf5,raw_Fs,Fsd,[wide_lcf wide_hcf]);
    lf8_f = get_lf_features(lf8,raw_Fs,Fsd,[wide_lcf wide_hcf]);
    wcv_f = get_lf_features(wcv_minus_spike,raw_Fs,Fsd,[wide_lcf wide_hcf]);
    if sess_data(d).thom_elec
        lf4_f = get_lf_features(lf4,raw_Fs,Fsd,[wide_lcf wide_hcf]);
        lf4_hf = get_hf_features(lf4,raw_Fs,Fsd,[high_lcf high_hcf],hf_smooth);
    end
    t_axis = (1:length(lf8_f))/Fsd;

    %% compute transition times
    load hsmm_state_seq_seg_lf_pert15
    load hsmm_state_seq8_seg_lf_pert15
    load hsmm_state_seq5_seg_lf_pert15
    mp_state_seq =  hsmm_bbstate_seq;
    mphf_state_seq =  hsmm_hfstate_seq;
    lfp_state_seq = hsmm_bbstate_seq8;
    lfphf_state_seq = hsmm_hfstate_seq8;
    lf5hf_state_seq = hsmm_hfstate_seq5;
    lf5_state_seq = hsmm_bbstate_seq5;
    if sess_data(d).thom_elec
        load hsmm_state_seq4_seg_lf_pert15
        lf4_state_seq = hsmm_bbstate_seq4;
        lf4hf_state_seq = hsmm_hfstate_seq4;
    end

    [lfp_up_state_ind,down_state_ind] = get_used_state_trans_seg(lfp_state_seq,Fsd,hmm8.Fs,hmm8.UDS_segs,2*Fsd);
    lfp_up_state_ind = down_state_ind;
    [lfphf_up_state_ind,down_state_ind] = get_used_state_trans_seg(lfphf_state_seq,Fsd,hmm8.Fs,hmm8.UDS_segs,2*Fsd);
    lfphf_up_state_ind = down_state_ind;
    [lf5hf_up_state_ind,down_state_ind] = get_used_state_trans_seg(lf5hf_state_seq,Fsd,hmm5.Fs,hmm5.UDS_segs,2*Fsd);
    lf5hf_up_state_ind = down_state_ind;
    [lf5_up_state_ind,down_state_ind] = get_used_state_trans_seg(lf5_state_seq,Fsd,hmm5.Fs,hmm5.UDS_segs,2*Fsd);
    lf5_up_state_ind = down_state_ind;
    [mp_up_state_ind,down_state_ind] = get_used_state_trans_seg(mp_state_seq,Fsd,hmm.Fs,hmm.UDS_segs,2*Fsd);
    mp_up_state_ind = down_state_ind;
    [mphf_up_state_ind,down_state_ind] = get_used_state_trans_seg(mphf_state_seq,Fsd,hmm.Fs,hmm.UDS_segs,2*Fsd);
    mphf_up_state_ind = down_state_ind;
    if sess_data(d).thom_elec
        [lf4_up_state_ind,down_state_ind] = get_used_state_trans_seg(lf4_state_seq,Fsd,hmm4.Fs,hmm4.UDS_segs,2*Fsd);
        lf4_up_state_ind = down_state_ind;
        [lf4hf_up_state_ind,down_state_ind] = get_used_state_trans_seg(lf4hf_state_seq,Fsd,hmm4.Fs,hmm4.UDS_segs,2*Fsd);
        lf4hf_up_state_ind = down_state_ind;
    end
    used_trans_ind = lfp_up_state_ind;
    n_trans = size(used_trans_ind,1);

    %% compute uptriggered matrices and distributions
    utrigm_lf8 = nan(n_trans,length(lags));
    utrigm_lf4 = nan(n_trans,length(lags));
    utrigm_lf5 = nan(n_trans,length(lags));
    utrigm_wcv = nan(n_trans,length(lags));
    utrigm_lf8h = nan(n_trans,length(lags));
    utrigm_lf4h = nan(n_trans,length(lags));
    utrigm_lf5h = nan(n_trans,length(lags));
    utrigm_wcvh = nan(n_trans,length(lags));
    delta_lfphfup = zeros(n_trans,1);
    delta_mpup = zeros(n_trans,1);
    delta_mphfup = zeros(n_trans,1);
    delta_lf4up = zeros(n_trans,1);
    delta_lf4hfup = zeros(n_trans,1);
    up_dur = used_trans_ind(:,2)-used_trans_ind(:,1);
    for i = 1:n_trans
        %find time diff between lfp and mp up trans
        [dummy,nearest_lfphf_up] = min(abs(lfp_up_state_ind(i)-lfphf_up_state_ind));
        delta_lfphfup(i) = (lfp_up_state_ind(i)-lfphf_up_state_ind(nearest_lfphf_up(1)))/Fsd;
        [dummy,nearest_lf5hf_up] = min(abs(lfp_up_state_ind(i)-lf5hf_up_state_ind));
        delta_lf5hfup(i) = (lfp_up_state_ind(i)-lf5hf_up_state_ind(nearest_lf5hf_up(1)))/Fsd;
        [dummy,nearest_lf5_up] = min(abs(lfp_up_state_ind(i)-lf5_up_state_ind));
        delta_lf5up(i) = (lfp_up_state_ind(i)-lf5_up_state_ind(nearest_lf5_up(1)))/Fsd;
        [dummy,nearest_mp_up] = min(abs(lfp_up_state_ind(i)-mp_up_state_ind));
        delta_mpup(i) = (lfp_up_state_ind(i)-mp_up_state_ind(nearest_mp_up(1)))/Fsd;
        [dummy,nearest_hfmp_up] = min(abs(lfp_up_state_ind(i)-mphf_up_state_ind));
        delta_mphfup(i) = (lfp_up_state_ind(i)-mphf_up_state_ind(nearest_hfmp_up(1)))/Fsd;
        cur_seg = used_trans_ind(i,1)-backlag:used_trans_ind(i,1)+forwardlag;
        utrigm_lf8(i,:) = lf8_f(cur_seg);
        utrigm_lf5(i,:) = lf5_f(cur_seg);
        utrigm_wcv(i,:) = wcv_f(cur_seg);
        utrigm_lf8h(i,:) = lf8_hf(cur_seg);
        utrigm_lf5h(i,:) = lf5_hf(cur_seg);
        utrigm_wcvh(i,:) = wcv_hf(cur_seg);
        if sess_data(d).thom_elec
            utrigm_lf4(i,:) = lf4_f(cur_seg);
            utrigm_lf4h(i,:) = lf4_hf(cur_seg);
            [dummy,nearest_lf4_up] = min(abs(lfp_up_state_ind(i)-lf4_up_state_ind));
            delta_lf4up(i) = (lfp_up_state_ind(i)-lf4_up_state_ind(nearest_lf4_up(1)))/Fsd;
            [dummy,nearest_lf4hf_up] = min(abs(lfp_up_state_ind(i)-lf4hf_up_state_ind));
            delta_lf4hfup(i) = (lfp_up_state_ind(i)-lf4hf_up_state_ind(nearest_lf4hf_up(1)))/Fsd;
        end

        %         next_trans = used_trans_ind(i,2)-used_trans_ind(i,1);
        %         if next_trans < forwardlag
        %             utrigm_lf8(i,backlag+next_trans+1:end) = nan;
        % %             utrigm_lf5(i,backlag+next_trans+1:end) = nan;
        %             utrigm_wcv(i,backlag+next_trans+1:end) = nan;
        %             utrigm_lf8h(i,backlag+next_trans+1:end) = nan;
        % %             utrigm_lf5h(i,backlag+next_trans+1:end) = nan;
        %             utrigm_wcvh(i,backlag+next_trans+1:end) = nan;
        %             if sess_data(d).thom_elec
        %                 utrigm_lf4(i,backlag+next_trans+1:end) = nan;
        %                 utrigm_lf4h(i,backlag+next_trans+1:end) = nan;
        %             end
        %         end
    end

        [sorted_delta_lf5s,delta_up_sort_lf5s] = sort(delta_lf5up);

    [sorted_delta_mp,delta_up_sort_mp] = sort(delta_mpup);
    sorted_delta_hflfp = delta_lfphfup(delta_up_sort_mp);
    sorted_delta_hflf5 = delta_lf5hfup(delta_up_sort_mp);
    sorted_delta_lf5 = delta_lf5up(delta_up_sort_mp);
    sorted_delta_hfmp = delta_mphfup(delta_up_sort_mp);
    if sess_data(d).thom_elec
        sorted_delta_lf4 = delta_lf4up(delta_up_sort_mp);
        sorted_delta_hflf4 = delta_lf4hfup(delta_up_sort_mp);
    end
%     


%     [r8,p8] = corrcoef(sorted_delta_mp,sorted_delta_hflfp);
%     deltacorr8(d) = r8(2,1);
%     deltap8(d) = p8(2,1);
%      [r5,p5] = corrcoef(sorted_delta_mp,sorted_delta_hflf5);
%     deltacorrh5(d) = r5(2,1);
%     deltaph5(d) = p5(2,1);
%      [r5,p5] = corrcoef(sorted_delta_mp,sorted_delta_lf5);
%     deltacorr5(d) = r5(2,1);
%     deltap5(d) = p5(2,1);
% if sess_data(d).thom_elec
%          [r5,p5] = corrcoef(sorted_delta_mp,sorted_delta_lf4);
%     deltacorr4(d) = r5(2,1);
%     deltap4(d) = p5(2,1);
%      [r5,p5] = corrcoef(sorted_delta_mp,sorted_delta_hflf4);
%     deltacorrh4(d) = r5(2,1);
%     deltaph4(d) = p5(2,1);
% end
    cname = strcat(sess_data(d).region,'_',sess_data(d).layer,'_',sess_data(d).name);

%     subplot(2,1,1)
%     pcolor(lags/Fsd,1:n_trans,utrigm_lf8(delta_up_sort_mp,:));shading flat, caxis([-1.5 2.5])
%     line([0 0],[1 n_trans],'color','w'), hold on
%     subplot(2,1,2)
%     pcolor(lags/Fsd,1:n_trans,utrigm_wcv(delta_up_sort_mp,:));shading flat, caxis([-1.5 2.5])
%     line([0 0],[1 n_trans],'color','w'), hold on
%     plot(-sorted_delta_mp,1:n_trans,'k')
%     t_names = ['G:\WC_Germany\parietal_cortical_2010\up_trig_lag_compare\dur_sort2_mp_lf8_' cname];
%     print('-dpng',t_names), close

%     subplot(2,1,1)
%     pcolor(lags/Fsd,1:n_trans,utrigm_lf8);shading flat, caxis([-1.5 2.5])
%     line([0 0],[1 n_trans],'color','w'), hold on
% %     plot(-sorted_delta_hflfp,1:n_trans,'r')
%     subplot(2,1,2)
%     pcolor(lags/Fsd,1:n_trans,utrigm_lf5);shading flat, caxis([-1.5 2.5])
%     line([0 0],[1 n_trans],'color','w'), hold on
% %     plot(-sorted_delta_mp,1:n_trans,'k')
%     t_names = ['G:\WC_Germany\parietal_cortical_2010\up_trig_lag_compare\lf8_lf5_' cname];
%     print('-dpng',t_names), close

%         subplot(2,1,1)
%     pcolor(lags/Fsd,1:n_trans,utrigm_lf8);shading flat, caxis([-1.5 2.5])
%     line([0 0],[1 n_trans],'color','w'), hold on
% %     plot(-sorted_delta_hflfp,1:n_trans,'r')
%     subplot(2,1,2)
%     pcolor(lags/Fsd,1:n_trans,utrigm_lf5h);shading flat, caxis([-1.5 2.5])
%     line([0 0],[1 n_trans],'color','w'), hold on
% %     plot(-sorted_delta_mp,1:n_trans,'k')
%     t_names = ['G:\WC_Germany\parietal_cortical_2010\up_trig_lag_compare\lf8_lf5h_' cname];
%     print('-dpng',t_names), close

%         subplot(2,1,1)
%     pcolor(lags/Fsd,1:n_trans,utrigm_wcv(delta_up_sort_mp,:));shading flat, caxis([-1.5 2.5])
%     line([0 0],[1 n_trans],'color','w'), hold on
% %     plot(-sorted_delta_hflfp,1:n_trans,'r')
%     subplot(2,1,2)
%     pcolor(lags/Fsd,1:n_trans,utrigm_lf8(delta_up_sort_mp,:));shading flat, caxis([-1.5 2.5])
%     line([0 0],[1 n_trans],'color','w'), hold on
%     plot(-sorted_delta_hflfp,1:n_trans,'r')
%     plot(-sorted_delta_mp,1:n_trans,'k')
%     t_names = ['G:\WC_Germany\parietal_cortical_2010\up_trig_lag_compare\ldelta_sort_mp_lf8_' cname];
%     print('-dpng',t_names), close
% 
            subplot(3,1,1)
    pcolor(lags/Fsd,1:n_trans,utrigm_lf8(delta_up_sort_mp,:));shading flat, caxis([-1.5 2.5])
    line([0 0],[1 n_trans],'color','w'), hold on
    plot(-sorted_delta_mp,1:n_trans,'k')
    subplot(3,1,2)
    pcolor(lags/Fsd,1:n_trans,utrigm_lf5(delta_up_sort_mp,:));shading flat, caxis([-1.5 2.5])
    line([0 0],[1 n_trans],'color','w'), hold on
    plot(-sorted_delta_mp,1:n_trans,'k')
    subplot(3,1,3)
    pcolor(lags/Fsd,1:n_trans,utrigm_lf5h(delta_up_sort_mp,:));shading flat, caxis([-1.5 2.5])
    line([0 0],[1 n_trans],'color','w'), hold on
    plot(-sorted_delta_mp,1:n_trans,'k')
    t_names = ['F:\WC_Germany\parietal_cortical_2010\up_trig_lag_compare\tripledown_lf8_lf5_lf5hf' cname];
    print('-dpng',t_names), close

%     if sess_data(d).thom_elec
%                 subplot(3,1,1)
%     pcolor(lags/Fsd,1:n_trans,utrigm_lf8(delta_up_sort_mp,:));shading flat, caxis([-1.5 2.5])
%     line([0 0],[1 n_trans],'color','w'), hold on
%     plot(-sorted_delta_mp,1:n_trans,'k')
%     subplot(3,1,2)
%     pcolor(lags/Fsd,1:n_trans,utrigm_lf4(delta_up_sort_mp,:));shading flat, caxis([-1.5 2.5])
%     line([0 0],[1 n_trans],'color','w'), hold on
%     plot(-sorted_delta_mp,1:n_trans,'k')
%     subplot(3,1,3)
%     pcolor(lags/Fsd,1:n_trans,utrigm_lf4h(delta_up_sort_mp,:));shading flat, caxis([-1.5 2.5])
%     line([0 0],[1 n_trans],'color','w'), hold on
%     plot(-sorted_delta_mp,1:n_trans,'k')
%     t_names = ['G:\WC_Germany\parietal_cortical_2010\up_trig_lag_compare\tripledown_lf8_lf4_lf4hf' cname];
%     print('-dpng',t_names), close
%     end
%  
%                 subplot(2,1,1)
%     pcolor(lags/Fsd,1:n_trans,utrigm_lf8(delta_up_sort_mp,:));shading flat, caxis([-1.5 2.5])
%     line([0 0],[1 n_trans],'color','w'), hold on
%     plot(-sorted_delta_mp,1:n_trans,'k')
%     subplot(2,1,2)
%     pcolor(lags/Fsd,1:n_trans,utrigm_wcv(delta_up_sort_mp,:));shading flat, caxis([-1.5 2.5])
%     line([0 0],[1 n_trans],'color','w'), hold on
%     plot(-sorted_delta_mp,1:n_trans,'k')
%     t_names = ['G:\WC_Germany\parietal_cortical_2010\up_trig_lag_compare\down_lf8_mp' cname];
%     print('-dpng',t_names), close

%                 subplot(2,1,1)
%     pcolor(lags/Fsd,1:n_trans,utrigm_lf8(delta_up_sort_mp,:));shading flat, caxis([-1.5 2.5])
%     line([0 0],[1 n_trans],'color','w'), hold on
% %     plot(-sorted_delta_mp,1:n_trans,'k')
%     subplot(2,1,2)
%     pcolor(lags/Fsd,1:n_trans,utrigm_lf5(delta_up_sort_mp,:));shading flat, caxis([-1.5 2.5])
%     line([0 0],[1 n_trans],'color','w'), hold on
% %     plot(-sorted_delta_lf5,1:n_trans,'r')
%     plot(-sorted_delta_lf5,1:n_trans,'k')
%     t_names = ['G:\WC_Germany\parietal_cortical_2010\up_trig_lag_compare\ldelta_sort_lf8_lf5_' cname];
%     print('-dpng',t_names), close
% 
%                 subplot(2,1,1)
%     pcolor(lags/Fsd,1:n_trans,utrigm_lf8(delta_up_sort_mp,:));shading flat, caxis([-1.5 2.5])
%     line([0 0],[1 n_trans],'color','w'), hold on
% %     plot(-sorted_delta_mp,1:n_trans,'k')
%     subplot(2,1,2)
%     pcolor(lags/Fsd,1:n_trans,utrigm_lf5h(delta_up_sort_mp,:));shading flat, caxis([-1.5 2.5])
%     line([0 0],[1 n_trans],'color','w'), hold on
% %     plot(-sorted_delta_hflf5,1:n_trans,'r')
%     plot(-sorted_delta_mp,1:n_trans,'k')
%     t_names = ['G:\WC_Germany\parietal_cortical_2010\up_trig_lag_compare\ldelta_sort_lf8_lf5hf_' cname];
%     print('-dpng',t_names), close

    %     subplot(2,1,1)
%     pcolor(lags/Fsd,1:n_trans,utrigm_lf8h(delta_up_sort_mp,:));shading flat, caxis([-1.5 2.5])
%     line([0 0],[1 n_trans],'color','w'), hold on
%     plot(-sorted_delta_hflfp,1:n_trans,'r')    
%     subplot(2,1,2)
%     pcolor(lags/Fsd,1:n_trans,utrigm_wcvh(delta_up_sort_mp,:));shading flat, caxis([-1.5 2.5])
%     line([0 0],[1 n_trans],'color','w'), hold on
%     plot(-sorted_delta_hfmp,1:n_trans,'r')        
%     t_names = ['G:\WC_Germany\parietal_cortical_2010\up_trig_lag_compare\dur_sort2_hfmp_hflf8_' cname];
%     print('-dpng',t_names), close

%     if sess_data(d).thom_elec
% 
%         subplot(2,1,1)
%         pcolor(lags/Fsd,1:n_trans,utrigm_lf4(delta_up_sort_mp,:));shading flat, caxis([-1.5 2.5])
%         line([0 0],[1 n_trans],'color','w'), hold on
%         plot(-sorted_delta_lf4,1:n_trans,'r')
%         subplot(2,1,2)
%         pcolor(lags/Fsd,1:n_trans,utrigm_wcv(delta_up_sort_mp,:));shading flat, caxis([-1.5 2.5])
%         line([0 0],[1 n_trans],'color','w'), hold on
%          plot(-sorted_delta_mp,1:n_trans,'k')       
%         t_names = ['G:\WC_Germany\parietal_cortical_2010\up_trig_lag_compare\dur_sort2_mp_lf4_' cname];
%         print('-dpng',t_names), close
% 
%         subplot(2,1,1)
%         pcolor(lags/Fsd,1:n_trans,utrigm_lf4h(delta_up_sort_mp,:));shading flat, caxis([-1.5 2.5])
%         line([0 0],[1 n_trans],'color','w'), hold on
%           plot(-sorted_delta_hflf4,1:n_trans,'r')      
%         subplot(2,1,2)
%         pcolor(lags/Fsd,1:n_trans,utrigm_wcv(delta_up_sort_mp,:));shading flat, caxis([-1.5 2.5])
%         line([0 0],[1 n_trans],'color','w'), hold on
%         plot(-sorted_delta_mp,1:n_trans,'k')
%         t_names = ['G:\WC_Germany\parietal_cortical_2010\up_trig_lag_compare\dur_sort2_mp_hf2lf4_' cname];
%         print('-dpng',t_names), close
% 
%         subplot(2,1,1)
%         pcolor(lags/Fsd,1:n_trans,utrigm_lf4(delta_up_sort_mp,:));shading flat, caxis([-1.5 2.5])
%         line([0 0],[1 n_trans],'color','w'), hold on
%           plot(-sorted_delta_lf4,1:n_trans,'r')              
%         subplot(2,1,2)
%         pcolor(lags/Fsd,1:n_trans,utrigm_lf8(delta_up_sort_mp,:));shading flat, caxis([-1.5 2.5])
%         line([0 0],[1 n_trans],'color','w')
%         t_names = ['G:\WC_Germany\parietal_cortical_2010\up_trig_lag_compare\dur_sort2_lf4_lf8_' cname];
%         print('-dpng',t_names), close
% 
%     end

clear delta_*

end

cd F:\WC_Germany\parietal_cortical_2010\
thomas_el = find_struct_field_vals(sess_data,'thom_elec',1);
thom_par = thomas_el(find(ismember(thomas_el,parietal)));
thom_pfc = setdiff(thomas_el,thom_par);
parietal = find_struct_field_vals(sess_data,'region','parietal');
frontal = setdiff(1:length(sess_data),parietal);
superficial = find_struct_field_vals(sess_data,'layer','23');
deep = setdiff(1:length(sess_data),superficial);


