clear all
close all

addpath('G:\WC_Germany\parietal_cortical_2010\')

cd G:\WC_Germany\parietal_cortical_2010\
load parietal_cortical_2010
load G:\WC_Germany\parietal_cortical_2010\desynch_times_mp_lf8

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
hf_smooth = 0.03;

backlag = round(0.3*Fsd);
forwardlag = round(0.6*Fsd);
lags = -backlag:forwardlag;

arange = linspace(-3,3,100);
utrigd_lf8 = zeros(n,length(lags),length(arange));
utrigd_wcv = zeros(n,length(lags),length(arange));
utrigd_lf4 = zeros(n,length(lags),length(arange));
utrigd_lf8hf = zeros(n,length(lags),length(arange));
utrigd_lf4hf = zeros(n,length(lags),length(arange));

for d = 1:n
    cdir = sess_data(d).directory;
    cdir(1) = 'G';
    cd(cdir)
    pwd

    %% load data
    load used_data lf3 lf8
%     lf5_hf = get_hf_features(lf5,raw_Fs,Fsd,[high_lcf high_hcf],hf_smooth);
    lf8_hf = get_hf_features(lf8,raw_Fs,Fsd,[high_lcf high_hcf],hf_smooth);
    lf3_hf = get_hf_features(lf3,raw_Fs,Fsd,[high_lcf high_hcf],hf_smooth);
    %     wcv_hf = get_hf_features(wcv_minus_spike,raw_Fs,Fsd,[high_lcf high_hcf],hf_smooth);
%     lf5_f = get_lf_features(lf5,raw_Fs,Fsd,[wide_lcf wide_hcf]);
    lf8_f = get_lf_features(lf8,raw_Fs,Fsd,[wide_lcf wide_hcf]);
    lf3_f = get_lf_features(lf3,raw_Fs,Fsd,[wide_lcf wide_hcf]);
    %     wcv_f = get_lf_features(wcv_minus_spike,raw_Fs,Fsd,[wide_lcf wide_hcf]);
%     if sess_data(d).thom_elec
%         lf4_f = get_lf_features(lf4,raw_Fs,Fsd,[wide_lcf wide_hcf]);
%         lf4_hf = get_hf_features(lf4,raw_Fs,Fsd,[high_lcf high_hcf],hf_smooth);
%     end
    t_axis = (1:length(lf8_f))/Fsd;

    %% compute transition times
    load hsmm_state_seq_seg_lf_pert15
    load hsmm_state_seq8_seg_lf_pert15
%     mp_state_seq =  hsmm_bbstate_seq;
    %     old_t = (1:length(hsmm_bbstate_seq))/Fs_bb;
    %     new_t = (1:length(lf8_f))/Fsd;
    %     mp_state_seq = round(interp1(old_t,mp_state_seq,new_t));
    lfp_state_seq = hsmm_bbstate_seq8;
    %     lfp_state_seq = round(interp1(old_t,lfp_state_seq,new_t));

    [lfp_up_state_ind,down_state_ind] = get_used_state_trans_seg(lfp_state_seq,Fsd,hmm8.Fs,hmm8.UDS_segs,2*Fsd);
    [mp_up_state_ind,down_state_ind] = get_used_state_trans_seg(mp_state_seq,Fsd,hmm.Fs,hmm.UDS_segs,2*Fsd);
    used_trans_ind = mp_up_state_ind;
    n_trans = size(used_trans_ind,1);

    %% compute uptriggered matrices and distributions
    utrigm_lf8 = nan(n_trans,length(lags));
%     utrigm_lf4 = nan(n_trans,length(lags));
    utrigm_lf3 = nan(n_trans,length(lags));
%     utrigm_wcv = nan(n_trans,length(lags));
    utrigm_lf8h = nan(n_trans,length(lags));
%     utrigm_lf4h = nan(n_trans,length(lags));
    utrigm_lf3h = nan(n_trans,length(lags));
%     utrigm_wcvh = nan(n_trans,length(lags));
    delta_up = zeros(n_trans,1);
    for i = 1:n_trans
        %find time diff between lfp and mp up trans
        [dummy,nearest_lfp_up] = min(abs(lfp_up_state_ind-mp_up_state_ind(i)));
        delta_up(i) = (mp_up_state_ind(i)-lfp_up_state_ind(nearest_lfp_up(1)))/Fsd;
        cur_seg = used_trans_ind(i)-backlag:used_trans_ind(i)+forwardlag;
        utrigm_lf8(i,:) = lf8_f(cur_seg);
        utrigm_lf3(i,:) = lf3_f(cur_seg);
%         utrigm_wcv(i,:) = wcv_f(cur_seg);
        utrigm_lf8h(i,:) = lf8_hf(cur_seg);
        utrigm_lf3h(i,:) = lf3_hf(cur_seg);
%         utrigm_wcvh(i,:) = wcv_hf(cur_seg);
%         if sess_data(d).thom_elec
%             utrigm_lf4(i,:) = lf4_f(cur_seg);
%             utrigm_lf4h(i,:) = lf4_hf(cur_seg);
%         end
        next_trans = used_trans_ind(i,2)-used_trans_ind(i,1);
        if next_trans < forwardlag
            utrigm_lf8(i,backlag+next_trans+1:end) = nan;
            utrigm_lf3(i,backlag+next_trans+1:end) = nan;
%             utrigm_wcv(i,backlag+next_trans+1:end) = nan;
            utrigm_lf8h(i,backlag+next_trans+1:end) = nan;
            utrigm_lf3h(i,backlag+next_trans+1:end) = nan;
%             utrigm_wcvh(i,backlag+next_trans+1:end) = nan;
%             if sess_data(d).thom_elec
%                 utrigm_lf4(i,backlag+next_trans+1:end) = nan;
%                 utrigm_lf4h(i,backlag+next_trans+1:end) = nan;
%             end
        end
    end

    %     [dummy,delta_up_sort] = sort(delta_up);

%     for i = 1:length(lags)
%         utrigd_lf8(d,i,:) = ksdensity(utrigm_lf8(:,i),arange);
%         utrigd_wcv(d,i,:) = ksdensity(utrigm_wcv(:,i),arange);
%         utrigd_lf8hf(d,i,:) = ksdensity(utrigm_lf8h(:,i),arange);
%         if sess_data(d).thom_elec
%             utrigd_lf4(d,i,:) = ksdensity(utrigm_lf4(:,i),arange);
%             utrigd_lf4hf(d,i,:) = ksdensity(utrigm_lf4h(:,i),arange);
%         end
%     end
%     utrig_meanlf8(d,:) = nanmean(utrigm_lf8);
%     utrig_stdlf8(d,:) = zscore(nanstd(utrigm_lf8));
%     utrig_meanwcv(d,:) = nanmean(utrigm_wcv);
%     utrig_stdwcv(d,:) = zscore(nanstd(utrigm_wcv));
%     utrig_meanlf8hf(d,:) = nanmean(utrigm_lf8h);
%     utrig_stdlf8hf(d,:) = zscore(nanstd(utrigm_lf8h));
% 
%     if sess_data(d).thom_elec
%         utrig_meanlf4(d,:) = nanmean(utrigm_lf4);
%         utrig_stdlf4(d,:) = zscore(nanstd(utrigm_lf4));
%         utrig_meanlf4hf(d,:) = nanmean(utrigm_lf4h);
%         utrig_stdlf4hf(d,:) = zscore(nanstd(utrigm_lf4h));
%     end

    cname = strcat(sess_data(d).region,'_',sess_data(d).layer,'_',sess_data(d).name);


    subplot(2,1,1)
    pcolor(lags/Fsd,arange,squeeze(utrigd_wcv(d,:,:))');shading flat, caxis([0 1])
    line([0 0],[arange(1) arange(end)],'color','w')
    hold on
    plot(lags/Fsd,utrig_meanwcv(d,:),'w')
    plot(lags/Fsd,utrig_stdwcv(d,:),'r')
    subplot(2,1,2)
    pcolor(lags/Fsd,arange,squeeze(utrigd_lf8(d,:,:))');shading flat, caxis([0 1])
    line([0 0],[arange(1) arange(end)],'color','w')
    hold on
    plot(lags/Fsd,utrig_meanlf8(d,:),'w')
    plot(lags/Fsd,utrig_stdlf8(d,:),'r')
    t_names = ['G:\WC_Germany\parietal_cortical_2010\up_trig_lag_compare\dist_mp_lf8_' cname];
    print('-dpng',t_names), close

    subplot(2,1,1)
    pcolor(lags/Fsd,arange,squeeze(utrigd_wcv(d,:,:))');shading flat, caxis([0 1])
    line([0 0],[arange(1) arange(end)],'color','w')
    hold on
    plot(lags/Fsd,utrig_meanwcv(d,:),'w')
    plot(lags/Fsd,utrig_stdwcv(d,:),'r')
    subplot(2,1,2)
    pcolor(lags/Fsd,arange,squeeze(utrigd_lf8hf(d,:,:))');shading flat, caxis([0 1])
    line([0 0],[arange(1) arange(end)],'color','w')
    hold on
    plot(lags/Fsd,utrig_meanlf8hf(d,:),'w')
    plot(lags/Fsd,utrig_stdlf8hf(d,:),'r')
    t_names = ['G:\WC_Germany\parietal_cortical_2010\up_trig_lag_compare\dist_mp_lf8hf_' cname];
    print('-dpng',t_names), close

    if sess_data(d).thom_elec


        subplot(2,1,1)
        pcolor(lags/Fsd,arange,squeeze(utrigd_wcv(d,:,:))');shading flat, caxis([0 1])
        line([0 0],[arange(1) arange(end)],'color','w')
        hold on
        plot(lags/Fsd,utrig_meanwcv(d,:),'w')
        plot(lags/Fsd,utrig_stdwcv(d,:),'r')
        subplot(2,1,2)
        pcolor(lags/Fsd,arange,squeeze(utrigd_lf4(d,:,:))');shading flat, caxis([0 1])
        line([0 0],[arange(1) arange(end)],'color','w')
        hold on
        plot(lags/Fsd,utrig_meanlf4(d,:),'w')
        plot(lags/Fsd,utrig_stdlf4(d,:),'r')
        t_names = ['G:\WC_Germany\parietal_cortical_2010\up_trig_lag_compare\dist_mp_lf4_' cname];
        print('-dpng',t_names), close

        subplot(2,1,1)
        pcolor(lags/Fsd,arange,squeeze(utrigd_wcv(d,:,:))');shading flat, caxis([0 1])
        line([0 0],[arange(1) arange(end)],'color','w')
        hold on
        plot(lags/Fsd,utrig_meanwcv(d,:),'w')
        plot(lags/Fsd,utrig_stdwcv(d,:),'r')
        subplot(2,1,2)
        pcolor(lags/Fsd,arange,squeeze(utrigd_lf4hf(d,:,:))');shading flat, caxis([0 1])
        line([0 0],[arange(1) arange(end)],'color','w')
        hold on
        plot(lags/Fsd,utrig_meanlf4hf(d,:),'w')
        plot(lags/Fsd,utrig_stdlf4hf(d,:),'r')
        t_names = ['G:\WC_Germany\parietal_cortical_2010\up_trig_lag_compare\dist_mp_lf4hf_' cname];
        print('-dpng',t_names), close

    end
end

cd G:\WC_Germany\parietal_cortical_2010\
thomas_el = find_struct_field_vals(sess_data,'thom_elec',1);
thom_par = thomas_el(find(ismember(thomas_el,parietal)));
thom_pfc = setdiff(thomas_el,thom_par);
parietal = find_struct_field_vals(sess_data,'region','parietal');
frontal = setdiff(1:length(sess_data),parietal);
superficial = find_struct_field_vals(sess_data,'layer','23');
deep = setdiff(1:length(sess_data),superficial);


