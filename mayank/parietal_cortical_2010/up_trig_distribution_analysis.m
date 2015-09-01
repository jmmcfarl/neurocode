clear all
close all

addpath('F:\WC_Germany\parietal_cortical_2010\')

cd G:\WC_Germany\parietal_cortical_2010\
load parietal_cortical_2010
load G:\WC_Germany\parietal_cortical_2010\desynch_times_mp_lf8_2_24_09

% %get rid of interneurons
% interneurons = find_struct_field_vals(sess_data,'cell_type','interneuron');
% sess_data(interneurons) = [];
% desynch_start_times(interneurons) = [];
% desynch_stop_times(interneurons) = [];

frontal = find_struct_field_vals(sess_data,'region','frontal');
prefrontal = find_struct_field_vals(sess_data,'region','prefrontal');
parietal = find_struct_field_vals(sess_data,'region','parietal');

% sess_data = sess_data(parietal);
% desynch_times = desynch_times(parietal);

%get rid of interneurons
interneurons = find_struct_field_vals(sess_data,'cell_type','interneuron');
sess_data(interneurons) = [];
desynch_times(interneurons) = [];

% %get rid of sessions without lfp
% bad_sessions = [];
% for i = 1:length(sess_data)
%     if ismember(0,sess_data(i).gains)
%         bad_sessions = [bad_sessions i];
%     end
% end
% sess_data(bad_sessions) = [];
% desynch_start_times(bad_sessions) = [];
% desynch_stop_times(bad_sessions) = [];

n = length(sess_data);

raw_Fs = 2016;
dsf = 8;
Fsd = raw_Fs/dsf;
niqf = raw_Fs/2;
hf_lcf = 10;
hf_hcf = 40;
lf_lcf = 0.05;
lf_hcf = 40;

backlag = round(0.4*Fsd);
forwardlag = round(0.75*Fsd);
lags = -backlag:forwardlag;

arange = linspace(-3,3,100);

for d = 1:n
    cdir = sess_data(d).directory;
    cdir(1) = 'G';
    cd(cdir)
    pwd

    %% load data
    load used_data lf5 lf4 lf8 lf3 wcv_minus_spike
    [lf5_hf,t_axis,Fs] = get_causal_hf_sig(lf5,raw_Fs,Fsd,[hf_lcf hf_hcf]);
    [lf3_hf,t_axis,Fs] = get_causal_hf_sig(lf3,raw_Fs,Fsd,[hf_lcf hf_hcf]);
    [lf8_hf,t_axis,Fs] = get_causal_hf_sig(lf8,raw_Fs,Fsd,[hf_lcf hf_hcf]);
    [wcv_hf,t_axis,Fs] = get_causal_hf_sig(wcv_minus_spike,raw_Fs,Fsd,[hf_lcf hf_hcf]);
    [lf5_f,t_axis,Fs] = get_lf_features(lf5,raw_Fs,Fsd,[lf_lcf lf_hcf]);
    [lf3_f,t_axis,Fs] = get_lf_features(lf3,raw_Fs,Fsd,[lf_lcf lf_hcf]);
    [lf8_f,t_axis,Fs] = get_lf_features(lf8,raw_Fs,Fsd,[lf_lcf lf_hcf]);
    [wcv_f,t_axis,Fs] = get_lf_features(wcv_minus_spike,raw_Fs,Fsd,[lf_lcf lf_hcf]);
    if sess_data(d).thom_elec
        [lf4_hf,t_axis,Fs] = get_causal_hf_sig(lf4,raw_Fs,Fsd,[hf_lcf hf_hcf]);
        [lf4_f,t_axis,Fs] = get_lf_features(lf4,raw_Fs,Fsd,[lf_lcf lf_hcf]);
    end

    %% compute transition times    
    load hsmm_state_seq_seg_lf_pert15
    load hsmm_state_seq8_seg_lf_pert15
    if sess_data(d).thom_elec
       load hsmm_state_seq4_seg_lf_pert15 
    end
    
    mp_state_seq =  hsmm_bbstate_seq;
    lf8_state_seq = hsmm_bbstate_seq8;
    if sess_data(d).thom_elec
        lf4_state_seq = hsmm_bbstate_seq4;
    end
    
    [up_state_ind,down_state_ind] = get_used_state_trans_seg(lf8_state_seq,Fsd,hmm.Fs,hmm.UDS_segs,2*Fsd);
    used_trans_ind = up_state_ind;
    n_trans = size(used_trans_ind,1);
    
    %% compute uptriggered matrices and distributions
    utrigm_lf8 = nan(n_trans,length(lags));
    utrigm_lf5 = nan(n_trans,length(lags));
    utrigm_lf3 = nan(n_trans,length(lags));
    utrigm_wcv = nan(n_trans,length(lags));
    utrigm_lf8h = nan(n_trans,length(lags));
    utrigm_lf5h = nan(n_trans,length(lags));
    utrigm_lf3h = nan(n_trans,length(lags));
    utrigm_wcvh = nan(n_trans,length(lags));
    if sess_data(d).thom_elec
        utrigm_lf4 = nan(n_trans,length(lags));
        utrigm_lf4h = nan(n_trans,length(lags));
    end
    for i = 1:n_trans
        cur_seg = used_trans_ind(i)-backlag:used_trans_ind(i)+forwardlag;
        utrigm_lf8(i,:) = lf8_f(cur_seg);
        utrigm_lf5(i,:) = lf5_f(cur_seg);
        utrigm_lf3(i,:) = lf3_f(cur_seg);
        utrigm_wcv(i,:) = wcv_f(cur_seg);
        utrigm_lf8h(i,:) = lf8_hf(cur_seg);
        utrigm_lf5h(i,:) = lf5_hf(cur_seg);
        utrigm_lf3h(i,:) = lf3_hf(cur_seg);
        utrigm_wcvh(i,:) = wcv_hf(cur_seg);
        if sess_data(d).thom_elec
            utrigm_lf4(i,:) = lf4_f(cur_seg);
            utrigm_lf4h(i,:) = lf4_hf(cur_seg);
        end
        next_trans = used_trans_ind(i,2)-used_trans_ind(i,1);
        if next_trans < forwardlag
            utrigm_lf8(i,backlag+next_trans+1:end) = nan;
            utrigm_lf5(i,backlag+next_trans+1:end) = nan;
            utrigm_lf3(i,backlag+next_trans+1:end) = nan;
            utrigm_wcv(i,backlag+next_trans+1:end) = nan;
            utrigm_lf8h(i,backlag+next_trans+1:end) = nan;
            utrigm_lf5h(i,backlag+next_trans+1:end) = nan;
            utrigm_lf3h(i,backlag+next_trans+1:end) = nan;
            utrigm_wcvh(i,backlag+next_trans+1:end) = nan;
            if sess_data(d).thom_elec
                utrigm_lf4(i,backlag+next_trans+1:end) = nan;
                utrigm_lf4h(i,backlag+next_trans+1:end) = nan;
            end
        end
    end
    

%     for i = 1:length(lags)
%         utrigd_lf8(d,:,i) = ksdensity(utrigm_lf8(:,i),arange);
%         utrigd_lf5(d,:,i) = ksdensity(utrigm_lf5(:,i),arange);
%         utrigd_lf4(d,:,i) = ksdensity(utrigm_lf4(:,i),arange);
%     end

    utrigmean_lf8(d,:) = nanmean(utrigm_lf8);
    utrigstd_lf8(d,:) = nanstd(utrigm_lf8h);
    utrigmean_lf5(d,:) = nanmean(utrigm_lf5);
    utrigstd_lf5(d,:) = nanstd(utrigm_lf5h);
    utrigmean_lf3(d,:) = nanmean(utrigm_lf3);
    utrigstd_lf3(d,:) = nanstd(utrigm_lf3h);
    utrigmean_wcv(d,:) = nanmean(utrigm_wcv);
    utrigstd_wcv(d,:) = nanstd(utrigm_wcvh);
    if sess_data(d).thom_elec
        utrigmean_lf4(d,:) = nanmean(utrigm_lf4);
        utrigstd_lf4(d,:) = nanstd(utrigm_lf4h);
    end
    
end

cd G:\WC_Germany\parietal_cortical_2010\
thom_el = find_struct_field_vals(sess_data,'thom_elec',1);
parietal = find_struct_field_vals(sess_data,'region','parietal');
thom_par = thom_el(find(ismember(thom_el,parietal)));
thom_pfc = setdiff(thom_el,thom_par);
parietal = find_struct_field_vals(sess_data,'region','parietal');
frontal = setdiff(1:length(sess_data),parietal);
superficial = find_struct_field_vals(sess_data,'layer','23');
deep = setdiff(1:length(sess_data),superficial);

figure
h = errorbar(lags/Fsd,nanmean(utrigmean_lf8),nanstd(utrigmean_lf8)/sqrt(length(sess_data)))
errorbar_tick(h,0.001,'units')
hold on
h=errorbar(lags/Fsd,nanmean(utrigmean_wcv(parietal,:)),nanstd(utrigmean_wcv(parietal,:))/sqrt(length(parietal)),'r')
errorbar_tick(h,0.001,'units')
h=errorbar(lags/Fsd,nanmean(utrigmean_wcv(frontal,:)),nanstd(utrigmean_wcv(frontal,:))/sqrt(length(frontal)),'c')
errorbar_tick(h,0.001,'units')
h=errorbar(lags/Fsd,nanmean(utrigmean_lf5),nanstd(utrigmean_lf5)/sqrt(length(sess_data)),'g')
errorbar_tick(h,0.001,'units')
h=errorbar(lags/Fsd,nanmean(utrigmean_lf4(thom_el,:)),nanstd(utrigmean_lf4(thom_el,:))/sqrt(length(thom_el)),'k')
errorbar_tick(h,0.001,'units')
h=errorbar(lags/Fsd,nanmean(utrigmean_lf3),nanstd(utrigmean_lf3)/sqrt(length(sess_data)),'color',[0.7 0.7 0.7])
errorbar_tick(h,0.001,'units')
xlim([-0.4 0.4])

figure, hold on
errorbar(lags/Fsd,nanmean(utrigmean_wcv(superficial,:)),nanstd(utrigmean_wcv(superficial,:))/sqrt(length(superficial)))
errorbar(lags/Fsd,nanmean(utrigmean_wcv(deep,:)),nanstd(utrigmean_wcv(deep,:))/sqrt(length(deep)),'r')

figure
h=errorbar(lags/Fsd,nanmean(utrigstd_lf8),nanstd(utrigstd_lf8)/sqrt(length(sess_data)))
errorbar_tick(h,0.001,'units')
hold on
h=errorbar(lags/Fsd,nanmean(utrigstd_wcv(parietal,:)),nanstd(utrigstd_wcv(parietal,:))/sqrt(length(parietal)),'r')
errorbar_tick(h,0.001,'units')
h=errorbar(lags/Fsd,nanmean(utrigstd_wcv(frontal,:)),nanstd(utrigstd_wcv(frontal,:))/sqrt(length(frontal)),'c')
errorbar_tick(h,0.001,'units')
h=errorbar(lags/Fsd,nanmean(utrigstd_lf5),nanstd(utrigstd_lf5)/sqrt(length(sess_data)),'g')
errorbar_tick(h,0.001,'units')
h=errorbar(lags/Fsd,nanmean(utrigstd_lf4(thom_el,:)),nanstd(utrigstd_lf4(thom_el,:))/sqrt(length(thom_el)),'k')
errorbar_tick(h,0.001,'units')
h=errorbar(lags/Fsd,nanmean(utrigstd_lf3),nanstd(utrigstd_lf3)/sqrt(length(sess_data)),'color',[0.7 0.7 0.7])
errorbar_tick(h,0.001,'units')
xlim([-0.4 0.5])

figure, hold on
errorbar(lags/Fsd,nanmean(utrigstd_wcv(superficial,:)),nanstd(utrigstd_wcv(superficial,:))/sqrt(length(superficial)))
errorbar(lags/Fsd,nanmean(utrigstd_wcv(deep,:)),nanstd(utrigstd_wcv(deep,:))/sqrt(length(deep)),'r')


