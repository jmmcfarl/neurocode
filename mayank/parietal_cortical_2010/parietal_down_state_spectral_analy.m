clear all
close all

load G:\WC_Germany\parietal_cortical_2010\parietal_cortical_2010
load G:\WC_Germany\parietal_cortical_2010\desynch_times_mp_lf8_2_24_09
addpath('G:\Code\Chronux\spectral_analysis\continuous\')
addpath('G:\WC_Germany\parietal_cortical_2010\')
addpath('G:\WC_Germany\hsmm_state_detection\')

raw_Fs = 2016;
dsf = 8;
Fsd = raw_Fs/dsf;
movingwin = [1 1];
niqf = 2016/2;
lcf = 0.05/niqf;
lcf2 = 1.5/niqf;
hcf = 45/niqf;
hcf2 = 8/niqf;
[b,a] = butter(2,[lcf hcf]);
[b2,a2] = butter(2,[lcf2 hcf2]);
back_wins = 0;
state_edges = [0 round(Fsd*0.15)];
maxlag = round(0.3*Fsd);


%get rid of interneurons
interneurons = find_struct_field_vals(sess_data,'cell_type','interneuron');
sess_data(interneurons) = [];
desynch_times(interneurons) = [];

frontal = find_struct_field_vals(sess_data,'region','frontal');
prefrontal = find_struct_field_vals(sess_data,'region','prefrontal');
parietal = find_struct_field_vals(sess_data,'region','parietal');

% thomas_el = find_struct_field_vals(sess_data,'thom_elec',1);
% sess_data = sess_data(thomas_el);
% desynch_times = desynch_times(thomas_el);

n = length(sess_data);

% for d = 1:n
    d=9
    to_dir = sess_data(d).directory;
    to_dir(1) = 'G';
disp(sess_data(d).name)
cd(to_dir);
    
    load used_data wcv_minus_spike lf8 lf7
%     lf8 = lf7;
    wcv_f = filtfilt(b,a,wcv_minus_spike);
    wcv_f = downsample(wcv_f,dsf)/sess_data(d).gains(1);
%     wcv_f2 = filtfilt(b2,a2,wcv_minus_spike);
%     wcv_f2 = downsample(wcv_f2,dsf)/sess_data(d).gains(1);
    lf8_f = filtfilt(b,a,lf8);
    lf8_f = downsample(lf8_f,dsf)/sess_data(d).gains(8);
%     lf8_f2 = filtfilt(b2,a2,lf8)/sess_data(d).gains(8);
%     lf8_f2 = downsample(lf8_f2,dsf);
%     lf3_f2 = filtfilt(b2,a2,lf2);
%     lf3_f2 = downsample(lf3_f2,dsf)/sess_data(d).gains(3);
    if sess_data(d).thom_elec
        lf4_f = filtfilt(b,a,lf4);
        lf4_f = downsample(lf4_f,dsf)/sess_data(d).gains(4);
    end
        
    t_axis = (1:length(lf8_f))/Fsd;
    
    plot(t_axis,lf8_f,'k'), hold on
    plot(t_axis,wcv_f/30)
%     plot(t_axis,lf8_f2*2,'r')
    
%     %% get used state transition times
% %     load hsmm_state_seq_seg_lf_pert15
%     load hsmm_state_seq8_seg_lf_pert15
% %     mp_state_seq =  hsmm_bbstate_seq;
%     lfp_state_seq = hsmm_bbstate_seq8;
%     
%     state_seq = lfp_state_seq;
%     hmm = hmm8;    
    
%     cur_down_inds = [];
%     UDS_segs = (hmm.UDS_segs-1)*5+1;
% % UDS_segs = hmm.UDS_segs;
%     for i = 1:hmm.Nsegs
%        cur_downs = UDS_segs(i,1) + find(state_seq{i}(1:end-1) == 2 & state_seq{i}(2:end) == 1);
%        cur_ups = UDS_segs(i,1) + find(state_seq{i}(1:end-1) == 1 & state_seq{i}(2:end) == 2);
%        cur_ups(cur_ups < cur_downs(1)) = [];
%        cur_downs(cur_downs > cur_ups(end)) = [];
%        cur_down_inds = [cur_down_inds; cur_downs(:) cur_ups(:)];
%     end
   
%     F = linspace(1,20,100);
%     n_downs = size(cur_down_inds,1);
%     lf8_Pxx = nan(n_downs,length(F));
%     wcv_Pxx = nan(n_downs,length(F));
%     down_durs = (cur_down_inds(:,2)-cur_down_inds(:,1))/Fsd;
%     for i = 1:n_downs
%         i
%         if down_durs(i) > 1
%             lf8_Pxx(i,:) = pwelch(detrend(lf8_f(cur_down_inds(i,1):cur_down_inds(i,2))),round(Fsd*1),[],F,Fsd);
%             wcv_Pxx(i,:) = pwelch(detrend(wcv_f(cur_down_inds(i,1):cur_down_inds(i,2))),round(Fsd*1),[],F,Fsd);
%         end
%     end
% 
% end

%%
% cd G:\WC_Germany\parietal_cortical_2010\
% save parietal_lf8_utrig_mtgrams_200ms_xc
% 
% n = length(sess_data);
% %%
% thom_el = find_struct_field_vals(sess_data,'thom_elec',1);
% thom_par = thom_el(find(ismember(thom_el,parietal)));
% thom_pfc = setdiff(thom_el,thom_par);
% parietal = find_struct_field_vals(sess_data,'region','parietal');
% frontal = setdiff(1:length(sess_data),parietal);
% superficial = find_struct_field_vals(sess_data,'layer','23');
% deep = setdiff(1:length(sess_data),superficial);
% sup_par = intersect(parietal,superficial);
% sup_deep = intersect(parietal,deep);
% all = 1:n-1;
% 
% cur_set = thom_el;
% 
% %%
