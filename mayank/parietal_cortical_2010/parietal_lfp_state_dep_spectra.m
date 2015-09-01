clear all
% close all

load G:\WC_Germany\parietal_cortical_2010\parietal_cortical_2010
load G:\WC_Germany\parietal_cortical_2010\desynch_times_mp_lf8
addpath('G:\Code\Chronux\spectral_analysis\continuous\')
addpath('G:\WC_Germany\parietal_cortical_2010\')

raw_Fs = 2016;
dsf = 8;
Fsd = raw_Fs/dsf;
params.Fs = Fsd;
params.tapers = [2 3];
params.pad = 1;
params.err = [2 0.05];
params.fpass = [5 100];
params.trialave = 1;
movingwin = [0.4 0.1];
niqf = 2016/2;
lcf = 5/niqf;
hcf = 110/niqf;
[b,a] = butter(2,[lcf hcf]);
maxlag = round(0.3*Fsd);
W = 3;

frontal = find_struct_field_vals(sess_data,'region','frontal');
prefrontal = find_struct_field_vals(sess_data,'region','prefrontal');
parietal = find_struct_field_vals(sess_data,'region','parietal');

sess_data = sess_data(parietal);
desynch_start_times = desynch_start_times(parietal);
desynch_stop_times = desynch_stop_times(parietal);

%get rid of interneurons
interneurons = find_struct_field_vals(sess_data,'cell_type','interneuron');
sess_data(interneurons) = [];
desynch_start_times(interneurons) = [];
desynch_stop_times(interneurons) = [];


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

frontal = find_struct_field_vals(sess_data,'region','frontal');
prefrontal = find_struct_field_vals(sess_data,'region','prefrontal');
parietal = find_struct_field_vals(sess_data,'region','parietal');

n = length(sess_data);
f_axis = linspace(params.fpass(1),params.fpass(2),1000);
for d = 1:n
    to_dir = sess_data(d).directory;
    to_dir(1) = 'G';
    disp(sprintf('session %d',d))
    cd(to_dir);
    
    load used_data wcv_minus_spike lf5 lf8
    wcv_f = filtfilt(b,a,wcv_minus_spike);
    wcv_f = downsample(wcv_f,dsf)/sess_data(d).gains(1);
    lf8_f = filtfilt(b,a,lf8);
    lf8_f = downsample(lf8_f,dsf)/sess_data(d).gains(8);
%     lf5_f = filtfilt(b,a,lf5);
%     lf5_f = downsample(lf5_f,dsf)/sess_data(d).gains(5);
    
    lfp_f = lf8_f;
    clear lf8* lf5*
        
    %% get used state transition times
    load hsmm_state_seq_lf_pert15
    load hsmm_state_seq8_lf_pert15
    mp_state_seq =  hsmm_bbstate_seq;
    old_t = (1:length(hsmm_bbstate_seq))/Fs_bb;
    new_t = (1:length(lfp_f))/Fsd;
    mp_state_seq = round(interp1(old_t,mp_state_seq,new_t));
    lfp_state_seq = hsmm_bbstate_seq8;
    lfp_state_seq = round(interp1(old_t,lfp_state_seq,new_t));
    
    %     mp_state_seq = fliplr(mp_state_seq);
    %     wcv_f = flipud(wcv_f);
    %     lf8_f = flipud(lf8_f);
    
    desynch_times = [desynch_start_times{d}' desynch_stop_times{d}'];
    
    [up_markers,down_markers] = get_used_state_trans(mp_state_seq,Fsd,desynch_times);
    
    [Cmn,Phimn,Smn,Smm,f_up] = coherencyc_unequal_length_trials_fixedW...
        ([wcv_f lfp_f], W, params, up_markers );
    aCmn_up(d,:) = atanh(interp1(f_up,Cmn,f_axis));   
    lSw_up(d,:) = log(interp1(f_up,Smm(:,1),f_axis));
    lS8_up(d,:) = log(interp1(f_up,Smm(:,2),f_axis));

    [Cmn,Phimn,Smn,Smm,f_down] = coherencyc_unequal_length_trials_fixedW...
        ([wcv_f lfp_f], W, params, down_markers );
    aCmn_down(d,:) = atanh(interp1(f_down,Cmn,f_axis));    
    lSw_down(d,:) = log(interp1(f_down,Smm(:,1),f_axis));
    lS8_down(d,:) = log(interp1(f_down,Smm(:,2),f_axis));
 
    
end

%%
cd G:\WC_Germany\parietal_cortical_2010\
n = length(sess_data);
save state_dep_spectra_lf8

%%
sup = find_struct_field_vals(sess_data,'layer','23');
deep = setdiff(1:n,sup);
pfc = setdiff(1:n,parietal);

sup_par = intersect(sup,parietal);
deep_par = intersect(deep,parietal);
sup_pfc = intersect(sup,pfc);
deep_pfc = intersect(deep,pfc);

cur_set = sup_par;

% figure
% h = errorbar(f_axis,mean(lSw_up(cur_set,:)),std(lSw_up(cur_set,:))/sqrt(length(cur_set))); hold on
% errorbar_tick(h,.01,'units')
% h=errorbar(f_axis,mean(lSw_down(cur_set,:)),std(lSw_down(cur_set,:))/sqrt(length(cur_set)),'r');
% errorbar_tick(h,.01,'units')

figure
h = errorbar(f_axis,mean(lS8_up(cur_set,:)),std(lS8_up(cur_set,:))/sqrt(length(cur_set))); hold on
errorbar_tick(h,.01,'units')
h=errorbar(f_axis,mean(lS8_down(cur_set,:)),std(lS8_down(cur_set,:))/sqrt(length(cur_set)),'r');
errorbar_tick(h,.01,'units')

Sw_ldiff = lSw_up - lSw_down;
S8_ldiff = lS8_up - lS8_down;

figure
% h = errorbar(f_axis,mean(Sw_ldiff(cur_set,:)),std(Sw_ldiff(cur_set,:))/sqrt(length(cur_set))); hold on
% errorbar_tick(h,.01,'units')
h=errorbar(f_axis,mean(S8_ldiff(cur_set,:)),std(S8_ldiff(cur_set,:))/sqrt(length(cur_set)),'r');
errorbar_tick(h,.01,'units')

figure
h = errorbar(f_axis,mean(aCmn_up(cur_set,:)),std(aCmn_up(cur_set,:))/sqrt(length(cur_set))); hold on
errorbar_tick(h,.01,'units')
h=errorbar(f_axis,mean(aCmn_down(cur_set,:)),std(aCmn_down(cur_set,:))/sqrt(length(cur_set)),'r');
errorbar_tick(h,.01,'units')