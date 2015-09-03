clear all
% close all

load F:\WC_Germany\parietal_cortical_2010\parietal_cortical_2010
load F:\WC_Germany\parietal_cortical_2010\desynch_times_mp_lf8
addpath('F:\Code\Chronux\spectral_analysis\continuous\')

raw_Fs = 2016;
dsf = 8;
Fsd = raw_Fs/dsf;
params.Fs = Fsd;
params.tapers = [2 3];
params.err = [2 0.05];
params.fpass = [3 95];
params.trialave = 1;
movingwin = [0.5 0.5];
niqf = 2016/2;
lcf = 2/niqf;
hcf = 110/niqf;
[b,a] = butter(4,[lcf hcf]);
delW = params.tapers(1)/movingwin(1);

frontal = find_struct_field_vals(sess_data,'region','frontal');
prefrontal = find_struct_field_vals(sess_data,'region','prefrontal');
parietal = find_struct_field_vals(sess_data,'region','parietal');

% sess_data = sess_data(parietal);
% desynch_start_times = desynch_start_times(parietal);
% desynch_stop_times = desynch_stop_times(parietal);

%get rid of interneurons
interneurons = find_struct_field_vals(sess_data,'cell_type','interneuron');
sess_data(interneurons) = [];
desynch_start_times(interneurons) = [];
desynch_stop_times(interneurons) = [];

%get rid of sessions without lfp
bad_sessions = [];
for i = 1:length(sess_data)
    if ismember(0,sess_data(i).gains)
        bad_sessions = [bad_sessions i];
    end
end
sess_data(bad_sessions) = [];
desynch_start_times(bad_sessions) = [];
desynch_stop_times(bad_sessions) = [];

frontal = find_struct_field_vals(sess_data,'region','frontal');
prefrontal = find_struct_field_vals(sess_data,'region','prefrontal');
parietal = find_struct_field_vals(sess_data,'region','parietal');

for d = 1:length(sess_data)
    to_dir = sess_data(d).directory;
    to_dir(1) = 'F';
    disp(sprintf('session %d',d))
    cd(to_dir);
    
    load used_data wcv_minus_spike lf8 lf7 lf6 lf5 lf4
    wcv_f = filtfilt(b,a,wcv_minus_spike);
    wcv_f = downsample(wcv_f,dsf)/sess_data(d).gains(1);
    lf8_f = filtfilt(b,a,lf8);
    lf8_f = downsample(lf8_f,dsf)/sess_data(d).gains(8);
    lf7_f = filtfilt(b,a,lf7);
    lf7_f = downsample(lf7_f,dsf)/sess_data(d).gains(7);
    lf6_f = filtfilt(b,a,lf6);
    lf6_f = downsample(lf6_f,dsf)/sess_data(d).gains(6);
    lf5_f = filtfilt(b,a,lf5);
    lf5_f = downsample(lf5_f,dsf)/sess_data(d).gains(5);
    lf4_f = filtfilt(b,a,lf4);
    lf4_f = downsample(lf4_f,dsf)/sess_data(d).gains(4);
    
    %% get used state transition times
    load hsmm_state_seq_lf_pert15
    %     load hsmm_state_seq8_lf_pert15
    mp_state_seq =  hsmm_bbstate_seq;
    old_t = (1:length(hsmm_bbstate_seq))/Fs_bb;
    new_t = (1:length(lf8_f))/Fsd;
    mp_state_seq = round(interp1(old_t,mp_state_seq,new_t));
    %     lfp_state_seq = hsmm_bbstate_seq8;
    %     lfp_state_seq = round(interp1(old_t,lfp_state_seq,new_t));
    
    desynch_times = [desynch_start_times{d}' desynch_stop_times{d}'];
    [up_segs,down_segs,up_nsegs,down_nsegs] = ...
        get_parsed_state_segments(mp_state_seq,Fsd,movingwin(1),0.1,desynch_times);
    
    up_markers1 = [up_segs(:) up_segs(:)+round(Fsd*movingwin(1))-1];
    down_markers1 = [down_segs(:) down_segs(:)+round(Fsd*movingwin(1))-1];
    % n_ups = size(up_markers1,1);
    % up_markers2 = up_markers1(randperm(n_ups),:);
    % n_downs = size(down_markers1,1);
    % down_markers2 = down_markers1(randperm(n_downs),:);
    up_markers2 = up_markers1;
    down_markers2 = down_markers1;
    
    %     [Cw8_wup(d,:),Phiw8_wup(d,:),Smn,S_up(d,:,:),f,confC_wup(d),dof_wup(d)] = ...
    %         coherencyc_unequal_length_trials_jmm_fixedsegs([wcv_f lf8_f], movingwin, params, up_markers1, up_markers2);
    %     [Cw8_wdown(d,:),Phiw8_wdown(d,:),Smn,S_down(d,:,:),f,confC_wdown(d),dof_wdown(d)] = ...
    %         coherencyc_unequal_length_trials_jmm_fixedsegs([wcv_f lf8_f], movingwin, params, down_markers1, down_markers2);
    
    [Smm,f] = dualspectra_unequal_length_trials_jmm_fixedsegs...
        ([lf8_f lf7_f], movingwin, params, up_markers1, up_markers2 );
    up_ldiff_87(d,:) = mean(log(Smm(:,:,1)) - log(Smm(:,:,2)));
    [Smm,f] = dualspectra_unequal_length_trials_jmm_fixedsegs...
        ([lf8_f lf7_f], movingwin, params, down_markers1, down_markers2 );
    down_ldiff_87(d,:) = mean(log(Smm(:,:,1)) - log(Smm(:,:,2)));
    
     [Smm,f] = dualspectra_unequal_length_trials_jmm_fixedsegs...
        ([lf8_f lf6_f], movingwin, params, up_markers1, up_markers2 );
    up_ldiff_86(d,:) = mean(log(Smm(:,:,1)) - log(Smm(:,:,2)));
    [Smm,f] = dualspectra_unequal_length_trials_jmm_fixedsegs...
        ([lf8_f lf6_f], movingwin, params, down_markers1, down_markers2 );
    down_ldiff_86(d,:) = mean(log(Smm(:,:,1)) - log(Smm(:,:,2)));

        [Smm,f] = dualspectra_unequal_length_trials_jmm_fixedsegs...
        ([lf8_f lf5_f], movingwin, params, up_markers1, up_markers2 );
    up_ldiff_85(d,:) = mean(log(Smm(:,:,1)) - log(Smm(:,:,2)));
    [Smm,f] = dualspectra_unequal_length_trials_jmm_fixedsegs...
        ([lf8_f lf5_f], movingwin, params, down_markers1, down_markers2 );
    down_ldiff_85(d,:) = mean(log(Smm(:,:,1)) - log(Smm(:,:,2)));

        [Smm,f] = dualspectra_unequal_length_trials_jmm_fixedsegs...
        ([lf8_f lf4_f], movingwin, params, up_markers1, up_markers2 );
    up_ldiff_84(d,:) = mean(log(Smm(:,:,1)) - log(Smm(:,:,2)));
    [Smm,f] = dualspectra_unequal_length_trials_jmm_fixedsegs...
        ([lf8_f lf4_f], movingwin, params, down_markers1, down_markers2 );
    down_ldiff_84(d,:) = mean(log(Smm(:,:,1)) - log(Smm(:,:,2)));

end

figure
h = errorbar(f,mean(up_ldiff_87),std(up_ldiff_87)/sqrt(n));
hold on
errorbar_tick(h,.01,'units')
h=errorbar(f,mean(up_ldiff_86),std(up_ldiff_86)/sqrt(n),'r');
errorbar_tick(h,.01,'units')
h=errorbar(f,mean(up_ldiff_85),std(up_ldiff_85)/sqrt(n),'g');
errorbar_tick(h,.01,'units')
h=errorbar(f,mean(up_ldiff_84),std(up_ldiff_84)/sqrt(n),'k');
errorbar_tick(h,.01,'units')
legend('87','86','85','84')

figure
h = errorbar(f,mean(down_ldiff_87),std(down_ldiff_87)/sqrt(n));
hold on
errorbar_tick(h,.01,'units')
h=errorbar(f,mean(down_ldiff_86),std(down_ldiff_86)/sqrt(n),'r');
errorbar_tick(h,.01,'units')
h=errorbar(f,mean(down_ldiff_85),std(down_ldiff_85)/sqrt(n),'g');
errorbar_tick(h,.01,'units')
h=errorbar(f,mean(down_ldiff_84),std(down_ldiff_84)/sqrt(n),'k');
errorbar_tick(h,.01,'units')
legend('87','86','85','84')

