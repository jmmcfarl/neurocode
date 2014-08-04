clear all
close all

ExptNum = 239;
load('~/Data/bruce/2_27_12/EyeCalibration/EyeCalParas.mat')
cd(sprintf('~/Data/bruce/M%d',ExptNum));
addpath('~/James_scripts/general_functions/');
addpath('~/James_scripts/bruce/2_27_scripts');

if ExptNum == 232
    bar_expts = [37 38 39 43 46 47]; %232
elseif ExptNum == 235
    bar_expts = [51:60]; %235
elseif ExptNum == 239
    bar_expts = [37 38 39 40 43 44 49]; %239
end

if ExptNum == 235
    bar_expts(bar_expts==51) = []; %this has different set of stim positions
end
if ExptNum == 239
    bar_expts(bar_expts==40) = []; %this has different set of stim positions
end


%%
blink_hcf = 10;
blink_peakthresh = 10;
blink_edgethresh = 5;
min_blink_dur = 0.075;
max_blink_dur = Inf;

sac_peakthresh = 10;
sac_edgethresh = 5;
min_sac_dur = 0.01;
max_sac_dur = 0.05;
for ee=1:length(bar_expts)
    load(sprintf('lemM%d.%d.em.mat',ExptNum,bar_expts(ee)))
    Eyedt = Expt.Header.CRrates(1); % resolution of eye signal (sec)
    
    niqf = 1/2/Eyedt;
    [blink_b,blink_a] = butter(2,blink_hcf/niqf,'low');
    
    n_trials(ee) = length(Expt.Trials);
    
    all_eye_ts{ee} = [];
    all_reye_vals{ee} = [];
    all_leye_vals{ee} = [];
    all_eyetrial_start_inds{ee} = [];
    all_blink_startinds{ee} = [];
    all_blink_stopinds{ee} = [];
    all_sac_startinds{ee} = [];
    all_sac_stopinds{ee} = [];
    all_sac_peakamps{ee} = [];
    all_eye_speed{ee} = [];
    all_deye_speed{ee} = [];
    all_dfeye_speed{ee} = [];
    
    cur_start_ts = [Expt.Trials(:).ftime]/1e4;
    % cur_start_ts = [Expt.Trials(:).bstime];
    cur_end_ts = [Expt.Trials(:).End]/1e4;
    cur_real_end_ts = [Expt.Trials(:).End]/1e4;
    for tt = 1:n_trials(ee)
        cur_nsamps = length(Expt.Trials(tt).Eyevals.rh);
        cur_real_end_ts(tt) = cur_start_ts(tt) + Eyedt*(cur_nsamps);
        cur_eyets = cur_start_ts(tt):Eyedt:cur_real_end_ts(tt);
        cur_reye_vals = [Expt.Trials(tt).Eyevals.rh Expt.Trials(tt).Eyevals.rv];
        cur_eyets(size(cur_reye_vals,1)+1:end) = [];
        
        r_len = size(cur_reye_vals,1);
        cur_leye_vals = nan(size(cur_reye_vals));
        cur_leye_vals(:,1) = Expt.Trials(tt).Eyevals.lh;
        cur_l_len = length(Expt.Trials(tt).Eyevals.lv);
        if cur_l_len < r_len
            cur_leye_vals(1:cur_l_len,2) = Expt.Trials(tt).Eyevals.lv;
        else
            cur_leye_vals(:,2) = Expt.Trials(tt).Eyevals.lv(1:r_len);
        end
        
        if ~isempty(all_eye_ts{ee})
            cur_sp = find(cur_eyets > max(all_eye_ts{ee}),1,'first');
        else
            cur_sp = 1;
        end
        cur_eyets = cur_eyets(cur_sp:end);
        cur_reye_vals = cur_reye_vals(cur_sp:end,:);
        cur_leye_vals = cur_leye_vals(cur_sp:end,:);
        cur(tt) = length(cur_eyets);
        
        if cur(tt) > 0
            
            avg_eye_vals = 0.5*cur_reye_vals+0.5*cur_leye_vals;
            avg_eye_vals(isnan(cur_leye_vals)) = cur_reye_vals(isnan(cur_leye_vals));
            
            eye_vels = [0 0; diff(avg_eye_vals)]/Eyedt;
            eye_speed = sqrt(sum(eye_vels.^2,2));
            
            diff_eye_vals = cur_reye_vals - cur_leye_vals;
            diff_eye_vels = [0 0; diff(diff_eye_vals)]/Eyedt;
            diff_eye_speed = max(abs(diff_eye_vels),[],2);
            to_filt = diff_eye_speed; to_filt(isnan(to_filt)) = 0;
            filt_eye_speed = filtfilt(blink_b,blink_a,to_filt);
            cur_blink_stats = locate_peaks_jmm(filt_eye_speed,blink_peakthresh,blink_edgethresh,...
                round(min_blink_dur/Eyedt),round(max_blink_dur/Eyedt),ones(size(filt_eye_speed)));
            
            cur_blink_starts = [cur_blink_stats(:).start_loc];
            cur_blink_stops = [cur_blink_stats(:).stop_loc];
            blink_vec = ones(size(eye_speed));
            for bb = 1:length(cur_blink_starts)
                blink_vec(cur_blink_starts(bb):cur_blink_stops(bb)) = 0;
            end
            
            cur_sac_stats = locate_peaks_jmm(eye_speed,sac_peakthresh,sac_edgethresh,...
                round(min_sac_dur/Eyedt),round(max_sac_dur/Eyedt),blink_vec);
            cur_sac_starts = [cur_sac_stats(:).start_loc];
            cur_sac_stops = [cur_sac_stats(:).stop_loc];
            cur_sac_amps = [cur_sac_stats(:).peak_amp];
            
            all_blink_startinds{ee} = [all_blink_startinds{ee} cur_blink_starts+length(all_eye_ts{ee})];
            all_blink_stopinds{ee} = [all_blink_stopinds{ee} cur_blink_stops+length(all_eye_ts{ee})];
            all_sac_startinds{ee} = [all_sac_startinds{ee} cur_sac_starts+length(all_eye_ts{ee})];
            all_sac_stopinds{ee} = [all_sac_stopinds{ee} cur_sac_stops+length(all_eye_ts{ee})];
            all_sac_peakamps{ee} = [all_sac_peakamps{ee} cur_sac_amps];
            
            all_eyetrial_start_inds{ee} = [all_eyetrial_start_inds{ee} cur_sp+length(all_eye_ts{ee})];
            all_eye_ts{ee} = [all_eye_ts{ee}; cur_eyets'];
            all_reye_vals{ee} = [all_reye_vals{ee}; cur_reye_vals];
            all_leye_vals{ee} = [all_leye_vals{ee}; cur_leye_vals];
            all_eye_speed{ee} = [all_eye_speed{ee}; eye_speed(:)];
            all_deye_speed{ee} = [all_deye_speed{ee}; diff_eye_speed(:)];
            all_dfeye_speed{ee} = [all_dfeye_speed{ee}; filt_eye_speed(:)];
        end
    end
end

for ee = 1:length(bar_expts)
    all_eye_vals_raw{ee} = [all_leye_vals{ee} all_reye_vals{ee}];
    all_eye_vals_cal{ee} = CalibrateEyeSignal(all_eye_vals_raw{ee}, EyePara.gain, EyePara.offset);
end

%%
min_isi = 0.05;

load(sprintf('./lemM%dExpts.mat',ExptNum));
all_delta_para = [];
all_delta_pos = [];
all_rel_times = [];
all_trial_num = [];
all_expt_num = [];
min_para_amp = 1;
first_range = [0.85 1.18]-0.1;
second_range = [1.55 1.75]-0.1;
for ee = 1:length(bar_expts)
    expt_trial_starts = [Expts{bar_expts(ee)}.Trials(:).TrialStart]/1e4;
    expt_trial_end = [Expts{bar_expts(ee)}.Trials(:).TrueEnd]/1e4;
    
    sac_start_times = all_eye_ts{ee}(all_sac_startinds{ee});
    sac_stop_times = all_eye_ts{ee}(all_sac_stopinds{ee});
    sac_amps = all_sac_peakamps{ee}';
    
    isis = [Inf; diff(sac_start_times)];
    bad_sacs = find(isis < min_isi);
    sac_start_times(bad_sacs) = [];
    sac_stop_times(bad_sacs) = [];
    sac_amps(bad_sacs) = [];
    all_sac_startinds{ee}(bad_sacs) = [];
    all_sac_stopinds{ee}(bad_sacs) = [];
    all_sac_peakamps{ee}(bad_sacs) = [];
    
    buff = round(0/Eyedt);
    avg_eye_xpos = 0.5*all_eye_vals_cal{ee}(:,1) + 0.5*all_eye_vals_cal{ee}(:,3);
    avg_eye_ypos = 0.5*all_eye_vals_cal{ee}(:,2) + 0.5*all_eye_vals_cal{ee}(:,4);
    sac_start_pos = [avg_eye_xpos(all_sac_startinds{ee}+buff) avg_eye_ypos(all_sac_startinds{ee}+buff)];
    sac_stop_pos = [avg_eye_xpos(all_sac_stopinds{ee}+buff) avg_eye_ypos(all_sac_stopinds{ee}+buff)];
    sac_delta_pos = sac_stop_pos - sac_start_pos;
    sac_delta_para = sac_delta_pos(:,1)*cos(Expts{bar_expts(1)}.Stimvals.or*pi/180) + sac_delta_pos(:,2)*sin( Expts{bar_expts(1)}.Stimvals.or*pi/180);
    
    intrial_sacs = zeros(size(sac_start_times));
    trial_num = nan(size(sac_start_times));
    rel_sac_times = zeros(size(sac_start_times));
    for tt = 1:length(expt_trial_starts)
        cur_sac_set = find(sac_start_times > expt_trial_starts(tt) & sac_start_times < expt_trial_end(tt));
        intrial_sacs(cur_sac_set) = 1;
        trial_num(cur_sac_set) = tt;
        rel_sac_times(cur_sac_set) = sac_start_times(cur_sac_set) - expt_trial_starts(tt);
    end
    all_delta_para = [all_delta_para; sac_delta_para(:)];
    all_delta_pos = [all_delta_pos; sac_delta_pos;];
    all_rel_times = [all_rel_times; rel_sac_times(:)];
    all_trial_num = [all_trial_num; trial_num(:)];
    all_expt_num = [all_expt_num; ones(length(trial_num),1)*ee];
end
all_isfirst = find(abs(all_delta_para) > min_para_amp & all_rel_times >= first_range(1) & all_rel_times < first_range(2));
all_issecond = find(abs(all_delta_para) > min_para_amp & all_rel_times >= second_range(1) & all_rel_times < second_range(2));

cnt = 1;
for ee = 1:length(bar_expts)
    expt_trial_starts = [Expts{bar_expts(ee)}.Trials(:).TrialStart]/1e4;
    for tt = 1:length(expt_trial_starts)
        cur_first(cnt) = sum(all_trial_num(all_isfirst)==tt & all_expt_num(all_isfirst)==ee);
        cur_second(cnt) = sum(all_trial_num(all_issecond)==tt & all_expt_num(all_issecond)==ee);
        cnt = cnt + 1;
    end
end
%%
cd(sprintf('~/Data/bruce/M%d',ExptNum));
save random_bar_eyedata_ftime all_* bar_expts
