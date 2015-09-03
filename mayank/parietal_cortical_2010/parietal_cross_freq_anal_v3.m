
clear all
close all
cd F:\WC_Germany\parietal_cortical_2010\
load parietal_cortical_2010
addpath('F:\WC_Germany\parietal_cortical_2010\')
addpath('F:\WC_Germany\persistent_2010\')
addpath('F:\WC_Germany\persistent_revised')
addpath('F:\Code\WC_anal\general\')
addpath('F:\Code\Information Breakdown ToolBox - v.1.0.5')
addpath('F:\WT_KO_analysis_4_29_09\positional_info')

raw_Fs = 2016;
niqf = raw_Fs/2;
dsf = 8;
Fsd = raw_Fs/dsf;

lcf = 0.05;
hcf = 20;

lcf2 = 1.5;
hcf2 = 10;

lcf3 = 2;
hcf3 = 5;

min_rsquared = 0.75;
min_dur = 1.5;

%get rid of interneurons
interneurons = find_struct_field_vals(sess_data,'cell_type','interneuron');
sess_data(interneurons) = [];


for d = 1:length(sess_data)

    cdir = sess_data(d).directory;
    cdir(1) = 'F';
    disp(sprintf('session %d',d))
    cd(cdir);

    load used_data lf8 wcv_minus_spike
    load hsmm_state_seq_seg_lf_4_10_10
    mp_state_seq = hsmm_bbstate_seq;
    load hsmm_state_seq8_seg_lf_4_10_10
    lfp_state_seq = hsmm_bbstate_seq8;

    used_state_seq = lfp_state_seq;

    lf8_f = get_lf_features(lf8,raw_Fs,Fsd,[lcf hcf]);
    lf8_f2 = get_lf_features(lf8,raw_Fs,Fsd,[lcf2 hcf2]);
    lf8_f3 = get_lf_features(lf8,raw_Fs,Fsd,[lcf3 hcf3]);
    t = (1:length(lf8_f))/Fsd;
    lf8_f3p = angle(hilbert(lf8_f3));

    %extract the index values of LF8 up transitions
    up_inds = [];
    down_inds = [];
    UDS_segs = (hmm8.UDS_segs-1)*2+1;
    for i = 1:length(lfp_state_seq)
        up_trans = UDS_segs(i,1)+find(used_state_seq{i}(1:end-1) == 1 & used_state_seq{i}(2:end) == 2);
        down_trans = UDS_segs(i,1)+find(used_state_seq{i}(1:end-1) == 2 & used_state_seq{i}(2:end) == 1);
        down_trans(down_trans < up_trans(1)) = [];
        up_trans(up_trans > down_trans(end)) = [];
        up_inds = [up_inds; up_trans];
        down_inds = [down_inds; down_trans];
    end
    %% Fit sigmoids to state transitions
    [urlid,urltime,urlamp,urlshift,urltau,urlerror,ut_90,ut_10] = ...
        get_lfp_wcv_sigmoid_fit_ut_12_4(up_inds,down_inds,lf8_f,t,Fsd);
    [drlid,drltime,drlamp,drlshift,drltau,drlerror,dt_90,dt_10] = ...
        get_lfp_wcv_sigmoid_fit_dt_12_4(down_inds,up_inds,-lf8_f,t,Fsd);


    % ignore state transitions where the sigmoid fit has too low of an R^2
    % value
    bad_lups = find(urlerror < min_rsquared);
    bad_ldowns = find(drlerror < min_rsquared);
    bad_states = union(bad_lups,bad_ldowns);
    up_inds(bad_states) = [];
    down_inds(bad_states) = [];
    ut_90(bad_states) = [];
    ut_10(bad_states) = [];
    dt_90(bad_states) = [];
    dt_10(bad_states) = [];
    ut_10 = ut_10(:);
    dt_90 = dt_90(:);

    % create randomized transition times
    jit_time = round(rand(size(ut_10(:)))*2*Fsd)-round(2*Fsd)/2;
    rut_10 = sort(ut_10(:) + jit_time);
    rdt_90 = sort(dt_90(:) + jit_time);
    % ut_10 = sort(round(rand(size(ut_10))*max(ut_10)));
    % dt_90 = sort(round(rand(size(dt_90))*max(dt_90)));
    % clear ut_90 dt_10
    rut_10(rut_10 > rdt_90(end)) = [];

    %% compute the up-triggered distribution of phase and amplitude of the theta-band signal
    % min_dur = 1.3;
    % min_pow = 0.5;
    backlag = round(Fsd*4);
    forwardlag = round(Fsd*2);
    or_states = find(ut_10 < backlag | dt_90 < backlag | ut_10 > length(lf8_f) - forwardlag | dt_90 > length(lf8_f) - forwardlag);
    ut_10(or_states) = [];
    dt_90(or_states) = [];
    ror_states = find(rut_10 < backlag | rdt_90 < backlag | rut_10 > length(lf8_f) - forwardlag | rdt_90 > length(lf8_f) - forwardlag);
    rut_10(or_states) = [];
    rdt_90(or_states) = [];

    n_ups = length(ut_10);

    lags = -backlag:forwardlag;
    phase_mat = nan(n_ups,length(lags));
    rphase_mat = nan(n_ups,length(lags));
    % amp_mat = nan(n_ups,length(lags));
    n_points = 3000;
    n_point = 1000;
    max_predown_amp = nan(n_ups,1);
    max_predown_freq = nan(n_ups,1);
    down_dur = nan(n_ups,1);
    for i = 1:n_ups

        cur_up_10 = ut_10(i);
        prev_down_90 = dt_90(find(dt_90 < cur_up_10,1,'last'));
        if isempty(prev_down_90)
            prev_down_90 = 1;
        end
        next_down_90 = dt_90(find(dt_90 > cur_up_10,1,'first'));

        %compute amp spectrum of previous down state
        prev_down_inds = prev_down_90:cur_up_10;
        next_up_length = next_down_90-cur_up_10;
        down_dur(i) = length(prev_down_inds)/Fsd;
        if down_dur(i) > min_dur
            L = length(prev_down_inds);
            NFFT = 2^nextpow2(L);
            Y = fft(lf8_f2(prev_down_inds),NFFT)/L;
            Y = 2*abs(Y(1:NFFT/2+1));
            f = Fsd/2*linspace(0,1,NFFT/2+1);
            [max_predown_amp(i),pkloc] = max(Y);
            max_predown_freq(i) = f(pkloc);

            used_inds = (cur_up_10 - backlag):(cur_up_10+forwardlag);
            out_of_range = find(used_inds > length(lf8_f) | used_inds < 1);
            map_to = 1:length(lags);
            map_to(out_of_range) = [];
            used_inds(out_of_range) = [];

            phase_mat(i,map_to) = lf8_f3p(used_inds);
%             amp_mat(i,map_to) = lf8_f3(used_inds);
            get_rid = backlag - length(prev_down_inds);
            if get_rid > 0
                phase_mat(i,1:get_rid) = nan;
%                 amp_mat(i,1:get_rid) = nan;
            end
            get_rid = forwardlag - next_up_length;
            if get_rid > 0
               phase_mat(i,end-get_rid+1:end) = nan; 
            end
        else
            down_dur(i) = nan;
            max_predown_freq(i) = nan;
            max_predown_amp(i) = nan;
%             max_nextup_amp(i) = nan;
%             max_nextup_freq(i) = nan;
        end
    end


    %% rescale time
    
    max_points = length(lags)*2;
    back_cycles = 7;
    forward_cycles = 5;
    ov_sc = max_points/length(lags);
    sc_backlag = ov_sc*backlag;
    sc_forwardlag = ov_sc*forwardlag;
    res_ax = linspace(-back_cycles,forward_cycles,max_points);
    used_axis = 1:max_points;
    sc_pre_down_phase = nan(n_ups,max_points);
    num_cycles = nan(n_ups,1);
    for i = 1:n_ups
        if max_predown_freq(i) > 1 && down_dur(i) > min_dur
            cur_data = phase_mat(i,:);
            period = round(Fsd/max_predown_freq(i));
            
            %for previous down state
            prev_down_data = cur_data(1:backlag);
            prev_down_data(isnan(prev_down_data)) = [];
            pd_cycles = length(prev_down_data)/period;
            needed_length = round(pd_cycles/back_cycles*sc_backlag);
            res_pd_data = interp1(linspace(0,pd_cycles,length(prev_down_data)),...
                prev_down_data,linspace(0,pd_cycles,needed_length));
            extra = length(res_pd_data)-sc_backlag;
            if extra > 0
                res_pd_data(1:extra) = [];
            end
            sp = sc_backlag - length(res_pd_data)+1;
            sc_pre_down_phase(i,sp:sc_backlag) = res_pd_data;
            
            %for next up state
            next_up_data = cur_data(backlag+1:end);
            next_up_data(isnan(next_up_data)) = [];
            ud_cycles = length(next_up_data)/period;
            needed_length = round(ud_cycles/forward_cycles*sc_forwardlag);
            res_ud_data = interp1(linspace(0,ud_cycles,length(next_up_data)),...
                next_up_data,linspace(0,ud_cycles,needed_length));
            extra = length(res_ud_data) - sc_forwardlag;
            if extra > 0
                res_ud_data(end-extra+1:end) = [];
            end
            ep = length(res_ud_data)+sc_backlag;
            sc_pre_down_phase(i,sc_backlag+1:ep) = res_ud_data;
            
            num_cycles(i) = pd_cycles;
        else
            num_cycles(i) = 0;
        end
    end

    good_cycles = find(num_cycles > 0);
   
    %%
    figure
    [dummy,cycle_sort] = sort(num_cycles(good_cycles));
    subplot(2,1,1)
    pcolor(res_ax,1:length(good_cycles),sc_pre_down_phase(good_cycles(cycle_sort),:));shading flat
    subplot(2,1,2)
    plot(res_ax,nanmean(sc_pre_down_phase)), hold on
    xlim([-back_cycles forward_cycles])
    ylim([-1 1])
    cname = strcat(sess_data(d).region,'_',sess_data(d).layer,'_',sess_data(d).name);
    t_names = ['F:\WC_Germany\parietal_cortical_2010\cross_freq\scale_downup_' cname];
    print(t_names,'-dpng')
    close


end