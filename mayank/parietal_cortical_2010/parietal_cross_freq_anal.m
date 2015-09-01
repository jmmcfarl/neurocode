
clear all
close all
cd G:\WC_Germany\parietal_cortical_2010\
load parietal_cortical_2010
addpath('G:\WC_Germany\parietal_cortical_2010\')
addpath('G:\WC_Germany\persistent_2010\')
addpath('G:\WC_Germany\persistent_revised')
addpath('G:\Code\WC_anal\general\')
addpath('G:\Code\Information Breakdown ToolBox - v.1.0.5')
addpath('G:\WT_KO_analysis_4_29_09\positional_info')

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

%get rid of interneurons
interneurons = find_struct_field_vals(sess_data,'cell_type','interneuron');
sess_data(interneurons) = [];


for d = 1:length(sess_data)

    cdir = sess_data(d).directory;
    cdir(1) = 'G';
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
    backlag = round(Fsd*3);
    forwardlag = round(Fsd*3);
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
    up_dur = nan(n_ups,1);
    rmax_predown_amp = nan(n_ups,1);
    rmax_predown_freq = nan(n_ups,1);
    rup_dur = nan(n_ups,1);
    for i = 1:n_ups

        cur_up_10 = ut_10(i);
        prev_down_90 = dt_90(find(dt_90 < cur_up_10,1,'last'));
        if isempty(prev_down_90)
            prev_down_90 = 1;
        end
        next_down_90 = dt_90(find(dt_90 > cur_up_10,1,'first'));

        %compute amp spectrum of previous down state
        prev_down_inds = prev_down_90:cur_up_10;
        prev_down_dur = length(prev_down_inds)/Fsd;
        if prev_down_dur > 0.1
            L = length(prev_down_inds);
            NFFT = 2^nextpow2(L);
            Y = fft(lf8_f2(prev_down_inds),NFFT)/L;
            Y = 2*abs(Y(1:NFFT/2+1));
            f = Fsd/2*linspace(0,1,NFFT/2+1);
            [max_predown_amp(i),pkloc] = max(Y);
            max_predown_freq(i) = f(pkloc);

%             %compute amp spectrum of current up state
%             cur_up_inds = cur_up_10:next_down_90;
%             up_dur(i) = length(cur_up_inds)/Fsd;
%             L = length(cur_up_inds);
%             NFFT = 2^nextpow2(L);
%             Y = fft(lf8_f2(cur_up_inds),NFFT)/L;
%             Y = 2*abs(Y(1:NFFT/2+1));
%             f = Fsd/2*linspace(0,1,NFFT/2+1);
%             [max_nextup_amp(i),pkloc] = max(Y);
%             max_nextup_freq(i) = f(pkloc);

            used_inds = (cur_up_10 - backlag):(cur_up_10+backlag);
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
        else
            up_dur(i) = nan;
            max_predown_freq(i) = nan;
            max_predown_amp(i) = nan;
%             max_nextup_amp(i) = nan;
%             max_nextup_freq(i) = nan;
        end
    end

    for i = 1:n_ups

        cur_up_10 = rut_10(i);
        prev_down_90 = rdt_90(find(rdt_90 < cur_up_10,1,'last'));
        if isempty(prev_down_90)
            prev_down_90 = 1;
        end
        next_down_90 = rdt_90(find(rdt_90 > cur_up_10,1,'first'));

        %compute amp spectrum of previous down state
        prev_down_inds = prev_down_90:cur_up_10;
        prev_down_dur = length(prev_down_inds)/Fsd;
        if prev_down_dur > 0.1
            L = length(prev_down_inds);
            NFFT = 2^nextpow2(L);
            Y = fft(lf8_f2(prev_down_inds),NFFT)/L;
            Y = 2*abs(Y(1:NFFT/2+1));
            f = Fsd/2*linspace(0,1,NFFT/2+1);
            [rmax_predown_amp(i),pkloc] = max(Y);
            rmax_predown_freq(i) = f(pkloc);

            used_inds = (cur_up_10 - backlag):(cur_up_10+backlag);
            out_of_range = find(used_inds > length(lf8_f) | used_inds < 1);
            map_to = 1:length(lags);
            map_to(out_of_range) = [];
            used_inds(out_of_range) = [];

            rphase_mat(i,map_to) = lf8_f3p(used_inds);
            get_rid = backlag - length(prev_down_inds);
            if get_rid > 0
                rphase_mat(i,1:get_rid) = nan;
            end
        else
            rup_dur(i) = nan;
            rmax_predown_freq(i) = nan;
            rmax_predown_amp(i) = nan;
        end
    end


    % figure
    % [dummy,up_dur_order] = sort(up_dur);
    % pcolor(lags/Fsd,1:n_ups,phase_mat(up_dur_order,:));shading flat
    % line([0 0],[1 410],'color','w')
    % hold on
    % plot(up_dur(up_dur_order),1:n_ups,'color','w')
    % shg
    %
    % figure
    % [dummy,up_dur_order] = sort(up_dur);
    % pcolor(lags/Fsd,1:n_ups,amp_mat(up_dur_order,:));shading flat
    % line([0 0],[1 410],'color','w')
    % hold on
    % plot(up_dur(up_dur_order),1:n_ups,'color','w')
    % shg

    % figure
    up_point = ceil(length(lags)/2);
    pre_down_phase = phase_mat(:,1:up_point);
    rpre_down_phase = rphase_mat(:,1:up_point);
    [dummy,predown_freorder] = sort(max_predown_freq);
    subplot(2,1,1)
    pcolor(lags(1:up_point)/Fsd,1:n_ups,pre_down_phase(predown_freorder,:));shading flat
    subplot(2,1,2)
    plot(lags(1:up_point)/Fsd,nanmean(pre_down_phase)), hold on
    plot(lags(1:up_point)/Fsd,nanmean(rpre_down_phase),'r')
    cname = strcat(sess_data(d).region,'_',sess_data(d).layer,'_',sess_data(d).name);
    t_names = ['G:\WC_Germany\parietal_cortical_2010\cross_freq\predown_' cname];
    print(t_names,'-dpng')
    close

    % figure
    % up_point = ceil(length(lags)/2);
    % pre_down_amp = amp_mat(:,1:up_point);
    % [dummy,predown_freorder] = sort(max_predown_freq);
    % pcolor(lags(1:up_point)/Fsd,1:n_ups,pre_down_amp(predown_freorder,:));shading flat


    %% rescale time
    max_points = 3000;
    max_cycles = 10;
    res_ax = linspace(0,max_cycles,max_points);
    used_axis = 1:max_points;
    pre_down_phase = fliplr(pre_down_phase);
    sc_pre_down_phase = nan(n_ups,max_points);
    rpre_down_phase = fliplr(rpre_down_phase);
    rsc_pre_down_phase = nan(n_ups,max_points);
    num_cycles = nan(n_ups,1);
    rnum_cycles = nan(n_ups,1);
    for i = 1:n_ups
        if max_predown_freq(i) > 1
            cur_data = pre_down_phase(i,:);
            cur_data(isnan(cur_data)) = [];
            period = round(Fsd/max_predown_freq(i));
            num_cycles(i) = length(cur_data)/period;
            needed_length = max_cycles*period;
            cur_data = [cur_data nan(1,needed_length-length(cur_data))];
            res_data = interp1(linspace(0,max_cycles,length(cur_data)),cur_data,res_ax);
            sc_pre_down_phase(i,:) = res_data;

            %         cur_data = pre_up_amp(i,:);
            %         cur_data(isnan(cur_data)) = [];
            %         period = round(Fsd/max_preup_freq(i));
            %         needed_length = max_cycles*period;
            %         cur_data = [cur_data nan(1,needed_length-length(cur_data))];
            %         res_data = interp1(linspace(0,max_cycles,length(cur_data)),cur_data,res_ax);
            %         sc_pre_up_amp(i,:) = res_data;
            %
        else
            num_cycles(i) = 0;
        end
        if rmax_predown_freq(i) > 1
            cur_data = rpre_down_phase(i,:);
            cur_data(isnan(cur_data)) = [];
            period = round(Fsd/rmax_predown_freq(i));
            rnum_cycles(i) = length(cur_data)/period;
            needed_length = max_cycles*period;
            cur_data = [cur_data nan(1,needed_length-length(cur_data))];
            res_data = interp1(linspace(0,max_cycles,length(cur_data)),cur_data,res_ax);
            rsc_pre_down_phase(i,:) = res_data;

        else
            rnum_cycles(i) = 0;
        end
    end

    figure
    [dummy,cycle_sort] = sort(num_cycles);
    subplot(2,1,1)
    pcolor(res_ax,1:n_ups,sc_pre_down_phase(cycle_sort,:));shading flat
    subplot(2,1,2)
    plot(res_ax,nanmean(sc_pre_down_phase)), hold on
    plot(res_ax,nanmean(rsc_pre_down_phase),'r')
    cname = strcat(sess_data(d).region,'_',sess_data(d).layer,'_',sess_data(d).name);
    t_names = ['G:\WC_Germany\parietal_cortical_2010\cross_freq\scale_predown_' cname];
    print(t_names,'-dpng')
    close

    % figure
    % [dummy,cycle_sort] = sort(num_cycles);
    % pcolor(res_ax,1:n_ups,sc_pre_down_amp(cycle_sort,:));shading flat
    %
    %%
    % min_n = 10;
    % for i = 1:max_points
    %     cur_set = sc_pre_down_phase(:,i);
    %     cur_set(isnan(cur_set)) = [];
    %     if length(cur_set) > min_n
    %         nt = length(cur_set);
    %         nb = round(nt/3);
    %         bin_edges = linspace(-2*pi,2*pi,nb);
    %         n_raw = histc(cur_set,bin_edges);
    %         p_dir = n_raw/sum(n_raw);
    %         ent_dir = -nansum(p_dir.*log(p_dir));
    %         max_ent = log(nt);
    %         cur_ent = estimate_pt_entropy(n_raw);
    %         phase_ent(i) = (max_ent-cur_ent)/max_ent;
    %     end
    % end
    % figure
    % plot(phase_ent)
    %
    % %%
    min_n = 10;
    n_pts = length(res_ax);
    phase_ent = nan(n_pts,1);
    rphase_ent = nan(n_pts,1);
    for i = 1:n_pts
        cur_set = sc_pre_down_phase(:,i);
        cur_set(isnan(cur_set)) = [];
        if length(cur_set) > min_n
            nt = length(cur_set);
            nb = round(nt/3);
            bin_edges = linspace(-2*pi,2*pi,nb);
            n_raw = histc(cur_set,bin_edges);
            p_dir = n_raw/sum(n_raw);
            ent_dir = -nansum(p_dir.*log(p_dir));
            max_ent = log(nt);
            cur_ent = estimate_pt_entropy(n_raw);
            phase_ent(i) = (max_ent-cur_ent)/max_ent;
        end
        cur_set = rsc_pre_down_phase(:,i);
        cur_set(isnan(cur_set)) = [];
        if length(cur_set) > min_n
            nt = length(cur_set);
            nb = round(nt/3);
            bin_edges = linspace(-2*pi,2*pi,nb);
            n_raw = histc(cur_set,bin_edges);
            p_dir = n_raw/sum(n_raw);
            ent_dir = -nansum(p_dir.*log(p_dir));
            max_ent = log(nt);
            cur_ent = estimate_pt_entropy(n_raw);
            rphase_ent(i) = (max_ent-cur_ent)/max_ent;
        end
    end
    
    
    figure
    plot(res_ax,phase_ent), hold on
    plot(res_ax,rphase_ent,'r')
    cname = strcat(sess_data(d).region,'_',sess_data(d).layer,'_',sess_data(d).name);
    t_names = ['G:\WC_Germany\parietal_cortical_2010\cross_freq\scale_predown_entropy_' cname];
    print(t_names,'-dpng')
    close

    %
    % %%
    % figure
    % plot(nanmean(sc_pre_down_phase),'r')

end