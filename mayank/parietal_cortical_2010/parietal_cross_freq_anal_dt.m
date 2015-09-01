
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


    %% compute the down-triggered distribution of phase and amplitude of the theta-band signal
    backlag = round(Fsd*3);
    forwardlag = round(Fsd*3);
    or_states = find(ut_10 < backlag | dt_90 < backlag | ut_10 > length(lf8_f) - forwardlag | dt_90 > length(lf8_f) - forwardlag);
    ut_10(or_states) = [];
    dt_90(or_states) = [];

    n_downs = length(dt_90);
    lags = -backlag:forwardlag;
    phase_mat = nan(n_downs,length(lags));
    n_points = 3000;
    n_point = 1000;
    max_nextdown_amp = nan(n_downs,1);
    max_nextdown_freq = nan(n_downs,1);
    up_dur = nan(n_downs,1);
    for i = 1:n_downs

        cur_dt_90 = dt_90(i);
        prev_up_10 = ut_10(find(ut_10 < cur_dt_90,1,'last'));
        if isempty(prev_up_10)
            prev_up_10 = 1;
        end
        next_up_10 = ut_10(find(ut_10 > cur_dt_90,1,'first'));

        %compute amp spectrum of following down state
        next_down_inds = cur_dt_90:next_up_10;
        next_down_dur = length(next_down_inds)/Fsd;
        if next_down_dur > 0.2
            L = length(next_down_inds);
            NFFT = 2^nextpow2(L);
            Y = fft(lf8_f2(next_down_inds),NFFT)/L;
            Y = 2*abs(Y(1:NFFT/2+1));
            f = Fsd/2*linspace(0,1,NFFT/2+1);
            [max_nextdown_amp(i),pkloc] = max(Y);
            max_nextdown_freq(i) = f(pkloc);

            extra = length(next_down_inds) - length(lags);
            if extra > 0
                next_down_inds(end-extra+1:end) = [];
            end
            map_to = 1:length(next_down_inds);
            phase_mat(i,map_to) = lf8_f3p(next_down_inds);
            down_dur(i) = length(next_down_inds)/Fsd;
        end
    end
    
    %% rescale time
    max_points = 3000;
    max_cycles = 10;
    Fsn = max_points/max_cycles;
    up_point = ceil(length(lags)/2);
    res_ax = linspace(0,max_cycles,max_points);
    fn = Fsn/2*linspace(0,1,L/2+1);
    Ln = length(res_ax);
    used_axis = 1:max_points;
    next_down_phase = phase_mat;
    num_cycles = nan(n_downs,1);
    sc_next_down_phase = nan(n_downs,length(res_ax));
    for i = 1:n_downs
        if max_nextdown_freq(i) > 1
            cur_data = next_down_phase(i,:);
            cur_data(isnan(cur_data)) = [];
            ndata = length(cur_data);
            period = round(Fsd/max_nextdown_freq(i));
            num_cycles(i) = length(cur_data)/period;
            needed_length = min(round(num_cycles(i)/max_cycles*max_points),max_points);
            res_data = interp1(linspace(0,num_cycles(i),ndata),cur_data,linspace(0,num_cycles(i),needed_length));
            sc_next_down_phase(i,1:length(res_data)) = res_data;
        else
            num_cycles(i) = 0;
        end
    end
avg_sc_next_down_phase = nanmean(sc_next_down_phase);
Y = fft(avg_sc_next_down_phase)/Ln;
Y = 2*abs(Y(1:Ln/2+1));
avg_spectrum = Y;
fn = Fsn/2*linspace(0,1,Ln/2+1);

min_n = 10;
n_pts = length(res_ax);
phase_ent = nan(n_pts,1);
for i = 1:n_pts
    cur_set = sc_next_down_phase(:,i);
    cur_set(isnan(cur_set)) = [];
    if length(cur_set) > min_n
        nt = length(cur_set);
        nb = round(nt/3);
        bin_edges = linspace(-2*pi,2*pi,nb);
        n_raw = histc(cur_set,bin_edges);
        max_ent = log(nt);
        cur_ent = estimate_pt_entropy(n_raw);
        phase_ent(i) = (max_ent-cur_ent)/max_ent;
    end
end

%% bootstrap resampling for CIs
    nreps = 100;
    rmax_predown_amp = nan(nreps,n_downs);
    rmax_predown_freq = nan(nreps,n_downs);
    rup_dur = nan(nreps,n_downs);
    rmean_phase = nan(nreps,length(lags));
    
    rmean_spectrum = nan(nreps,length(avg_spectrum));  
    rnum_cycles = nan(nreps,n_downs);
    rsc_mean_phase = nan(nreps,length(res_ax));
    rphase_ent = nan(nreps,n_pts);
    for r = 1:nreps
        % create randomized transition times
        jit_time = round(rand(size(ut_10(:)))*2*Fsd)-round(2*Fsd)/2;
        rut_10 = sort(ut_10(:) + jit_time);
        rdt_90 = sort(dt_90(:) + jit_time);
        rut_10(rut_10 > rdt_90(end)) = [];
        ror_states = find(rut_10 < backlag | rdt_90 < backlag);
        rut_10(ror_states) = [];
        rdt_90(ror_states) = [];
        if length(rut_10) > n_downs
            rut_10(n_ups+1:end) = [];
            rdt_90(n_ups+1:end) = [];
        end
        rphase_mat = nan(n_downs,length(lags));
        rn_ups = length(rut_10);
        for i = 1:rn_ups
            
            cur_up_10 = rdt_90(i);
            prev_down_90 = rut_10(find(rut_10 < cur_up_10,1,'last'));
            if isempty(prev_down_90)
                prev_down_90 = 1;
            end
            next_down_90 = rut_10(find(rut_10 > cur_up_10,1,'first'));
            
            %compute amp spectrum of previous down state
            prev_down_inds = prev_down_90:cur_up_10;
            prev_down_dur = length(prev_down_inds)/Fsd;
            if prev_down_dur > 0.2
                L = length(prev_down_inds);
                NFFT = 2^nextpow2(L);
                Y = fft(lf8_f2(prev_down_inds),NFFT)/L;
                Y = 2*abs(Y(1:NFFT/2+1));
                f = Fsd/2*linspace(0,1,NFFT/2+1);
                [rmax_predown_amp(r,i),pkloc] = max(Y);
                rmax_predown_freq(r,i) = f(pkloc);
                
                extra = length(prev_down_inds) - length(lags);
                if extra > 0
                    prev_down_inds(end-extra+1:end) = [];
                end
                map_to = 1:length(prev_down_inds);
                rphase_mat(i,map_to) = lf8_f3p(prev_down_inds);
            end
        end
        rmean_phase(r,:) = nanmean(rphase_mat);
        rpre_down_phase = rphase_mat;
        rsc_pre_down_phase = nan(rn_ups,max_points);
        for i = 1:rn_ups           
            if rmax_predown_freq(r,i) > 1
                cur_data = rpre_down_phase(i,:);
                cur_data(isnan(cur_data)) = [];
                ndata = length(cur_data);
                period = round(Fsd/rmax_predown_freq(r,i));
                rnum_cycles(r,i) = ndata/period;
                needed_length = min(round(rnum_cycles(r,i)/max_cycles*max_points),max_points);
                res_data = interp1(linspace(0,rnum_cycles(r,i),ndata),cur_data,linspace(0,rnum_cycles(r,i),needed_length));
                rsc_pre_down_phase(i,1:needed_length) = res_data;     
            else
                rnum_cycles(r,i) = 0;
            end          
        end
        rsc_mean_phase(r,:) = nanmean(rsc_pre_down_phase);
        Y = fft(rsc_mean_phase(r,:))/Ln;
        Y = 2*abs(Y(1:Ln/2+1));
        rmean_spectrum(r,:) = Y;
        
        for i = 1:n_pts
            cur_set = rsc_pre_down_phase(:,i);
            cur_set(isnan(cur_set)) = [];
            if length(cur_set) > min_n
                nt = length(cur_set);
                nb = round(nt/3);
                bin_edges = linspace(-2*pi,2*pi,nb);
                n_raw = histc(cur_set,bin_edges);
                max_ent = log(nt);
                cur_ent = estimate_pt_entropy(n_raw);
                rphase_ent(r,i) = (max_ent-cur_ent)/max_ent;
            end
        end

    end

    uci_rsc_phase = prctile(rsc_mean_phase,97.5);
    lci_rsc_phase = prctile(rsc_mean_phase,2.5);
    uci_rsc_spec = prctile(rmean_spectrum,95);
    uci_phase_ent = prctile(rphase_ent,95);
    
    %%
    f1 = find(fn > 0.999,1,'first');
    r_coh_ind_95(d) = uci_rsc_spec(f1);
    coh_ind(d) = avg_spectrum(f1);
    
    %%
    [dummy,cycle_sort] = sort(num_cycles);
    figure;
    set(gcf,'PaperUnits','inches')
    set(gcf,'PaperSize',[10 15])
    set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
    subplot(3,1,1)
    pcolor(res_ax,1:n_downs,sc_next_down_phase(cycle_sort,:));shading flat
    subplot(3,1,2)
    plot(res_ax,uci_rsc_phase,'r'), hold on
    plot(res_ax,lci_rsc_phase,'r')
    plot(res_ax,avg_sc_next_down_phase)
    title(sprintf('coh: %.3f  rcoh: %.3f',coh_ind(d),r_coh_ind_95(d)));
   cname = strcat(sess_data(d).region,'_',sess_data(d).layer,'_',sess_data(d).name);
   subplot(3,1,3)
   plot(res_ax,phase_ent), hold on
   plot(res_ax,uci_phase_ent,'r')
  t_names = ['F:\WC_Germany\parietal_cortical_2010\cross_freq\scale_nextdown_spec' cname];
    print(t_names,'-dpng')
    close


end