
clear all
data_dir = cell(0);
data_type = cell(0);
data_exp = cell(0);
data_ep = [];
data_dp = [];
data_hpc_lfp = [];

cd C:\WC_Germany\persistent_downs\
load ./overall_EC_dir

addpath('C:\WC_Germany\parietal_cortical_2010\');
addpath('C:\WC_Germany\persistent_9_27_2010\');
addpath('C:\WC_Germany\new_mec\');
addpath('C:\WC_Germany\overall_EC');
addpath('C:\WC_Germany\hsmm_state_detection');
addpath('C:\WC_Germany\sven_thomas_combined');
addpath('C:\Code\general_functions\');
%%

raw_Fs = 2016;
dsf = 8;
Fsd = raw_Fs/dsf;
niqf = raw_Fs/2;
hcf = 10;
lcf = 0.05;
lcf_hf = 15;
hcf_hf = 80;
hcf_sm = 0.025;

maxlag = round(4*Fsd);
backlag = 4*Fsd;
forwardlag = 4*Fsd;

params.Fs = raw_Fs;
params.fpass = [0 100];
params.tapers = [3 5];
win = 25;

%robust persistence requires MP skips states longer than this
thresh_lf8_updur = 0.5;
thresh_lf8_downdur = 0.5;

%for state duration distribution calculations
up_range = [0.1 15];
down_range = [0.1 15];
numBins = 60;
log_dur_grid = logspace(log10(up_range(1)),log10(up_range(2)),numBins+1);

%%
min_rec_dur = 500; %in sec
used_dirs = find(data_ep > min_rec_dur);
for d = 1:length(sess_data)
    cd(sess_data(d).directory)
    pwd
    
    load ./used_data lf7 wcv_minus_spike
    [lf8_lf,t_axis] = get_lf_features(lf7,raw_Fs,Fsd,[lcf hcf]);
    wcv_lf = get_lf_features(wcv_minus_spike,raw_Fs,Fsd,[lcf hcf]);
    
    %     if data_hpc_lfp(d) == 3
    %         load ./used_data lf3 lf5
    %         if ismember(d,old_data_inds)
    %             lf3 = lf3 + lf5; %redefine LF3 wrt gnd
    %         end
    %         hpc_hf = get_hf_features(lf3,raw_Fs,Fsd,[lcf_hf hcf_hf],hcf_sm);
    %     else
%     load ./used_data lf2
%     hpc_hf = get_hf_features(lf2,raw_Fs,Fsd,[lcf_hf hcf_hf],hcf_sm);
    %     end
    
    if sess_ep(d) ~= 0 & sess_dp(d) ~= 0
        end_time = min(sess_ep(d),sess_dp(d));
        ep = find(t_axis >= end_time,1);
        if ~isempty(ep)
            lf8_lf(ep+1:end) = []; wcv_lf(ep+1:end) = []; t_axis(ep+1:end) = []; 
        else
            ep = length(t_axis);
        end
    else
        ep = length(t_axis);
    end
    rec_dur(d) = range(t_axis);
    
    if exist('./all_combined_mp_uds.mat','file')
        load ./all_combined_mp_uds.mat
        load ./all_combined_lf7_uds.mat
    else
        load ./pa_hsmm_state_seq7_combined_fin_nd.mat
        load ./pa_hsmm_state_seq_combined_fin_nd.mat
    end
    hmm8 = hmm7;
    lfp_state_seq = hmm_bbstate_seq7;
    mp_state_seq = hmm_bbstate_seq;
    
    %         if you want to flip the roles of the MP and LFP
    %         temp = hsmm_bbstate_seq8;
    %         hsmm_bbstate_seq8 = hsmm_bbstate_seq;
    %         hsmm_bbstate_seq = temp;
    
    if isempty(hmm) || isempty(hmm8)
        use_data(d) = false;
    else
        use_data(d) = true;
    end
    
    if exist('./allEC_ctx_period_data.mat','file')
        load ./allEC_ctx_period_data.mat
    elseif exist('./combined_lf7_period_data_fin_nd.mat','file')
        load ./combined_lf7_period_data_fin_nd
    else
        use_data(d) = false;
    end
    
    if use_data(d)
        [new_seg_inds] = resample_uds_seg_inds_v2(hmm.UDS_segs,hmm.Fs,Fsd,mp_state_seq);
        dur_uds = sum(diff(new_seg_inds,[],2))/Fsd;
        [up_trans_inds,down_trans_inds] = compute_state_transitions_seg(new_seg_inds,mp_state_seq);
        [up_trans_inds8,down_trans_inds8] = compute_state_transitions_seg(new_seg_inds,lfp_state_seq);
        
        bad_mp_states = find(up_trans_inds > ep | down_trans_inds > ep);
        up_trans_inds(bad_mp_states) = []; down_trans_inds(bad_mp_states) = [];
        bad_lfp_states = find(up_trans_inds8 > ep | down_trans_inds8 > ep);
        up_trans_inds8(bad_lfp_states) = []; down_trans_inds8(bad_lfp_states) = [];
        
        fract_uds_dur(d) = dur_uds/t_axis(end);
        
        %%
        lf8_period_vec = nan(size(wcv_lf));
        for i = 1:size(new_seg_inds,1)
            cur_inds = new_seg_inds(i,1):new_seg_inds(i,2);
            cur_inds_used = find(cur_inds <= ep);
            lf8_period_vec(cur_inds(cur_inds_used)) = lf8_period_f{i}(cur_inds_used);
        end
        
        %% compute state duration distributions
        [mp_state_durations{d}] = compute_state_durations_seg(mp_state_seq,Fsd);
        [lfp_state_durations{d}] = compute_state_durations_seg(lfp_state_seq,Fsd);
        mp_state_durations{d}{1}(bad_mp_states) = []; mp_state_durations{d}{2}(bad_mp_states) = [];
        lfp_state_durations{d}{1}(bad_lfp_states) = []; lfp_state_durations{d}{2}(bad_lfp_states) = [];
        
        [up_amps,down_amps] = get_state_amplitudes(lf8_lf,up_trans_inds8,down_trans_inds8);
        
        up_dur_loghist(d,:) = histc(mp_state_durations{d}{2},log_dur_grid);
        up_dur_loghist8(d,:) = histc(lfp_state_durations{d}{2},log_dur_grid);
        down_dur_loghist(d,:) = histc(mp_state_durations{d}{1},log_dur_grid);
        down_dur_loghist8(d,:) = histc(lfp_state_durations{d}{1},log_dur_grid);
        
        mp_updurs = mp_state_durations{d}{2};
        mp_downdurs = mp_state_durations{d}{1};
        lf8_updurs= lfp_state_durations{d}{2};
        lf8_downdurs = lfp_state_durations{d}{1};
        
        mp_dutycycs = mp_updurs./(mp_updurs+mp_downdurs);
        lf8_dutycycs = lf8_updurs./(lf8_updurs+lf8_downdurs);
        median_dc(d) = nanmedian(mp_dutycycs);
        median_dc8(d) = nanmedian(lf8_dutycycs);
        mean_dc(d) = nanmean(mp_dutycycs);
        mean_dc8(d) = nanmean(lf8_dutycycs);
        
        mean_updur8(d) = nanmean(lfp_state_durations{d}{2});
        mean_downdur8(d) = nanmean(lfp_state_durations{d}{1});
        median_updur8(d) = nanmedian(lfp_state_durations{d}{2});
        median_downdur8(d) = nanmedian(lfp_state_durations{d}{1});
        mean_updur(d) = nanmean(mp_state_durations{d}{2});
        mean_downdur(d) = nanmean(mp_state_durations{d}{1});
        median_updur(d) = nanmedian(mp_state_durations{d}{2});
        median_downdur(d) = nanmedian(mp_state_durations{d}{1});
        max_updur8(d) = nanmax(lfp_state_durations{d}{2});
        max_updur(d) = nanmax(mp_state_durations{d}{2});
        max_downdur8(d) = nanmax(lfp_state_durations{d}{1});
        max_downdur(d) = nanmax(mp_state_durations{d}{1});
        
        %% compute corresponding state transitions and transition lags
        [corresp_lf8_upinds{d},corresp_lf8_downinds{d}] = greedy_find_corresponding_ncx_state_transitions_simp(...
            up_trans_inds,down_trans_inds,up_trans_inds8,down_trans_inds8);
        
        %find non-skipped states
        non_skipped_mp_up_states{d} = find(~isnan(corresp_lf8_upinds{d}));
        non_skipped_mp_down_states{d} = find(~isnan(corresp_lf8_downinds{d}));
        
        %compute transition lags for non-skipped states
        mp_uplags{d} = nan(size(up_trans_inds));
        mp_downlags{d} = nan(size(up_trans_inds));
        mp_uplags{d}(non_skipped_mp_up_states{d}) = up_trans_inds(non_skipped_mp_up_states{d}) - up_trans_inds8(corresp_lf8_upinds{d}(non_skipped_mp_up_states{d}));
        mp_downlags{d}(non_skipped_mp_down_states{d}) = down_trans_inds(non_skipped_mp_down_states{d}) - down_trans_inds8(corresp_lf8_downinds{d}(non_skipped_mp_down_states{d}));
        
        %compute transition lags relative to ctx state durations
        mp_reldownlags{d} = nan(size(up_trans_inds));
        mp_reluplags{d} = nan(size(up_trans_inds));
        for i = 1:length(up_trans_inds)
            if ~isnan(corresp_lf8_downinds{d}(i))
                mp_reldownlags{d}(i) = mp_downlags{d}(i)/lf8_downdurs(corresp_lf8_downinds{d}(i))/Fsd;
            end
            if ~isnan(corresp_lf8_upinds{d}(i))
                mp_reluplags{d}(i) = mp_uplags{d}(i)/lf8_updurs(corresp_lf8_upinds{d}(i))/Fsd;
            end
        end
        
        %% compute persistence
        [mp_upskipped{d},mp_downskipped{d}] = greedy_find_skipped_ncx_states(...
            corresp_lf8_upinds{d},corresp_lf8_downinds{d},lfp_state_durations{d}{2},lfp_state_durations{d}{1},thresh_lf8_downdur,thresh_lf8_updur);
        
        rt2_ups{d} = find(mp_upskipped{d}.rnum_skipped > 0);
        t2_ups{d} = find(mp_upskipped{d}.num_skipped > 0);
        nrt2_ups{d} = find(mp_upskipped{d}.num_skipped == 0);
        fract_rt2_ups(d) = length(rt2_ups{d})/(length(rt2_ups{d}) + length(nrt2_ups{d}));
        fract_t2_ups(d) = length(t2_ups{d})/(length(t2_ups{d}) + length(nrt2_ups{d}));
        
        rt2_downs{d} = find(mp_downskipped{d}.rnum_skipped > 0);
        nrt2_downs{d} = find(mp_downskipped{d}.num_skipped==0); %number of down states that didn't skip any LFP up states
        fract_rt2_downs(d) = length(rt2_downs{d})/(length(rt2_downs{d}) + length(nrt2_downs{d}));
        
        robust_non_skipped_mp_ups{d} = [rt2_ups{d}; nrt2_ups{d}];
        robust_non_skipped_mp_downs{d} = [rt2_downs{d}; nrt2_downs{d}];
        
        mp_rel_updurs{d} = nan(size(up_trans_inds));
        mp_rel_downdurs{d} = nan(size(up_trans_inds));
        mp_rel_updurs{d}(robust_non_skipped_mp_ups{d}) = mp_state_durations{d}{2}(robust_non_skipped_mp_ups{d}) - ...
            lfp_state_durations{d}{2}(corresp_lf8_upinds{d}(robust_non_skipped_mp_ups{d}));
        mp_rel_downdurs{d}(robust_non_skipped_mp_downs{d}) = mp_state_durations{d}{1}(robust_non_skipped_mp_downs{d}) - ...
            lfp_state_durations{d}{1}(corresp_lf8_downinds{d}(robust_non_skipped_mp_downs{d}));
        mp_rel_updurs_t2{d} = mp_rel_updurs{d}(rt2_ups{d});
        mp_rel_updurs_nt2{d} = mp_rel_updurs{d}(nrt2_ups{d});
        mp_rel_downdurs_t2{d} = mp_rel_downdurs{d}(rt2_downs{d});
        mp_rel_downdurs_nt2{d} = mp_rel_downdurs{d}(nrt2_downs{d});
        
        fract_rt1_ups(d) = nansum(mp_rel_updurs{d} > 0)/sum(~isnan(mp_rel_updurs{d}));
        fract_rt1_ups_nt2(d) = nansum(mp_rel_updurs_nt2{d} > 0)/sum(~isnan(mp_rel_updurs_nt2{d}));
        fract_rt1_downs(d) = nansum(mp_rel_downdurs{d} > 0)/sum(~isnan(mp_rel_downdurs{d}));
        
        %% compute durations in units of Ncx UDS cycles
        [mp_updurs_lfpc{d},mp_downdurs_lfpc{d}] = find_duration_ncx_uds_cycles(up_trans_inds,down_trans_inds,...
            mp_uplags{d},mp_downlags{d},lf8_period_vec);
        mp_corresp_lf8_mindur{d} = [mp_upskipped{d}.min_dur];
        mp_corresp_lf8_mindur{d}(nrt2_ups{d}) = lf8_updurs(corresp_lf8_upinds{d}(nrt2_ups{d}));
        rmp_updurs_lfpc{d} = mp_updurs_lfpc{d};
        rmp_updurs_lfpc{d}(mp_upskipped{d}.min_dur < thresh_lf8_downdur) = [];
        rmp_downdurs_lfpc{d} = mp_downdurs_lfpc{d};
        rmp_downdurs_lfpc{d}(mp_downskipped{d}.min_dur < thresh_lf8_updur) = [];
        
        %%
        median_uplag(d) = nanmedian(mp_uplags{d});
        median_uplag_t2(d) = nanmedian(mp_uplags{d}(rt2_ups{d}));
        median_uplag_nt2(d) = nanmedian(mp_uplags{d}(nrt2_ups{d}));
        median_reluplag_t2(d) = nanmedian(mp_reluplags{d}(rt2_ups{d}));
        median_reluplag_nt2(d) = nanmedian(mp_reluplags{d}(nrt2_ups{d}));
        median_downlag(d) = nanmedian(mp_downlags{d});
        median_downlag_t2(d) = nanmedian(mp_downlags{d}(rt2_ups{d}));
        median_downlag_nt2(d) = nanmedian(mp_downlags{d}(nrt2_ups{d}));
        median_reldownlag_t2(d) = nanmedian(mp_reldownlags{d}(rt2_ups{d}));
        median_reldownlag_nt2(d) = nanmedian(mp_reldownlags{d}(nrt2_ups{d}));
        mean_uplag(d) = nanmean(mp_uplags{d});
        mean_downlag(d) = nanmean(mp_downlags{d});
        median_reluplag(d) = nanmedian(mp_reluplags{d});
        median_reldownlag(d) = nanmedian(mp_reldownlags{d});
        mean_reluplag(d) = nanmean(mp_reluplags{d});
        mean_reldownlag(d) = nanmean(mp_reldownlags{d});
        
        
    end
end

cd C:\WC_Germany\persistent_downs\
save new_down_core_analysis_allEC


%%
mec = find_struct_field_vals(sess_data,'region','MEC');
lec = find_struct_field_vals(sess_data,'region','LEC');
layer3 = find_struct_field_vals(sess_data,'layer','3');
layer2 = find_struct_field_vals(sess_data,'layer','2');
layer23 = find_struct_field_vals(sess_data,'layer','2/3');
layer56 = find_struct_field_vals(sess_data,'layer','5/6');
layer35 = find_struct_field_vals(sess_data,'layer','3/5');
superficial = sort(unique([layer2 layer3 layer23]));
pyramidal = find_struct_field_vals(sess_data,'cell_type','pyramidal');
interneuron = find_struct_field_vals(sess_data,'cell_type','interneuron');
multipolar = find_struct_field_vals(sess_data,'cell_type','multipolar');
fan = find_struct_field_vals(sess_data,'cell_type','fan');
stellate = find_struct_field_vals(sess_data,'cell_type','stellate');
CA1 = find_struct_field_vals(sess_data,'region','CA1');
CA3 = find_struct_field_vals(sess_data,'region','CA3');
prefrontal = find_struct_field_vals(sess_data,'region','prefrontal');
frontal = find_struct_field_vals(sess_data,'region','frontal');
parietal = find_struct_field_vals(sess_data,'region','parietal');
V1 = find_struct_field_vals(sess_data,'region','V1');

DG = find_struct_field_vals(sess_data,'region','DG');
l3mec = intersect(mec,layer3);
l35mec = intersect(mec,layer35);
l56mec = intersect(mec,layer56);
deepmec = unique([l35mec l56mec]);
l3lec = intersect(lec,layer3);
l2mec = intersect(mec,layer2);
l2lec = intersect(lec,layer2);
CA1_pyr = intersect(CA1,pyramidal);
CA1_int = intersect(CA1,interneuron);
neocortex = unique([prefrontal(:); frontal(:); parietal(:); V1(:)]);

%%
cd C:\persDowns_paper\Figs\
group = nan(length(fract_rt2_ups),1);
group(l3mec) = 1;
group(l3lec) = 2;
group(prefrontal) = 3;
group(frontal) = 4;
group(parietal) = 5;

close all
names = {'L3-MEC','L3-LEC','Prefrontal','Frontal','Parietal'};
uset = find(~isnan(group));
figure
boxplot(fract_rt2_ups(uset)',names(group(uset)),'plotstyle','compact');
set(gca,'fontsize',12,'fontname','arial');
ylabel('Fraction persistent up states','fontsize',14);
box off
fillPage(gcf,'papersize',[5 5]);

figure
boxplot(fract_rt2_downs(uset)',names(group(uset)),'plotstyle','compact');
set(gca,'fontsize',12,'fontname','arial');
ylabel('Fraction persistent down states','fontsize',14);
box off
fillPage(gcf,'papersize',[5 5]);

figure
boxplot(median_downlag(uset)'/Fsd,names(group(uset)),'plotstyle','compact');
set(gca,'fontsize',12,'fontname','arial');
ylabel('Down-transition lag (s)','fontsize',14);
box off
fillPage(gcf,'papersize',[5 5]);

figure
boxplot(median_uplag(uset)'/Fsd,names(group(uset)),'plotstyle','compact');
set(gca,'fontsize',12,'fontname','arial');
ylabel('Up-transition lag (s)','fontsize',14);
box off
fillPage(gcf,'papersize',[5 5]);

%%
cd C:\WC_Germany\persistent_downs\
group = nan(length(fract_rt2_ups),1);
group(l2mec) = 1;
group(l3mec) = 2;
group(deepmec) = 3;
group(l2lec) = 4;
group(l3lec) = 5;
group(CA1_pyr) = 6;
group(CA1_int) = 7;
group(CA3) = 8;
group(DG) = 9;
group(neocortex) = 10;

close all
names = {'L2-MEC','L3-MEC','Deep-MEC','L2-LEC','L3-LEC','CA1-pyr','CA1-int','CA3','DG','Neocortex'};
uset = find(~isnan(group));
figure
boxplot(fract_rt2_ups(uset)',names(group(uset)),'plotstyle','compact');
set(gca,'fontsize',12,'fontname','arial');
fillPage(gcf,'papersize',[5 5]);
ylabel('Fraction persistent up states','fontsize',14);
box off

figure
boxplot(fract_rt2_downs(uset)',names(group(uset)),'plotstyle','compact');
set(gca,'fontsize',12,'fontname','arial');
fillPage(gcf,'papersize',[5 5]);
ylabel('Fraction persistent down states','fontsize',14);
box off
%%
uset = [l2mec l3mec deepmec l2lec l3lec CA1_int DG neocortex'];
figure
% plot(median_downlag(uset)/Fsd,fract_rt2_ups(uset),'.','markersize',8);
plot(mean_downlag(uset)/Fsd,fract_rt2_ups(uset),'.','markersize',8);
xlabel('Average down-transition lag (s)','fontsize',14);
ylabel('Fraction persistent ups','fontsize',14);
box off;
set(gca,'fontsize',12,'fontname','arial');
fillPage(gcf,'papersize',[5 5]);

figure
% plot(median_uplag(uset)/Fsd,fract_rt2_downs(uset),'.','markersize',8);
plot(mean_uplag(uset)/Fsd,fract_rt2_downs(uset),'.','markersize',8);
xlabel('Average up-transition lag (s)','fontsize',14);
ylabel('Fraction persistent downs','fontsize',14);
box off;
set(gca,'fontsize',12,'fontname','arial');
fillPage(gcf,'papersize',[5 5]);

%%
l3mec = find(strcmp(data_type(used_dirs),'L3MEC'));
l3lec = find(strcmp(data_type(used_dirs),'L3LEC'));
figure;
subplot(2,1,1)
plot(median_downlag(l3mec)/Fsd,fract_rt2_ups(l3mec),'o')
hold on
plot(median_downlag(l3lec)/Fsd,fract_rt2_ups(l3lec),'ro')

subplot(2,1,2)
plot(median_uplag(l3mec)/Fsd,fract_rt2_downs(l3mec),'o')
hold on
plot(median_uplag(l3lec)/Fsd,fract_rt2_downs(l3lec),'ro')