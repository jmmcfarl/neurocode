clear all
close all

cd ~/Data/bruce/2_27_12/M232
load ./random_bar_eyedata_ftime.mat

load ./lemM232Expts.mat

load ./CellList.mat
good_sus = find(all(CellList(bar_expts,:,1) > 0));
use_sus = 1:24;

%%
Fs = 1000; %assume sample rate exactly 1 kHz for LFP!
dsf = 3;
Fsd = Fs/dsf;
niqf = Fsd/2;

maxlag = round(Fsd*0.5); %time range around stimulus to compute triggered average
maxlag_spk = round(Fsd*0.5); %time range around stimulus to compute triggered average

%bandpass filter for LFP
bb_lcf = 1;
[b_bb,a_bb] = butter(2,[bb_lcf]/niqf,'high');

all_lfps = [];
all_lfp_t_axis = [];
all_expt = [];
all_binned_spks = [];
for ee = 1:length(bar_expts);
    fprintf('Expt %d of %d\n',ee,length(bar_expts));
    fname = sprintf('lemM232A.%d.lfp.mat',bar_expts(ee));
    load(fname);
    fname = sprintf('Expt%dClusterTimes.mat',bar_expts(ee));
    load(fname);
    
    %     Fs = 1/LFP.Header.CRsamplerate;
    
    n_trials(ee) = length(LFP.Trials);
    %     lfp_trial_starts = [LFP.Trials(:).Start]/1e4;
    lfp_trial_starts = [LFP.Trials(:).ftime]/1e4;
    lfp_trial_ends = [LFP.Trials(:).End]/1e4;
    
    expt_lfp_t_axis = [];
    expt_lfps = [];
    expt_binned_spks = [];
    cur_npts = nan(n_trials(ee),1);
    cur_t_end = nan(n_trials(ee),1);
    for tt = 1:n_trials(ee)
        
        cur_npts(tt) = size(LFP.Trials(tt).LFP,1);
        cur_t_end(tt) = lfp_trial_starts(tt)+(cur_npts(tt)-1)/Fs;
        cur_t_axis = lfp_trial_starts(tt):1/Fs:cur_t_end(tt);
        
        if ~isempty(expt_lfp_t_axis)
            cur_sp = find(cur_t_axis > max(expt_lfp_t_axis),1,'first');
        else
            cur_sp = 1;
        end
        cur_t_axis = cur_t_axis(cur_sp:end);
        cur_LFP = [LFP.Trials(tt).LFP];
        cur_LFP = cur_LFP(cur_sp:end,:);
        cur_LFP_d = [];
        for ch = 1:24
            cur_LFP_d = [cur_LFP_d; decimate(cur_LFP(:,ch),dsf)'];
        end
        cur_LFP_d = filtfilt(b_bb,a_bb,cur_LFP_d);
        
        cur_t_axis_d = downsample(cur_t_axis,dsf);
        expt_lfp_t_axis = [expt_lfp_t_axis; cur_t_axis_d(:)];
        expt_lfps = [expt_lfps; cur_LFP_d'];
        
        cur_binned_spks = nan(length(cur_t_axis_d),length(use_sus));
        for cc = 1:length(use_sus)
            cur_spikes = Clusters{use_sus(cc)}.times;
            cur_spikes = cur_spikes(cur_spikes >= cur_t_axis_d(1) & cur_spikes < cur_t_axis_d(end));
            cur_binned_spks(:,cc) = hist(cur_spikes,cur_t_axis_d);
        end
        expt_binned_spks = [expt_binned_spks; cur_binned_spks];
        
    end
    all_lfps = [all_lfps; expt_lfps];
    all_lfp_t_axis = [all_lfp_t_axis; expt_lfp_t_axis];
    all_expt = [all_expt; ee*ones(length(expt_lfp_t_axis),1)];
    all_binned_spks = [all_binned_spks; expt_binned_spks];
end

%%
load ./un_bar_pos.mat
un_bar_pos(1:2) = [];
flen = 20;

all_trial_start = [];
all_trial_stop = [];
all_trial_expt = [];
all_stim_starts = [];
all_stim_start_times = [];
all_stim_trials = [];
all_bar_locs = [];
all_stim_class =[];
all_rel_start_times = [];
all_bar_Xmat = [];
all_stim_expts = [];
for ee = 1:length(bar_expts)
    expt_trial_starts = [Expts{bar_expts(ee)}.Trials(:).TrialStart]/1e4;
    expt_trial_end = [Expts{bar_expts(ee)}.Trials(:).TrueEnd]/1e4;
    
    %     buff = round(0.2*100);
    buff = 1;
    cur_expt_start_times = [];
    cur_rel_times = [];
    cur_expt_bar_locs = [];
    cur_bar_Xmat = [];
    cur_stim_trial_num = [];
    for t = 1:length(expt_trial_starts)
        cur_expt_start_times = [cur_expt_start_times; Expts{bar_expts(ee)}.Trials(t).Start];
        cur_rel_times = [cur_rel_times; Expts{bar_expts(ee)}.Trials(t).Start/1e4 - expt_trial_starts(t)];
        temp_bar_locs = Expts{bar_expts(ee)}.Trials(t).Op;
        cur_expt_bar_locs = [cur_expt_bar_locs; temp_bar_locs];
        cur_stim_trial_num = [cur_stim_trial_num; ones(length(temp_bar_locs),1)*t];
        
        cur_bar_mat = zeros(length(temp_bar_locs),length(un_bar_pos));
        for bb = 1:length(un_bar_pos)
            cur_set = find(temp_bar_locs==un_bar_pos(bb));
            cur_bar_mat(cur_set,bb) = 1;
        end
        cur_bar_Xmat = [cur_bar_Xmat; makeStimRows(cur_bar_mat,flen)];
    end
    cur_expt_start_times = cur_expt_start_times/1e4;
    
    %     sac_start_times = all_eye_ts{ee}(all_sac_startinds{ee});
    %     sac_stop_times = all_eye_ts{ee}(all_sac_stopinds{ee});
    %
    %     cur_sac_set = find(all_expt_num==ee);
    %     if length(cur_sac_set) ~= length(sac_start_times)
    %         error('saccade mismatch')
    %     end
    %     intrial_sac_set = find(~isnan(all_trial_num(cur_sac_set)));
    %     infirst_sac_set = find(ismember(cur_sac_set,all_isfirst));
    %     insecond_sac_set = find(ismember(cur_sac_set,all_issecond));
    
    cur_lfp_set = find(all_expt == ee);
    
    stim_start_inds = round(interp1(all_lfp_t_axis(cur_lfp_set),1:length(cur_lfp_set),cur_expt_start_times));
    rid_stims = find(stim_start_inds < maxlag | stim_start_inds > length(cur_lfp_set) - maxlag);
    rid_stims = [rid_stims 1+find(diff(cur_expt_start_times) == 0)];
    stim_start_inds(rid_stims) = []; cur_expt_bar_locs(rid_stims) = []; 
    cur_rel_times(rid_stims) = []; cur_expt_start_times(rid_stims) = [];
    cur_stim_trial_num(rid_stims) = [];
    cur_bar_Xmat(rid_stims,:) = [];
    
    stim_class = zeros(size(cur_expt_start_times));
    stim_class(cur_rel_times < 0.8) = 1;
    stim_class(cur_rel_times > 0.8 & cur_rel_times < 1.5) = 2;
    
    all_stim_class = [all_stim_class(:); stim_class];
    all_stim_starts = [all_stim_starts(:); cur_lfp_set(stim_start_inds)];
    all_stim_start_times = [all_stim_start_times; cur_expt_start_times];
    all_stim_trials = [all_stim_trials; cur_stim_trial_num];
    all_stim_expts = [all_stim_expts(:); ones(length(stim_start_inds),1)*ee];
    all_bar_locs = [all_bar_locs(:); cur_expt_bar_locs];
    all_trial_start = [all_trial_start expt_trial_starts];
    all_trial_stop = [all_trial_stop expt_trial_end];
    all_trial_expt = [all_trial_expt ones(1,length(expt_trial_starts))*ee];
    all_rel_start_times = [all_rel_start_times; cur_rel_times];
    all_bar_Xmat = [all_bar_Xmat; cur_bar_Xmat];
end

%%
spk_sm = round(Fsd*0.01);
sm_binned_spks = zeros(size(all_binned_spks));
for i = 1:length(use_sus)
    sm_binned_spks(:,i) = jmm_smooth_1d_cor(all_binned_spks(:,i),spk_sm);
    
end

%%
load ./random_bar_unit_models_ftime.mat
all_pred_rate = [];
for ee = 1:length(bar_expts)
    cur_set_lfp = find(all_expt==ee);
    cur_set = find(all_stim_expts==ee);
    pred_rate_set = nan(length(cur_set_lfp),length(use_sus));
    for i = 1:length(use_sus)
        [~, ~, ~, prate] = getLL_GNM(glm_fit(i),all_bar_Xmat(cur_set,:),[],'none');
        interp_rate = interp1(all_stim_start_times(cur_set),prate,all_lfp_t_axis(cur_set_lfp));
        pred_rate_set(:,i) = interp_rate;
    end
    all_pred_rate = [all_pred_rate; pred_rate_set];
end
%%
lags = (-maxlag:maxlag);
lags_spk = (-maxlag_spk:maxlag_spk);

st_bounds = [0.3 1.7];

[un_bar_pos,ia,ic] = unique(all_bar_locs);

stim_trig_avgs = nan(length(un_bar_pos),length(lags),24);
stim_trig_avg_spks = nan(length(un_bar_pos),length(lags_spk),length(use_sus));
stim_trig_avg_pred = nan(length(un_bar_pos),length(lags_spk),length(use_sus));
for i = 1:length(un_bar_pos)
    fprintf('Bar pos %d of %d\n',i,length(un_bar_pos));
    cur_barset = find(all_bar_locs==un_bar_pos(i) & all_rel_start_times >= st_bounds(1) & all_rel_start_times <= st_bounds(2));
    for j = 1:length(lags)
        stim_trig_avgs(i,j,:) =  mean(all_lfps(all_stim_starts(cur_barset)+lags(j),:));
    end
    for j = 1:length(lags_spk)
        stim_trig_avg_spks(i,j,:) =  mean(sm_binned_spks(all_stim_starts(cur_barset)+lags_spk(j),:));
        stim_trig_avg_pred(i,j,:) =  nanmean(all_pred_rate(all_stim_starts(cur_barset)+lags_spk(j),:));
    end
end
ms_stim_trig_avgs = bsxfun(@minus,stim_trig_avgs,mean(stim_trig_avgs));


%%
% for i = 1:24
%     i
%     subplot(2,1,1)
%     imagesc(lags/Fsd,1:14,squeeze(stim_trig_avgs(5:end,:,i)));
%     xlim([-0.1 0.3])
%     ca = caxis();
%     ca = max(abs(ca));
%     caxis([-ca ca])
%     subplot(2,1,2)
%     imagesc(lags_spk/Fsd,1:14,squeeze(stim_trig_avg_pred(5:end,:,i)))
%     xlim([-0.1 0.3])
%     pause
%     clf
% end

%%
% lags = (-maxlag:maxlag);
%
% [un_bar_pos,ia,ic] = unique(all_bar_locs);
%
% stim_trig_avgs0 = nan(length(un_bar_pos),length(lags),24);
% stim_trig_avgs1 = nan(length(un_bar_pos),length(lags),24);
% stim_trig_avgs2 = nan(length(un_bar_pos),length(lags),24);
% for i = 1:length(un_bar_pos)
%     cur_barset = find(all_bar_locs==un_bar_pos(i));
%     cur_set0 = cur_barset(find(all_stim_class(cur_barset) == 0));
%     cur_set1 = cur_barset(find(all_stim_class(cur_barset) == 1));
%     cur_set2 = cur_barset(find(all_stim_class(cur_barset) == 2));
%     for j = 1:length(lags)
%         stim_trig_avgs0(i,j,:) =  mean(all_lfps(all_stim_starts(cur_set0)+lags(j),:));
%         stim_trig_avgs1(i,j,:) =  mean(all_lfps(all_stim_starts(cur_set1)+lags(j),:));
%         stim_trig_avgs2(i,j,:) =  mean(all_lfps(all_stim_starts(cur_set2)+lags(j),:));
%     end
% end
%%
% cur_ch = 14;
% i = 10;
%    cur_barset = find(all_bar_locs==un_bar_pos(i));
% nstim_trig_avgs = zeros(length(cur_barset),length(lags));
% for j = 1:length(cur_barset)
%    cur = (all_stim_starts(cur_barset(j)) - maxlag):(all_stim_starts(cur_barset(j)) + maxlag);
%    nstim_trig_avgs(j,:) = all_lfps(cur,cur_ch);
% end

%%
clear lin_pred
for cur_ch = 1:24;
    lin_pred(cur_ch,:) = zeros(1,length(all_lfp_t_axis));
    for i = 1:length(un_bar_pos)
        cur_lin_filt = squeeze(ms_stim_trig_avgs(i,:,cur_ch));
        %                    cur_lin_filt = cur_lin_filt - mean(cur_lin_filt);
        cur_barset = find(all_bar_locs==un_bar_pos(i));
        null_sig = zeros(size(all_lfp_t_axis));
        null_sig(all_stim_starts(cur_barset)) = 1;
        cur_pred = conv(null_sig,cur_lin_filt,'same');
        lin_pred(cur_ch,:) = lin_pred(cur_ch,:) + cur_pred';
    end
end
lin_pred = lin_pred';
lin_res = all_lfps-lin_pred;
%%
% params.tapers = [3 5];
% params.Fs = Fsd;
% movingwin = [1.3 1.3];
% params.fpass = [1 150];
% buffer_w = round(Fsd*0.3);
% for ee = 1:length(bar_expts)
%     cur_trial_set = find(all_trial_expt==ee);
%     cur_lfp_set = find(all_expt == ee);
%     trial_start_inds = round(interp1(all_lfp_t_axis(cur_lfp_set),1:length(cur_lfp_set),all_trial_start(cur_trial_set)));
%     trial_stop_inds = round(interp1(all_lfp_t_axis(cur_lfp_set),1:length(cur_lfp_set),all_trial_stop(cur_trial_set)));
%     sMarkers = [trial_start_inds(:)+buffer_w trial_stop_inds(:)-buffer_w];
%     for cur_ch = 1:24
%         cur_data = [all_lfps(cur_lfp_set,cur_ch) lin_pred(cur_lfp_set,cur_ch)];
%         [Cmn(ee,cur_ch,:),Phimn,Smn,Smm,f] = coherencyc_unequal_length_trials( cur_data, movingwin, params, sMarkers );
%     end
% end
% avg_Cmn = squeeze(mean(Cmn));

%%
% cur_ch = 14;
% cur_ch2 = 2;
% cur_set = find(all_expt==2);
% figure
% plot(all_lfp_t_axis(cur_set),all_lfps(cur_set,cur_ch));
% hold on
% plot(all_lfp_t_axis(cur_set),lin_pred(cur_ch,cur_set),'r');
% plot(all_lfp_t_axis(cur_set),all_lfps(cur_set,cur_ch2),'k');
% plot(all_lfp_t_axis(cur_set),lin_pred(cur_ch2,cur_set),'g');

%%
% close all
% use_expt = 2;
% test_ch = 14;
% 
% cur_trial_set = find(all_trial_expt==use_expt);
% cur_lfp_set = find(all_expt==use_expt);
% trial_start_inds = round(interp1(all_lfp_t_axis(cur_lfp_set),1:length(cur_lfp_set),all_trial_start(cur_trial_set)));
% trial_stop_inds = round(interp1(all_lfp_t_axis(cur_lfp_set),1:length(cur_lfp_set),all_trial_stop(cur_trial_set)));
% 
% buffer = round(Fs*0.0);
% for i = 1:length(cur_trial_set)
%     cur_inds = (trial_start_inds(i)-buffer):(trial_stop_inds(i)+buffer);
% 
%     subplot(2,2,1)
%     imagesc((cur_inds-cur_inds(1))/Fs,1:24,all_lfps(cur_inds,:)');
%     caxis([-2 2])
%     subplot(2,2,3)
%     imagesc((cur_inds-cur_inds(1))/Fs,1:24,lin_pred(cur_inds,:)');
%     caxis([-2 2])
%     subplot(2,2,2)
%     plot(all_lfp_t_axis(cur_inds),lin_pred(cur_inds,test_ch),'b')
%     hold on
%     plot(all_lfp_t_axis(cur_inds),all_lfps(cur_inds,test_ch),'r')
%     plot(all_lfp_t_axis(cur_inds),all_lfps(cur_inds,test_ch)-lin_pred(cur_inds,test_ch),'k')
%     axis tight
%     subplot(2,2,4)
%     imagesc((cur_inds-cur_inds(1))/Fs,1:24,all_lfps(cur_inds,:)'-lin_pred(cur_inds,:)');
%     caxis([-2 2])
%     pause
%     clf
% 
% end


%% COMPUTE OBSERVED AND PREDICTED PHASE AND AMPLITUDE SPECTRA
scales = logspace(log10(5.5),log10(80),25);
% scales = [scales];
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);
nwfreqs = length(wfreqs);

all_res_phasegram = [];
all_res_ampgram = [];
all_pred_phasegram = [];
all_pred_ampgram = [];
for ee = 1:length(bar_expts)
    fprintf('Computing CWT Expt %d of %d\n',ee,length(bar_expts));
    cur_lfp_set = find(all_expt==ee);
    
    cur_phasegram = nan(length(cur_lfp_set),length(wfreqs),24);
    cur_ampgram = nan(length(cur_lfp_set),length(wfreqs),24);
    pred_phasegram = nan(length(cur_lfp_set),length(wfreqs),24);
    pred_ampgram = nan(length(cur_lfp_set),length(wfreqs),24);
    for ll = 1:24
        temp = cwt(lin_res(cur_lfp_set,ll),scales,'cmor1-1');
        cur_phasegram(:,:,ll) = angle(temp)';
        cur_ampgram(:,:,ll) = abs(temp)';
        temp = cwt(lin_pred(cur_lfp_set,ll),scales,'cmor1-1');
        pred_phasegram(:,:,ll) = angle(temp)';
        pred_ampgram(:,:,ll) = abs(temp)';
    end
    
    all_res_phasegram = cat(1,all_res_phasegram,cur_phasegram);
    all_res_ampgram = cat(1,all_res_ampgram,cur_ampgram);
    all_pred_phasegram = cat(1,all_pred_phasegram,pred_phasegram);
    all_pred_ampgram = cat(1,all_pred_ampgram,pred_ampgram);
    
end
all_res_n_ampgram = bsxfun(@rdivide,all_res_ampgram,std(all_res_ampgram));
all_pred_n_ampgram = bsxfun(@rdivide,all_pred_ampgram,std(all_pred_ampgram));

%% COMPUTE USED INDS
buffer = round(0.25*Fsd);
use_inds = [];
interp_inds = [];
for ee = 1:length(bar_expts)
   
    cur_lfp_set = find(all_expt==ee);
    cur_trial_set = find(all_trial_expt==ee);
    trial_start_inds = round(interp1(all_lfp_t_axis(cur_lfp_set),1:length(cur_lfp_set),all_trial_start(cur_trial_set)));
    trial_stop_inds = round(interp1(all_lfp_t_axis(cur_lfp_set),1:length(cur_lfp_set),all_trial_stop(cur_trial_set)));
    bad = find(isnan(trial_start_inds) | isnan(trial_stop_inds));
    trial_start_inds(bad) = []; trial_stop_inds(bad) = [];
    cur_use_inds = [];
    for i = 1:length(trial_start_inds)
        cur_use_inds = [cur_use_inds (trial_start_inds(i)+buffer):(trial_stop_inds(i)-buffer)];
    end
    use_inds = [use_inds; cur_lfp_set(cur_use_inds(:))];
    
    cur_set = find(all_stim_expts==ee);
    cur_int = round(interp1(all_stim_start_times(cur_set),1:length(cur_set),all_lfp_t_axis(cur_lfp_set(cur_use_inds))));
    %     cur_int = round(interp1(all_lfp_t_axis(cur_lfp_set),1:length(cur_lfp_set),all_stim_start_times(cur_set)));
    uu = find(~isnan(cur_int));
    temp = nan(size(cur_int));
    temp(uu) = cur_set(cur_int(uu));
    interp_inds = [interp_inds; temp(:)];
end

use_inds(isnan(interp_inds)) = [];
interp_inds(isnan(interp_inds)) = [];

[c,ia,ic] = unique([all_stim_trials all_stim_expts],'rows');
n_trials = length(ia);

xv_frac = 0.2;
n_xv_trials = round(n_trials*xv_frac);
xv_set = randperm(n_trials);
xv_set(n_xv_trials+1:end) = [];
xv_inds = find(ismember(ic,xv_set));
tr_inds = find(~ismember(ic,xv_set))';

xv_inds_new = find(ismember(interp_inds,xv_inds));
tr_inds_new = find(ismember(interp_inds,tr_inds));

%%
NL_type = 1;
reg_params.dl1_ind = 0;
reg_params.dl1_dep = 0;
reg_params.dl2_ind = 10000;
reg_params.dl2_dep = 10000;
reg_params.dl2_freq_ind = 5000;
reg_params.dl2_freq_dep = 5000;
reg_params.dl2_time_ind = 0;
reg_params.dl2_time_dep = 0;
reg_params.l1_ind = 0;
reg_params.l1_dep = 0;
reg_params.l1t_ind = 0;
reg_params.l1t_dep = 0;
reg_params.is_phase = 0;


res_phase_set = [reshape(cos(all_res_phasegram(use_inds,:,:)),length(use_inds),length(wfreqs)*24) reshape(sin(all_res_phasegram(use_inds,:,:)),length(use_inds),length(wfreqs)*24)];
pred_phase_set = [reshape(cos(all_pred_phasegram(use_inds,:,:)),length(use_inds),length(wfreqs)*24) reshape(sin(all_pred_phasegram(use_inds,:,:)),length(use_inds),length(wfreqs)*24)];

res_ampphase_set = [reshape(all_res_n_ampgram(use_inds,:,:),length(use_inds),length(wfreqs)*24).*reshape(cos(all_res_phasegram(use_inds,:,:)),length(use_inds),length(wfreqs)*24) ...
    reshape(all_res_n_ampgram(use_inds,:,:),length(use_inds),length(wfreqs)*24).*reshape(sin(all_res_phasegram(use_inds,:,:)),length(use_inds),length(wfreqs)*24)];
pred_ampphase_set = [reshape(all_pred_n_ampgram(use_inds,:,:),length(use_inds),length(wfreqs)*24).*reshape(cos(all_pred_phasegram(use_inds,:,:)),length(use_inds),length(wfreqs)*24) ...
    reshape(all_pred_n_ampgram(use_inds,:,:),length(use_inds),length(wfreqs)*24).*reshape(sin(all_pred_phasegram(use_inds,:,:)),length(use_inds),length(wfreqs)*24)];

phase_elec_set = [repmat(1:24,1,length(wfreqs)) repmat(1:24,1,length(wfreqs))];
phase_elec_set = phase_elec_set(:);

%%
silent = 1;
for cc = 1:24
    fprintf('Cell %d of %d\n',cc,24);
    Robs = all_binned_spks(use_inds(tr_inds_new),cc);
    tr_spkbns = convert_to_spikebins(Robs);
   
    glm_kern = get_k_mat(glm_fit(cc));
    stim_out = all_bar_Xmat*glm_kern;
    stim_out_interp = stim_out(interp_inds(tr_inds_new));

    Xmat = [stim_out_interp];
    klen = size(Xmat,2);
    K0 = zeros(klen+1,1);
    [fitp_so,grad] = GLMsolve_jmm(Xmat, tr_spkbns, K0, 1, [], [],[], [], [], [], NL_type);

    use_elecs = 1:24;
%     use_elecs(cc) = [];
    use_set = find(ismember(phase_elec_set,use_elecs));
    
    Xmat = [res_phase_set(tr_inds_new,use_set) stim_out_interp];
    stim_params = [length(wfreqs),length(use_elecs)];
    klen = size(Xmat,2);
    K0 = zeros(klen+1,1);
    [fitp_rphase] = fit_GLM_phase_model(Xmat, Robs, K0,silent, stim_params, reg_params,0,NL_type);
    res_phase_cfilt(cc,:) = fitp_rphase.k(1:length(use_elecs)*length(wfreqs));
    res_phase_sfilt(cc,:) = fitp_rphase.k((length(use_elecs)*length(wfreqs)+1):length(use_elecs)*length(wfreqs)*2);
    res_phase_g(cc) = fitp_rphase.k(end-1);
    
    
    Xmat = [pred_phase_set(tr_inds_new,use_set) stim_out_interp];
    stim_params = [length(wfreqs),length(use_elecs)];
    klen = size(Xmat,2);
    K0 = zeros(klen+1,1);
    [fitp_pphase] = fit_GLM_phase_model(Xmat, Robs, K0,silent, stim_params, reg_params,0,NL_type);
    pred_phase_cfilt(cc,:) = fitp_pphase.k(1:length(use_elecs)*length(wfreqs));
    pred_phase_sfilt(cc,:) = fitp_pphase.k((length(use_elecs)*length(wfreqs)+1):length(use_elecs)*length(wfreqs)*2);
    pred_phase_g(cc) = fitp_pphase.k(end-1);

    
    Xmat = [res_ampphase_set(tr_inds_new,use_set) stim_out_interp];
    stim_params = [length(wfreqs),length(use_elecs)];
    klen = size(Xmat,2);
    K0 = zeros(klen+1,1);
    [fitp_raphase] = fit_GLM_phase_model(Xmat, Robs, K0,silent, stim_params, reg_params,0,NL_type);
    res_ampphase_cfilt(cc,:) = fitp_raphase.k(1:length(use_elecs)*length(wfreqs));
    res_ampphase_sfilt(cc,:) = fitp_raphase.k((length(use_elecs)*length(wfreqs)+1):length(use_elecs)*length(wfreqs)*2);
    res_ampphase_g(cc) = fitp_raphase.k(end-1);
    
    Xmat = [pred_ampphase_set(tr_inds_new,use_set) stim_out_interp];
    stim_params = [length(wfreqs),length(use_elecs)];
    klen = size(Xmat,2);
    K0 = zeros(klen+1,1);
    [fitp_paphase] = fit_GLM_phase_model(Xmat, Robs, K0,silent, stim_params, reg_params,0,NL_type);
    pred_ampphase_cfilt(cc,:) = fitp_paphase.k(1:length(use_elecs)*length(wfreqs));
    pred_ampphase_sfilt(cc,:) = fitp_paphase.k((length(use_elecs)*length(wfreqs)+1):length(use_elecs)*length(wfreqs)*2);
    pred_ampphase_g(cc) = fitp_paphase.k(end-1);
    
    Xmat = [res_ampphase_set(tr_inds_new,use_set)];
    stim_params = [length(wfreqs),length(use_elecs)];
    klen = size(Xmat,2);
    K0 = zeros(klen+1,1);
    [fitp_raphase_ns] = fit_GLM_phase_model(Xmat, Robs, K0,silent, stim_params, reg_params,0,NL_type);
    res_ampphase_ns_cfilt(cc,:) = fitp_raphase_ns.k(1:length(use_elecs)*length(wfreqs));
    res_ampphase_ns_sfilt(cc,:) = fitp_raphase_ns.k((length(use_elecs)*length(wfreqs)+1):length(use_elecs)*length(wfreqs)*2);

    Xmat = [pred_ampphase_set(tr_inds_new,use_set)];
    stim_params = [length(wfreqs),length(use_elecs)];
    klen = size(Xmat,2);
    K0 = zeros(klen+1,1);
    [fitp_paphase_ns] = fit_GLM_phase_model(Xmat, Robs, K0,silent, stim_params, reg_params,0,NL_type);
    pred_ampphase_ns_cfilt(cc,:) = fitp_paphase_ns.k(1:length(use_elecs)*length(wfreqs));
    pred_ampphase_ns_sfilt(cc,:) = fitp_paphase_ns.k((length(use_elecs)*length(wfreqs)+1):length(use_elecs)*length(wfreqs)*2);

    %%
    xv_Robs = all_binned_spks(use_inds(xv_inds_new),cc);
    stim_out_interp = stim_out(interp_inds(xv_inds_new));
    
    Xmat = [stim_out_interp];
    xv_so_pred_rate = Xmat*fitp_so.k(1:end-1) + fitp_so.k(end);
    if NL_type == 1
        xv_so_pred_rate = exp(xv_so_pred_rate);
    else
        xv_so_pred_rate = log(1+exp(xv_so_pred_rate));
    end
    xv_so_LL(cc) = -sum(xv_Robs.*log(xv_so_pred_rate)-xv_so_pred_rate)/sum(xv_Robs);
            
    Xmat = [res_phase_set(xv_inds_new,use_set) stim_out_interp];
    xv_rphase_pred_rate = Xmat*fitp_rphase.k(1:end-1) + fitp_rphase.k(end);
    if NL_type == 0
        xv_rphase_pred_rate = log(1+exp(xv_rphase_pred_rate));
    else
        xv_rphase_pred_rate = exp(xv_rphase_pred_rate);
    end
    xv_rphase_LL(cc) = -sum(xv_Robs.*log(xv_rphase_pred_rate)-xv_rphase_pred_rate)/sum(xv_Robs);
 
    Xmat = [pred_phase_set(xv_inds_new,use_set) stim_out_interp];
    xv_pphase_pred_rate = Xmat*fitp_pphase.k(1:end-1) + fitp_pphase.k(end);
    if NL_type == 0
        xv_pphase_pred_rate = log(1+exp(xv_pphase_pred_rate));
    else
        xv_pphase_pred_rate = exp(xv_pphase_pred_rate);
    end
    xv_pphase_LL(cc) = -sum(xv_Robs.*log(xv_pphase_pred_rate)-xv_pphase_pred_rate)/sum(xv_Robs);
    
    Xmat = [res_ampphase_set(xv_inds_new,use_set) stim_out_interp];
    xv_raphase_pred_rate = Xmat*fitp_raphase.k(1:end-1) + fitp_raphase.k(end);
    if NL_type == 0
        xv_raphase_pred_rate = log(1+exp(xv_raphase_pred_rate));
    else
        xv_raphase_pred_rate = exp(xv_raphase_pred_rate);
    end
    xv_raphase_LL(cc) = -sum(xv_Robs.*log(xv_raphase_pred_rate)-xv_raphase_pred_rate)/sum(xv_Robs);
 
    Xmat = [pred_ampphase_set(xv_inds_new,use_set) stim_out_interp];
    xv_paphase_pred_rate = Xmat*fitp_paphase.k(1:end-1) + fitp_paphase.k(end);
    if NL_type == 0
        xv_paphase_pred_rate = log(1+exp(xv_paphase_pred_rate));
    else
        xv_paphase_pred_rate = exp(xv_paphase_pred_rate);
    end
    xv_paphase_LL(cc) = -sum(xv_Robs.*log(xv_paphase_pred_rate)-xv_paphase_pred_rate)/sum(xv_Robs);
    
    Xmat = [res_ampphase_set(xv_inds_new,use_set)];
    xv_raphase_ns_pred_rate = Xmat*fitp_raphase_ns.k(1:end-1) + fitp_raphase_ns.k(end);
    if NL_type == 0
        xv_raphase_ns_pred_rate = log(1+exp(xv_raphase_ns_pred_rate));
    else
        xv_raphase_ns_pred_rate = exp(xv_raphase_ns_pred_rate);
    end
    xv_raphase_ns_LL(cc) = -sum(xv_Robs.*log(xv_raphase_ns_pred_rate)-xv_raphase_ns_pred_rate)/sum(xv_Robs);

    Xmat = [pred_ampphase_set(xv_inds_new,use_set)];
    xv_paphase_ns_pred_rate = Xmat*fitp_paphase_ns.k(1:end-1) + fitp_paphase_ns.k(end);
    if NL_type == 0
        xv_paphase_ns_pred_rate = log(1+exp(xv_paphase_ns_pred_rate));
    else
        xv_paphase_ns_pred_rate = exp(xv_paphase_ns_pred_rate);
    end
    xv_paphase_ns_LL(cc) = -sum(xv_Robs.*log(xv_paphase_ns_pred_rate)-xv_paphase_ns_pred_rate)/sum(xv_Robs);
    
    avg_rate = mean(Robs);
    null_pred = avg_rate*ones(size(xv_Robs));
    xv_null_LL(cc) = -sum(xv_Robs.*log(null_pred)-null_pred)/sum(xv_Robs);


end
%%
use_c = 1:24;
xv_so_imp = xv_null_LL(use_c) - xv_so_LL(use_c);
% xv_phase_imp = xv_so_LL(use_c) - xv_phase_LL(use_c);
% xv_rphase_imp = xv_so_LL(use_c) - xv_rphase_LL(use_c);
xv_paphase_imp = xv_so_LL(use_c) - xv_paphase_LL(use_c);
xv_paphase_ns_imp = xv_null_LL(use_c) - xv_paphase_ns_LL(use_c);
xv_raphase_imp = xv_so_LL(use_c) - xv_raphase_LL(use_c);
xv_raphase_ns_imp = xv_null_LL(use_c) - xv_raphase_ns_LL(use_c);

% res_phase_phasekern = -atan2(res_phase_cfilt,res_phase_sfilt)+pi/2;
% pred_phase_phasekern = -atan2(pred_phase_cfilt,pred_phase_sfilt)+pi/2;
% pred_phase_ampkern = sqrt(pred_phase_cfilt.^2+pred_phase_sfilt.^2);
% res_phase_ampkern = sqrt(res_phase_cfilt.^2+res_phase_sfilt.^2);

res_ampphase_ampkern = sqrt(res_ampphase_cfilt.^2+res_ampphase_sfilt.^2);
res_ampphase_phasekern = -atan2(res_ampphase_cfilt,res_ampphase_sfilt)+pi/2;
pred_ampphase_ampkern = sqrt(pred_ampphase_cfilt.^2+pred_ampphase_sfilt.^2);
pred_ampphase_phasekern = -atan2(pred_ampphase_cfilt,pred_ampphase_sfilt)+pi/2;

%%
close all
for cc = 1:24
    subplot(2,2,1)
    pcolor(wfreqs,1:24,reshape(res_ampphase_ampkern(cc,:),length(wfreqs),24)');shading flat
%     set(gca,'xscale','log')
    xlabel('Frequency (Hz)','fontsize',16)
    ylabel('Channel','fontsize',16)
%     title('Cosine kernel')
    subplot(2,2,3)
    pcolor(wfreqs,1:24,reshape(pred_ampphase_ampkern(cc,:),length(wfreqs),24)');shading flat
%     set(gca,'xscale','log')
    xlabel('Frequency (Hz)','fontsize',16)
    ylabel('Channel','fontsize',16)
%     title('Sin kernel')

    subplot(2,2,2)
    pcolor(wfreqs,1:24,reshape(res_ampphase_phasekern(cc,:),length(wfreqs),24)');shading flat
%     set(gca,'xscale','log')
    xlabel('Frequency (Hz)','fontsize',16)
    ylabel('Channel','fontsize',16)
%     title('Cosine kernel')
    subplot(2,2,4)
    pcolor(wfreqs,1:24,reshape(pred_ampphase_phasekern(cc,:),length(wfreqs),24)');shading flat
%     set(gca,'xscale','log')
    xlabel('Frequency (Hz)','fontsize',16)
    ylabel('Channel','fontsize',16)
%     title('Sin kernel')

cc
    pause
    clf
end
