clear all
close all

addpath(genpath('~/James_scripts/iCSD/'))

ExptNum = 235;
cd(['~/Data/bruce/M' num2str(ExptNum)]);
load ./random_bar_eyedata_ftime.mat

load(['lemM' num2str(ExptNum) 'Expts.mat']);
load ./bar_params.mat

if ExptNum == 235
    bar_expts(bar_expts==51) = []; %this has different set of stim positions
end
if ExptNum == 239
    bar_expts(bar_expts==40) = []; %this has different set of stim positions
end

%%
Fs = 1000; %assume sample rate exactly 1 kHz for LFP!
niqf = Fs/2;
dsf = 4;
Fsd = Fs/dsf;

maxlag = round(Fsd*0.4); %time range around stimulus to compute triggered average

%bandpass filter for LFP
bb_lcf = 1;
bb_hcf = 50;
[b_bb,a_bb] = butter(2,[bb_lcf bb_hcf]/niqf);

all_lfps = [];
all_lfp_t_axis = [];
all_expt = [];
for ee = 1:length(bar_expts);
    fprintf('Expt %d of %d\n',ee,length(bar_expts));
    fname = sprintf('lemM%dA.%d.lfp.mat',ExptNum,bar_expts(ee));
    load(fname);
    
%     Fs = 1/LFP.Header.CRsamplerate;
    
    n_trials(ee) = length(LFP.Trials);
%     lfp_trial_starts = [LFP.Trials(:).Start]/1e4;
    lfp_trial_starts = [LFP.Trials(:).ftime]/1e4;
    lfp_trial_ends = [LFP.Trials(:).End]/1e4;
    
    expt_lfp_t_axis = [];
    expt_lfps = [];
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
        cur_LFP = filtfilt(b_bb,a_bb,cur_LFP);
        
        cur_LFP = downsample(cur_LFP,dsf);
        cur_t_axis = downsample(cur_t_axis,dsf);
        
        if length(cur_LFP) > 10
        %compute CSD
        vars.Fs = Fsd;
        vars.BrainBound = 1;
        vars.ChanSep = 0.05;
        vars.diam = 2; %0.5
        CSD = PettersenCSD(cur_LFP','spline',vars)';
        cur_LFP = CSD;
        
        expt_lfp_t_axis = [expt_lfp_t_axis; cur_t_axis(:)];
        expt_lfps = [expt_lfps; cur_LFP];
        end 
    end
    all_lfps = [all_lfps; expt_lfps];
    all_lfp_t_axis = [all_lfp_t_axis; expt_lfp_t_axis];
    all_expt = [all_expt; ee*ones(length(expt_lfp_t_axis),1)];
end

%%
all_trial_start = [];
all_trial_stop = [];
all_trial_expt = [];
all_stim_starts = [];
all_bar_locs = [];
for ee = 1:length(bar_expts)
    expt_trial_starts = [Expts{bar_expts(ee)}.Trials(:).TrialStart]/1e4;
    expt_trial_end = [Expts{bar_expts(ee)}.Trials(:).TrueEnd]/1e4;
    
    cur_expt_start_times = [];
    cur_expt_bar_locs = [];
    for t = 1:length(expt_trial_starts)
       cur_expt_start_times = [cur_expt_start_times; Expts{bar_expts(ee)}.Trials(t).Start];
       cur_expt_bar_locs = [cur_expt_bar_locs; Expts{bar_expts(ee)}.Trials(t).Op];
    end
    cur_expt_start_times = cur_expt_start_times/1e4;
    
    
    cur_lfp_set = find(all_expt == ee);
    
    stim_start_inds = round(interp1(all_lfp_t_axis(cur_lfp_set),1:length(cur_lfp_set),cur_expt_start_times));
    rid_stims = find(stim_start_inds < maxlag | stim_start_inds > length(cur_lfp_set) - maxlag);
    stim_start_inds(rid_stims) = []; cur_expt_bar_locs(rid_stims) = [];

    all_stim_starts = [all_stim_starts(:); cur_lfp_set(stim_start_inds)];
    all_bar_locs = [all_bar_locs(:); cur_expt_bar_locs];
    all_trial_start = [all_trial_start expt_trial_starts];
    all_trial_stop = [all_trial_stop expt_trial_end];
    all_trial_expt = [all_trial_expt ones(1,length(expt_trial_starts))*ee];
end
%%
lags = (-maxlag:maxlag);

[un_bar_pos,ia,ic] = unique(all_bar_locs);

stim_trig_avgs = nan(length(un_bar_pos),length(lags),24);
for i = 1:length(un_bar_pos)
   cur_barset = find(all_bar_locs==un_bar_pos(i));
   for j = 1:24
       stim_trig_avgs(i,:,j) = get_event_trig_avg(all_lfps(:,j),all_stim_starts(cur_barset),maxlag,maxlag);
   end
   
%    for j = 1:length(lags)
%       stim_trig_avgs(i,j,:) =  mean(all_lfps(all_stim_starts(cur_barset)+lags(j),:));
%    end
end
ms_stim_trig_avgs = bsxfun(@minus,stim_trig_avgs,mean(stim_trig_avgs));

%%
% ca_ov = max(abs(ms_stim_trig_avgs(:)));
close all
for ii = 1:length(un_bar_pos)
    temp_mat = squeeze(ms_stim_trig_avgs(ii,:,:))';
    imagesc(lags/Fsd,1:24,temp_mat);
    ca = max(abs(temp_mat(:)));
    caxis([-ca ca]*0.95);
% caxis([-ca_ov ca_ov]*0.95);
    colorbar
    xlim([-0.1 0.3]);
    title(sprintf('Bar pos %.3f',un_bar_pos(ii)));
    pause
    clf
end
    
%%
for cc = 1:24
     temp_mat = squeeze(ms_stim_trig_avgs(:,:,cc));
    imagesc(lags/Fsd,1:length(un_bar_pos),temp_mat);
    ca = max(abs(temp_mat(:)));
    caxis([-ca ca]*0.95);
    xlim([-0.1 0.3]);
    title(sprintf('CH %d',cc));
    pause
    clf
end

%%
% clear lin_pred
% for cur_ch = 1:24;
%     lin_pred(cur_ch,:) = zeros(1,length(all_lfp_t_axis));
%     for i = 1:length(un_bar_pos)
%         cur_lin_filt = squeeze(ms_stim_trig_avgs(i,:,cur_ch));
%         %    cur_lin_filt = cur_lin_filt - mean(cur_lin_filt);
%         cur_barset = find(all_bar_locs==un_bar_pos(i));
%         null_sig = zeros(size(all_lfp_t_axis));
%         null_sig(all_stim_starts(cur_barset)) = 1;
%         cur_pred = conv(null_sig,cur_lin_filt,'same');
%         lin_pred(cur_ch,:) = lin_pred(cur_ch,:) + cur_pred';
%     end
% end
% 
% %%
% params.tapers = [4 7];
% params.Fs = Fs;
% movingwin = [1 1];
% params.fpass = [1 150];
% for ee = 1:length(bar_expts)
%     cur_trial_set = find(all_trial_expt==ee);
%     cur_lfp_set = find(all_expt == ee);
%     trial_start_inds = round(interp1(all_lfp_t_axis(cur_lfp_set),1:length(cur_lfp_set),all_trial_start(cur_trial_set)));
%     trial_stop_inds = round(interp1(all_lfp_t_axis(cur_lfp_set),1:length(cur_lfp_set),all_trial_stop(cur_trial_set)));
%     sMarkers = [trial_start_inds(:) trial_stop_inds(:)];
%     for cur_ch = 1:24
%         cur_data = [all_lfps(cur_lfp_set,cur_ch) lin_pred(cur_ch,cur_lfp_set)'];
%         [Cmn(ee,cur_ch,:),Phimn,Smn,Smm,f] = coherencyc_unequal_length_trials( cur_data, movingwin, params, sMarkers );
%     end
% end
% avg_Cmn = squeeze(mean(Cmn));
% 
% %%
% % cur_ch = 16;
% % cur_set = find(all_expt==1);
% % plot(all_lfp_t_axis(cur_set),all_lfps(cur_set,cur_ch));
% % hold on
% % plot(all_lfp_t_axis(cur_set),lin_pred(cur_ch,cur_set),'r');
% 
% 
