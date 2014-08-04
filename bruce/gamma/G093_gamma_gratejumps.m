% close all
clear all


dir_prefix = '~';
% dir_prefix = '/Volumes/james';
Expt_name = 'G093';
data_dir = [dir_prefix '/Data/bruce/' Expt_name];
addpath([dir_prefix '/James_scripts/bruce/temporary_scripts/'])

cd(data_dir);

load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct

%%
dsf = 1; %lfps originally sampled at 400Hz
use_lfps = 1:2:96;
lcf = 0.5;

min_trial_dur = 0.5;
cur_expt_set = [41 43 44 54:57];

%%
fprintf('Computing prep data\n');
trial_cnt = 0;

all_t_axis = [];
all_tsince_start = [];
all_exptvec = [];
all_trialvec = [];
all_trial_durs = [];
all_trial_start_times = [];
all_trial_end_times = [];
all_V = [];
all_trial_exptvec = [];
all_stim_times = [];
all_stim_ph = [];
all_stim_co = [];
for ee = 1:length(cur_expt_set);
    fprintf('Expt %d of %d\n',ee,length(cur_expt_set));
    cur_expt = cur_expt_set(ee);
    
    trial_start_times = [Expts{cur_expt}.Trials(:).TrialStart]/1e4;
    trial_end_times = [Expts{cur_expt}.Trials(:).TrueEnd]/1e4;
    trial_durs = [Expts{cur_expt}.Trials(:).dur]/1e4;
    trial_ids = [Expts{cur_expt}.Trials(:).id];
    
%     [un_ids,id_inds] = unique(trial_ids);
    
    all_trial_start_times = cat(1,all_trial_start_times,trial_start_times');
    all_trial_end_times = cat(1,all_trial_end_times,trial_end_times');
    all_trial_durs = cat(1,all_trial_durs,trial_durs');
    all_trial_exptvec = cat(1,all_trial_exptvec,ee*ones(length(trial_ids),1));
    
    n_trials = length(trial_ids);
    
    for tt = 1:n_trials
        cur_stim_times = [Expts{cur_expt}.Trials(tt).Start]/1e4;
        cur_stim_ph = [Expts{cur_expt}.Trials(tt).ph];
        cur_stim_co = [Expts{cur_expt}.Trials(tt).co];
  
        all_stim_times = [all_stim_times;  cur_stim_times(:)];
        all_stim_ph = [all_stim_ph;  cur_stim_ph(:)];
        all_stim_co = [all_stim_co;  cur_stim_co(:)];
    end
    
    %load lfps
    lfp_fname = sprintf('Expt%d_LFP.mat',cur_expt);
    load(lfp_fname);
    Fs = lfp_params.Fsd;
    cur_lfps = bsxfun(@times,double(lfp_mat(:,use_lfps)),lfp_int2V(use_lfps)');
    if lcf > 0
        [filt_b,filt_a] = butter(2,lcf/(Fs/2),'high');
        for ii = 1:length(use_lfps)
            cur_lfps(:,ii) = filtfilt(filt_b,filt_a,cur_lfps(:,ii));
        end
    end
    if dsf > 1
        [filt_b,filt_a] = butter(4,0.8/dsf,'low');
        for ii = 1:length(use_lfps)
            cur_lfps(:,ii) = filtfilt(filt_b,filt_a,cur_lfps(:,ii));
        end
        cur_lfps = downsample(cur_lfps,dsf);
        cur_lfp_t = downsample(lfp_t_ax,dsf);
    else
        cur_lfp_t = lfp_t_ax;
    end
    all_V = cat(1,all_V,cur_lfps);
    all_t_axis = [all_t_axis; cur_lfp_t'];
    
    expt_tsince_start = nan(length(cur_lfp_t),1);
    expt_exptvec = nan(length(cur_lfp_t),1);
    expt_trialvec = nan(length(cur_lfp_t),1);
    for tt = 1:n_trials
        cur_samples = find(cur_lfp_t >= trial_start_times(tt) & cur_lfp_t <= trial_end_times(tt));
        
        expt_tsince_start(cur_samples) = cur_lfp_t(cur_samples) - trial_start_times(tt);
        expt_exptvec(cur_samples) = ee;
        expt_trialvec(cur_samples) = tt + trial_cnt;
    end
    trial_cnt = trial_cnt + n_trials;
    
    all_exptvec = cat(1,all_exptvec,expt_exptvec);
    all_tsince_start = cat(1,all_tsince_start,expt_tsince_start);
    all_trialvec = cat(1,all_trialvec,expt_trialvec);
end
Fsd = Fs/dsf;


all_stim_co(all_stim_ph == -1009) = 0;
all_stim_ph(all_stim_ph == -1009) = 0;

%%
event_inds = round(interp1(all_t_axis,1:length(all_t_axis),all_stim_times));
event_types = zeros(size(event_inds));
bad_events = find(isnan(event_inds));
event_inds(isnan(event_inds)) = 1; %just to avoid errors

trial_starts = find(all_tsince_start(event_inds) < 0.1);
gray_stims = find(all_stim_co == 0);
cont_dec = 1+find(all_stim_co(2:end) < all_stim_co(1:end-1)); %contrast dec-step
cont_inc = 1+find(all_stim_co(2:end) > all_stim_co(1:end-1)); %contrast inc-step
cont_change = 1+find(all_stim_co(1:end-1) ~= all_stim_co(2:end));
cont_on = 1+find(all_stim_co(2:end) ~= 0 & all_stim_co(1:end-1) == 0); %contrast inc-step
cont_off = 1+find(all_stim_co(2:end) == 0 & all_stim_co(1:end-1) ~= 0); %contrast inc-step
ph_change = 1+find(all_stim_co(2:end) ~= 0 & all_stim_co(1:end-1) ~= 0 & all_stim_ph(1:end-1) ~= all_stim_ph(2:end));
same_ph = 1+find(all_stim_ph(2:end) == all_stim_ph(1:end-1) & all_stim_co(1:end-1) ~= 0 & all_stim_co(2:end) ~= 0); %no phase change
same_co = 1+find(all_stim_co(2:end) == all_stim_co(1:end-1) & all_stim_co(2:end) ~= 0); %no phase change
half_co = find(all_stim_co == 0.5); 
full_co = find(all_stim_co == 1); 
half_co_const = 1+find(all_stim_co(1:end-1) == 0.5 & all_stim_co(2:end) == 0.5); 
full_co_const = 1+find(all_stim_co(1:end-1) == 1 & all_stim_co(2:end) == 1); 
neither_gray = 1 + find(all_stim_co(1:end-1) ~= 0 & all_stim_co(2:end) ~= 0);
either_gray = 1 + find(all_stim_co(1:end-1) == 0 | all_stim_co(2:end) == 0);

clear event_types
event_types{1} = (setdiff(trial_starts,gray_stims)); %trial onsets (no gray stims)
event_types{2} = (setdiff(ph_change,trial_starts)); %phase change (no gray)
event_types{3} = (setdiff(intersect(ph_change,cont_change),trial_starts)); %phase change and cont change
event_types{4} = (setdiff(intersect(ph_change,cont_inc),trial_starts)); %phase change and cont inc
event_types{5} = (setdiff(intersect(ph_change,cont_dec),trial_starts)); %phase change and cont dec
event_types{6} = (setdiff(intersect(ph_change,half_co_const),trial_starts)); 
event_types{7} = (setdiff(intersect(ph_change,full_co_const),trial_starts)); 
event_types{8} = (setdiff(cont_change,[either_gray; ph_change; trial_starts])); %cont change and same ph
event_types{9} = (setdiff(cont_inc,[ph_change; either_gray; trial_starts])); %cont inc and same ph
event_types{10} = (setdiff(cont_dec,[ph_change; either_gray; trial_starts])); %cont dec and same ph
event_types{11} = (setdiff(cont_off,[ph_change; trial_starts])); %cont on
event_types{12} = (setdiff(cont_on,[ph_change; trial_starts])); %cont off
event_types{13} = (setdiff(intersect(same_ph,same_co),trial_starts)); 
event_types{14} = (setdiff(intersect(same_ph,half_co_const),trial_starts)); 
event_types{15} = (setdiff(intersect(same_ph,full_co_const),trial_starts)); 
event_types{16} = (intersect(trial_starts,gray_stims)); %trial onsets (only gray stims)
event_types{17} = setdiff(intersect(cont_on,half_co),trial_starts);
event_types{18} = setdiff(intersect(cont_on,full_co),trial_starts);

event_names{1} = 'trial onset';
event_names{2} = 'phase change';
event_names{3} = 'phase change and contrast change';
event_names{4} = 'phase change and contrast inc';
event_names{5} = 'phase change and contrast dec';
event_names{6} = 'phase change, half contrast const';
event_names{7} = 'phase change, full contrast const';
event_names{8} = 'same phase and contrast change';
event_names{9} = 'same phase and contrast inc';
event_names{10} = 'same phase and contrast dec';
event_names{11} = 'contrast on-to-off';
event_names{12} = 'contrast off-to-on';
event_names{13} = 'constant stim';
event_names{14} = 'constant stim, half contrast';
event_names{15} = 'constant stim, full contrast';
event_names{16} = 'trial start with gray stim';
event_names{17} = 'contrast 0-50';
event_names{18} = 'contrast 0-100';

n_events = length(event_types);
%to check
is_tstart = zeros(length(event_inds),1);
is_tstart(trial_starts) = 1;
pre_posts = [is_tstart(2:end) all_stim_co(2:end) all_stim_ph(2:end) all_stim_co(1:end-1) all_stim_ph(1:end-1)];
[C,IA,IC] = unique(pre_posts,'rows');
% for nn = 2:n_events
%     disp(event_names{nn})
%     cur_set = IC(event_types{nn}-1);
%     cur_set = unique(cur_set);
%     C(cur_set,:)
%     pause
% end
    
%%
%WAVELET SCALES this gives log-freq spacing ~(2-100) hz
nwfreqs = 40;
min_freq = 2; max_freq = 120;
min_scale = 1/max_freq*Fsd;
max_scale = 1/min_freq*Fsd;
scales = logspace(log10(min_scale),log10(max_scale),nwfreqs);
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);

forlag = round(Fsd*0.75);
backlag = round(0.2*Fsd);
inc_prev_trial = 1;
for ll = 1:length(use_lfps);
    fprintf('LFP %d of %d\n',ll,length(use_lfps));
    temp = cwt(all_V(:,ll),scales,'cmor1-1');
    cur_ampgram = abs(temp)';
    ov_avg_ampgram(ll,:) = mean(cur_ampgram);
    ov_std_ampgram(ll,:) = std(cur_ampgram);

    
%     un_wfreqs = linspace(wfreqs(1),wfreqs(end),100);
%     avg_pow_int = interp1(wfreqs,avg_pow,un_wfreqs);
%     
%     powfun = @(a,x)(a(1)*x.^a(2));
%     BETA = nlinfit(un_wfreqs,avg_pow_int,powfun,[max(avg_pow) -2]);
%     powfit = un_wfreqs.^BETA(2)*BETA(1);
%     interp_powfit = interp1(un_wfreqs,powfit,wfreqs);
%     
%     cur_ampgram_norm = bsxfun(@rdivide,cur_ampgram,interp_powfit);
    
    
    %%
    for cc = 1:n_events
        fprintf('Event type %d of %d\n',cc,n_events);
        use_events = event_types{cc};
        use_events(ismember(use_events,bad_events)) = [];
        [cond_trig_spec(ll,cc,:,:),lags] = get_event_trig_avg(cur_ampgram,event_inds(use_events),backlag,forlag,0,all_trialvec,inc_prev_trial);
        [cond_trig_lfp(ll,cc,:),lags] = get_event_trig_avg(all_V(:,ll),event_inds(use_events),backlag,forlag,0,all_trialvec,inc_prev_trial);
    end
end

%%
save_dir = ['~/Analysis/bruce/' Expt_name];
cd(save_dir);
if ~exist('gamma','dir')
    system('mkdir gamma');
end
cd('gamma');
save grateev_cond_trig_avgs use_lfps wfreqs lags Fsd cond_* C ov_* event_names
 
%%
% avg_trig_spec_cond = squeeze(mean(trig_spec_cond));
% avg_trig_sig_cond = squeeze(mean(trig_sig_cond));
% sem_trig_sig_cond = squeeze(std(trig_sig_cond))/sqrt(length(use_lfps));
% %%
% close all
% to_print = 1;
% figname = ['~/Analysis/bruce/G093/grating_etas'];
% for cc = 1:n_events
%     
%     fprintf('Event type %d\n',cc);
%     disp(event_names{cc})
%     subplot(2,1,1)
%     pcolor(lags/Fsd,wfreqs,squeeze(avg_trig_spec_cond(cc,:,:))');shading flat
%     caxis([0.5 1.75])
%     xlim([-0.1 0.6]);
%     xlabel('Time (s)'); ylabel('Frequency (Hz)');
%     title(sprintf('Event type: %s',event_names{cc}));
%     subplot(2,1,2)
%     shadedErrorBar(lags/Fsd,avg_trig_sig_cond(cc,:),sem_trig_sig_cond(cc,:));
%     xlim([-0.1 0.6]);
%     ylim([-12 4]*1e-5);
%     xlabel('Time (s)'); ylabel('Amplitude (V)');
%     xl = xlim();
%     line(xl,[0 0],'color','k','linestyle','--');
%     
%     if to_print == 1
%         fillPage(gcf,'Papersize',[6 8]);
%         if cc == 1
%             print(figname,'-dpsc');
%         else
%             print(figname,'-dpsc','-append');
%         end
%         close;
%     else
%         pause
%         clf
%     end
% end
% 
