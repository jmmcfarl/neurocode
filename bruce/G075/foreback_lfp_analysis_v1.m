clear all
% close all
cd ~/Data/bruce/G075/
load jbeG075Expts.mat

fore_x0_val = 0.35;
load ~/Data/bruce/foreground_data_freq_fin.mat
freq_x0_seq = x0_sequence;
freq_y0_seq = y0_sequence;
freq_fore_seq = fore_freq_sequence;

load ~/Data/bruce/foreground_data_fin.mat
or_x0_seq = x0_sequence;
or_y0_seq = y0_sequence;
or_fore_seq = fore_orientation_sequence;
or_back_seq = back_orientation_sequence;
or_infore = (or_x0_seq==fore_x0_val);
or_stim = nan(size(or_x0_seq));
or_stim(or_infore==1) = or_fore_seq(or_infore==1);
or_stim(or_infore==0) = or_back_seq(or_infore==0);
%%
stim_dur = 0.47;
foreback_blocks = [9 16 21 30 35 45 47] - 6;
stim_fs = 1e4/117.5;

new_stim_inds = [20:40:300];
stim_start_inds = [0:40:280];

% repeat_inds = [1e3 2e3 3e3 2.01e5 2.02e5 2.03e5 4.01e5 4.02e5 4.03e5];
all_stim_expt_ids = [];
all_stim_ids = [];
all_stim_start_times = [];
for bb = 1:length(foreback_blocks)
    n_trials(bb) = length(Expts{foreback_blocks(bb)}.Trials);
    trial_start_times{bb} = [Expts{foreback_blocks(bb)}.Trials(:).Start]/1e4;
    trial_stop_times{bb} = [Expts{foreback_blocks(bb)}.Trials(:).End]/1e4;
    trial_completed{bb} = [Expts{foreback_blocks(bb)}.Trials(:).Result];
    trial_durs{bb} = trial_stop_times{bb} - trial_start_times{bb};
    
    %     cur_used = find(trial_completed{bb}==1);
    cur_used_trials = 1:length(trial_completed{bb});
    cur_stim_ids = [];
    cur_stim_start_times = [];
    cur_used_stims = [];
    for i = 1:length(cur_used_trials)
        cur_n_frames(i) = (length(Expts{foreback_blocks(bb)}.Trials(cur_used_trials(i)).Seedseq));
        cur_n_stims(i) = floor(cur_n_frames(i)/40);
        if cur_n_stims(i) > 1 %throw out first stim of each trial
            cur_seedseq = Expts{foreback_blocks(bb)}.Trials(cur_used_trials(i)).Seedseq;
            cur_stim_ids = [cur_stim_ids; cur_seedseq(new_stim_inds(2:cur_n_stims(i)))];
            cur_stim_start_times = [cur_stim_start_times; trial_start_times{bb}(cur_used_trials(i))+(stim_start_inds(2:cur_n_stims(i)))'/stim_fs];
        end
    end
    
    all_stim_ids = [all_stim_ids; cur_stim_ids];
    all_stim_start_times = [all_stim_start_times; cur_stim_start_times];
    all_stim_expt_ids = [all_stim_expt_ids; ones(length(cur_stim_ids),1)*foreback_blocks(bb)];
end
all_stim_stop_times = all_stim_start_times + stim_dur;

all_is_freq = zeros(length(all_stim_ids),1);
all_x0_seq = zeros(length(all_stim_ids),1);
all_freq_val = nan(length(all_stim_ids),1);
all_or_val = nan(length(all_stim_ids),1);

cur_freq_set = find(all_stim_ids > 1000);
freq_ids = all_stim_ids(cur_freq_set) - 1000;
all_is_freq(cur_freq_set) = 1;
all_x0_seq(cur_freq_set) = freq_x0_seq(freq_ids);
cur_fore_set = find(freq_x0_seq(freq_ids)==fore_x0_val);
all_freq_val(cur_freq_set(cur_fore_set)) = freq_fore_seq(freq_ids(cur_fore_set));
cur_back_set = find(freq_x0_seq(freq_ids)~=fore_x0_val);
all_freq_val(cur_freq_set(cur_back_set)) = mod(freq_fore_seq(freq_ids(cur_back_set)),2)+1;

cur_or_set = find(all_stim_ids < 1000);
or_ids = all_stim_ids(cur_or_set);
all_x0_seq(cur_or_set) = or_x0_seq(or_ids);
all_or_val(cur_or_set) = or_stim(or_ids);

all_infore = all_x0_seq==fore_x0_val;
all_noback = isnan(all_x0_seq);
all_back = ~isnan(all_x0_seq) & all_x0_seq~=fore_x0_val;
%%
Fs = 3e4;
dsf = 80;Fsd = Fs/dsf;
% dsf = 60;Fsd = Fs/dsf;
niqf = Fsd/2;
% [filt_b,filt_a] = butter(2,[2 30]/niqf);
[filt_b,filt_a] = butter(2,[1 20]/niqf);
% [gfilt_b,gfilt_a] = butter(2,[35 70]/niqf);
use_lfps = [1:48:96];
backlag = round(0*Fsd);
forwardlag = round(4*Fsd);
lags = (-backlag:forwardlag);

% all_Vmat = [];
all_Vmatf = [];
% all_Vmatg = [];
% all_Vmatga = [];
all_t_ax = [];
all_stim_start_inds = [];
all_expt_inds = [];
for bb = 1:length(foreback_blocks)
    fprintf('Analyzing block %d of %d\n',bb,length(foreback_blocks));
    % bb = 1;
    
    filename = sprintf('Expt%dFullVmean.mat',foreback_blocks(bb));
    load(filename);
    
    Vmatf = [];
    for ll = 1:length(use_lfps)
        fprintf('Electrode %d of %d\n',ll,length(use_lfps));
        filename = sprintf('Expt%d.p%dFullV.mat',foreback_blocks(bb),use_lfps(ll));
        load(filename);
        V = double(FullV.V);
        V = V + sumv*FullV.sumscale;
        V = V*FullV.intscale(1)/FullV.intscale(2);
        nparts = length(FullV.blklen);
        dVf = [];
        cur_pt = 1;
        for pp = 1:nparts
            cur_range = cur_pt:(cur_pt + FullV.blklen(pp)-1);
            cur_range(cur_range > length(V)) = [];
            dVf = [dVf filtfilt(filt_b,filt_a,decimate(V(cur_range),dsf))]; %do some high-pass filtering
            cur_pt = cur_pt + FullV.blklen(pp);
        end
        Vmatf(ll,:) = dVf;
   end
        
    t_ax = [];
    for pp = 1:nparts
        cur_t_ax = linspace(FullV.blkstart(pp),FullV.blkstart(pp)+FullV.blklen(pp)/Fs,FullV.blklen(pp));
        t_ax = [t_ax downsample(cur_t_ax,dsf)];
    end
    t_ax(size(Vmatf,2)+1:end) = [];   
    
    cur_use_trials = find(all_stim_expt_ids==foreback_blocks(bb));    
    use_stim_start_inds = round(interp1(t_ax,1:length(t_ax),all_stim_start_times(cur_use_trials)));
    all_stim_start_inds = [all_stim_start_inds; use_stim_start_inds + length(all_t_ax)];
    
%     all_Vmat = [all_Vmat Vmat];
    all_Vmatf = [all_Vmatf Vmatf];
    all_t_ax = [all_t_ax t_ax];
    all_expt_inds = [all_expt_inds ones(1,length(t_ax))*bb];
    
end

%%
desired_dt = .005;
stim_t_edges = 0:desired_dt:stim_dur;
stim_t_cents = 0.5*stim_t_edges(1:end-1) + 0.5*stim_t_edges(2:end);

all_interp_t = [];
all_interp_stimid = [];
all_interp_exptid = [];
all_interp_lfps = [];
all_interp_spks = [];
for bb = 1:length(foreback_blocks)
    fprintf('Analyzing Expt %d of %d\n',bb,length(foreback_blocks));
    load(sprintf('Expt%dClusterTimes.mat',foreback_blocks(bb)));
    cur_tset = find(all_stim_expt_ids==foreback_blocks(bb));
    
    expt_t_ax = [];
    expt_stim_ids = [];
    expt_binned_spks = [];
    for tt = 1:length(cur_tset)
        
        cur_t_edges = all_stim_start_times(cur_tset(tt)) + stim_t_edges;
        cur_t_cents = all_stim_start_times(cur_tset(tt)) + stim_t_cents;
        
        cur_binned_spks = nan(96,length(cur_t_cents));
        for j = 1:96
            temp = histc(Clusters{j}.times,cur_t_edges);
            cur_binned_spks(j,:) = temp(1:end-1);
        end
        
        expt_t_ax = [expt_t_ax cur_t_cents];
        expt_stim_ids = [expt_stim_ids ones(size(cur_t_cents))*all_stim_ids(cur_tset(tt))];
        expt_binned_spks = [expt_binned_spks; cur_binned_spks'];
    end
    
    all_interp_t = [all_interp_t; expt_t_ax(:)];
    all_interp_stimid = [all_interp_stimid; expt_stim_ids(:)];
    all_interp_exptid = [all_interp_exptid; ones(length(expt_t_ax),1)*foreback_blocks(bb)];
    all_interp_spks = [all_interp_spks; expt_binned_spks];
   
    interp_lfps = interp1(all_t_ax,all_Vmatf',expt_t_ax(:));
    all_interp_lfps = [all_interp_lfps; interp_lfps];
    
end
clear all_sm_spks
for i = 1:96
    all_sm_spks(:,i) = jmm_smooth_1d_cor(all_interp_spks(:,i),1);
end
all_stim_start_inds = round(interp1(all_interp_t,1:length(all_interp_t),all_stim_start_times));
norm_interp_spks = bsxfun(@rdivide,all_interp_spks,mean(all_interp_spks));
norm_sm_spks = bsxfun(@rdivide,all_sm_spks,mean(all_interp_spks));
all_sm_spks = all_sm_spks/desired_dt;
%%
backlag = round(0.1/desired_dt);
forwardlag = round(0.47/desired_dt);
lags = (-backlag:forwardlag)*desired_dt;
or_stims = find(~isnan(all_or_val) & ~isnan(all_stim_start_inds));
poss_ors = unique(all_or_val(or_stims));
freq_stims = find(~isnan(all_freq_val) & ~isnan(all_stim_start_inds));
n_poss_ors = length(poss_ors);
n_poss_freqs = 2;
clear or_trg*
for i = 1:length(use_lfps)
    
    for j = 1:n_poss_ors
        cur_set = or_stims(all_or_val(or_stims)==poss_ors(j) & all_infore(or_stims) == 1);
        [cur_mat_Vf,tlags] = get_event_trig_mat(all_interp_lfps(:,i),all_stim_start_inds(cur_set),backlag,forwardlag);
        or_trg_avg(i,j,:) = mean(cur_mat_Vf);
        or_trg_sem(i,j,:) = std(cur_mat_Vf)/sqrt(length(cur_set));
    end
    for j = 1:2
        cur_set = freq_stims(all_freq_val(freq_stims)==j & all_infore(freq_stims) == 1);
        [cur_mat_Vf,tlags] = get_event_trig_mat(all_interp_lfps(:,i),all_stim_start_inds(cur_set),backlag,forwardlag);
        freq_trg_avg(i,j,:) = mean(cur_mat_Vf);
        freq_trg_sem(i,j,:) = std(cur_mat_Vf)/sqrt(length(cur_set));
        cur_set = freq_stims(all_freq_val(freq_stims)==j & all_back(freq_stims) == 1);
        [cur_mat_Vf,tlags] = get_event_trig_mat(all_interp_lfps(:,i),all_stim_start_inds(cur_set),backlag,forwardlag);
        freq_btrg_avg(i,j,:) = mean(cur_mat_Vf);
        freq_btrg_sem(i,j,:) = std(cur_mat_Vf)/sqrt(length(cur_set));
        cur_set = freq_stims(all_freq_val(freq_stims)==j & all_noback(freq_stims) == 1);
        [cur_mat_Vf,tlags] = get_event_trig_mat(all_interp_lfps(:,i),all_stim_start_inds(cur_set),backlag,forwardlag);
        freq_b2trg_avg(i,j,:) = mean(cur_mat_Vf);
        freq_b2trg_sem(i,j,:) = std(cur_mat_Vf)/sqrt(length(cur_set));
    end
    ov_oravg_avg(i,1,:) = get_event_trig_avg(all_interp_lfps(:,i),all_stim_start_inds(or_stims),backlag,forwardlag);
    ov_favg_avg(i,1,:) = get_event_trig_avg(all_interp_lfps(:,i),all_stim_start_inds(freq_stims),backlag,forwardlag);
end
for i = 1:96
    for j = 1:n_poss_ors
       cur_set = or_stims(all_or_val(or_stims)==poss_ors(j) & all_infore(or_stims) == 1);
       [or_trg_spk_avg(i,j,:),tlags] = get_event_trig_avg(norm_sm_spks(:,i),all_stim_start_inds(cur_set),backlag,forwardlag);
    end
     for j = 1:2
       cur_set = freq_stims(all_freq_val(freq_stims)==j & all_infore(freq_stims) == 1);
       [freq_trg_spk_avg(i,j,:),tlags] = get_event_trig_avg(all_sm_spks(:,i),all_stim_start_inds(cur_set),backlag,forwardlag);
       cur_set = freq_stims(all_freq_val(freq_stims)==j & all_back(freq_stims) == 1);
       [freq_btrg_spk_avg(i,j,:),tlags] = get_event_trig_avg(norm_sm_spks(:,i),all_stim_start_inds(cur_set),backlag,forwardlag);
    end
end

or_trg_avg_ms = bsxfun(@minus,or_trg_avg,ov_oravg_avg);
freq_trg_avg_ms = bsxfun(@minus,freq_trg_avg,ov_favg_avg);
freq_btrg_avg_ms = bsxfun(@minus,freq_btrg_avg,ov_favg_avg);

or_stim_avgs = squeeze(mean(or_trg_spk_avg,3));
for i = 1:96
    [~,or_sort] = sort(or_stim_avgs(i,:));
    sort_or_trg_spk_avg(i,:,:) = or_trg_spk_avg(i,or_sort,:);
end

close all
ll = 1;
% subplot(2,1,1)
cmap = jet(n_poss_ors);
for j = 1:n_poss_ors
    subplot(2,1,1)
    shadedErrorBar(lags,squeeze(or_trg_avg_ms(ll,j,:)),squeeze(or_trg_sem(ll,j,:)),{'color',cmap(j,:)});
    hold on
    subplot(2,1,2)
    shadedErrorBar(lags,squeeze(or_trg_avg(ll,j,:)),squeeze(or_trg_sem(ll,j,:)),{'color',cmap(j,:)});
    hold on
end
% figure
% subplot(2,1,1)
% shadedErrorBar(lags,squeeze(freq_trg_avg_ms(ll,1,:)),squeeze(freq_trg_sem(ll,1,:)),{'color','b'});
% hold on
% shadedErrorBar(lags,squeeze(freq_trg_avg_ms(ll,2,:)),squeeze(freq_trg_sem(ll,2,:)),{'color','r'});
% shadedErrorBar(lags,squeeze(freq_btrg_avg_ms(ll,1,:)),squeeze(freq_btrg_sem(ll,1,:)),{'color','g'});
% shadedErrorBar(lags,squeeze(freq_btrg_avg_ms(ll,2,:)),squeeze(freq_btrg_sem(ll,2,:)),{'color','k'});
% subplot(2,1,2)
% shadedErrorBar(lags,squeeze(freq_trg_avg(ll,1,:)),squeeze(freq_trg_sem(ll,1,:)),{'color','b'});
% hold on
% shadedErrorBar(lags,squeeze(freq_trg_avg(ll,2,:)),squeeze(freq_trg_sem(ll,2,:)),{'color','r'});
% shadedErrorBar(lags,squeeze(freq_btrg_avg(ll,1,:)),squeeze(freq_btrg_sem(ll,1,:)),{'color','g'});
% shadedErrorBar(lags,squeeze(freq_btrg_avg(ll,2,:)),squeeze(freq_btrg_sem(ll,2,:)),{'color','k'});
% % shadedErrorBar(lags,squeeze(freq_b2trg_avg(ll,1,:)),squeeze(freq_b2trg_sem(ll,1,:)),{'color','c'});
% % shadedErrorBar(lags,squeeze(freq_b2trg_avg(ll,2,:)),squeeze(freq_b2trg_sem(ll,2,:)),{'color','y'});

figure
for i = 1:96
    for j = 1:n_poss_ors
        plot(lags,squeeze(or_trg_spk_avg(i,j,:)),'color',cmap(j,:));
        hold on
    end
    i
    pause
    clf
end
figure
for i = 1:96
    plot(lags,squeeze(freq_trg_spk_avg(i,1,:)),'r')
    hold on
    plot(lags,squeeze(freq_trg_spk_avg(i,2,:)),'b')
%      plot(lags,squeeze(freq_btrg_spk_avg(i,1,:)),'r--')
%     plot(lags,squeeze(freq_btrg_spk_avg(i,2,:)),'b--')
   i
    pause
    clf
end

% figure
% shadedErrorBar(lags,squeeze(freq_trg_avg(ll,1,:)),squeeze(freq_trg_sem(ll,1,:)),{'color','r'});
% hold on
% shadedErrorBar(lags,squeeze(freq_trg_avg(ll,2,:)),squeeze(freq_trg_sem(ll,2,:)),{'color','b'});

%%



backlag = round(Fsd*0.1);
forwardlag = round(Fsd*0.5);
fore_lf = find(all_infore==1 & all_freq_val==1);
fore_hf = find(all_infore==1 & all_freq_val==2);
% back_lf = find(all_infore==0 & all_freq_val==1);
% back_hf = find(all_infore==0 & all_freq_val==2); 
% % back_lf = find(all_back==1 & all_freq_val==1);
% % back_hf = find(all_back==1 & all_freq_val==2); 
% back_lf = find(abs(all_x0_seq)>1&abs(all_x0_seq)<2 & all_freq_val==1);
% back_hf = find(abs(all_x0_seq)>1&abs(all_x0_seq)<2 & all_freq_val==2);
% back2_lf = find(abs(all_x0_seq)>2& all_freq_val==1);
% back2_hf = find(abs(all_x0_seq)>2& all_freq_val==2);
back2_lf = find(all_noback==1 & all_freq_val==1);
back2_hf = find(all_noback==1 & all_freq_val==2); 
for i = 1:length(use_lfps)
    [fore_lf_mat_Vf,tlags] = get_event_trig_mat(all_Vmatf(i,:),all_stim_start_inds(fore_lf),backlag,forwardlag);
    fore_lf_avg(i,:) = mean(fore_lf_mat_Vf);
    fore_lf_sem(i,:) = std(fore_lf_mat_Vf)/sqrt(length(fore_lf));
    [fore_hf_mat_Vf,tlags] = get_event_trig_mat(all_Vmatf(i,:),all_stim_start_inds(fore_hf),backlag,forwardlag);
    fore_hf_avg(i,:) = mean(fore_hf_mat_Vf);
    fore_hf_sem(i,:) = std(fore_hf_mat_Vf)/sqrt(length(fore_hf));
    
    [back_lf_mat_Vf,tlags] = get_event_trig_mat(all_Vmatf(i,:),all_stim_start_inds(back_lf),backlag,forwardlag);
    back_lf_avg(i,:) = mean(back_lf_mat_Vf);
    back_lf_sem(i,:) = std(back_lf_mat_Vf)/sqrt(length(back_lf));
    [back_hf_mat_Vf,tlags] = get_event_trig_mat(all_Vmatf(i,:),all_stim_start_inds(back_hf),backlag,forwardlag);
    back_hf_avg(i,:) = mean(back_hf_mat_Vf);
    back_hf_sem(i,:) = std(back_hf_mat_Vf)/sqrt(length(back_hf));

        
    [back2_lf_mat_Vf,tlags] = get_event_trig_mat(all_Vmatf(i,:),all_stim_start_inds(back2_lf),backlag,forwardlag);
    back2_lf_avg(i,:) = mean(back2_lf_mat_Vf);
    back2_lf_sem(i,:) = std(back2_lf_mat_Vf)/sqrt(length(back2_lf));
    [back2_hf_mat_Vf,tlags] = get_event_trig_mat(all_Vmatf(i,:),all_stim_start_inds(back2_hf),backlag,forwardlag);
    back2_hf_avg(i,:) = mean(back2_hf_mat_Vf);
    back2_hf_sem(i,:) = std(back2_hf_mat_Vf)/sqrt(length(back2_hf));

    ov_avg_avg(i,:) = get_event_trig_avg(all_Vmatf(i,:),all_stim_start_inds,backlag,forwardlag);
end
close all
figure
for i = 1:length(use_lfps)
    subplot(2,1,1)
    hold on
%     shadedErrorBar(tlags/Fsd,fore_lf_avg(i,:)-ov_avg_avg(i,:),fore_lf_sem(i,:),{'color','b'});
%     hold on
%     shadedErrorBar(tlags/Fsd,fore_hf_avg(i,:)-ov_avg_avg(i,:),fore_hf_sem(i,:),{'color','r'});
%      shadedErrorBar(tlags/Fsd,back2_lf_avg(i,:)-ov_avg_avg(i,:),back2_lf_sem(i,:),{'color','k'});
    shadedErrorBar(tlags/Fsd,fore_lf_avg(i,:),fore_lf_sem(i,:),{'color','b'});
     shadedErrorBar(tlags/Fsd,back_lf_avg(i,:),back_lf_sem(i,:),{'color','r'});
   xlim([-0.1 0.5])
    subplot(2,1,2)
    hold on
%     shadedErrorBar(tlags/Fsd,back2_hf_avg(i,:)-ov_avg_avg(i,:),back2_hf_sem(i,:),{'color','k'});
    shadedErrorBar(tlags/Fsd,fore_hf_avg(i,:),fore_hf_sem(i,:),{'color','b'});
    shadedErrorBar(tlags/Fsd,back_hf_avg(i,:),back_hf_sem(i,:),{'color','r'});
    xlim([-0.1 0.5])
    pause
    clf
end

figure
subplot(2,1,1)
hold on
plot(tlags/Fsd,mean(fore_lf_avg),'b')
plot(tlags/Fsd,mean(back_lf_avg),'r')
plot(tlags/Fsd,mean(ov_avg_avg),'k')
xlim([-0.1 0.5])
legend('Foreground','Background','Avg')
title('Low freq')
subplot(2,1,2)
hold on
plot(tlags/Fsd,mean(fore_hf_avg),'b')
plot(tlags/Fsd,mean(back_hf_avg),'r')
plot(tlags/Fsd,mean(ov_avg_avg),'k')
legend('Foreground','Background','Avg')
title('High freq')
    xlim([-0.1 0.5])

    figure
subplot(2,1,1)
hold on
plot(tlags/Fsd,mean(fore_lf_avg-ov_avg_avg),'b')
plot(tlags/Fsd,mean(back_lf_avg-ov_avg_avg),'r')
legend('Foreground','Background')
title('Low freq')
xlim([-0.1 0.5])
subplot(2,1,2)
hold on
plot(tlags/Fsd,mean(fore_hf_avg-ov_avg_avg),'b')
plot(tlags/Fsd,mean(back_hf_avg-ov_avg_avg),'r')
    xlim([-0.1 0.5])
legend('Foreground','Background')
title('High freq')

%%
backlag = round(Fsd*0.1);
forwardlag = round(Fsd*0.5);
or_stims = find(~isnan(all_or_val));
poss_ors = unique(all_or_val(or_stims));
n_poss_ors = length(poss_ors);
cmap = jet(n_poss_ors);
use_or = 2;
for i = 1:length(use_lfps)
    
%     for j = 1:n_poss_ors
%        cur_set = or_stims(all_or_val(or_stims)==poss_ors(j));
%        [cur_mat_Vf,tlags] = get_event_trig_mat(all_Vmatg(i,:),all_stim_start_inds(cur_set),backlag,forwardlag);
%        or_avg(j,:) = mean(cur_mat_Vf);
%        or_sem(j,:) = std(cur_mat_Vf)/sqrt(length(cur_set));
%        shadedErrorBar(tlags/Fsd,or_avg(j,:),or_sem(j,:),{'color',cmap(j,:)});
%        hold on
%     end
    

    fore_set = or_stims(all_infore(or_stims)==1 & all_or_val(or_stims)==poss_ors(use_or));
    [cur_mat_Vf,tlags] = get_event_trig_mat(all_Vmatg(i,:),all_stim_start_inds(fore_set),backlag,forwardlag);
    fore_avg = mean(cur_mat_Vf);
    fore_sem = std(cur_mat_Vf)/sqrt(length(fore_set));
    back_set = or_stims(all_back(or_stims)==1 & all_or_val(or_stims)==poss_ors(use_or));
    [cur_mat_Vf,tlags] = get_event_trig_mat(all_Vmatg(i,:),all_stim_start_inds(back_set),backlag,forwardlag);
    back_avg = mean(cur_mat_Vf);
    back_sem = std(cur_mat_Vf)/sqrt(length(back_set));
    nback_set = or_stims(all_noback(or_stims)==1 & all_or_val(or_stims)==poss_ors(use_or));
    [cur_mat_Vf,tlags] = get_event_trig_mat(all_Vmatg(i,:),all_stim_start_inds(nback_set),backlag,forwardlag);
    nback_avg = mean(cur_mat_Vf);
    nback_sem = std(cur_mat_Vf)/sqrt(length(nback_set));
    
    shadedErrorBar(tlags/Fsd,fore_avg,fore_sem,{'color','b'});
    hold on
    shadedErrorBar(tlags/Fsd,back_avg,back_sem,{'color','r'});
%     shadedErrorBar(tlags/Fsd,nback_avg,nback_sem,{'color','k'});

    pause
    clf
end

%%
backlag = round(Fsd*0.1);
forwardlag = round(Fsd*0.5);
or_stims = find(~isnan(all_or_val));
poss_ors = unique(all_or_val(or_stims));
n_poss_ors = length(poss_ors);
cmap = jet(n_poss_ors);

close all
for j = 1:length(use_lfps)
     [ov_mat,tlags] = get_event_trig_mat(all_Vmatf(j,:),all_stim_start_inds,backlag,forwardlag);
     ov_avg = mean(ov_mat);
for i = 1:n_poss_ors
    subplot(2,3,i)
fprintf('ori %d of %d\n',i,n_poss_ors);
    cur_set = or_stims(all_or_val(or_stims)==poss_ors(i) & all_infore(or_stims)==1);
    cur_or_in(i) = length(cur_set);
    [cur_mat_Vf,tlags] = get_event_trig_mat(all_Vmatf(j,:),all_stim_start_inds(cur_set),backlag,forwardlag);
    fore_avg = mean(cur_mat_Vf);
    fore_sem = std(cur_mat_Vf)/sqrt(length(cur_set));
%     shadedErrorBar(tlags/Fsd,fore_avg-ov_avg,fore_sem,{'color',cmap(i,:)});
%     hold on

    hold on

%     cur_set = or_stims(all_or_val(or_stims)==poss_ors(i) & all_back(or_stims)==1);
%     [cur_mat_Vf,tlags] = get_event_trig_mat(all_Vmatf(j,:),all_stim_start_inds(cur_set),backlag,forwardlag);
%     back_avg = mean(cur_mat_Vf);
%     back_sem = std(cur_mat_Vf)/sqrt(length(cur_set));
%     
%     shadedErrorBar(tlags/Fsd,back_avg-ov_avg,back_sem,{'color','r'});
    
    cur_set = or_stims(all_or_val(or_stims)==poss_ors(i) & all_back(or_stims)==1);
    cur_or_back(i) = length(cur_set);
    [cur_mat_Vf,tlags] = get_event_trig_mat(all_Vmatf(j,:),all_stim_start_inds(cur_set),backlag,forwardlag);
    noback_avg = mean(cur_mat_Vf);
    noback_sem = std(cur_mat_Vf)/sqrt(length(cur_set));
%     shadedErrorBar(tlags/Fsd,noback_avg-ov_avg,noback_sem,{'color','k'});
    shadedErrorBar(tlags/Fsd,noback_avg,noback_sem,{'color','k'});
    
% shadedErrorBar(tlags/Fsd,fore_avg-ov_avg,fore_sem,{'color','b'});
shadedErrorBar(tlags/Fsd,fore_avg,fore_sem,{'color','b'});
    xlim([-0.1 0.5])
%     pause
%     clf
end
pause
clf
end
% 
% figure
% for i = 1:length(use_lfps)
%     shadedErrorBar(tlags/Fsd,fore_lf_avg(i,:),fore_lf_sem(i,:),{'color','b'});
%     hold on
%     shadedErrorBar(tlags/Fsd,fore_hf_avg(i,:),fore_hf_sem(i,:),{'color','r'});
%     shadedErrorBar(tlags/Fsd,back_lf_avg(i,:),back_lf_sem(i,:),{'color','k'});
%     shadedErrorBar(tlags/Fsd,back_hf_avg(i,:),back_hf_sem(i,:),{'color','g'});
%     pause
%     clf
% end
% 
