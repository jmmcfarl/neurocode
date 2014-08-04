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
foreback_blocks = [9 16 21 30 35 45 47 51] - 6;
stim_fs = 1e4/117.5;

new_stim_inds = [20:40:300];

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
    
    cur_used = find(trial_completed{bb}==1);
    cur_stim_ids = [];
    cur_stim_start_times = [];
    for i = 1:length(cur_used)
       cur_n_stims(i) = ceil(length(Expts{foreback_blocks(bb)}.Trials(cur_used(i)).Seedseq)); 
       cur_seedseq = Expts{foreback_blocks(bb)}.Trials(cur_used(i)).Seedseq;
       cur_stim_ids = [cur_stim_ids; cur_seedseq(new_stim_inds)];
       cur_stim_start_times = [cur_stim_start_times; trial_start_times{bb}(cur_used(i))+(0:40:280)'/stim_fs];
    end
    
    all_stim_ids = [all_stim_ids; cur_stim_ids];
    all_stim_start_times = [all_stim_start_times; cur_stim_start_times];
    all_stim_expt_ids = [all_stim_expt_ids; ones(length(cur_stim_ids),1)*foreback_blocks(bb)];
end


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
dsf = 120;Fsd = Fs/dsf;
niqf = Fsd/2;
[filt_b,filt_a] = butter(2,[1 30]/niqf);
use_lfps = [1:48:96];
backlag = round(0*Fsd);
forwardlag = round(4*Fsd);
lags = (-backlag:forwardlag);

all_Vmat = [];
all_t_ax = [];
all_stim_start_inds = [];
all_expt_inds = [];
for bb = 1:length(foreback_blocks)
    fprintf('Analyzing block %d of %d\n',bb,length(foreback_blocks));
    % bb = 1;
    
    filename = sprintf('Expt%dFullVmean.mat',foreback_blocks(bb));
    load(filename);
    
    Vmat = [];
    for ll = 1:length(use_lfps)
        fprintf('Electrode %d of %d\n',ll,length(use_lfps));
        filename = sprintf('Expt%d.p%dFullV.mat',foreback_blocks(bb),use_lfps(ll));
        load(filename);
        V = double(FullV.V);
        V = V + sumv*FullV.sumscale;
        V = V*FullV.intscale(1)/FullV.intscale(2);
        nparts = length(FullV.blklen);
        dV = [];
        cur_pt = 1;
        for pp = 1:nparts
            cur_range = cur_pt:(cur_pt + FullV.blklen(pp)-1);
            cur_range(cur_range > length(V)) = [];
            dV = [dV filtfilt(filt_b,filt_a,decimate(V(cur_range),dsf))]; %do some high-pass filtering
            cur_pt = cur_pt + FullV.blklen(pp);
        end
        Vmat(ll,:) = dV;
   end
        
    t_ax = [];
    for pp = 1:nparts
        cur_t_ax = linspace(FullV.blkstart(pp),FullV.blkstart(pp)+FullV.blklen(pp)/Fs,FullV.blklen(pp));
        t_ax = [t_ax downsample(cur_t_ax,dsf)];
    end
    t_ax(size(Vmat,2)+1:end) = [];   
    
    cur_use_trials = find(all_stim_expt_ids==foreback_blocks(bb));    
    use_stim_start_inds = round(interp1(t_ax,1:length(t_ax),all_stim_start_times(cur_use_trials)));
    all_stim_start_inds = [all_stim_start_inds; use_stim_start_inds + length(all_t_ax)];
    
    all_Vmat = [all_Vmat Vmat];
    all_t_ax = [all_t_ax t_ax];
    all_expt_inds = [all_expt_inds ones(1,length(t_ax))*bb];
    
end

%%
backlag = round(Fsd*0.1);
forwardlag = round(Fsd*0.5);
fore_lf = find(all_infore==1 & all_freq_val==1);
fore_hf = find(all_infore==1 & all_freq_val==2);
back_lf = find(all_infore==0 & all_freq_val==1);
back_hf = find(all_infore==0 & all_freq_val==2); 
% % back_lf = find(all_back==1 & all_freq_val==1);
% % back_hf = find(all_back==1 & all_freq_val==2); 
% back_lf = find(abs(all_x0_seq)>1&abs(all_x0_seq)<2 & all_freq_val==1);
% back_hf = find(abs(all_x0_seq)>1&abs(all_x0_seq)<2 & all_freq_val==2);
% back2_lf = find(abs(all_x0_seq)>2& all_freq_val==1);
% back2_hf = find(abs(all_x0_seq)>2& all_freq_val==2);
back2_lf = find(all_noback==1 & all_freq_val==1);
back2_hf = find(all_noback==1 & all_freq_val==2); 
for i = 1:length(use_lfps)
    [fore_lf_mat_Vf,tlags] = get_event_trig_mat(all_Vmat(i,:),all_stim_start_inds(fore_lf),backlag,forwardlag);
    fore_lf_avg(i,:) = mean(fore_lf_mat_Vf);
    fore_lf_sem(i,:) = std(fore_lf_mat_Vf)/sqrt(length(fore_lf));
    [fore_hf_mat_Vf,tlags] = get_event_trig_mat(all_Vmat(i,:),all_stim_start_inds(fore_hf),backlag,forwardlag);
    fore_hf_avg(i,:) = mean(fore_hf_mat_Vf);
    fore_hf_sem(i,:) = std(fore_hf_mat_Vf)/sqrt(length(fore_hf));
    
    [back_lf_mat_Vf,tlags] = get_event_trig_mat(all_Vmat(i,:),all_stim_start_inds(back_lf),backlag,forwardlag);
    back_lf_avg(i,:) = mean(back_lf_mat_Vf);
    back_lf_sem(i,:) = std(back_lf_mat_Vf)/sqrt(length(back_lf));
    [back_hf_mat_Vf,tlags] = get_event_trig_mat(all_Vmat(i,:),all_stim_start_inds(back_hf),backlag,forwardlag);
    back_hf_avg(i,:) = mean(back_hf_mat_Vf);
    back_hf_sem(i,:) = std(back_hf_mat_Vf)/sqrt(length(back_hf));

        
    [back2_lf_mat_Vf,tlags] = get_event_trig_mat(all_Vmat(i,:),all_stim_start_inds(back2_lf),backlag,forwardlag);
    back2_lf_avg(i,:) = mean(back2_lf_mat_Vf);
    back2_lf_sem(i,:) = std(back2_lf_mat_Vf)/sqrt(length(back2_lf));
    [back2_hf_mat_Vf,tlags] = get_event_trig_mat(all_Vmat(i,:),all_stim_start_inds(back2_hf),backlag,forwardlag);
    back2_hf_avg(i,:) = mean(back2_hf_mat_Vf);
    back2_hf_sem(i,:) = std(back2_hf_mat_Vf)/sqrt(length(back2_hf));

    ov_avg_avg(i,:) = get_event_trig_avg(all_Vmat(i,:),all_stim_start_inds,backlag,forwardlag);
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
      shadedErrorBar(tlags/Fsd,back2_lf_avg(i,:),back2_lf_sem(i,:),{'color','k'});
  xlim([-0.1 0.5])
    subplot(2,1,2)
    hold on
%     shadedErrorBar(tlags/Fsd,back2_hf_avg(i,:)-ov_avg_avg(i,:),back2_hf_sem(i,:),{'color','k'});
    shadedErrorBar(tlags/Fsd,fore_hf_avg(i,:),fore_hf_sem(i,:),{'color','b'});
    shadedErrorBar(tlags/Fsd,back_hf_avg(i,:),back_hf_sem(i,:),{'color','r'});
    shadedErrorBar(tlags/Fsd,back2_hf_avg(i,:),back2_hf_sem(i,:),{'color','k'});
    xlim([-0.1 0.5])
    pause
    clf

%     subplot(3,1,1)
%     hold on
%     shadedErrorBar(tlags/Fsd,fore_lf_avg(i,:),fore_lf_sem(i,:),{'color','b'});
%      shadedErrorBar(tlags/Fsd,fore_hf_avg(i,:),fore_hf_sem(i,:),{'color','r'});
%    xlim([-0.1 0.5])
%     subplot(3,1,2)
%     hold on
%     shadedErrorBar(tlags/Fsd,back_lf_avg(i,:),back_lf_sem(i,:),{'color','b'});
%     shadedErrorBar(tlags/Fsd,back_hf_avg(i,:),back_hf_sem(i,:),{'color','r'});
%     xlim([-0.1 0.5])
%     subplot(3,1,3)
%     hold on
%     shadedErrorBar(tlags/Fsd,back2_lf_avg(i,:),back2_lf_sem(i,:),{'color','b'});
%     shadedErrorBar(tlags/Fsd,back2_hf_avg(i,:),back2_hf_sem(i,:),{'color','r'});
% %     shadedErrorBar(tlags/Fsd,back2_hf_avg(i,:),back2_hf_sem(i,:),{'color','k'});
%     xlim([-0.1 0.5])
%     pause
%     clf

end

% figure
% subplot(2,1,1)
% hold on
% plot(tlags/Fsd,mean(fore_lf_avg),'b')
% plot(tlags/Fsd,mean(back_lf_avg),'r')
% plot(tlags/Fsd,mean(ov_avg_avg),'k')
% xlim([-0.1 0.5])
% legend('Foreground','Background','Avg')
% title('Low freq')
% subplot(2,1,2)
% hold on
% plot(tlags/Fsd,mean(fore_hf_avg),'b')
% plot(tlags/Fsd,mean(back_hf_avg),'r')
% plot(tlags/Fsd,mean(ov_avg_avg),'k')
% legend('Foreground','Background','Avg')
% title('High freq')
%     xlim([-0.1 0.5])
% 
%     figure
% subplot(2,1,1)
% hold on
% plot(tlags/Fsd,mean(fore_lf_avg-ov_avg_avg),'b')
% plot(tlags/Fsd,mean(back_lf_avg-ov_avg_avg),'r')
% legend('Foreground','Background')
% title('Low freq')
% xlim([-0.1 0.5])
% subplot(2,1,2)
% hold on
% plot(tlags/Fsd,mean(fore_hf_avg-ov_avg_avg),'b')
% plot(tlags/Fsd,mean(back_hf_avg-ov_avg_avg),'r')
%     xlim([-0.1 0.5])
% legend('Foreground','Background')
% title('High freq')

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
    [ov_mat,tlags] = get_event_trig_mat(all_Vmat(j,:),all_stim_start_inds,backlag,forwardlag);
    ov_avg = mean(ov_mat);
    for i = 1:n_poss_ors
        subplot(2,3,i)
        fprintf('ori %d of %d\n',i,n_poss_ors);
        cur_set = or_stims(all_or_val(or_stims)==poss_ors(i) & all_infore(or_stims)==1);
        cur_or_in(i) = length(cur_set);
        [cur_mat_Vf,tlags] = get_event_trig_mat(all_Vmat(j,:),all_stim_start_inds(cur_set),backlag,forwardlag);
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
        
%         cur_set = or_stims(all_or_val(or_stims)==poss_ors(i) & all_back(or_stims)==1);
        cur_set = or_stims(all_or_val(or_stims)==poss_ors(i) & all_noback(or_stims)==1);
        cur_or_back(i) = length(cur_set);
        [cur_mat_Vf,tlags] = get_event_trig_mat(all_Vmat(j,:),all_stim_start_inds(cur_set),backlag,forwardlag);
        noback_avg = mean(cur_mat_Vf);
        noback_sem = std(cur_mat_Vf)/sqrt(length(cur_set));
%         shadedErrorBar(tlags/Fsd,noback_avg-ov_avg,noback_sem,{'color','k'});
            shadedErrorBar(tlags/Fsd,noback_avg,noback_sem,{'color','k'});
        
%         shadedErrorBar(tlags/Fsd,fore_avg-ov_avg,fore_sem,{'color','b'});
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
