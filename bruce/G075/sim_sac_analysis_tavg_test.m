clear all
% close all
cd ~/Data/bruce/G075/
load jbeG075Expts.mat

%%
sim_sac_blocks = [10 14 18 22 24 28 32 37 40 42 50] - 6;
% sim_sac_blocks = [18 28 40] - 6; %sim sac
% sim_sac_blocks = [10 22 32 42] - 6; %sim sac a
% sim_sac_blocks = [14 24 37 50] - 6; %sim sac b
% sim_sac_blocks = [14 18 24 28 37 40 50] -6; %sim sac and sim sac a

% repeat_inds = [1e3 2e3 3e3 2.01e5 2.02e5 2.03e5 4.01e5 4.02e5 4.03e5];
all_expt_id = [];
for bb = 1:length(sim_sac_blocks)
    n_trials(bb) = length(Expts{sim_sac_blocks(bb)}.Trials);
    trial_start_times{bb} = [Expts{sim_sac_blocks(bb)}.Trials(:).Start]/1e4;
    trial_stop_times{bb} = [Expts{sim_sac_blocks(bb)}.Trials(:).End]/1e4;
    trial_seof{bb} = [Expts{sim_sac_blocks(bb)}.Trials(:).seof];
    trial_completed{bb} = [Expts{sim_sac_blocks(bb)}.Trials(:).Result];
    trial_durs{bb} = trial_stop_times{bb} - trial_start_times{bb};
    all_expt_id = [all_expt_id sim_sac_blocks(bb)*ones(1,n_trials(bb))];
end
all_trial_start_times = cell2mat(trial_start_times);
all_trial_stop_times = cell2mat(trial_stop_times);
all_trial_seof = cell2mat(trial_seof);
all_trial_completed = cell2mat(trial_completed);
all_trial_durs = cell2mat(trial_durs);

%%
stim_fs = 1e4/117.5;
Fs = 3e4;
% dsf = 120;Fsd = Fs/dsf;
dsf = 60;Fsd = Fs/dsf;
niqf = Fsd/2;
[filt_b,filt_a] = butter(2,[1 100]/niqf);
use_lfps = [1:8:96];
backlag = round(0*Fsd);
forwardlag = round(4*Fsd);
lags = (-backlag:forwardlag);

all_Vmat = [];
all_Vmatf = [];
all_t_ax = [];
all_expt_inds = [];
all_trial_start_inds = [];
all_stim_start_inds = [];
all_trial_start_expts = [];
all_stim_start_expts = [];
for bb = 1:length(sim_sac_blocks)
    fprintf('Analyzing block %d of %d\n',bb,length(sim_sac_blocks));
    % bb = 1;
    
    filename = sprintf('Expt%dFullVmean.mat',sim_sac_blocks(bb));
    load(filename);
    
    Vmat = [];
    Vmatf = [];
    for ll = 1:length(use_lfps)
        fprintf('Electrode %d of %d\n',ll,length(use_lfps));
        filename = sprintf('Expt%d.p%dFullV.mat',sim_sac_blocks(bb),use_lfps(ll));
        load(filename);
        V = double(FullV.V);
        V = V + sumv*FullV.sumscale;
        V = V*FullV.intscale(1)/FullV.intscale(2);
        nparts = length(FullV.blklen);
        dV = [];
        dVf = [];
        cur_pt = 1;
        for pp = 1:nparts
            cur_range = cur_pt:(cur_pt + FullV.blklen(pp)-1);
            cur_range(cur_range > length(V)) = [];
            dV = [dV decimate(V(cur_range),dsf)];
            dVf = [dVf filtfilt(filt_b,filt_a,decimate(V(cur_range),dsf))]; %do some high-pass filtering
            cur_pt = cur_pt + FullV.blklen(pp);
        end
        Vmat(ll,:) = dV;
        Vmatf(ll,:) = dVf;
    end
        
    t_ax = [];
    for pp = 1:nparts
        cur_t_ax = linspace(FullV.blkstart(pp),FullV.blkstart(pp)+FullV.blklen(pp)/Fs,FullV.blklen(pp));
        t_ax = [t_ax downsample(cur_t_ax,dsf)];
    end
    t_ax(size(Vmat,2)+1:end) = [];   
    
    cur_use_trials = find(all_expt_id==sim_sac_blocks(bb) & all_trial_completed==1);
    use_trial_start_inds = round(interp1(t_ax,1:length(t_ax),all_trial_start_times(cur_use_trials)));
    use_trial_stop_inds = round(interp1(t_ax,1:length(t_ax),all_trial_stop_times(cur_use_trials)));
    
    cur_trial_start_times = all_trial_start_times(cur_use_trials);
    cur_stim_start_times = bsxfun(@plus,cur_trial_start_times,(0:40:240)'/stim_fs);
    cur_stim_start_times = sort(cur_stim_start_times(:));
    use_stim_start_inds = round(interp1(t_ax,1:length(t_ax),cur_stim_start_times));

    all_trial_start_inds = [all_trial_start_inds use_trial_start_inds + length(all_t_ax)];
    all_trial_start_expts = [all_trial_start_expts; ones(length(use_trial_start_inds),1)*sim_sac_blocks(bb)];    
    all_stim_start_inds = [all_stim_start_inds; use_stim_start_inds + length(all_t_ax)];
    all_stim_start_expts = [all_stim_start_expts; ones(length(use_stim_start_inds),1)*sim_sac_blocks(bb)];    
    
    all_Vmat = [all_Vmat Vmat];
    all_Vmatf = [all_Vmatf Vmatf];
    all_t_ax = [all_t_ax t_ax];
    all_expt_inds = [all_expt_inds ones(1,length(t_ax))*bb];
    
end
bad_starts = find(isnan(all_stim_start_inds));
all_stim_start_inds(bad_starts) = [];
all_stim_start_expts(bad_starts) = [];
bad_starts = find(isnan(all_trial_start_inds));
all_trial_start_inds(bad_starts) = [];
all_trial_start_expts(bad_starts) = [];
%%
% backlag = round(Fsd*0.1);
% forwardlag = round(Fsd*0.5);
% for i = 1:length(use_lfps)
%     [stim_trig_avg_Vf(i,:),tlags] = get_event_trig_avg(all_Vmatf(i,:),all_stim_start_inds,backlag,forwardlag);
%     [stim_trig_avg_V(i,:),tlags] = get_event_trig_avg(all_Vmat(i,:),all_stim_start_inds,backlag,forwardlag);
% end

%%
backlag = round(Fsd*0.1);
forwardlag = round(Fsd*0.5);
sim_sac_ii = [18 28 40] - 6; %sim sac
sim_noise_ii = [10 22 32 42] - 6; %sim sac a
sim_flip_ii = [14 24 37 50] - 6; %sim sac b

sim_sacs = find(ismember(all_stim_start_expts,sim_sac_ii));
sim_noise = find(ismember(all_stim_start_expts,sim_noise_ii));
sim_flips = find(ismember(all_stim_start_expts,sim_flip_ii));
for i = 1:length(use_lfps)
    [stim_trig_avg_Vf(i,:),tlags] = get_event_trig_avg(all_Vmatf(i,:),all_stim_start_inds,backlag,forwardlag);
    [stim_trig_avg_V(i,:),tlags] = get_event_trig_avg(all_Vmat(i,:),all_stim_start_inds,backlag,forwardlag);
    [noise_trig_avg_Vf(i,:),tlags] = get_event_trig_avg(all_Vmatf(i,:),all_stim_start_inds(sim_noise),backlag,forwardlag);
    [noise_trig_avg_V(i,:),tlags] = get_event_trig_avg(all_Vmat(i,:),all_stim_start_inds(sim_noise),backlag,forwardlag);
    [sac_trig_avg_Vf(i,:),tlags] = get_event_trig_avg(all_Vmatf(i,:),all_stim_start_inds(sim_sacs),backlag,forwardlag);
    [sac_trig_avg_V(i,:),tlags] = get_event_trig_avg(all_Vmat(i,:),all_stim_start_inds(sim_sacs),backlag,forwardlag);
    [flip_trig_avg_Vf(i,:),tlags] = get_event_trig_avg(all_Vmatf(i,:),all_stim_start_inds(sim_flips),backlag,forwardlag);
    [flip_trig_avg_V(i,:),tlags] = get_event_trig_avg(all_Vmat(i,:),all_stim_start_inds(sim_flips),backlag,forwardlag);
end

%%
 for i = 1:length(use_lfps)
plot(tlags/Fsd,noise_trig_avg_Vf(i,:))
hold on
plot(tlags/Fsd,sac_trig_avg_Vf(i,:),'r')
plot(tlags/Fsd,flip_trig_avg_Vf(i,:),'k')
pause
clf
end

%%
figure
shadedErrorBar(tlags/Fsd,mean(noise_trig_avg_Vf),std(noise_trig_avg_Vf),{'color','b'});
hold on
shadedErrorBar(tlags/Fsd,mean(sac_trig_avg_Vf),std(sac_trig_avg_Vf),{'color','r'});
shadedErrorBar(tlags/Fsd,mean(flip_trig_avg_Vf),std(flip_trig_avg_Vf),{'color','k'});
xlim([0 0.5])
xlabel('Time (s)','fontsize',16)
ylabel('Amplitude (V)','fontsize',16)