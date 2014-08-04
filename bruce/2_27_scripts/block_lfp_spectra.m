clear all
close all

cd ~/Data/bruce/2_27_12
    load Blocks

for blockid = 1:5;
    fprintf('Block %d of %d\n',blockid,5);
    load(sprintf('lemM232A.5%d.lfp.mat',blockid));
    
    %%
    %get start times of each LFP trial
    n_lfp_trials = length(LFP.Trials);
    lfp_trial_start = nan(n_lfp_trials,1);
    for i = 1:n_lfp_trials
        lfp_trial_start(i) = LFP.Trials(i).Start/1e4; %originally in tenths of ms
        lfp_trial_stop(i) = LFP.Trials(i).End/1e4;
        lfp_dur(i) = size(LFP.Trials(i).LFP,1)/1000;
    end
    % lfp_trial_stop = lfp_trial_start+lfp_dur';
    
    %get block start times
    block_trial_times = Blocks{blockid}.blocktimes;
    
    %%
    ov_t = block_trial_times(1,1):.001:block_trial_times(2,end);
    
    lfp_time = [];
    lfp_samps = [];
    for i = 1:n_lfp_trials
        lfp_time = [lfp_time linspace(lfp_trial_start(i),lfp_trial_start(i)+lfp_dur(i),size(LFP.Trials(i).LFP,1))];
        lfp_samps = [lfp_samps; LFP.Trials(i).LFP];
    end
    % lfp_samps = zscore(lfp_samps);
    
    Fs = 1000;
    niqf = Fs/2;
    [b,a] = butter(2,[1 200]/niqf);
    lfp_samps = filtfilt(b,a,lfp_samps);
    
    dsf = 2;
    lfp_sampsd = downsample(lfp_samps,dsf);
    Fsd = Fs/dsf;
    
    %%
    win = 20;
    params.Fs = Fsd;
    params.tapers = [4 7];
    n_lfps = size(lfp_samps,2);
    for i = 1:n_lfps
        [S{blockid}(i,:),f]=mtspectrumsegc(lfp_sampsd(:,i),win,params);
    end
    
end
%%
cd ~/James_scripts/bruce/block_lfp_spectra/
cmap = jet(5);
for i = 1:n_lfps
    avg_spec = zeros(5,length(f));
    for j = 1:5
        avg_spec(j,:) = log10(S{j}(i,:));
    end
    for j = 1:5
        cur_diff = bsxfun(@minus,avg_spec(j,:),mean(avg_spec));
%         plot(f,cur_diff,'color',cmap(j,:))
        hold on
        plot(f,smooth(cur_diff,40),'color',cmap(j,:),'linewidth',2)
    end
    set(gca,'xscale','log')
    print('-dpng',sprintf('LFP_%d',i));
    close
end