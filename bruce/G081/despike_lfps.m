clear all
close all
addpath('~/James_scripts/bruce/G081/');
addpath(genpath('~/James_scripts/despikingtoolbox/'));

cd ~/Data/bruce/G081/
load jbeG081Expts.mat
load ./CellList.mat

Fs = 3e4;
dsf = 30;
Fsd = Fs/dsf;

% target_expts = [8:60];
target_expts = [23 24 25 35 36 37 47 52 53];
target_lfps = 1:96;
%%
for cur_expt = target_expts
    
    trial_start_times = [Expts{cur_expt}.Trials(:).TrialStart]/1e4;
    trial_end_times = [Expts{cur_expt}.Trials(:).TrueEnd]/1e4;
    fname = sprintf('Expt%dClusterTimes.mat',cur_expt);
    load(fname);
    
    %%
    for cur_lfp = target_lfps;
        filename = sprintf('Expt%d.p%dFullV.mat',cur_expt,cur_lfp);
        load(filename);
        V = double(FullV.V);
        V = V*FullV.intscale(1)/FullV.intscale(2);
        nparts = length(FullV.blklen);
        dV = [];
        %splice together multiple blocks
        cur_pt = 1;
        for pp = 1:nparts
            cur_range = cur_pt:(cur_pt + FullV.blklen(pp)-1);
            cur_range(cur_range > length(V)) = [];
            dV = [dV V(cur_range)];
            cur_pt = cur_pt + FullV.blklen(pp);
        end
        V = V';
        
        t_ax = [];
        for pp = 1:nparts
            cur_t_ax = linspace(FullV.blkstart(pp),FullV.blkstart(pp)+FullV.blklen(pp)/Fs,FullV.blklen(pp));
            t_ax = [t_ax cur_t_ax];
        end
        
        %%
        spike_times = Clusters{cur_lfp}.times;
        spike_inds = round(interp1(t_ax,1:length(t_ax),spike_times));
        
        %%
        spk_win = round(Fs*0.003);
        spk_backlag = round(Fs*0.001);
        n_trials = length(trial_start_times);
        despiked_lfp = V;
        in_trial_inds = zeros(size(V));
        for nn = 1:n_trials
            fprintf('Processing Trial %d of %d\n',nn,n_trials);
            
            cur_inds = find(t_ax >= trial_start_times(nn) & t_ax <= trial_end_times(nn));
            trial_len = length(cur_inds);
            if mod(trial_len,2) ~= 0
                cur_inds(end) = [];
            end
            cur_durs(nn) = trial_len;
            
            cur_V = V(cur_inds);
            g = fitLFPpowerSpectrum(cur_V,0.1,250,Fs);
            
            cur_spike_inds = find(ismember(cur_inds,spike_inds));
            cur_spike_inds(cur_spike_inds <= spk_backlag) = [];
            Bs = eye(spk_win);
            S = zeros(length(cur_inds),1);
            S(cur_spike_inds - spk_backlag) = 1;
            
            opts.displaylevel = 0;
            tic;
            results = despikeLFP(cur_V,S,Bs,g,opts);
            toc;
            %    plot(results.phi);
            
            despiked_lfp(cur_inds) = results.z;
            in_trial_inds(cur_inds) = 1;
        end
        
        %%
        t_axd = downsample(t_ax,dsf)';
        despiked_lfp = decimate(despiked_lfp,dsf);
        
        sname = sprintf('Expt%d_p%d_despkLFP',cur_expt,cur_lfp);
        cd /home/james/Data/bruce/G081
        save(sname,'t_axd','despiked_lfp');
        
    end
end