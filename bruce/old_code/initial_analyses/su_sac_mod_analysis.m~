clear all
close all

cd ~/Data/bruce/2_27_12

load Blocks
cur_block = 1;
cur_cell = 10;

ov_sac_times = [];
ov_sac_amps = [];
ov_sac_durs = [];
ov_sac_vel = [];

ov_sac_trig_avg = [];
for cur_block = 1:4
    fprintf('Block %d of %d\n',cur_block,5);
    block_times = Blocks{cur_block}.blocktimes;
    block_durs = diff(block_times);
    n_trials = size(block_times,2);
    
    stim_times = Blocks{cur_block}.stimtime;
    %%
    cd ~/Data/bruce/2_27_12/saccades/
    eval(sprintf('load lemM232.5%d.em.sac.mat',cur_block));
    
    sac_times =  Expt.Trials.EyeMsac.sacT;
    sac_amps = Expt.Trials.EyeMsac.amplitude;
    sac_vels = Expt.Trials.EyeMsac.peakV;
    sac_durs = Expt.Trials.EyeMsac.duration;
    
    %% Find used saccades
    buffer = 0.5;
    
    is_sac_used = zeros(size(sac_times));
    for i = 1:length(sac_times)
        prev_block_start = find(block_times(1,:) < sac_times(i),1,'last');
        if sac_times(i) < block_times(2,prev_block_start) - buffer & ...
                sac_times(i) > block_times(1,prev_block_start) + buffer
            is_sac_used(i) = 1;
        end
    end
    
    used_sac_times = sac_times(is_sac_used==1);
    used_sac_amps = sac_amps(is_sac_used==1);
    used_sac_durs = sac_durs(is_sac_used==1);
    used_sac_vels = sac_vels(is_sac_used==1);
    
    ov_sac_times = [ov_sac_times used_sac_times];
    ov_sac_amps = [ov_sac_amps used_sac_amps];
    ov_sac_durs = [ov_sac_durs used_sac_durs];
    ov_sac_vel = [ov_sac_vel used_sac_vels];
    
    %% Cycle through single units and compute saccade and image trig averages
    n_units = length(Blocks{cur_block}.spktimes);
    
    dt = 0.005;
    %create time axis for binning
    bin_axis = [];
    for n = 1:n_trials
        bin_axis = [bin_axis block_times(1,n):dt:block_times(2,n)];
    end
    
    binned_su = zeros(n_units,length(bin_axis));
    for c = 1:n_units
        binned_su(c,:) = histc(Blocks{cur_block}.spktimes{c},bin_axis);
    end
    
    mean_rates(cur_block,:) = mean(binned_su,2);
    std_rates(cur_block,:) = std(binned_su,[],2);
    
    %%
    
    
    sactrig_lag = round(1/dt);
    sac_trig_avg = zeros(length(used_sac_times),2*sactrig_lag+1);
    sac_axis = (-sactrig_lag:sactrig_lag)*dt;
    for s = 1:length(used_sac_times)
        cur_bin_centers = used_sac_times(s) + sac_axis;
        cur_dist = hist(Blocks{cur_block}.spktimes{cur_cell},cur_bin_centers);
        cur_dist([1 end]) = 0;
        sac_trig_avg(s,:) = cur_dist;
    end
    sac_trig_avg = bsxfun(@minus,sac_trig_avg,mean_rates(cur_block,cur_cell)');
    sac_trig_avg = bsxfun(@rdivide,sac_trig_avg,std_rates(cur_block,cur_cell));
    
    ov_sac_trig_avg = [ov_sac_trig_avg; sac_trig_avg];
end

[~,amp_ord] = sort(ov_sac_amps);
[~,vel_ord] = sort(ov_sac_vel);
[~,dur_ord] = sort(ov_sac_durs);

%%
n_bins = 5;
criterion = ov_sac_amps;
sort_edges = prctile(criterion,linspace(0,100,n_bins+1));
for i = 1:n_bins
    cur_set = find(criterion >= sort_edges(i) & criterion < sort_edges(i+1));
    prc_avg(i,:) = mean(ov_sac_trig_avg(cur_set,:));
    prc_sem(i,:) = std(ov_sac_trig_avg(cur_set,:))/sqrt(length(cur_set));
end

figure; hold on
cmap = colormap(jet(n_bins));
for i = 1:n_bins
    shadedErrorBar(sac_axis,prc_avg(i,:),prc_sem(i,:),{'color',cmap(i,:)},1);
end
xlim([-0.1 0.2])