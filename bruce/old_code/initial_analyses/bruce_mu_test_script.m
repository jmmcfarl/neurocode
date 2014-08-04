clear all
close all

cd ~/Data/bruce/2_27_12

load Blocks
% cur_block = 1;

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
    
    %% Cycle through single units and compute saccade and image trig averages
    n_units = length(Blocks{cur_block}.mutimes);
    
    dt = 0.005;
    %create time axis for binning
    bin_axis = [];
    for n = 1:n_trials
        bin_axis = [bin_axis block_times(1,n):dt:block_times(2,n)];
    end
    
    binned_su = zeros(n_units,length(bin_axis));
    for c = 1:n_units
        binned_su(c,:) = histc(Blocks{cur_block}.mutimes{c},bin_axis);
    end
    
    mean_rates(cur_block,:) = mean(binned_su,2);
    std_rates(cur_block,:) = std(binned_su,[],2);
    
    %%
    sactrig_lag = round(1/dt);
    sac_trig_avgs = zeros(n_units,2*sactrig_lag+1);
    sac_axis = (-sactrig_lag:sactrig_lag)*dt;
    for s = 1:length(used_sac_times)
        cur_bin_centers = used_sac_times(s) + sac_axis;
        for c = 1:n_units
            cur_dist = hist(Blocks{cur_block}.mutimes{c},cur_bin_centers);
            cur_dist([1 end]) = 0;
            sac_trig_avgs(c,:) = sac_trig_avgs(c,:) + cur_dist;
        end
    end
    sac_trig_avgs = sac_trig_avgs/length(used_sac_times);
    sac_trig_avgsz(cur_block,:,:) = bsxfun(@minus,sac_trig_avgs,mean_rates(cur_block,:)');
    sac_trig_avgsz(cur_block,:,:) = bsxfun(@rdivide,sac_trig_avgsz(cur_block,:,:),std_rates(cur_block,:));
    
    %%
    imagetrig_lag = round(1/dt);
    im_trig_avgs = zeros(n_units,2*imagetrig_lag+1);
    im_axis = (-imagetrig_lag:imagetrig_lag)*dt;
    for s = 1:length(stim_times)
        cur_bin_centers = stim_times(s) + im_axis;
        for c = 1:n_units
            cur_dist = hist(Blocks{cur_block}.mutimes{c},cur_bin_centers);
            cur_dist([1 end]) = 0;
            im_trig_avgs(c,:) = im_trig_avgs(c,:) + cur_dist;
        end
    end
    im_trig_avgs = im_trig_avgs/length(stim_times);
    im_trig_avgsz(cur_block,:,:) = bsxfun(@minus,im_trig_avgs,mean_rates(cur_block,:)');
    im_trig_avgsz(cur_block,:,:) = bsxfun(@rdivide,im_trig_avgsz(cur_block,:,:),std_rates(cur_block,:));
    
end


%%
figure
cmap = colormap(jet(n_units));
hold on
for i = 1:n_units
    cur_avg = squeeze(mean(im_trig_avgsz(:,i,:)));
    cur_sem = squeeze(std(im_trig_avgsz(:,i,:)))/2;
    shadedErrorBar(im_axis,cur_avg+i/3,cur_sem,{'color',cmap(i,:)})
end
yl = ylim();
line([0 0],yl,'color','k')
% xlim([-0.2 0.5])
xlabel('Relative Time (s)','fontsize',14)
ylabel('Relative Firing rate (z)','fontsize',14)


%%
span = 10;
sac_trig_avgs_sm = zeros(size(sac_trig_avgsz));
for b = 1:4
    for c = 1:n_units
        sac_trig_avgs_sm(b,c,:) = smooth(squeeze(sac_trig_avgsz(b,c,:)),span,'lowess');
    end
end

[a,b] = max(sac_trig_avgs_sm,[],3);
peak_delays = sac_axis(b);
avg_delays = mean(peak_delays);
avg_delay_inds = round(mean(b));

xx = linspace(0,1400,100);
p = polyfit(Blocks{1}.muprobes*50,avg_delays,1);

figure
errorbar(Blocks{1}.muprobes*50,avg_delays,std(peak_delays)/2,'o','linewidth',1)
xlabel('Depth from superficial probe (um)','fontsize',14)
ylabel('Saccade response delay (s)','fontsize',14)
hold on
plot(xx,polyval(p,xx),'k','linewidth',1)

%%
figure
cmap = colormap(jet(n_units));
hold on
for i = 1:n_units
    cur_avg = squeeze(mean(sac_trig_avgsz(:,i,:)));
    cur_sem = squeeze(std(sac_trig_avgsz(:,i,:)))/2;
    shadedErrorBar(sac_axis,cur_avg+i/5,cur_sem,{'color',cmap(i,:)})
    plot(sac_axis(avg_delay_inds(i)),cur_avg(avg_delay_inds(i))+i/5,'ko','markersize',10,'linewidth',2)
end
% ylim([0 2.7])
yl = ylim();
line([0 0],yl,'color','k')
xlim([-0.2 0.5])
xlabel('Relative Time (s)','fontsize',14)
ylabel('Relative Firing rate (z)','fontsize',14)
