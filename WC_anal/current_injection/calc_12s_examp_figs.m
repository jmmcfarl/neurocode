clear all


mp_min = -95;
mp_max = 10;
nbins = 500;
mp_range = linspace(mp_min,mp_max,nbins);

Fs = 2e4;
dsf = 40;
Fsd = Fs/dsf;

minspk = 0.7;
maxspk = 2.5;

niqf = Fs/2;
hcf = 2/niqf;
[b,a] = butter(2,hcf,'low');

%% FOR NEW DATA
stim_dur = 12;
pause_dur = 18;
first_stim = 9;


depol_array{1} = 'C:\WC_Germany\april_09_data\2009-04-05_CWC_LFP\2009-4-13-20-41-27_12s_CI';
depol_array{2} = 'C:\WC_Germany\april_09_data\2009-04-13_A\2009-4-13-18\2009-04-13_CWC_LFP_A_12s_CI';
depol_array{3} = 'C:\WC_Germany\april_09_data\2009-04-13_B\2009-04-13_CWC_LFP_B_12s_CI_smaller';

prior_array{1} = 'C:\WC_Germany\april_09_data\2009-04-05_CWC_LFP\2009-4-13-20-41-27_spontaneous';
prior_array{2} = 'C:\WC_Germany\april_09_data\2009-04-13_A\2009-4-13-18\2009-04-13_CWC_LFP_A_spontaneous';
prior_array{3} = 'C:\WC_Germany\april_09_data\2009-04-13_B\2009-04-13_CWC_LFP_B_spontaneous';
    
d = 1
    
    load(prior_array{d});
    new_prior_dist(d,:) = gpkde(data,2,[mp_min;mp_max;nbins]);
    new_mean_prior(d) = mean(data);
    new_std_prior(d) = std(data);
    
        [pr_peak_amps,pr_peak_locs] = findpeaks(new_prior_dist(d,:),'minpeakheight',0.01);

    if length(pr_peak_locs) == 2
        n_pr_down_mp(d) = mp_range(pr_peak_locs(1));
        n_pr_up_mp(d) = mp_range(pr_peak_locs(2));
    elseif length(pr_peak_locs) > 2
        [dummy, peakorder] = sort(pr_peak_amps,'descend');
        pr_peak_locs(peakorder(3:end)) = [];
        n_pr_down_mp(d) = mp_range(pr_peak_locs(1));
        n_pr_up_mp(d) = mp_range(pr_peak_locs(2));
    else
        n_pr_down_mp(d) = mp_range(pr_peak_locs(1));
        n_pr_up_mp(d) = nan;
    end

    load(depol_array{d});
    diff_data = [0; diff(data)];
    spklocs = find(diff_data > minspk & diff_data < maxspk);
    spklocs = round(spklocs/dsf);
    
    data = filtfilt(b,a,data);
    time = (1:length(data))/Fs;
    data = downsample(data,dsf);
    time = downsample(time,dsf);

    num_cis = floor((max(time)-2*first_stim)/(stim_dur+pause_dur))+1;

    sweep_time = (stim_dur+pause_dur)*Fsd;
    tvec = (1:sweep_time)/Fsd;

    sweep_mat = zeros(num_cis,sweep_time);

    for i = 1:num_cis

        begpt = (i-1)*30*Fsd+1;
        endpt = begpt+30*Fsd;
        if endpt > length(data)
            endpt = length(data);
            cur_length = endpt-begpt+1;
            sweep_mat(i,1:cur_length) = data(begpt:endpt);
            sweep_mat(i,cur_length+1:end) = nan;
        else
            sweep_mat(i,:) = data(begpt:endpt-1);
        end
        
    end

    figure
   bar(mp_range,new_prior_dist(d,:))
   line([n_pr_down_mp(d) n_pr_down_mp(d)],[0 0.02],'Color','k')
   line([n_pr_up_mp(d) n_pr_up_mp(d)],[0 0.02],'Color','k')
   xlim([-90 10])
   
   figure
    i=1
       plot(tvec,sweep_mat(i,:))
       line([0 30],[n_pr_down_mp(d) n_pr_down_mp(d)],'Color','k')
       line([0 30],[n_pr_up_mp(d) n_pr_up_mp(d)],'Color','k')
       loc = 0;
       dur = 30;
used_spikes = find(spklocs/Fsd > loc & spklocs/Fsd < loc+dur);
% plot(t(spike_ids)-66,ones(size(spike_ids))*3,'ro')
line_bottom = -30;
line_top = -27;
for i = 1:length(used_spikes)
   line([tvec(spklocs(used_spikes(i)))-loc tvec(spklocs(used_spikes(i)))-loc],...
       [line_bottom line_top],'Color','k','linewidth',0.5) 
end
ylim([-90 10])
shg
