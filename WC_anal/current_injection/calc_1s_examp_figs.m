clear all

depol_array{1} = 'C:\WC_Germany\april_09_data\2009-04-07\2009-4-7-19\2009-04-07_CWC_LFP_1s_NEG_CI';
depol_array{2} = 'C:\WC_Germany\april_09_data\2009-04-13_A\2009-4-13-18\2009-04-13_CWC_LFP_A_1s_NEG';

prior_array{1} = 'C:\WC_Germany\april_09_data\2009-04-07\2009-4-7-19\2009-04-07_CWC_LFP_spontaneous';
prior_array{2} = 'C:\WC_Germany\april_09_data\2009-04-13_A\2009-4-13-18\2009-04-13_CWC_LFP_A_spontaneous';


dep_time = 1;
wait_time = 2;
first_pulse = 1;

Fs = 2e4;
dsf = 20;
Fsd = Fs/dsf;

mp_min = -95;
mp_max = 10;
nbins = 500;
mp_range = linspace(mp_min,mp_max,nbins);

d=2

    %% get prior distribution
    load(prior_array{d})

    prior_dist(d,:) = gpkde(data,2,[mp_min;mp_max;nbins]);
    mean_prior(d) = mean(data);
    std_prior(d) = std(data);

    %% get depolarizing CI
    load(depol_array{d})

    data = downsample(data,dsf);
    time = downsample(time,dsf);

    num_cis = floor((max(time)-2*first_pulse)/(dep_time+wait_time))

    sweep_time = (dep_time+wait_time)*Fsd;

    sweep_mat = zeros(num_cis,sweep_time);

    for i = 1:num_cis

        begpt = (i-1)*3*Fsd+1;
        endpt = begpt+3*Fsd;
        sweep_mat(i,:) = data(begpt:endpt-1);

    end

    tvec = (1:sweep_time)/Fsd;

    [pr_peak_amps,pr_peak_locs] = findpeaks(prior_dist(d,:),'minpeakheight',0.01);

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
        n_pr_midpt(d) = length(mp_range);
    end

%    bar(mp_range,prior_dist)
%    line([n_pr_down_mp(d) n_pr_down_mp(d)],[0 0.02],'Color','k')
%    line([n_pr_up_mp(d) n_pr_up_mp(d)],[0 0.02],'Color','k')
%    xlim([-90 10])
%    figure
   
   
   
   for i=1:58
       plot(tvec,sweep_mat(i,:))
       line([0 3],[n_pr_down_mp(d) n_pr_down_mp(d)],'Color','k')
       line([0 3],[n_pr_up_mp(d) n_pr_up_mp(d)],'Color','k')
       ylim([-90 20])
       i
       pause
       clf
   end