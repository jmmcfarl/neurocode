clear all

hypol_array{1} = 'C:\WC_Germany\april_09_data\2009-04-05_CWC_LFP\2009-4-13-20-41-27_1s_NEG_CI';
hypol_array{2} = 'C:\WC_Germany\april_09_data\2009-04-13_A\2009-4-13-18\2009-04-13_CWC_LFP_A_1s_NEG';

prior_array{1} = 'C:\WC_Germany\april_09_data\2009-04-05_CWC_LFP\2009-4-13-20-41-27_spontaneous';
prior_array{2} = 'C:\WC_Germany\april_09_data\2009-04-13_A\2009-4-13-18\2009-04-13_CWC_LFP_A_spontaneous';


hyp_time = 1;
wait_time = 2;
first_pulse = 1;

Fs = 2e4;
dsf = 20;
Fsd = Fs/dsf;

mp_min = -110;
mp_max = 10;
nbins = 500;
mp_range = linspace(mp_min,mp_max,nbins);

for d = 1:2

    %% get prior distribution
    load(prior_array{d})

    prior_dist(d,:) = gpkde(data,2,[mp_min;mp_max;nbins]);
    mean_prior(d) = mean(data);
    std_prior(d) = std(data);

    %% get hypolarizing CI
    load(hypol_array{d})

    time = (1:length(data))/Fs;
    data = downsample(data,dsf);
    time = downsample(time,dsf);

    num_cis = floor((max(time)-2*first_pulse)/(hyp_time+wait_time))

    sweep_time = (hyp_time+wait_time)*Fsd;

    sweep_mat = zeros(num_cis,sweep_time);

    for i = 1:num_cis

        begpt = (i-1)*3*Fsd+1;
        endpt = begpt+3*Fsd;
        sweep_mat(i,:) = data(begpt:endpt-1);

    end

    tvec = (1:sweep_time)/Fsd;

    thyp_dist = zeros(length(tvec),length(mp_range));

    for t = 1:length(tvec)
        thyp_dist(t,:) = gpkde(sweep_mat(:,t),2,[mp_min;mp_max;nbins]);
    end

    norm_sweep_mat = (sweep_mat - mean_prior(d))/std_prior(d);

    hyp_end = find(tvec > 2.0,1,'first');
    hyp_start = find(tvec > 1.0,1,'first');
    mean_afterhyp_mp(d,:) = mean(mean(sweep_mat(:,hyp_end+100:hyp_end+round(0.5*Fsd))));
    mean_afterhyp_traj(d,:) = mean(norm_sweep_mat);
    mean_afterhyp(d,:) = mean(thyp_dist(hyp_end+100:hyp_end+round(0.5*Fsd),:));
    mean_beforehyp(d,:) = mean(thyp_dist(round(0.5*Fsd):hyp_start-100,:));
    mean_duringhyp(d,:) = mean(thyp_dist(hyp_start:hyp_end,:));

end

used_tvec = tvec(round(0.5*Fsd)+1:end-round(0.5*Fsd));
used_mean_afterhyp_traj = mean_afterhyp_traj(:,round(0.5*Fsd)+1:end-round(0.5*Fsd));

cd C:\WC_Germany\current_injection\pyramidal\L3MEC_Pyramids_CIs_Pulsedata
save april_09_1s_hyp_traj_data used*


%% comparison
for d = 1:2

    [ad_peak_amps,ad_peak_locs] = findpeaks(mean_afterhyp(d,:),'minpeakheight',0.01);
    [pr_peak_amps,pr_peak_locs] = findpeaks(prior_dist(d,:),'minpeakheight',0.01);

    if length(ad_peak_locs) == 2
        n_ad_down_mp(d) = mp_range(ad_peak_locs(1));
        n_ad_up_mp(d) = mp_range(ad_peak_locs(2));
        [dummy,minloc] = min(mean_afterhyp(d,ad_peak_locs(1):ad_peak_locs(2)));
        n_ad_midpt(d) = minloc+ad_peak_locs(1);
    elseif length(ad_peak_locs) > 2
        [dummy, peakorder] = sort(ad_peak_amps,'descend');
        ad_peak_locs(peakorder(3:end)) = [];
        n_ad_down_mp(d) = mp_range(ad_peak_locs(1));
        n_ad_up_mp(d) = mp_range(ad_peak_locs(2));
        [dummy,minloc] = min(mean_afterhyp(d,ad_peak_locs(1):ad_peak_locs(2)));
        n_ad_midpt(d) = minloc+ad_peak_locs(1);
    else
        n_ad_down_mp(d) = mp_range(ad_peak_locs(1));
        n_ad_up_mp(d) = nan;
        n_ad_midpt(d) = length(mp_range);
    end

    if length(pr_peak_locs) == 2
        n_pr_down_mp(d) = mp_range(pr_peak_locs(1));
        n_pr_up_mp(d) = mp_range(pr_peak_locs(2));
        [dummy,minloc] = min(prior_dist(d,pr_peak_locs(1):pr_peak_locs(2)));
        n_pr_midpt(d) = minloc+pr_peak_locs(1);
    elseif length(pr_peak_locs) > 2
        [dummy, peakorder] = sort(pr_peak_amps,'descend');
        pr_peak_locs(peakorder(3:end)) = [];
        n_pr_down_mp(d) = mp_range(pr_peak_locs(1));
        n_pr_up_mp(d) = mp_range(pr_peak_locs(2));
        [dummy,minloc] = min(prior_dist(d,pr_peak_locs(1):pr_peak_locs(2)));
        n_pr_midpt(d) = minloc+pr_peak_locs(1);
    else
        n_pr_down_mp(d) = mp_range(pr_peak_locs(1));
        n_pr_up_mp(d) = nan;
        n_pr_midpt(d) = length(mp_range);
    end

    n_ad_down_fract(d) = trapz(mean_afterhyp(d,1:n_ad_midpt(d)))/trapz(mean_afterhyp(d,:));
    n_pr_down_fract(d) = trapz(prior_dist(d,1:n_pr_midpt(d)))/trapz(prior_dist(d,:));

    d
end

cd C:\WC_Germany\current_injection\pyramidal\L3MEC_Pyramids_CIs_Pulsedata
save april_09_1s_hyp_data *fract *down_mp *up_mp
