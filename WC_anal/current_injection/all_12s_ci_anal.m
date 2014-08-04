clear all


mp_min = -1.0;
mp_max = 0.1;
nbins = 500;
mp_range = linspace(mp_min,mp_max,nbins);

Fs = 2e4;
dsf = 40;
Fsd = Fs/dsf;


%% FOR OLD DATA
old_depol_array{1} = 'C:\WC_Germany\longdepol\2007-05-31_CWC_LFP\2007-05-31_CWC_LFP_Hekadata\sweep_data';
old_depol_array{2} = 'C:\WC_Germany\longdepol\2007-06-04_CWC_LFP_B\2007-06-04_CWC_LFP_B_Hekadata\sweep_data';

old_prior_array{1} = 'C:\WC_Germany\EC_MPasci\A2007_05_31_CWC_LFP';
old_prior_array{2} = 'C:\WC_Germany\EC_MPasci\A2007_06_04_CWC_LFP_B';

% for d = 1:2
%     
%     load(old_prior_array{d})
%     data = data_minus_spike{d};
%     data_string = old_prior_array{d}(25:end);
%     data_string = strcat(data_string,'_MP');
%     eval(['cur_data = ' data_string ';'])
%     prior_dist(d,:) = gpkde(cur_data,.02,[mp_min;mp_max;nbins]);
%     mean_prior(d) = mean(cur_data);
%     std_prior(d) = std(cur_data);
%     
%    d
% end

load old_12s_dep_data_minus_spike
for d = 1:2
   
    load(old_depol_array{d})
    sweep_data = data_minus_spike{d};
    
%     n_sweep_data = (sweep_data - mean_prior(d))/std_prior(d);
%     mean_norm_sweep(d,:) = mean(n_sweep_data);
    mean_sweep(d,:) = mean(sweep_data);
    num_cis(d) = size(sweep_data,1);
end


%% FOR NEW DATA
mp_min = -100;
mp_max = 10;
nbins = 500;
mp_range = linspace(mp_min,mp_max,nbins);


stim_dur = 12;
pause_dur = 18;
first_stim = 9;


% depol_array{1} = 'C:\WC_Germany\april_09_data\2009-04-05_CWC_LFP\2009-4-13-20-41-27_12s_CI';
depol_array{1} = 'C:\WC_Germany\april_09_data\2009-04-07\2009-4-7-19\2009-04-07_CWC_LFP_12s_CI';
depol_array{2} = 'C:\WC_Germany\april_09_data\2009-04-13_A\2009-4-13-18\2009-04-13_CWC_LFP_A_12s_CI';
depol_array{3} = 'C:\WC_Germany\april_09_data\2009-04-13_B\2009-04-13_CWC_LFP_B_12s_CI_smaller';

% prior_array{1} = 'C:\WC_Germany\april_09_data\2009-04-05_CWC_LFP\2009-4-13-20-41-27_spontaneous';
prior_array{1} = 'C:\WC_Germany\april_09_data\2009-04-07\2009-4-7-19\2009-04-07_CWC_LFP_spontaneous';
prior_array{2} = 'C:\WC_Germany\april_09_data\2009-04-13_A\2009-4-13-18\2009-04-13_CWC_LFP_A_spontaneous';
prior_array{3} = 'C:\WC_Germany\april_09_data\2009-04-13_B\2009-04-13_CWC_LFP_B_spontaneous';
 load new_12s_dep_data_minus_spike   
% for d = 1:3
%     
%     load(prior_array{d})
%     new_prior_dist(d,:) = gpkde(data,2,[mp_min;mp_max;nbins]);
%     new_mean_prior(d) = mean(data);
%     new_std_prior(d) = std(data);
%     
%    d
% end

for d = 1:3
    
load(depol_array{d})
data = data_minus_spike{d};
    time = (1:length(data))/Fs;
    data = downsample(data,dsf);
    time = downsample(time,dsf);

    num_cis(d+2) = floor((max(time)-2*first_stim)/(stim_dur+pause_dur))+1;

    sweep_time = (stim_dur+pause_dur)*Fsd;

    sweep_mat = zeros(num_cis,sweep_time);

    for i = 1:num_cis(d+2)

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

%     n_sweep_mat = (sweep_mat - new_mean_prior(d))/new_std_prior(d);
%     new_mean_norm_sweep(d,:) = nanmean(n_sweep_mat);
    new_mean_sweep(d,:) = nanmean(sweep_mat);
    
end

%  save all_12s_ci_data *sweep time t_axis