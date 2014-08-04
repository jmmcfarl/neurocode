clear all
% close all
cd 
addpath('/Users/James/Data/bruce/7_15_12/')

cd /Users/James/Data/bruce/7_15_12/G029/  

load ./jbeG029.em.mat
em_data = Expt; clear Expt
load ./CellList.mat
load ./G029Expts.mat

%%
Expt_nu = [30 31]; %these are the grating expts
n_allunits = 96;
for ee = 1:length(Expt_nu)
    fprintf('Analyzing Expt %d of %d\n',ee,length(Expt_nu));
    
    single_units{ee} = find(CellList(Expt_nu(ee),:,1) > 0);
    n_sus = length(single_units);
    multi_units{ee} = setdiff(1:n_allunits,single_units{ee});
    
    load(sprintf('Expt%dClusterTimes.mat',Expt_nu(ee)));
    
    Trial_starts = [Expts{Expt_nu(ee)}.Trials(:).Start]/1e4;
    Trial_ends = [Expts{Expt_nu(ee)}.Trials(:).End]/1e4;
    endEvents = [Expts{Expt_nu(ee)}.Trials(:).endevent];
    Trial_durs = (Trial_ends-Trial_starts);
    line_oris = [Expts{Expt_nu(ee)}.Trials(:).or];
    xo = [Expts{Expt_nu(ee)}.Trials(:).xo];
    yo = [Expts{Expt_nu(ee)}.Trials(:).yo];
    se = [Expts{Expt_nu(ee)}.Trials(:).se]; %changes
    rw = [Expts{Expt_nu(ee)}.Trials(:).rw]; %changes
    st = [Expts{Expt_nu(ee)}.Trials(:).st]; %changes
    sO = [Expts{Expt_nu(ee)}.Trials(:).sO]; %changes
    Op = [Expts{Expt_nu(ee)}.Trials(:).Op]; %changes
    ar = [Expts{Expt_nu(ee)}.Trials(:).ar]; %changes
    Pp = [Expts{Expt_nu(ee)}.Trials(:).Pp]; %changes
    completed_trials = find(Trial_durs > 0.4);
    
    %% GET NEEDED EYE TRACKING DATA
    [eye_vals,eye_ts,eye_insample] = get_eye_tracking_data(em_data,[Trial_starts(1) Trial_ends(end)]);
    
    %%
    n_used_trials(ee) = length(completed_trials);
    spk_counts = nan(n_allunits,n_used_trials(ee));
    for su = 1:n_allunits
        %     cur_spk_times = Clusters{single_units(su)}.times;
        cur_spk_times = Clusters{su}.times;
        for i = 1:n_used_trials(ee)
            spk_counts(su,i) = sum(cur_spk_times > Trial_starts(completed_trials(i)) & ...
                cur_spk_times < Trial_ends(completed_trials(i)));
        end
    end
    spk_rates = bsxfun(@rdivide,spk_counts,Trial_durs(completed_trials));
    
    %%
    unique_oris = unique(line_oris);
    unique_ars = unique(ar);
    tot_perms = length(unique_oris)*length(unique_ars);
    n_pos = 15;
    
    avg_profile(ee,:,:,:) = nan(n_allunits,tot_perms,n_pos);
    n_used_stims(ee,:) = nan(1,tot_perms);
    x_vals(ee,:,:) = nan(tot_perms,n_pos);
    y_vals(ee,:,:) = nan(tot_perms,n_pos);
    so_vals(ee,:,:) = nan(tot_perms,n_pos);
    or_vals(ee,:) = nan(1,tot_perms);
    ar_vals(ee,:) = nan(1,tot_perms);
    cnt = 1;
    for i = 1:length(unique_oris)
        for k = 1:length(unique_ars)
            avg_line_oris(cnt) = unique_oris(i);
            avg_bar_ars(cnt) = unique_ars(k);
            cur_set = find(line_oris(completed_trials) == unique_oris(i) & ar(completed_trials) == unique_ars(k));
            csO = sO(completed_trials(cur_set));
            n_used_stims(ee,cnt) = length(cur_set);
            poss_sO = unique(csO);
            for j = 1:n_pos
                cur_pos_set = cur_set(sO(completed_trials(cur_set)) == poss_sO(j));
                avg_profile(ee,:,cnt,j) = mean(spk_rates(:,cur_pos_set),2);
                x_vals(ee,cnt,j) = unique(xo(completed_trials(cur_pos_set)));
                y_vals(ee,cnt,j) = unique(yo(completed_trials(cur_pos_set)));
                so_vals(ee,cnt,j) = unique(sO(completed_trials(cur_pos_set)));
            end
            cnt = cnt + 1;
        end
    end
end

%%
ov_x_vals = squeeze(mean(x_vals));
ov_y_vals = squeeze(mean(y_vals));
ov_so_vals = squeeze(mean(so_vals));
ov_avg_profile = squeeze(mean(avg_profile));
max_profiles = max(ov_avg_profile,[],3);
min_profiles = min(ov_avg_profile,[],3);
avg_bar_oris = avg_line_oris;
unique_bar_oris = unique(avg_bar_oris);
for i = 1:length(unique_bar_oris)
   used_set = find(avg_bar_oris == unique_bar_oris(i)); 
   mod_index(i,:) = (max(max_profiles(:,used_set),[],2) - min(min_profiles(:,used_set),[],2));
end
[bar_ori_maxrates,pref_bar_oris] = max(mod_index);
pref_bar_oris = unique_bar_oris(pref_bar_oris);
for i = 1:96
    cur_set = find(avg_bar_ars < 1 & avg_bar_oris == pref_bar_oris(i));
    cur_par_x = ov_x_vals(cur_set,:);
    cur_par_y = ov_y_vals(cur_set,:);
    cur_par_so = ov_so_vals(cur_set,:);
    parallel_profile(i,:) =  squeeze(ov_avg_profile(i,cur_set,:));
    
    cur_set = find(avg_bar_ars > 1 & avg_bar_oris == pref_bar_oris(i));
    cur_orth_x = ov_x_vals(cur_set,:);
    cur_orth_y = ov_y_vals(cur_set,:);
    cur_orth_so = ov_so_vals(cur_set,:);
    orth_profile(i,:) = squeeze(ov_avg_profile(i,cur_set,:));
    
    par_prof_norm = parallel_profile(i,:)/sum(parallel_profile(i,:));
    orth_prof_norm = orth_profile(i,:)/sum(orth_profile(i,:));
    
    mean_x(i) = 0.5*cur_par_x*par_prof_norm' + 0.5*cur_orth_x*orth_prof_norm';
    mean_y(i) = 0.5*cur_par_y*par_prof_norm' + 0.5*cur_orth_y*orth_prof_norm';
    
    mean_par_s(i) = cur_par_so*par_prof_norm';
    mean_orth_s(i) = cur_orth_so*orth_prof_norm';
    
    std_par(i) = sqrt((cur_par_so - mean_par_s(i)).^2*par_prof_norm');
    std_orth(i) = sqrt((cur_orth_so - mean_orth_s(i)).^2*orth_prof_norm');
    
end

ov_peak_rates = max(max(parallel_profile,[],2),max(orth_profile,[],2));
norm_parallel_profile = bsxfun(@rdivide,parallel_profile,ov_peak_rates);
norm_orth_profile = bsxfun(@rdivide,orth_profile,ov_peak_rates);

%%
used_set = single_units{1};
figure; hold on
cmap = jet(length(used_set));
for i = 1:length(used_set)
    
    ellipse(std_orth(used_set(i)),std_par(used_set(i)),degtorad(pref_bar_oris(used_set(i))),...
        mean_x(used_set(i)),mean_y(used_set(i)),cmap(i,:))
    plot(mean_x(used_set(i)),mean_y(used_set(i)),'o','color',cmap(i,:))
end


