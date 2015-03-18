% close all
clear all

dir_prefix = '~';
% dir_prefix = '/Volumes/james';
Expt_name = 'G096';
data_dir = [dir_prefix '/Data/bruce/' Expt_name];
addpath([dir_prefix '/James_scripts/bruce/temporary_scripts/'])

cd(data_dir);

load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct

for ii = 1:27
    expt_name{ii} = Expts{ii}.Header.expname;
end
%%
dsf = 1; %lfps originally sampled at 400Hz
use_lfps = 1:4:96;
lcf = 0.5;

min_trial_dur = 1;
% cur_expt_set = [6:19];
cur_expt_set = [1:27];
cur_expt_set([1 2 5]) = [];
%%
fprintf('Computing prep data\n');
trial_cnt = 0;

all_t_axis = [];
all_tsince_start = [];
all_exptvec = [];
all_trial_result = [];
all_trialvec = [];
all_trial_durs = [];
all_trial_start_times = [];
all_trial_end_times = [];
all_V = [];
all_trial_exptvec = [];
all_trial_sf = [];
all_trial_tf = [];
all_trial_or = [];
all_trial_sz = [];
all_trial_jv = [];
all_trial_nsf = [];
all_trial_st = [];
all_trial_co = [];
all_trial_ph = [];
all_trial_fh = [];
all_trial_Bc = [];
all_trial_hi = [];
for ee = 1:length(cur_expt_set);
%     fprintf('Expt %d of %d\n',ee,length(cur_expt_set));
    cur_expt = cur_expt_set(ee);
    
    trial_start_times = [Expts{cur_expt}.Trials(:).TrialStart]/1e4;
    trial_end_times = [Expts{cur_expt}.Trials(:).TrueEnd]/1e4;
    trial_durs = [Expts{cur_expt}.Trials(:).dur]/1e4;
    trial_ids = [Expts{cur_expt}.Trials(:).id];
    trial_result = [Expts{cur_expt}.Trials(:).Result];
    
    fprintf('Expt %d, %d Trials\n',ee,length(trial_ids));
    
    if isfield(Expts{cur_expt}.Trials,'sf')
        trial_sf = [Expts{cur_expt}.Trials(:).sf];
    else
%         trial_sf = ones(size(trial_start_times))*Expts{cur_expt}.Stimvals.sf;
        trial_sf = nan(size(trial_start_times));
    end
    if isfield(Expts{cur_expt}.Trials,'tf')
        trial_tf = [Expts{cur_expt}.Trials(:).tf];
    else
%         trial_tf = ones(size(trial_start_times))*Expts{cur_expt}.Stimvals.tf;
        trial_tf = nan(size(trial_start_times));
    end
    if isfield(Expts{cur_expt}.Trials,'or')
        trial_or = [Expts{cur_expt}.Trials(:).or];
    else
        trial_or = ones(size(trial_start_times))*Expts{cur_expt}.Stimvals.or;
    end
    trial_sz = [Expts{cur_expt}.Trials(:).sz];
    if isfield(Expts{cur_expt}.Trials,'jv')
        trial_jv = [Expts{cur_expt}.Trials(:).jv];
    else
%         trial_jv = ones(size(trial_start_times))*Expts{cur_expt}.Stimvals.jv;
        trial_jv = nan(size(trial_start_times));
    end
    if isfield(Expts{cur_expt}.Trials,'nsf')
        trial_nsf = reshape([Expts{cur_expt}.Trials(:).nsf],10,length(trial_start_times));
    else
        fprintf('No nsf in expt %d\n',ee);
        trial_nsf = nan(10,length(trial_jv));
%         trial_nsf = [];
    end
    if isfield(Expts{cur_expt}.Trials,'st')
        trial_st = [Expts{cur_expt}.Trials(:).st];
    else
        trial_st = ones(size(trial_start_times))*Expts{cur_expt}.Stimvals.st;
    end
    if isfield(Expts{cur_expt}.Trials,'co')
        trial_co = [Expts{cur_expt}.Trials(:).co];
    else
%         trial_co = ones(size(trial_start_times))*Expts{cur_expt}.Stimvals.co;
        trial_co = nan(size(trial_start_times));
    end
    if isfield(Expts{cur_expt}.Trials,'ph')
        trial_ph = [Expts{cur_expt}.Trials(:).ph];
    else
%         trial_ph = ones(size(trial_start_times))*Expts{cur_expt}.Stimvals.ph;
        trial_ph = nan(size(trial_start_times));
    end
    if isfield(Expts{cur_expt}.Trials,'fh')
        trial_fh = [Expts{cur_expt}.Trials(:).fh];
    else
%         trial_fh = ones(size(trial_start_times))*Expts{cur_expt}.Stimvals.fh;
        trial_fh = nan(size(trial_start_times));
    end
     if isfield(Expts{cur_expt}.Trials,'Bc')
        trial_Bc = [Expts{cur_expt}.Trials(:).Bc];
    else
%         trial_Bc = ones(size(trial_start_times))*Expts{cur_expt}.Stimvals.Bc;
        trial_Bc = nan(size(trial_start_times));
    end
     if isfield(Expts{cur_expt}.Trials,'hi')
        trial_hi = [Expts{cur_expt}.Trials(:).hi];
    else
%         trial_hi = ones(size(trial_start_times))*Expts{cur_expt}.Stimvals.hi;
        trial_hi = nan(size(trial_start_times));
    end
   
    
%     [un_ids,id_inds] = unique(trial_ids);
    
    all_trial_start_times = cat(1,all_trial_start_times,trial_start_times');
    all_trial_end_times = cat(1,all_trial_end_times,trial_end_times');
    all_trial_durs = cat(1,all_trial_durs,trial_durs');
    all_trial_result = cat(1,all_trial_result,trial_result');
    all_trial_exptvec = cat(1,all_trial_exptvec,ee*ones(length(trial_ids),1));
    
    all_trial_sf = cat(1,all_trial_sf,trial_sf');
    all_trial_tf = cat(1,all_trial_tf,trial_tf');
    all_trial_or = cat(1,all_trial_or,trial_or');
    all_trial_sz = cat(1,all_trial_sz,trial_sz');
    all_trial_jv = cat(1,all_trial_jv,trial_jv');
    all_trial_nsf = cat(1,all_trial_nsf,trial_nsf');
    all_trial_st = cat(1,all_trial_st,trial_st');
    all_trial_co = cat(1,all_trial_co,trial_co');
    all_trial_ph = cat(1,all_trial_ph,trial_ph');
    all_trial_fh = cat(1,all_trial_fh,trial_fh');
    all_trial_Bc = cat(1,all_trial_Bc,trial_Bc');
    all_trial_hi = cat(1,all_trial_hi,trial_hi');
    
    n_trials = length(trial_ids);

    
    %load lfps
    lfp_fname = sprintf('Expt%d_LFP.mat',cur_expt);
    load(lfp_fname);
    Fs = lfp_params.Fsd;
    cur_lfps = bsxfun(@times,double(lfp_mat(:,use_lfps)),lfp_int2V(use_lfps)');
    if lcf > 0
        [filt_b,filt_a] = butter(2,lcf/(Fs/2),'high');
        for ii = 1:length(use_lfps)
            cur_lfps(:,ii) = filtfilt(filt_b,filt_a,cur_lfps(:,ii));
        end
    end
    if dsf > 1
        [filt_b,filt_a] = butter(4,0.8/dsf,'low');
        for ii = 1:length(use_lfps)
            cur_lfps(:,ii) = filtfilt(filt_b,filt_a,cur_lfps(:,ii));
        end
        cur_lfps = downsample(cur_lfps,dsf);
        cur_lfp_t = downsample(lfp_t_ax,dsf);
    else
        cur_lfp_t = lfp_t_ax;
    end
    all_V = cat(1,all_V,cur_lfps);
    all_t_axis = [all_t_axis; cur_lfp_t'];
    
    expt_tsince_start = nan(length(cur_lfp_t),1);
    expt_exptvec = nan(length(cur_lfp_t),1);
    expt_trialvec = nan(length(cur_lfp_t),1);
    for tt = 1:n_trials
        cur_samples = find(cur_lfp_t >= trial_start_times(tt) & cur_lfp_t <= trial_end_times(tt));
        
        expt_tsince_start(cur_samples) = cur_lfp_t(cur_samples) - trial_start_times(tt);
        expt_exptvec(cur_samples) = ee;
        expt_trialvec(cur_samples) = tt + trial_cnt;
    end
    trial_cnt = trial_cnt + n_trials;
    
    all_exptvec = cat(1,all_exptvec,expt_exptvec);
    all_tsince_start = cat(1,all_tsince_start,expt_tsince_start);
    all_trialvec = cat(1,all_trialvec,expt_trialvec);
end
Fsd = Fs/dsf;


%%
addpath(genpath('~/James_scripts/chronux/spectral_analysis/'));
start_buffer = round(Fsd*0.5);
all_trial_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),all_trial_start_times));
all_trial_end_inds = round(interp1(all_t_axis,1:length(all_t_axis),all_trial_end_times));

%% SINGLE GRATINGS
base_gratings = find(all_trial_st == 3 & isnan(all_trial_Bc) & all_trial_sz == 6.0057);
stim_values = [all_trial_tf(base_gratings) all_trial_sf(base_gratings)];
[C,IA,IC] = unique(stim_values,'rows');
n_types = size(C,1);

params.Fs = Fsd;
params.tapers = [2 3];
params.fpass = [0 125];
movingwin = [2 2];
sMarkers = [all_trial_start_inds(base_gratings)+start_buffer all_trial_end_inds(base_gratings)];
used_trial_durs = (sMarkers(:,2) - sMarkers(:,1))/Fsd;
params.err = [0];

clear S
for ii = 1:n_types
    fprintf('Condition %d of %d\n',ii,n_types);
    cur_trials = find(IC == ii & used_trial_durs > movingwin(1));
    
    for ll = 1:length(use_lfps)
        [S(ii,ll,:), f]= mtspectrumc_unequal_length_trials( all_V(:,ll), movingwin, params, sMarkers(cur_trials,:));
    end
        
end

avg_pow = squeeze(mean(log10(S),2));
%%
% close all
% for ii = 1:n_types
%     C(ii,:)
%     plot(f,avg_pow,'k');
%     hold on
%     plot(f,avg_pow(ii,:),'r','linewidth',2);
%     pause
%     clf
% end

%%
poss_tf = [0 1 2 4 8];
poss_sf = [1 2 4];
lwds = [1 2 3];
cmap = jet(length(poss_tf));

cmap(1,:) = [0 0 1];
cmap(end-1,:) = [1 0.5 0];
cmap(end,:) = [1 0 0];

close all
% for ii = 1:length(poss_sf)
for ii = 3
    for jj = 1:length(poss_tf)
        
        cur_type = find(C(:,1) == poss_tf(jj) & C(:,2) == poss_sf(ii));
        plot(f,avg_pow(cur_type,:),'color',cmap(jj,:),'linewidth',lwds(ii));
        hold on
        
    end
    if ii == 1
        legend('0 Hz','1 Hz','2 Hz','4 Hz','8 Hz');
    end
end
axis tight
xlabel('Frequency (Hz)');
ylabel('Log Power');
figufy(gcf);
% cur_fname = '/home/james/Desktop/gamma_sf_tf';
% print(gcf,'-dpdf',cur_fname);

% cur_fname = '/home/james/Desktop/gamma_sf4_tf';
% print(gcf,'-dpdf',cur_fname);

%% GRATINGS VS RLS
base_rls = find(all_trial_st == 15);
stim_values = [all_trial_jv(base_rls) all_trial_tf(base_rls)];
[C_rls,IA,IC] = unique(stim_values,'rows');
n_types_rls = size(C_rls,1);

params.Fs = Fsd;
params.tapers = [3 5];
params.fpass = [0 125];
movingwin = [1.5 1.5];
sMarkers = [all_trial_start_inds(base_rls)+start_buffer all_trial_end_inds(base_rls)];
used_trial_durs = (sMarkers(:,2) - sMarkers(:,1))/Fsd;
params.err = [0];

clear S
for ii = 1:n_types_rls
    fprintf('Condition %d of %d\n',ii,n_types_rls);
    cur_trials = find(IC == ii & used_trial_durs > movingwin(1));
    
    for ll = 1:length(use_lfps)
        [S(ii,ll,:), f]= mtspectrumc_unequal_length_trials( all_V(:,ll), movingwin, params, sMarkers(cur_trials,:));
    end
        
end

avg_pow_rls = squeeze(mean(log10(S),2));
%%
close all
for ii = 1:n_types_rls
    C_rls(ii,:)
    plot(f,avg_pow_rls,'k');
    hold on
    plot(f,avg_pow_rls(ii,:),'r','linewidth',2);
    pause
    clf
end


%% COMPOUND GRATINGS
compound_gratings = find(all_trial_st == 18 & isnan(all_trial_Bc) & all_trial_sz == 6.0057 & ~isnan(all_trial_nsf(:,1)));
stim_values = [all_trial_nsf(compound_gratings,:) all_trial_jv(compound_gratings) all_trial_tf(compound_gratings)];
[C_comp,IA,IC] = unique(stim_values,'rows');
n_types_comp = size(C_comp,1);

params.Fs = Fsd;
params.tapers = [3 5];
params.fpass = [0 125];
movingwin = [1.5 1.5];
sMarkers = [all_trial_start_inds(compound_gratings)+start_buffer all_trial_end_inds(compound_gratings)];
used_trial_durs = (sMarkers(:,2) - sMarkers(:,1))/Fsd;
params.err = [0];

clear S
for ii = 1:n_types_comp
    fprintf('Condition %d of %d\n',ii,n_types_comp);
    cur_trials = find(IC == ii & used_trial_durs > movingwin(1));
    
    for ll = 1:length(use_lfps)
        [S(ii,ll,:), f]= mtspectrumc_unequal_length_trials( all_V(:,ll), movingwin, params, sMarkers(cur_trials,:));
    end
        
end

avg_pow_comp = squeeze(mean(log10(S),2));

% %only have 1 instance of this stim type so eliminate
% avg_pow_comp(1,:) = [];
% n_types_comp = n_types_comp - 1;

%%
jv_two = find(C_comp(:,11) == 2);
jv_one = find(C_comp(:,11) == 1);
jv_p5 = find(C_comp(:,11) == 0.5);
f1 = figure(); hold
plot(f,avg_pow_comp(jv_two(1),:),'k','linewidth',3);
plot(f,avg_pow_comp(jv_two(2),:),'r','linewidth',3);
plot(f,avg_pow_comp(jv_two(3),:),'g','linewidth',3);
plot(f,avg_pow_comp(jv_two(4),:),'b','linewidth',3);
legend('10-grating','1-2-4','2-4','2-4-8');

plot(f,avg_pow_comp(jv_one(1),:),'r','linewidth',2);
plot(f,avg_pow_comp(jv_one(2),:),'g','linewidth',2);
plot(f,avg_pow_comp(jv_one(3),:),'b','linewidth',2);

plot(f,avg_pow_comp(jv_p5(1),:),'r','linewidth',1);
plot(f,avg_pow_comp(jv_p5(2),:),'g','linewidth',1);
plot(f,avg_pow_comp(jv_p5(3),:),'b','linewidth',1);
axis tight
xlabel('Frequency (Hz)');
ylabel('Log Power');
figufy(gcf);
cur_fname = '/home/james/Desktop/gamma_compound_drift';
print(gcf,'-dpdf',cur_fname);

close(gcf);
f1 = figure(); hold
plot(f,avg_pow_comp(jv_two(1),:),'k','linewidth',3);
plot(f,avg_pow_comp(jv_two(2),:),'r','linewidth',3);
plot(f,avg_pow_comp(jv_two(3),:),'g','linewidth',3);
plot(f,avg_pow_comp(jv_two(4),:),'b','linewidth',3);
legend('10-grating','1-2-4','2-4','2-4-8');
axis tight
xlabel('Frequency (Hz)');
ylabel('Log Power');
figufy(gcf);
cur_fname = '/home/james/Desktop/gamma_compound_drift_jv2';
print(gcf,'-dpdf',cur_fname);

close(gcf);
f1 = figure(); hold
plot(f,avg_pow_comp(jv_one(1),:),'r','linewidth',2);
plot(f,avg_pow_comp(jv_one(2),:),'g','linewidth',2);
plot(f,avg_pow_comp(jv_one(3),:),'b','linewidth',2);
legend('1-2-4','2-4','2-4-8');
axis tight
xlabel('Frequency (Hz)');
ylabel('Log Power');
figufy(gcf);
cur_fname = '/home/james/Desktop/gamma_compound_drift_jv1';
print(gcf,'-dpdf',cur_fname);

close(gcf);
f1 = figure(); hold
plot(f,avg_pow_comp(jv_p5(1),:),'r','linewidth',1);
plot(f,avg_pow_comp(jv_p5(2),:),'g','linewidth',1);
plot(f,avg_pow_comp(jv_p5(3),:),'b','linewidth',1);
legend('1-2-4','2-4','2-4-8');
axis tight
xlabel('Frequency (Hz)');
ylabel('Log Power');
figufy(gcf);
cur_fname = '/home/james/Desktop/gamma_compound_drift_jvp5';
print(gcf,'-dpdf',cur_fname);

jv_one = find(C_comp(:,11) == 1);
tf_4 = find(C_comp(:,11) == 0 & C_comp(:,12) == 4);
f2 = figure(); hold 
plot(f,avg_pow_comp(jv_one(1),:),'r','linewidth',2);
plot(f,avg_pow_comp(jv_one(2),:),'g','linewidth',2);
plot(f,avg_pow_comp(jv_one(3),:),'b','linewidth',2);
legend('1-2-4','2-4','2-4-8');

plot(f,avg_pow_comp(tf_4(1),:),'r','linewidth',1);
plot(f,avg_pow_comp(tf_4(2),:),'g','linewidth',1);
plot(f,avg_pow_comp(tf_4(3),:),'b','linewidth',1);
axis tight
xlabel('Frequency (Hz)');
ylabel('Log Power');
figufy(gcf);
cur_fname = '/home/james/Desktop/gamma_compound_tf_jv';
print(gcf,'-dpdf',cur_fname);

%%
poss_tf = [0 1 2 4 8];
poss_sf = [1 2 4];
lwds = [1 2 3];
cmap = jet(length(poss_tf));

cmap(1,:) = [0 0 1];
cmap(end-1,:) = [1 0.5 0];
cmap(end,:) = [1 0 0];

close all
for zz = 1:n_types_comp
    for ii = 1:length(poss_sf)
        for jj = 1:length(poss_tf)
            
            cur_type = find(C(:,1) == poss_tf(jj) & C(:,2) == poss_sf(ii));
            plot(f,avg_pow(cur_type,:),'color',cmap(jj,:),'linewidth',lwds(ii));
            hold on
            
        end
    end
    C_comp(zz,:)
    plot(f,avg_pow_comp(zz,:),'k','linewidth',4);
    pause
    clf
end

%%
close all
for ii = 1:n_types_comp
    C(ii,:)
    plot(f,avg_pow_comp,'k');
    hold on
    plot(f,avg_pow_comp(ii,:),'r','linewidth',2);
    pause
    clf
end

%% DRIFTING COMPOUND GRATINGS 
compound_gratings = find(all_trial_st == 18 & isnan(all_trial_Bc) & all_trial_sz == 6.0057 & ~isnan(all_trial_nsf(:,1)) & all_trial_jv > 0);
stim_values = [all_trial_nsf(compound_gratings,:) all_trial_jv(compound_gratings) all_trial_tf(compound_gratings)];
[C_comp,IA,IC] = unique(stim_values,'rows');
n_types_comp = size(C_comp,1);

params.Fs = Fsd;
params.tapers = [3 5];
params.fpass = [0 125];
movingwin = [1.5 1.5];
sMarkers = [all_trial_start_inds(compound_gratings)+start_buffer all_trial_end_inds(compound_gratings)];
used_trial_durs = (sMarkers(:,2) - sMarkers(:,1))/Fsd;
params.err = [0];

clear S
for ii = 1:n_types_comp
    fprintf('Condition %d of %d\n',ii,n_types_comp);
    cur_trials = find(IC == ii & used_trial_durs > movingwin(1));
    
    for ll = 1:length(use_lfps)
        [S(ii,ll,:), f]= mtspectrumc_unequal_length_trials( all_V(:,ll), movingwin, params, sMarkers(cur_trials,:));
    end
        
end

avg_pow_comp = squeeze(mean(log10(S),2));

%% SIZE/BACKGROUND
size_gratings = find(all_trial_st == 3 & ~isnan(all_trial_Bc) & isnan(all_trial_co));
stim_values = [all_trial_sz(size_gratings,:) all_trial_Bc(size_gratings)];
[C_size,IA,IC] = unique(stim_values,'rows');
n_types_size = size(C_size,1);

params.Fs = Fsd;
params.tapers = [3 5];
params.fpass = [0 125];
movingwin = [1.5 1.5];
sMarkers = [all_trial_start_inds(size_gratings)+start_buffer all_trial_end_inds(size_gratings)];
used_trial_durs = (sMarkers(:,2) - sMarkers(:,1))/Fsd;
params.err = [0];

clear S
for ii = 1:n_types_size
    fprintf('Condition %d of %d\n',ii,n_types_size);
    cur_trials = find(IC == ii & used_trial_durs > movingwin(1));
    
    for ll = 1:length(use_lfps)
        [S(ii,ll,:), f]= mtspectrumc_unequal_length_trials( all_V(:,ll), movingwin, params, sMarkers(cur_trials,:));
    end
        
end

avg_pow_size = squeeze(mean(log10(S),2));

%%
poss_sz = unique(C_size(:,1));
poss_Bc = [0 1];
lwds = [1 2];
cmap = jet(length(poss_sz));

cmap(1,:) = [0 0 1];
cmap(end-1,:) = [1 0.5 0];
cmap(end,:) = [1 0 0];


base_compare = find(C(:,1) == 4 & C(:,2) == 4);
close all
plot(f,avg_pow(base_compare,:),'k','linewidth',2); hold on
for ii = 1:length(poss_Bc)
% for ii = 1
    for jj = 1:length(poss_sz)
        
        cur_type = find(C_size(:,1) == poss_sz(jj) & C_size(:,2) == poss_Bc(ii));
        plot(f,avg_pow_size(cur_type,:),'color',cmap(jj,:),'linewidth',lwds(ii));
        hold on
        
    end
    if ii == 1
       legend('0.5','1','2','4','6'); 
    end
end

axis tight
xlabel('Frequency (Hz)');
ylabel('Log Power');
figufy(gcf);
cur_fname = '/home/james/Desktop/gamma_cent_surround';
print(gcf,'-dpdf',cur_fname);

%% SIZE/BACKGROUND HOLE
size_gratings = find(all_trial_exptvec == 24);
stim_values = [all_trial_sz(size_gratings,:) all_trial_co(size_gratings)];
[C_size,IA,IC] = unique(stim_values,'rows');
n_types_size = size(C_size,1);

params.Fs = Fsd;
params.tapers = [3 5];
params.fpass = [0 125];
movingwin = [1.5 1.5];
sMarkers = [all_trial_start_inds(size_gratings)+start_buffer all_trial_end_inds(size_gratings)];
used_trial_durs = (sMarkers(:,2) - sMarkers(:,1))/Fsd;
params.err = [0];

clear S
for ii = 1:n_types_size
    fprintf('Condition %d of %d\n',ii,n_types_size);
    cur_trials = find(IC == ii & used_trial_durs > movingwin(1));
    
    for ll = 1:length(use_lfps)
        [S(ii,ll,:), f]= mtspectrumc_unequal_length_trials( all_V(:,ll), movingwin, params, sMarkers(cur_trials,:));
    end
        
end

avg_pow_size = squeeze(mean(log10(S),2));

%%
poss_sz = unique(C_size(:,1));
poss_co = [0 1];
lwds = [1 2];
cmap = jet(length(poss_sz));

cmap(1,:) = [0 0 1];
cmap(end-1,:) = [1 0.5 0];
cmap(end,:) = [1 0 0];


base_compare = find(C(:,1) == 4 & C(:,2) == 4);
close all
for ii = 1:length(poss_co)
% for ii = 1
    for jj = 1:length(poss_sz)
        
        cur_type = find(C_size(:,1) == poss_sz(jj) & C_size(:,2) == poss_co(ii));
        plot(f,avg_pow_size(cur_type,:),'color',cmap(jj,:),'linewidth',lwds(ii));
        hold on
        
    end
    if ii == 1
       legend('0.5','1','2','4','6'); 
    end
end
% plot(f,avg_pow(base_compare,:),'k','linewidth',2); hold on

axis tight
xlabel('Frequency (Hz)');
ylabel('Log Power');
figufy(gcf);
cur_fname = '/home/james/Desktop/gamma_cent_hole';
print(gcf,'-dpdf',cur_fname);

%%
close all
for ii = 1:n_types_size
    C_size(ii,:)
    plot(f,avg_pow_size,'k');
    hold on
    plot(f,avg_pow_size(ii,:),'r','linewidth',2);
    pause
    clf
end

%% SIZE/BACKGROUND
% size_gratings = find(all_trial_st == 3 & (~isnan(all_trial_Bc) | ~isnan(all_trial_co)));
% all_trial_Bc(size_gratings(isnan(all_trial_Bc(size_gratings)))) = 1;
% all_trial_co(size_gratings(isnan(all_trial_co(size_gratings)))) = 1;
% 
% stim_values = [all_trial_sz(size_gratings,:) all_trial_Bc(size_gratings) all_trial_co(size_gratings)];
% [C_size,IA,IC] = unique(stim_values,'rows');
% n_types_size = size(C_size,1);
% 
% params.Fs = Fsd;
% params.tapers = [3 5];
% params.fpass = [0 125];
% movingwin = [1.5 1.5];
% sMarkers = [all_trial_start_inds(size_gratings)+start_buffer all_trial_end_inds(size_gratings)];
% used_trial_durs = (sMarkers(:,2) - sMarkers(:,1))/Fsd;
% params.err = [0];
% 
% clear S
% for ii = 1:n_types_size
%     fprintf('Condition %d of %d\n',ii,n_types_size);
%     cur_trials = find(IC == ii & used_trial_durs > movingwin(1));
%     
%     for ll = 1:length(use_lfps)
%         [S(ii,ll,:), f]= mtspectrumc_unequal_length_trials( all_V(:,ll), movingwin, params, sMarkers(cur_trials,:));
%     end
%         
% end
% 
% avg_pow_size = squeeze(mean(log10(S),2));
% 
% %%
% close all
% for ii = 1:n_types_size
%     C_size(ii,:)
%     plot(f,avg_pow_size,'k');
%     hold on
%     plot(f,avg_pow_size(ii,:),'r','linewidth',2);
%     pause
%     clf
% end
