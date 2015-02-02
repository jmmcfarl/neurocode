clear all

data_dir = '~/Data/Sam/MEA_01_28_15/';
cd(data_dir);

%%
load spikeSorting.mat
load Images1_visStimInfo.mat

all_spike_times = t;
clear t
%%
n_units = size(ic,2);
all_spike_cellids = nan(size(all_spike_times));
for cc = 1:n_units
   cur_ind_range = ic(3,cc):ic(4,cc);
   all_spike_cellids(cur_ind_range) = cc;
   spike_times{cc} = all_spike_times(cur_ind_range);
end

nStims = 2500;
trial_bin_edges = vs.stimTimes(1:nStims);
[trial_totspks,all_spike_trialids] = histc(all_spike_times,trial_bin_edges);

%%
count_win = [0 3];
trial_stimids = vs.imgSequence(vs.order);
nImages = length(unique(trial_stimids));
nRpts = 50;

spike_count_mat = nan(n_units,nRpts,nImages);
for ii = 1:nImages
    fprintf('Binning spikes for image number %d/%d\n',ii,nImages);
    cur_trial_set = find(trial_stimids==ii);
    for jj = 1:length(cur_trial_set)
        cur_spike_set = find(all_spike_trialids == cur_trial_set(jj));
        cur_spike_reltimes = (all_spike_times(cur_spike_set) - vs.stimTimes(cur_trial_set(jj)))/1e3;
        cur_spike_cellids = all_spike_cellids(cur_spike_set);
        
        inWindow = cur_spike_reltimes >= count_win(1) & cur_spike_reltimes <= count_win(2);
        spike_count_mat(:,jj,ii) = hist(cur_spike_cellids(inWindow),1:n_units);
    end    
end
avg_spk_counts = squeeze(nanmean(spike_count_mat,2));
norm_spk_counts = nanzscore(avg_spk_counts')';
image_norm_spkcounts = nanmean(norm_spk_counts);
[~,best_image] = nanmax(image_norm_spkcounts);
neuron_avg_spk_counts = mean(avg_spk_counts,2);
%%
% close all
% line_height = 0.9;
% line_width = 2;
% xl = [0 2];
% f1 = figure();
% 
% for trialid = 1:2500;
%     
%     cur_spike_set = find(all_spike_trialids == trialid);
%     cur_spike_reltimes = (all_spike_times(cur_spike_set) - vs.stimTimes(trialid))/1e3;
%     cur_spike_cellids = all_spike_cellids(cur_spike_set);
%     for ii = 1:length(cur_spike_set)
%         line(cur_spike_reltimes(ii) + [0 0],cur_spike_cellids(ii)+[0 line_height],'color','k','linewidth',line_width);
%     end
%     xlim(xl);
%     pause
%     clf
% end
% 
%%
close all
fig_dir = '/home/james/Desktop/SamFigs/';

min_rate = 0.2;
use_units = find(neuron_avg_spk_counts >= min_rate);
n_use_units = length(use_units);

line_height = 1;
line_width = 0.5;
line_length = 0.01;
xl = [0 2];
trial_offset = 0.1;
to_print = true;

fig_width = 3; rel_height = 1.5;
n_use_images = 5;
f1 = figure();

cur_trial_set = find(trial_stimids == best_image);
[~,trialord] = sort(trial_totspks(cur_trial_set),'descend');
trialord = trialord(1:n_use_images);
for ii = 1:length(trialord)
    trialid = cur_trial_set(trialord(ii));
    
    cur_spike_set = find(all_spike_trialids == trialid);
    cur_spike_cellids = all_spike_cellids(cur_spike_set);
    cur_spike_set = cur_spike_set(ismember(cur_spike_cellids,use_units));
    cur_spike_cellids = all_spike_cellids(cur_spike_set);
    cur_spike_reltimes = (all_spike_times(cur_spike_set) - vs.stimTimes(trialid))/1e3;
    for jj = 1:length(cur_spike_set)
        line(cur_spike_reltimes(jj) + [0 line_length]-trial_offset,cur_spike_cellids(jj)+[0 0],'color','r','linewidth',line_width);
    end
%     for jj = 1:n_use_units
%         line(cur_spike_reltimes(jj) + [0 0]-trial_offset,jj+[0 line_height],'color','r','linewidth',line_width);
%     end
    xlim(xl);
    ylim([0 n_use_units]);
    
    if to_print
        fig_name = [fig_dir sprintf('Allcellraster_trial%d',ii)];
        figufy(f1);
        exportfig(f1,fig_name,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
    clf
    else
    pause
    clf
    end
end

%%
trial_dur = 5;
dt = 0.025;
trial_tbin_edges = (0:dt:trial_dur)*1e3;
trial_tbin_cents = 0.5*trial_tbin_edges(1:end-1) + 0.5*trial_tbin_edges(2:end);

n_trials = length(trialStartT);
tbt_spikes = cell(n_trials,n_units);
psths = zeros(n_units,length(trial_tbin_cents));
for ii = 1:n_trials
    for cc = 1:n_units
       cur_trial_spikes = spike_times{cc}(spike_times{cc} >= trialStartTMS(ii) & ...
           spike_times{cc} <= trialEndTMS(ii));
       tbt_spikes{ii,cc} = cur_trial_spikes - trialStartTMS(ii);
       
       cur_hist = histc(tbt_spikes{ii,cc},trial_tbin_edges);
       psths(cc,:) = psths(cc,:) + cur_hist(1:end-1);
    end
end
psths = psths/n_trials;

trial_spike_counts = cellfun(@(x) length(x),tbt_spikes);
ov_trial_avgs = mean(trial_spike_counts);

%%
sm_sig = round(0.1/dt);
norm_psths = bsxfun(@rdivide,psths/dt,ov_trial_avgs');
for ii = 1:n_units
   norm_psths(ii,:) = jmm_smooth_1d_cor(norm_psths(ii,:),sm_sig); 
end

uset = find(ov_trial_avgs >= 0.01);
f1 = figure();
shadedErrorBar(trial_tbin_cents/1e3,mean(norm_psths(uset,:)),std(norm_psths(uset,:))/sqrt(length(uset)));
xlabel('Time (s)');
ylabel('Rate (Hz)');

%%
un_oris = unique(angleOrder);
ori_rates = nan(length(un_oris),n_units);
for ii = 1:length(un_oris)
   cur_trials = find(angleOrder == un_oris(ii));
   ori_rates(ii,:) = mean(trial_spike_counts(cur_trials,:));
end
ori_nrates = bsxfun(@rdivide,ori_rates,ov_trial_avgs);

%%
un_sfs = unique(spatialFreq);
sf_rates = nan(length(un_sfs),n_units);
for ii = 1:length(un_sfs)
   cur_trials = find(spatialFreq == un_sfs(ii));
   sf_rates(ii,:) = mean(trial_spike_counts(cur_trials,:));
end
sf_nrates = bsxfun(@rdivide,sf_rates,ov_trial_avgs);

%%
un_tfs = unique(tempFreq);
tf_rates = nan(length(un_tfs),n_units);
for ii = 1:length(un_tfs)
   cur_trials = find(tempFreq == un_tfs(ii));
   tf_rates(ii,:) = mean(trial_spike_counts(cur_trials,:));
end
tf_nrates = bsxfun(@rdivide,tf_rates,ov_trial_avgs);

%% LOAD LFP DATA 
use_chs = 1:1:252;
lfp_dsf = 6;
load('hdr.mat')
LFP_Fs = hdr.filteredFS;
% LFP_Fs = 1e3;
LFP_Fsd = LFP_Fs/lfp_dsf;
[bb,aa] = butter(2,[0.1 0.8*LFP_Fsd]/LFP_Fs);

lfp_files = dir('chunk*.mat');
n_lfp_files = length(lfp_files);
n_use_files = 4;
all_lfps = [];
for ii = 1:n_use_files
    cur_fname = sprintf('chunk%d.mat',ii);
    fprintf('Loading %s\n',cur_fname);
    load(cur_fname);
    lfp = squeeze(lfp);
    cur_lfps = double(lfp(use_chs,:)');
    cur_lfps = filtfilt(bb,aa,cur_lfps);
    all_lfps = cat(1,all_lfps,downsample(cur_lfps,lfp_dsf));
end
all_lfps = zscore(all_lfps);

lfp_taxis = (1:size(all_lfps,1))/LFP_Fsd;

%%
close all

fig_dir = '/home/james/Desktop/SamFigs/';

MEA_nums = hdr.MEAlayoutNums;

trial_toffset = 0.1;
trial_dispid =2;
trial_dur = 10;
trial_onset_time = vs.stimTimes(trial_dispid)/1e3;
use_inds = find(lfp_taxis >= trial_onset_time & lfp_taxis <= trial_onset_time+trial_dur);
xl = [0 1];

[~,sort_ord] = sort(probe_Y);

cur_lfp_T = lfp_taxis(use_inds) - trial_onset_time - trial_toffset;
f1 = figure();
imagesc(cur_lfp_T,1:252,all_lfps(use_inds,:)');
% imagesc(1:length(use_inds),1:252,all_lfps(use_inds,:)');
xlim(xl);
caxis([-13 13]);
xlabel('Time (s)');
ylabel('Channel Number');

   
[Xo,Yo] = meshgrid(1:16);
[Xq,Yq] = meshgrid(1:0.25:16);

ca = [-14 12];

targ_ind1 = 41;
f2 = figure();
cur_mat = nan(16,16);
cur_mat(~isnan(MEA_nums)) = all_lfps(use_inds(targ_ind1),MEA_nums(~isnan(MEA_nums)));
% subplot(2,1,1);
% imagescnan(cur_mat);
uset = ~isnan(cur_mat);
F = TriScatteredInterp(Xo(uset),Yo(uset),cur_mat(uset));
cur_mat_interp = F(Xq,Yq);
%    caxis(ca);
% subplot(2,1,2);
imagescnan(cur_mat_interp);
%    caxis(ca);
set(gca,'xtick',[],'ytick',[]);

targ_ind2 = 49;
f3 = figure();
cur_mat = nan(16,16);
cur_mat(~isnan(MEA_nums)) = all_lfps(use_inds(targ_ind2),MEA_nums(~isnan(MEA_nums)));
% subplot(2,1,1);
% imagescnan(cur_mat);
uset = ~isnan(cur_mat);
F = TriScatteredInterp(Xo(uset),Yo(uset),cur_mat(uset));
cur_mat_interp = F(Xq,Yq);
%    caxis(ca);
% subplot(2,1,2);
imagescnan(cur_mat_interp);
%    caxis(ca);
set(gca,'xtick',[],'ytick',[]);

targ_ind3 = 71;
f4 = figure();
cur_mat = nan(16,16);
cur_mat(~isnan(MEA_nums)) = all_lfps(use_inds(targ_ind3),MEA_nums(~isnan(MEA_nums)));
% subplot(2,1,1);
% imagescnan(cur_mat);
uset = ~isnan(cur_mat);
F = TriScatteredInterp(Xo(uset),Yo(uset),cur_mat(uset));
cur_mat_interp = F(Xq,Yq);
%    caxis(ca);
% subplot(2,1,2);
imagescnan(cur_mat_interp);
%    caxis(ca);
set(gca,'xtick',[],'ytick',[]);

% targ_ind = 79;
% f5 = figure();
% cur_mat = nan(16,16);
% cur_mat(~isnan(MEA_nums)) = all_lfps(use_inds(targ_ind),MEA_nums(~isnan(MEA_nums)));
% subplot(2,1,1);
% imagescnan(cur_mat);
% uset = ~isnan(cur_mat);
% F = TriScatteredInterp(Xo(uset),Yo(uset),cur_mat(uset));
% cur_mat_interp = F(Xq,Yq);
%    caxis(ca);
% subplot(2,1,2);
% imagescnan(cur_mat_interp);
%    caxis(ca);

figure(f1);
yl = ylim();
line(cur_lfp_T(41) + [0 0],yl,'color','w');
% line(cur_lfp_T(79) + [0 0],yl,'color','w');
% line(cur_lfp_T(66) + [0 0],yl,'color','w');
line(cur_lfp_T(71) + [0 0],yl,'color','w');
line(cur_lfp_T(49) + [0 0],yl,'color','w');


fig_width = 4; rel_height = 1.2;
fig_name = [fig_dir 'Example_LFP_stack.pdf'];
figufy(f1);
exportfig(f1,fig_name,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f1);

fig_width = 3; rel_height = 1;
fig_name = [fig_dir sprintf('LFP_slice_%d.pdf',targ_ind1)];
figufy(f2);
exportfig(f2,fig_name,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f2);

fig_name = [fig_dir sprintf('LFP_slice_%d.pdf',targ_ind2)];
figufy(f3);
exportfig(f3,fig_name,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f3)

fig_name = [fig_dir sprintf('LFP_slice_%d.pdf',targ_ind3)];
figufy(f4);
exportfig(f4,fig_name,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f4)

%%
use_frames = find(lfp_taxis(use_inds) - trial_onset_time > 0.2 & lfp_taxis(use_inds) - trial_onset_time < 1);
movie_X = nan(length(use_frames),size(Xq,1),size(Xq,2));
for ii = 1:length(use_frames)
    cur_mat = nan(16,16);
    cur_mat(~isnan(MEA_nums)) = all_lfps(use_inds(ii),MEA_nums(~isnan(MEA_nums)));
    uset = ~isnan(cur_mat);
    F = TriScatteredInterp(Xo(uset),Yo(uset),cur_mat(uset));
    movie_X(ii,:,:) = F(Xq,Yq);
end

%%
close all
f1 = figure();
for ii = 1:length(use_frames)
   imagescnan(squeeze(movie_X(ii,:,:)));
   caxis([-12 12]);
   pause(0.3);   
end
%%

clear sig_var
tstart_inds = round(interp1(lfp_taxis,1:length(lfp_taxis),trialStartTMS/1e3));
tstart_inds(tstart_inds > length(lfp_taxis) | isnan(tstart_inds)) = [];

backlag = round(0*LFP_Fsd);
forlag = round(3*LFP_Fsd);
[temp,templags,~,~,tempmat] = get_event_trig_avg_v3(all_lfps,tstart_inds,backlag,forlag,2);
tempmat = tempmat';
tempmat = reshape(tempmat,[forlag+backlag+1, size(all_lfps,2), size(tempmat,2)]);

%%
tempFreq_bin_edges = prctile(tempFreq,[0:8:100]);
[~,tfreq_bin] = histc(tempFreq,tempFreq_bin_edges);
%%
params.Fs = LFP_Fsd;
params.tapers = [5 9];
params.trialave = 0;
for cc = 1:size(tempmat,2)
    cc
cur_data = squeeze(tempmat(:,cc,:));
[S,f] = mtspectrumc(cur_data,params);

%%
tot_trials = size(S,2);
un_oris = unique(angleOrder);
ori_Savg = nan(length(un_oris),length(f));
for oo = 1:length(un_oris)
    tset = find(angleOrder == un_oris(oo));
    tset(tset > tot_trials) = [];
    ori_Savg(oo,:) = nanmean(log10(S(:,tset)),2);
end
ori_Savg = bsxfun(@minus,ori_Savg,nanmean(log10(S),2)');
ori_Savg = bsxfun(@rdivide,ori_Savg,nanstd(log10(S),[],2)');
% f1 = figure();
% imagesc(f,un_oris,ori_Savg);
%%
tot_trials = size(S,2);
un_sfs = unique(spatialFreq);
sf_Savg = nan(length(un_sfs),length(f));
for oo = 1:length(un_sfs)
    tset = find(spatialFreq == un_sfs(oo));
    tset(tset > tot_trials) = [];
    sf_Savg(oo,:) = nanmean(log10(S(:,tset)),2);
end
sf_Savg = bsxfun(@minus,sf_Savg,nanmean(log10(S),2)');
sf_Savg = bsxfun(@rdivide,sf_Savg,nanstd(log10(S),[],2)');
% f1 = figure();
% imagesc(f,un_sfs,ori_Savg);

%%
tot_trials = size(S,2);
un_tfs = unique(tfreq_bin);
tf_Savg = nan(length(un_tfs),length(f));
for oo = 1:length(un_tfs)
    tset = find(tfreq_bin == un_tfs(oo));
    tset(tset > tot_trials) = [];
    tf_Savg(oo,:) = nanmean(log10(S(:,tset)),2);
end
tf_Savg = bsxfun(@minus,tf_Savg,nanmean(log10(S),2)');
tf_Savg = bsxfun(@rdivide,tf_Savg,nanstd(log10(S),[],2)');

%%
all_ori_Savg(cc,:,:) = ori_Savg;
all_sf_Savg(cc,:,:) = sf_Savg;
all_tf_Savg(cc,:,:) = tf_Savg;
end

%%
for cc = 1:size(tempmat,2)
%    imagesc(f,un_oris,squeeze(all_ori_Savg(cc,:,:)))
   imagesc(f,un_sfs,squeeze(all_sf_Savg(cc,:,:)))
%    imagesc(f,un_tfs,squeeze(all_tf_Savg(cc,:,:)))
   pause
   clf
   
    
end

%%
pow_range = find(f > 60);
ov_ori_pow = nanmean(all_ori_Savg(:,:,pow_range),3);
ov_sf_pow = nanmean(all_sf_Savg(:,:,pow_range),3);
ov_tf_pow = nanmean(all_tf_Savg(:,:,pow_range),3);
f1 = figure();
subplot(3,1,1);
shadedErrorBar(un_oris,nanmean(ov_ori_pow),nanstd(ov_ori_pow)/sqrt(size(ov_ori_pow,1)));
xlabel('Orientation');
ylabel('Avg pow (z)');
axis tight
subplot(3,1,2);
shadedErrorBar(un_sfs,nanmean(ov_sf_pow),nanstd(ov_sf_pow)/sqrt(size(ov_sf_pow,1)));
xlabel('spatial freq');
ylabel('Avg pow (z)');
axis tight
subplot(3,1,3);
shadedErrorBar(un_tfs,nanmean(ov_tf_pow),nanstd(ov_tf_pow)/sqrt(size(ov_tf_pow,1)));
xlabel('temporal freq');
ylabel('Avg pow (z)');
axis tight
%%
MEA_nums = hdr.MEAlayoutNums;
[mea_vals,mea_ord] = sort(MEA_nums(:));
[XX,YY] = meshgrid(1:12);
xvals = XX(mea_ord);
yvals = YY(mea_ord);
xvals(isnan(mea_vals)) = [];
yvals(isnan(mea_vals)) = [];

[~,xsort] = sort(xvals);
[~,ysort] = sort(yvals);
[~,dsort] = sort(xvals.^2 + yvals.^2);

%%
% for ii = 1:length(templags)
ii = 168
    templags(ii)/LFP_Fsd
Vq = griddata(xvals,yvals,temp(ii,:)',XX,YY);
imagescnan(Vq); caxis([-0.35 0.35])
xlabel('X'); ylabel('Y'); 
% pause(0.2)
% clf
% end

%%
params.Fs = LFP_Fsd;
params.tapers = [2 3];
params.trialave = 1;
movingwin = [0.25 0.1];
win = [-0.5 5];
E = trialStartTMS/1e3;
E(E > lfp_taxis(end)) = [];
E(end) = [];
for cc = 1:size(all_lfps,2)
[S,t,f]=mtspecgramtrigc(all_lfps(:,cc),E,win,movingwin,params);
all_Sgram(cc,:,:) = log10(S);
end
