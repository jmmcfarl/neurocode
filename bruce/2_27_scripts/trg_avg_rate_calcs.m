clear all
% close all
addpath(genpath('~/James_scripts'));

cd ~/Data/bruce/2_27_12/
load ./Blocks.mat

cd /Users/James/Data/bruce/2_27_12/stimrecon
load ./fixation_data_v3

dt = .005;

cellids = 1:10;
muaids = 1:14;

spk_cnts = [spk_cnts(:,cellids) mua_cnts(:,muaids)];

n_used_cells = size(spk_cnts,2);

%%
all_time_axis = [];
all_blockids = [];
all_fixids = [];
spikes_binned = [];
time_since_fix = [];
for blockid = 1:4;
    
    blocktimes = Blocks{blockid}.blocktimes;
    
    cur_all_taxis = blocktimes(1,1):dt:blocktimes(2,end);
    all_time_axis = [all_time_axis cur_all_taxis];
    all_blockids = [all_blockids blockid*ones(size(cur_all_taxis))];
    
    temp_binned = nan(n_used_cells,length(cur_all_taxis));
    for i = 1:size(blocktimes,2)
        cur_taxis = blocktimes(1,i):dt:blocktimes(2,i);
        cur_tinds = round(interp1(cur_all_taxis,1:length(cur_all_taxis),cur_taxis));
        for c = 1:n_used_cells
            if c <= 10
                cur_spk_set = find(Blocks{blockid}.spktimes{cellids(c)} >= cur_taxis(1) & ...
                    Blocks{blockid}.spktimes{cellids(c)} < cur_taxis(end));
                temp = hist(Blocks{blockid}.spktimes{cellids(c)}(cur_spk_set),cur_taxis);
            else
                cur_spk_set = find(Blocks{blockid}.mutimes{muaids(c-10)} >= cur_taxis(1) & ...
                    Blocks{blockid}.mutimes{muaids(c-10)} < cur_taxis(end));
                temp = hist(Blocks{blockid}.mutimes{muaids(c-10)}(cur_spk_set),cur_taxis);
            end
            temp(end) = 0;
            temp_binned(c,cur_tinds) = temp;
        end
    end
    spikes_binned = [spikes_binned temp_binned];
end

%%
all_eye_speed = [];
for blockid = 1:4
    
    fprintf('Processing block %d of %d...\n',blockid,4);
    
    block_times = Blocks{blockid}.blocktimes;
    stim_times = Blocks{blockid}.stimtime;
    stimID = Blocks{blockid}.stimids;
    
    n_stims = length(stim_times);
    
    cd ~/Data/bruce/2_27_12/saccades/
        load(sprintf('lemM232.5%d.em.sac.mat',blockid))    
    
        % identify saccade start and stop times
    EyeStartT = Expt.Trials.Start/10000; % time of first eye sample
    EyeEndT = Expt.Trials.End/10000; % time of last eye sample
    Eyedt = Expt.Header.CRrates(1); % resolution of eye signal (sec)
    eyets = EyeStartT:Eyedt:(EyeEndT-Eyedt); %eye tracking time axis (sec)
    
    reye_pos = [Expt.Trials.Eyevals.rh Expt.Trials.Eyevals.rv];
    leye_pos = [Expt.Trials.Eyevals.lh Expt.Trials.Eyevals.lv];
    
    avg_eyepos = (reye_pos + leye_pos)/2;
    clear sm_avg_eyepos eye_vel
    sm_avg_eyepos(:,1) = smooth(avg_eyepos(:,1),3);
    sm_avg_eyepos(:,2) = smooth(avg_eyepos(:,2),3);
    eye_vel(:,1) = [0; diff(sm_avg_eyepos(:,1))]/Eyedt;
    eye_vel(:,2) = [0; diff(sm_avg_eyepos(:,2))]/Eyedt;
    eye_speed = sqrt(eye_vel(:,1).^2+eye_vel(:,2).^2);

    cur_set = find(all_blockids==blockid);
    interp_eye_speed = interp1(eyets,eye_speed,all_time_axis(cur_set));
    all_eye_speed = [all_eye_speed interp_eye_speed];
end
    
%%

all_fix_start_inds = [];
all_fix_stop_inds = [];
all_stim_inds = [];
for blockid = 1:4
    
   cur_set = find(all_blockids==blockid);
   cur_fix_set = find(blockids==blockid);
   cur_taxis = all_time_axis(cur_set);
   
   cur_fix_start_inds = round(interp1(cur_taxis,1:length(cur_taxis),all_fix_start_times(cur_fix_set)));
   cur_fix_stop_inds = round(interp1(cur_taxis,1:length(cur_taxis),all_fix_stop_times(cur_fix_set)));
    
   cur_stim_start_inds = round(interp1(cur_taxis,1:length(cur_taxis),Blocks{blockid}.stimtime));
   
   all_fix_start_inds = [all_fix_start_inds; cur_set(cur_fix_start_inds(:))'];
   all_fix_stop_inds = [all_fix_stop_inds; cur_set(cur_fix_stop_inds(:))'];
   all_stim_inds = [all_stim_inds; cur_set(cur_stim_start_inds(:))'];
   
end

%%
backwin = round(1/dt);
forwin = round(1/dt);
lags = -backwin:forwin;
n_fixs = length(all_fix_start_inds);
n_stims = length(all_stim_inds);

fix_trg_mat = nan(n_used_cells,n_fixs,length(lags));
for i = 1:n_fixs
    cur_set = (all_fix_start_inds(i)-backwin):all_fix_stop_inds(i);
    cur_set(length(lags)+1:end) = [];
    
    fix_trg_mat(:,i,1:length(cur_set)) = spikes_binned(:,cur_set);    
end

fix_trg_eye = nan(n_fixs,length(lags));
for i = 1:n_fixs
    cur_set = (all_fix_start_inds(i)-backwin):all_fix_stop_inds(i);
    cur_set(length(lags)+1:end) = [];
    
    fix_trg_eye(i,1:length(cur_set)) = all_eye_speed(cur_set);    
end

n_sacs = n_fixs;
sac_trg_mat = nan(n_used_cells,n_fixs,length(lags));
for i = 1:n_sacs
    cur_set = (all_fix_stop_inds(i)-backwin):(all_fix_stop_inds(i)+forwin);
    badset = find(cur_set < all_fix_start_inds(i));
    cur_set(cur_set > size(spikes_binned,2)) = [];
    sac_trg_mat(:,i,1:length(cur_set)) = spikes_binned(:,cur_set); 
    sac_trg_mat(:,i,badset) = nan;
end

stim_trg_mat = nan(n_used_cells,n_stims,length(lags));
for i = 1:n_stims
    cur_set = (all_stim_inds(i)-backwin):(all_stim_inds(i)+forwin);
    cur_set(length(lags)+1:end) = [];
    
    stim_trg_mat(:,i,1:length(cur_set)) = spikes_binned(:,cur_set);    
end

avg_rates = nanmean(spikes_binned,2)/dt;
fix_trg_avgs = squeeze(nanmean(fix_trg_mat,2))/dt;
stim_trg_avgs = squeeze(nanmean(stim_trg_mat,2))/dt;
sac_trg_avgs = squeeze(nanmean(sac_trg_mat,2))/dt;

fix_trg_navgs = bsxfun(@rdivide,fix_trg_avgs,avg_rates);
sac_trg_navgs = bsxfun(@rdivide,sac_trg_avgs,avg_rates);
stim_trg_navgs = bsxfun(@rdivide,stim_trg_avgs,avg_rates);

sm_fix_trg_navgs = fix_trg_navgs;
sm_sac_trg_navgs = sac_trg_navgs;
sm_stim_trg_navgs = stim_trg_navgs;
for i = 1:n_used_cells
   sm_fix_trg_navgs(i,:) = smooth(sm_fix_trg_navgs(i,:),15,'lowess'); 
   sm_sac_trg_navgs(i,:) = smooth(sm_sac_trg_navgs(i,:),15,'lowess'); 
   sm_stim_trg_navgs(i,:) = smooth(sm_stim_trg_navgs(i,:),15,'lowess'); 
end

%%
suprobes = Blocks{1}.suprobes;
muprobes = Blocks{1}.muprobes;
allprobes = [suprobes muprobes];
[dord,ord] = sort(allprobes);

%%
close all
figure
set(gca,'fontname','arial')
imagesc(lags*dt,1:24,sm_fix_trg_navgs(ord,:))
xlim([-0.3 0.4])
set(gca,'ydir','normal')
caxis([0 3])
xlim([-0.1 0.3])
yl = ylim();
line([0 0],yl,'color','w','linewidth',2)
xlabel('Time since fixation onset','fontsize',16)
ylabel('Probe','fontsize',16)
colorbar

figure
imagesc(lags*dt,1:24,sm_stim_trg_navgs(ord,:))
xlim([-0.3 0.4])
set(gca,'ydir','normal')
caxis([0.5 4])
xlim([-0.1 0.3])
caxis([0.5 4])
yl = ylim();
line([0 0],yl,'color','w','linewidth',2)
xlabel('Time since stimulus onset','fontsize',16)
ylabel('Probe','fontsize',16)
colorbar

figure
imagesc(lags*dt,1:24,sm_sac_trg_navgs(ord,:))
xlim([-0.3 0.4])
set(gca,'ydir','normal')
caxis([0.35 2.75])
xlim([-0.15 0.4])
yl = ylim();
line([0 0],yl,'color','w','linewidth',2)
xlabel('Time since stimulus onset','fontsize',16)
ylabel('Probe','fontsize',16)
colorbar

%%
shadedErrorBar(lags*dt,nanmean(sm_fix_trg_navgs),nanstd(sm_fix_trg_navgs)/sqrt(24),{'b'});
hold on
shadedErrorBar(lags*dt,nanmean(sm_stim_trg_navgs),nanstd(sm_stim_trg_navgs)/sqrt(24),{'r'});
% shadedErrorBar(lags*dt,nanmean(sm_sac_trg_navgs),nanstd(sm_sac_trg_navgs)/sqrt(24),{'k'});
xlabel('Time since onset (s)','fontsize',20)
ylabel('Relative rate','fontsize',20)
xlim([-0.5 0.5])
ylim([0 3])
%%
close all
c = 2;
plot(lags*dt,fix_trg_avgs(c,:),'o-')
hold on
plot(lags*dt,smooth(stim_trg_avgs(c,:),5,'lowess'),'ro-')
plot(lags*dt,sac_trg_avgs(c,:),'ko-')
xlim([-0.1 0.3])
set(gca,'fontsize',14)
xlabel('Time (s)','fontsize',16)
ylabel('Firing rate (Hz)','fontsize',16)
