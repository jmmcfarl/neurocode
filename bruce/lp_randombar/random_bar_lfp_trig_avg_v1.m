clear all
close all

ExptNum = 235;
cd(['~/Data/bruce/M' num2str(ExptNum)]);
load ./random_bar_eyedata_ftime.mat

load(['lemM' num2str(ExptNum) 'Expts.mat']);
%%
Fs = 1000;
niqf = Fs/2;
bb_lcf = 1;
bb_hcf = 100;
[b_bb,a_bb] = butter(2,[bb_lcf bb_hcf]/niqf);

dsf = 4;
Fsd = Fs/dsf;

maxlag = round(Fsd*0.5);

all_lfps = [];
all_gam = [];
all_lfp_t_axis = [];
all_expt = [];
for ee = 1:length(bar_expts);
    fprintf('Expt %d of %d\n',ee,length(bar_expts));
    fname = sprintf('lemM%dA.%d.lfp.mat',ExptNum,bar_expts(ee));
    load(fname);
    
    Fs = 1/LFP.Header.CRsamplerate;
    
    n_trials(ee) = length(LFP.Trials);
    lfp_trial_starts = [LFP.Trials(:).ftime]/1e4;
    lfp_trial_ends = [LFP.Trials(:).End]/1e4;
    expt_lfp_t_axis = [];
    expt_lfps = [];
    %     expt_gam = [];
    for tt = 1:n_trials(ee)
        
        cur_npts = size(LFP.Trials(tt).LFP,1);
        cur_t_end(tt) = lfp_trial_starts(tt)+(cur_npts-1)/Fs;
        cur_t_axis = lfp_trial_starts(tt):1/Fs:cur_t_end(tt);
        
        
        if ~isempty(expt_lfp_t_axis)
            cur_sp = find(cur_t_axis > max(expt_lfp_t_axis),1,'first');
        else
            cur_sp = 1;
        end
        cur_t_axis = cur_t_axis(cur_sp:end);
        if length(cur_t_axis) >= 50
            cur_LFP = [LFP.Trials(tt).LFP];
            cur_LFP = cur_LFP(cur_sp:end,:);
            
            cur_LFP = filtfilt(b_bb,a_bb,cur_LFP);
            cur_LFP = downsample(cur_LFP,dsf);
            cur_t_axis = downsample(cur_t_axis,dsf);
            
            expt_lfp_t_axis = [expt_lfp_t_axis; cur_t_axis(:)];
            expt_lfps = [expt_lfps; cur_LFP];
        end
    end
    all_lfps = [all_lfps; expt_lfps];
    all_lfp_t_axis = [all_lfp_t_axis; expt_lfp_t_axis];
    all_expt = [all_expt; ee*ones(length(expt_lfp_t_axis),1)];
end
all_lfp_std = std(all_lfps);

%%
all_fix_inds = [];
all_first_inds = [];
all_second_inds = [];
all_micro_inds = [];
all_trial_start = [];
all_trial_stop = [];
for ee = 1:length(bar_expts)
    expt_trial_starts = [Expts{bar_expts(ee)}.Trials(:).TrialStart]/1e4;
    expt_trial_end = [Expts{bar_expts(ee)}.Trials(:).TrueEnd]/1e4;
    
    sac_start_times = all_eye_ts{ee}(all_sac_startinds{ee});
    sac_stop_times = all_eye_ts{ee}(all_sac_stopinds{ee});
    
    cur_sac_set = find(all_expt_num==ee);
    if length(cur_sac_set) ~= length(sac_start_times)
        error('saccade mismatch')
    end
    intrial_sac_set = find(~isnan(all_trial_num(cur_sac_set)));
    infirst_sac_set = find(ismember(cur_sac_set,all_isfirst));
    insecond_sac_set = find(ismember(cur_sac_set,all_issecond));
    
    micro_sac_set = find(abs(all_delta_pos(cur_sac_set)) < 1 & ~isnan(all_trial_num(cur_sac_set)));
    rid = find(ismember(micro_sac_set,infirst_sac_set) | ismember(micro_sac_set,insecond_sac_set));
    micro_sac_set(rid) = [];
    
    fix_times = sac_start_times(intrial_sac_set);
    first_fix_times = sac_start_times(infirst_sac_set);
    second_fix_times = sac_start_times(insecond_sac_set);
    micro_sac_times = sac_start_times(micro_sac_set);
    
    cur_lfp_set = find(all_expt == ee);
    
    fix_inds = round(interp1(all_lfp_t_axis(cur_lfp_set),1:length(cur_lfp_set),fix_times));
    fix_inds(fix_inds < maxlag | fix_inds > length(cur_lfp_set) - maxlag) = [];
    first_fix_inds = round(interp1(all_lfp_t_axis(cur_lfp_set),1:length(cur_lfp_set),first_fix_times));
    first_fix_inds(first_fix_inds < maxlag | first_fix_inds > length(cur_lfp_set) - maxlag) = [];
    second_fix_inds = round(interp1(all_lfp_t_axis(cur_lfp_set),1:length(cur_lfp_set),second_fix_times));
    second_fix_inds(second_fix_inds < maxlag | second_fix_inds > length(cur_lfp_set) - maxlag) = [];
    
    micro_sac_inds = round(interp1(all_lfp_t_axis(cur_lfp_set),1:length(cur_lfp_set),micro_sac_times));
    micro_sac_inds(micro_sac_inds < maxlag | micro_sac_inds > length(cur_lfp_set) - maxlag) = [];
    
    trial_start_inds = round(interp1(all_lfp_t_axis(cur_lfp_set),1:length(cur_lfp_set),expt_trial_starts));
    trial_start_inds(trial_start_inds < maxlag | trial_start_inds > length(cur_lfp_set) - maxlag) = [];
    
    all_fix_inds = [all_fix_inds; cur_lfp_set(fix_inds(:))];
    all_first_inds = [all_first_inds; cur_lfp_set(first_fix_inds(:))];
    all_second_inds = [all_second_inds; cur_lfp_set(second_fix_inds(:))];
    all_micro_inds = [all_micro_inds; cur_lfp_set(micro_sac_inds(:))];
    all_trial_start = [all_trial_start; cur_lfp_set(trial_start_inds(:))];
end
%%
lags = (-maxlag:maxlag)/Fsd;
tmaxlag = round(Fsd*2);
tlags = (0:tmaxlag)/Fsd;

all_sac_inds = unique([all_first_inds; all_second_inds]);
all_trial_start(all_trial_start > size(all_lfps,1) - tmaxlag);

n_fixs = length(all_trial_start);
trial_trig_mat = nan(n_fixs,length(tlags),24);
for i = 1:n_fixs
    cur_set = (all_trial_start(i)):(all_trial_start(i)+tmaxlag);
    trial_trig_mat(i,:,:) = all_lfps(cur_set,:);
end

% n_fixs = length(all_fix_inds);
% fix_trig_mat = nan(n_fixs,length(lags),24);
% fix_trig_gmat = nan(n_fixs,length(lags),24);
% for i = 1:n_fixs
%    cur_set = (all_fix_inds(i)-maxlag):(all_fix_inds(i)+maxlag);
%    fix_trig_mat(i,:,:) = all_lfps(cur_set,:);
%    fix_trig_gmat(i,:,:) = all_gam(cur_set,:);
% end

n_micros = length(all_micro_inds);
micro_trig_mat = nan(n_micros,length(lags),24);
% fix_trig_gmat = nan(n_micros,length(lags),24);
for i = 1:n_micros
    cur_set = (all_micro_inds(i)-maxlag):(all_micro_inds(i)+maxlag);
    micro_trig_mat(i,:,:) = all_lfps(cur_set,:);
    %    micro_trig_gmat(i,:,:) = all_gam(cur_set,:);
end

n_ffixs = length(all_first_inds);
first_trig_mat = nan(n_ffixs,length(lags),24);
% first_trig_gmat = nan(n_ffixs,length(lags),24);
for i = 1:n_ffixs
    cur_set = (all_first_inds(i)-maxlag):(all_first_inds(i)+maxlag);
    first_trig_mat(i,:,:) = all_lfps(cur_set,:);
    %    first_trig_gmat(i,:,:) = all_gam(cur_set,:);
end

n_sfixs = length(all_second_inds);
second_trig_mat = nan(n_sfixs,length(lags),24);
% second_trig_gmat = nan(n_sfixs,length(lags),24);
for i = 1:n_sfixs
    cur_set = (all_second_inds(i)-maxlag):(all_second_inds(i)+maxlag);
    second_trig_mat(i,:,:) = all_lfps(cur_set,:);
    %    second_trig_gmat(i,:,:) = all_gam(cur_set,:);
end

n_fixs = length(all_sac_inds);
sac_trig_mat = nan(n_fixs,length(lags),24);
% first_trig_gmat = nan(n_ffixs,length(lags),24);
for i = 1:n_fixs
    cur_set = (all_sac_inds(i)-maxlag):(all_sac_inds(i)+maxlag);
    sac_trig_mat(i,:,:) = all_lfps(cur_set,:);
    %    first_trig_gmat(i,:,:) = all_gam(cur_set,:);
end

% n_fixs = length(all_fix_inds);
% n_micros = length(all_micro_inds);
% n_ffixs = length(all_first_inds);
% n_sfixs = length(all_second_inds);
%% COMPUTE TRIG AVG LFPs
sac_trig_avg = squeeze(mean(sac_trig_mat));
msac_trig_avg = squeeze(mean(micro_trig_mat));
fsac_trig_avg = squeeze(mean(first_trig_mat));
ssac_trig_avg = squeeze(mean(second_trig_mat));
trial_trig_avg = squeeze(mean(trial_trig_mat));
%%
to_print = 1;

cd(sprintf('~/Analysis/bruce/M%d',ExptNum));

close all
ca = [-0.6 0.6];

figure
imagesc(lags,1:24,sac_trig_avg');colorbar
xlim([-0.2 0.4])
caxis(ca);
if to_print
    fname = 'Sac_trig_LFP';
    print('-dpng',fname);
    close
end

figure
imagesc(lags,1:24,msac_trig_avg');colorbar
xlim([-0.2 0.4])
caxis(ca);
if to_print
    fname = 'Msac_trig_LFP';
    print('-dpng',fname);
    close
end

figure
subplot(2,1,1)
imagesc(lags,1:24,fsac_trig_avg');colorbar
xlim([-0.2 0.4])
caxis(ca);
subplot(2,1,2)
imagesc(lags,1:24,ssac_trig_avg');colorbar
xlim([-0.2 0.4])
caxis(ca);
if to_print
    fname = 'FSsac_trig_LFP';
    print('-dpng',fname);
    close
end

figure
imagesc(tlags,1:24,trial_trig_avg');colorbar
caxis(ca);
if to_print
    fname = 'Trial_trig_LFP';
    print('-dpng',fname);
    close
end

%%
addpath(genpath('~/James_scripts/iCSD/'));
%base parameters
vars.Fs = Fs; %sample freq
vars.BrainBound = 1; %first channel that is in the brain
vars.ChanSep = 0.05; %channel sep in mm

%parameters for spline and step methods
vars.diam = 2; %current disc diameter (in mm)

%parameters for standard method
vars.useVaknin = 'true'; %method for handling boundary conditions
vars.useHamming = 'true'; %use hamming smoothing
method = 'spline';

Data = permute(micro_trig_mat,[3 2 1]);
mCSD = PettersenCSD(Data,method,vars);
micro_CSD_avg = squeeze(mean(mCSD,3));

Data = permute(sac_trig_mat,[3 2 1]);
CSD = PettersenCSD(Data,method,vars);
sac_CSD_avg = squeeze(mean(CSD,3));

Data = permute(trial_trig_mat,[3 2 1]);
CSD = PettersenCSD(Data,method,vars);
trial_CSD_avg = squeeze(mean(CSD,3));

Data = permute(first_trig_mat,[3 2 1]);
fCSD = PettersenCSD(Data,method,vars);
first_CSD_avg = squeeze(mean(fCSD,3));

Data = permute(second_trig_mat,[3 2 1]);
sCSD = PettersenCSD(Data,method,vars);
second_CSD_avg = squeeze(mean(sCSD,3));

%%
% % close all
% ca = 2e6;
%
% ch1 = 3;
% ch2 = 12;
% ch3 = 21;
% figure
% subplot(3,1,1)
% hold on
% shadedErrorBar(lags,squeeze(mean(boot_fCSD(ch1,:,:),2)),squeeze(std(boot_fCSD(ch1,:,:),[],2)));
% shadedErrorBar(lags,squeeze(mean(boot_sCSD(ch1,:,:),2)),squeeze(std(boot_sCSD(ch1,:,:),[],2)),{'r'});
% shadedErrorBar(lags,squeeze(mean(boot_mCSD(ch1,:,:),2)),squeeze(std(boot_mCSD(ch1,:,:),[],2)),{'b'});
% xlim([-0.4 0.4])
% % ca = max(abs(ylim()));
% ylim([-ca ca])
% line([0 0],[-ca ca],'color','k')
% xl = xlim();
% line(xl,[0 0],'color','k')
% xlabel('Time since saccade onset (s)','fontsize',16)
% ylabel('Current source','fontsize',16)
%
% subplot(3,1,2)
% hold on
% shadedErrorBar(lags,squeeze(mean(boot_fCSD(ch2,:,:),2)),squeeze(std(boot_fCSD(ch2,:,:),[],2)));
% shadedErrorBar(lags,squeeze(mean(boot_sCSD(ch2,:,:),2)),squeeze(std(boot_sCSD(ch2,:,:),[],2)),{'r'});
% shadedErrorBar(lags,squeeze(mean(boot_mCSD(ch2,:,:),2)),squeeze(std(boot_mCSD(ch2,:,:),[],2)),{'b'});
% xlim([-0.4 0.4])
% % ca = max(abs(ylim()));
% ylim([-ca ca])
% line([0 0],[-ca ca],'color','k')
% xl = xlim();
% line(xl,[0 0],'color','k')
% xlabel('Time since saccade onset (s)','fontsize',16)
% ylabel('Current source','fontsize',16)
%
% subplot(3,1,3)
% hold on
% shadedErrorBar(lags,squeeze(mean(boot_fCSD(ch3,:,:),2)),squeeze(std(boot_fCSD(ch3,:,:),[],2)));
% shadedErrorBar(lags,squeeze(mean(boot_sCSD(ch3,:,:),2)),squeeze(std(boot_sCSD(ch3,:,:),[],2)),{'r'});
% shadedErrorBar(lags,squeeze(mean(boot_mCSD(ch3,:,:),2)),squeeze(std(boot_mCSD(ch3,:,:),[],2)),{'b'});
% xlim([-0.4 0.4])
% % ca = max(abs(ylim()));
% ylim([-ca ca])
% xl = xlim();
% line(xl,[0 0],'color','k')
% line([0 0],[-ca ca],'color','k')
% xlabel('Time since saccade onset (s)','fontsize',16)
% ylabel('Current source','fontsize',16)

%%
close all
to_print = 0;

cd(sprintf('~/Analysis/bruce/M%d',ExptNum));

figure
imagesc(lags,(1:24),micro_CSD_avg);
% set(gca,'ydir','normal');
xlim([-0.2 0.4])
xlabel('Time lag (s)','fontsize',16)
ylabel('Channel','fontsize',16)
colorbar
colormap(colormap_redblackblue);
ca = max(abs(caxis()));
caxis([-ca ca]*0.75);
yl = ylim();
line([0 0],yl,'color','w');
if to_print
    fname = 'Micro_trig_CSD';
    print('-dpng',fname);
    close
end


figure
imagesc(tlags,(1:24),trial_CSD_avg);
% set(gca,'ydir','normal');
xlim([0 0.4])
xlabel('Time lag (s)','fontsize',16)
ylabel('Channel','fontsize',16)
colorbar
colormap(colormap_redblackblue);
ca = max(abs(caxis()));
caxis([-ca ca]*0.75);
yl = ylim();
line([0 0],yl,'color','w');
if to_print
    fname = 'Trial_trig_CSD';
    print('-dpng',fname);
    close
end


figure
imagesc(lags,(1:24),sac_CSD_avg);
% set(gca,'ydir','normal');
xlim([-0.2 0.4])
xlabel('Time lag (s)','fontsize',16)
ylabel('Channel','fontsize',16)
colorbar
colormap(colormap_redblackblue);
ca = max(abs(caxis()));
caxis([-ca ca]*0.75);
yl = ylim();
line([0 0],yl,'color','w');
if to_print
    fname = 'Sac_trig_CSD';
    print('-dpng',fname);
    close
end

figure
imagesc(lags,(1:24),first_CSD_avg);
% set(gca,'ydir','normal');
xlim([-0.2 0.4])
xlabel('Time lag (s)','fontsize',16)
ylabel('Channel','fontsize',16)
colorbar
colormap(colormap_redblackblue);
ca = max(abs(caxis()));
caxis([-ca ca]*0.75);
yl = ylim();
line([0 0],yl,'color','w');
if to_print
    fname = 'Firstsac_trig_CSD';
    print('-dpng',fname);
    close
end

figure
imagesc(lags,(1:24),second_CSD_avg);
% set(gca,'ydir','normal');
xlim([-0.2 0.4])
xlabel('Time lag (s)','fontsize',16)
ylabel('Channel','fontsize',16)
colorbar
colormap(colormap_redblackblue);
ca = max(abs(caxis()));
caxis([-ca ca]*0.75);
yl = ylim();
line([0 0],yl,'color','w');
if to_print
    fname = 'Secondsac_trig_CSD';
    print('-dpng',fname);
    close
end

%%
% % Data = permute(micro_trig_mat,[3 2 1]);
% % CSD = PettersenCSD(Data,'spline',vars);
%
% for i = 1:n_micros
%     subplot(2,1,1)
%     imagesc(lags,1:24,squeeze(micro_trig_mat(i,:,:))')
%     ca = max(abs(caxis()));
%     caxis([-ca ca]);
%     subplot(2,1,2)
%     imagesc(lags,1:24,squeeze(CSD(:,:,i)))
%     ca = max(abs(caxis()));
%     caxis([-ca ca]);
%     pause
%     clf
% end
%%
% close all
%
% cur_el = 1;
% figure
% % shadedErrorBar(lags,squeeze(mean(fix_trig_mat(:,:,cur_el))),squeeze(std(fix_trig_mat(:,:,cur_el)))/sqrt(n_fixs),{'color','k'})
% shadedErrorBar(lags,squeeze(mean(micro_trig_mat(:,:,cur_el))),squeeze(std(micro_trig_mat(:,:,cur_el)))/sqrt(n_micros),{'color','k'})
% hold on
% shadedErrorBar(lags,squeeze(mean(first_trig_mat(:,:,cur_el))),squeeze(std(first_trig_mat(:,:,cur_el)))/sqrt(n_ffixs),{'color','r'})
% shadedErrorBar(lags,squeeze(mean(second_trig_mat(:,:,cur_el))),squeeze(std(second_trig_mat(:,:,cur_el)))/sqrt(n_sfixs),{'color','b'})
%
% cur_el = 14;
% figure
% % shadedErrorBar(lags,squeeze(mean(fix_trig_mat(:,:,cur_el))),squeeze(std(fix_trig_mat(:,:,cur_el)))/sqrt(n_fixs),{'color','k'})
% shadedErrorBar(lags,squeeze(mean(micro_trig_mat(:,:,cur_el))),squeeze(std(micro_trig_mat(:,:,cur_el)))/sqrt(n_micros),{'color','k'})
% hold on
% shadedErrorBar(lags,squeeze(mean(first_trig_mat(:,:,cur_el))),squeeze(std(first_trig_mat(:,:,cur_el)))/sqrt(n_ffixs),{'color','r'})
% shadedErrorBar(lags,squeeze(mean(second_trig_mat(:,:,cur_el))),squeeze(std(second_trig_mat(:,:,cur_el)))/sqrt(n_sfixs),{'color','b'})
%
% cur_el = 24;
% figure
% % shadedErrorBar(lags,squeeze(mean(fix_trig_mat(:,:,cur_el))),squeeze(std(fix_trig_mat(:,:,cur_el)))/sqrt(n_fixs),{'color','k'})
% shadedErrorBar(lags,squeeze(mean(micro_trig_mat(:,:,cur_el))),squeeze(std(micro_trig_mat(:,:,cur_el)))/sqrt(n_micros),{'color','k'})
% hold on
% shadedErrorBar(lags,squeeze(mean(first_trig_mat(:,:,cur_el))),squeeze(std(first_trig_mat(:,:,cur_el)))/sqrt(n_ffixs),{'color','r'})
% shadedErrorBar(lags,squeeze(mean(second_trig_mat(:,:,cur_el))),squeeze(std(second_trig_mat(:,:,cur_el)))/sqrt(n_sfixs),{'color','b'})

% %%
% micro_trig_avgs = squeeze(mean(micro_trig_mat))';
% micro_trig_avgs_n = bsxfun(@times,micro_trig_avgs,all_lfp_std');
% csd = 2*micro_trig_avgs_n(2:end-1,:) - micro_trig_avgs_n(1:end-2,:) - micro_trig_avgs_n(3:end,:);
%
% %%
% micro_trig_mat = bsxfun(@times,micro_trig_mat,reshape(all_lfp_std,[1,1,24]));
% csd_mat = 2*micro_trig_mat(:,:,2:end-1) - micro_trig_mat(:,:,1:end-2) - micro_trig_mat(:,:,3:end);
%%
% close all
%
% cur_el = 1;
% figure
% shadedErrorBar(lags,squeeze(mean(fix_trig_gmat(:,:,cur_el))),squeeze(std(fix_trig_gmat(:,:,cur_el)))/sqrt(n_fixs),{'color','k'})
% hold on
% shadedErrorBar(lags,squeeze(mean(first_trig_gmat(:,:,cur_el))),squeeze(std(first_trig_gmat(:,:,cur_el)))/sqrt(n_ffixs),{'color','r'})
% shadedErrorBar(lags,squeeze(mean(second_trig_gmat(:,:,cur_el))),squeeze(std(second_trig_gmat(:,:,cur_el)))/sqrt(n_sfixs),{'color','b'})
%
% cur_el = 14;
% figure
% shadedErrorBar(lags,squeeze(mean(fix_trig_gmat(:,:,cur_el))),squeeze(std(fix_trig_gmat(:,:,cur_el)))/sqrt(n_fixs),{'color','k'})
% hold on
% shadedErrorBar(lags,squeeze(mean(first_trig_gmat(:,:,cur_el))),squeeze(std(first_trig_gmat(:,:,cur_el)))/sqrt(n_ffixs),{'color','r'})
% shadedErrorBar(lags,squeeze(mean(second_trig_gmat(:,:,cur_el))),squeeze(std(second_trig_gmat(:,:,cur_el)))/sqrt(n_sfixs),{'color','b'})
%
% cur_el = 24;
% figure
% shadedErrorBar(lags,squeeze(mean(fix_trig_gmat(:,:,cur_el))),squeeze(std(fix_trig_gmat(:,:,cur_el)))/sqrt(n_fixs),{'color','k'})
% hold on
% shadedErrorBar(lags,squeeze(mean(first_trig_gmat(:,:,cur_el))),squeeze(std(first_trig_gmat(:,:,cur_el)))/sqrt(n_ffixs),{'color','r'})
% shadedErrorBar(lags,squeeze(mean(second_trig_gmat(:,:,cur_el))),squeeze(std(second_trig_gmat(:,:,cur_el)))/sqrt(n_sfixs),{'color','b'})
