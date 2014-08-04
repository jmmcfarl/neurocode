clear all
% close all

data_dir_base = '~/Analysis/bruce';
Expt_nums =  [93];
%% LOAD OVERALL SU DATA
load ~/Analysis/bruce/summary_analysis/su_data.mat
mahal_thresh = su_data.mah_thresh;

all_sua_data = [];
all_mua_data = [];
all_mu_probes = [];
all_mu_exnums = [];
all_su_probes = [];
all_su_exnums = [];
all_su_avgspks = [];
all_su_stdspks = [];
all_su_avgmahal = [];
for ex = 1:length(Expt_nums)
    fprintf('Loading data from expt %d\n',Expt_nums(ex));
    fname = [data_dir_base sprintf('/G0%d',Expt_nums(ex)) '/sac_mod/full_wi_sacmod_data.mat'];
    load(fname);
    
    all_mua_data = cat(1,all_mua_data,mua_data');
    all_mu_probes = cat(1,all_mu_probes,(1:96)');
    all_mu_exnums = cat(1,all_mu_exnums,Expt_nums(ex)*ones(96,1));

    all_sua_data = cat(1,all_sua_data,sua_data');
    all_su_probes = cat(1,all_su_probes,su_probes');
    all_su_exnums = cat(1,all_su_exnums,Expt_nums(ex)*ones(length(su_probes),1));
    
    cur_su_avgmahal = su_data.avg_mahal(ex,su_probes);
    all_su_avgmahal = [all_su_avgmahal; cur_su_avgmahal'];
end
lags = anal_params.lags*anal_params.dt;
dt = anal_params.dt;

%% view all expt-avg mu sacmod
cd ~/Analysis/bruce/summary_analysis/

bad_probes = [16; 92];
close all
n_un_wi = 3;
% f1 = figure();
% f2 = figure();
% f3 = figure();
ex = 1;
% for ex = 1:length(Expt_nums)
    fprintf('Expt %d\n',Expt_nums(ex));
    cur_mu_set = find(all_mu_exnums == Expt_nums(ex));
    cur_su_set = find(all_su_exnums == Expt_nums(ex));
    if length(cur_mu_set) ~= 96
        error('Wrong mu Num');
    end
    exclude = [bad_probes; all_su_probes(cur_su_set)];
    cur_mu_set(exclude) = [];
    
    mu_gsac_data = get_struct_data(all_mua_data,cur_mu_set,'wi_gsac_avg');
    mu_msac_data = get_struct_data(all_mua_data,cur_mu_set,'wi_msac_avg');
    mu_simsac_data = get_struct_data(all_mua_data,cur_mu_set,'wi_simsac_avg');

    mu_gsac_data = reshape(mu_gsac_data,[length(cur_mu_set) n_un_wi length(lags)]);
    mu_msac_data = reshape(mu_msac_data,[length(cur_mu_set) n_un_wi length(lags)]);
    mu_simsac_data = reshape(mu_simsac_data,[length(cur_mu_set) n_un_wi length(lags)]);
    
    figure(f1);
    subplot(3,1,1);hold on; grid on %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h1 = shadedErrorBar(lags,squeeze(nanmean(mu_gsac_data(:,1,:))),squeeze(nanstd(mu_gsac_data(:,1,:)))/sqrt(length(cur_mu_set)),{'color','r'});
    h2 = shadedErrorBar(lags,squeeze(nanmean(mu_gsac_data(:,2,:))),squeeze(nanstd(mu_gsac_data(:,2,:)))/sqrt(length(cur_mu_set)),{'color','b'});
    h3 = shadedErrorBar(lags,squeeze(nanmean(mu_gsac_data(:,3,:))),squeeze(nanstd(mu_gsac_data(:,3,:)))/sqrt(length(cur_mu_set)),{'color','k'});
    legend([h1.mainLine h2.mainLine h3.mainLine],{'1 deg','2 deg','4 deg'});legend('boxoff');
    xlim([-0.2 0.5]);
    ylim([0.7 1.4]);
    xl = xlim(); yl = ylim();
    line(xl,[1 1],'color','k','linestyle','--');
    line([0 0],yl,'color','k','linestyle','--');
    xlabel('Time since saccade onset (s)','fontsize',14)
    ylabel('Relative firing rate','fontsize',14)
    set(gca,'fontsize',12); box off
    title('Guided Saccades','fontsize',16);
    subplot(3,1,2);hold on; grid on %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h1 = shadedErrorBar(lags,squeeze(nanmean(mu_msac_data(:,1,:))),squeeze(nanstd(mu_msac_data(:,1,:)))/sqrt(length(cur_mu_set)),{'color','r'});
    h2 = shadedErrorBar(lags,squeeze(nanmean(mu_msac_data(:,2,:))),squeeze(nanstd(mu_msac_data(:,2,:)))/sqrt(length(cur_mu_set)),{'color','b'});
    h3 = shadedErrorBar(lags,squeeze(nanmean(mu_msac_data(:,3,:))),squeeze(nanstd(mu_msac_data(:,3,:)))/sqrt(length(cur_mu_set)),{'color','k'});
    legend([h1.mainLine h2.mainLine h3.mainLine],{'1 deg','2 deg','4 deg'});legend('boxoff');
    xlim([-0.2 0.5]);
    ylim([0.7 1.3]);
    xl = xlim(); yl = ylim();
    line(xl,[1 1],'color','k','linestyle','--');
    line([0 0],yl,'color','k','linestyle','--');
    xlabel('Time since saccade onset (s)','fontsize',14)
    ylabel('Relative firing rate','fontsize',14)
    set(gca,'fontsize',12); box off
    title('Micro-saccades','fontsize',16);
    subplot(3,1,3);hold on; grid on %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h1 = shadedErrorBar(lags,squeeze(nanmean(mu_simsac_data(:,1,:))),squeeze(nanstd(mu_simsac_data(:,1,:)))/sqrt(length(cur_mu_set)),{'color','r'});
    h2 = shadedErrorBar(lags,squeeze(nanmean(mu_simsac_data(:,2,:))),squeeze(nanstd(mu_simsac_data(:,2,:)))/sqrt(length(cur_mu_set)),{'color','b'});
    h3 = shadedErrorBar(lags,squeeze(nanmean(mu_simsac_data(:,3,:))),squeeze(nanstd(mu_simsac_data(:,3,:)))/sqrt(length(cur_mu_set)),{'color','k'});
    legend([h1.mainLine h2.mainLine h3.mainLine],{'1 deg','2 deg','4 deg'});legend('boxoff');
    xlim([-0.2 0.5]);
    ylim([0.85 1.15]);
    xl = xlim(); yl = ylim();
    line(xl,[1 1],'color','k','linestyle','--');
    line([0 0],yl,'color','k','linestyle','--');
    xlabel('Time since saccade onset (s)','fontsize',14)
    ylabel('Relative firing rate','fontsize',14)
    set(gca,'fontsize',12); box off
    title('Simulated Saccades','fontsize',16);
    
    fname = 'mu_stripewidth_trigavgs';
    fillPage(gcf,'papersize',[6 15]);
    print(fname,'-dpdf','-painters');close
    

%     pause
%     figure(f1); clf;
%     figure(f2); clf;
%     figure(f3); clf;
% end

%% view all expt-avg mu sacmod
bad_probes = [16; 92];
close all
n_un_wi = 3;
% f1 = figure();
% f2 = figure();
% f3 = figure();
% f4 = figure();
% f5 = figure();
% for ex = 1:length(Expt_nums)
    fprintf('Expt %d\n',Expt_nums(ex));
    cur_mu_set = find(all_mu_exnums == Expt_nums(ex));
    cur_su_set = find(all_su_exnums == Expt_nums(ex));
    if length(cur_mu_set) ~= 96
        error('Wrong mu Num');
    end
    exclude = [bad_probes; all_su_probes(cur_su_set)];
    cur_mu_set(exclude) = [];
    
    mu_gsac_data = get_struct_data(all_mua_data,cur_mu_set,'wi_gsac_avg');
    mu_msac_data = get_struct_data(all_mua_data,cur_mu_set,'wi_msac_avg');
    mu_simsac_data = get_struct_data(all_mua_data,cur_mu_set,'wi_simsac_avg');
    mu_gsac_im_data = get_struct_data(all_mua_data,cur_mu_set,'wi_gsac_im_avg');
    mu_gsac_gray_data = get_struct_data(all_mua_data,cur_mu_set,'wi_gsac_gray_avg');
    mu_msac_im_data = get_struct_data(all_mua_data,cur_mu_set,'wi_msac_im_avg');
    mu_msac_gray_data = get_struct_data(all_mua_data,cur_mu_set,'wi_msac_gray_avg');

    mu_gsac_data = reshape(mu_gsac_data,[length(cur_mu_set) n_un_wi length(lags)]);
    mu_msac_data = reshape(mu_msac_data,[length(cur_mu_set) n_un_wi length(lags)]);
    mu_simsac_data = reshape(mu_simsac_data,[length(cur_mu_set) n_un_wi length(lags)]);
    mu_gsac_im_data = reshape(mu_gsac_im_data,[length(cur_mu_set) n_un_wi length(lags)]);
    mu_gsac_gray_data = reshape(mu_gsac_gray_data,[length(cur_mu_set) n_un_wi length(lags)]);
    mu_msac_im_data = reshape(mu_msac_im_data,[length(cur_mu_set) n_un_wi length(lags)]);
    mu_msac_gray_data = reshape(mu_msac_gray_data,[length(cur_mu_set) n_un_wi length(lags)]);
    
    all_cond_avgs = get_struct_data(all_mua_data,cur_mu_set,'cond_avg_rate');
    mu_gsac_data = bsxfun(@rdivide,mu_gsac_data,all_cond_avgs);
    mu_msac_data = bsxfun(@rdivide,mu_msac_data,all_cond_avgs);
    mu_simsac_data = bsxfun(@rdivide,mu_simsac_data,all_cond_avgs);
    mu_gsac_im_data = bsxfun(@rdivide,mu_gsac_im_data,all_cond_avgs);
    mu_gsac_gray_data = bsxfun(@rdivide,mu_gsac_gray_data,all_cond_avgs);
    mu_msac_im_data = bsxfun(@rdivide,mu_msac_im_data,all_cond_avgs);
    mu_msac_gray_data = bsxfun(@rdivide,mu_msac_gray_data,all_cond_avgs);
    
    
    figure(f1);
    subplot(3,1,1);hold on; grid on %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h1 = shadedErrorBar(lags,squeeze(nanmean(mu_gsac_data(:,1,:))),squeeze(nanstd(mu_gsac_data(:,1,:)))/sqrt(length(cur_mu_set)),{'color','r'});
    h2 = shadedErrorBar(lags,squeeze(nanmean(mu_gsac_data(:,2,:))),squeeze(nanstd(mu_gsac_data(:,2,:)))/sqrt(length(cur_mu_set)),{'color','b'});
    h3 = shadedErrorBar(lags,squeeze(nanmean(mu_gsac_data(:,3,:))),squeeze(nanstd(mu_gsac_data(:,3,:)))/sqrt(length(cur_mu_set)),{'color','k'});
    legend([h1.mainLine h2.mainLine h3.mainLine],{'1 deg','2 deg','4 deg'});legend('boxoff');
    xlim([-0.2 0.5]);
    ylim([0.65 1.4]);
    xl = xlim(); yl = ylim();
    line(xl,[1 1],'color','k','linestyle','--');
    line([0 0],yl,'color','k','linestyle','--');
    xlabel('Time since saccade onset (s)','fontsize',14)
    ylabel('Relative firing rate','fontsize',14)
    set(gca,'fontsize',12); box off
    title('Guided Saccades','fontsize',16);
    subplot(3,1,2);hold on; grid on %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h1 = shadedErrorBar(lags,squeeze(nanmean(mu_msac_data(:,1,:))),squeeze(nanstd(mu_msac_data(:,1,:)))/sqrt(length(cur_mu_set)),{'color','r'});
    h2 = shadedErrorBar(lags,squeeze(nanmean(mu_msac_data(:,2,:))),squeeze(nanstd(mu_msac_data(:,2,:)))/sqrt(length(cur_mu_set)),{'color','b'});
    h3 = shadedErrorBar(lags,squeeze(nanmean(mu_msac_data(:,3,:))),squeeze(nanstd(mu_msac_data(:,3,:)))/sqrt(length(cur_mu_set)),{'color','k'});
    legend([h1.mainLine h2.mainLine h3.mainLine],{'1 deg','2 deg','4 deg'});legend('boxoff');
    xlim([-0.2 0.5]);
    ylim([0.7 1.3]);
    xl = xlim(); yl = ylim();
    line(xl,[1 1],'color','k','linestyle','--');
    line([0 0],yl,'color','k','linestyle','--');
    xlabel('Time since saccade onset (s)','fontsize',14)
    ylabel('Relative firing rate','fontsize',14)
    set(gca,'fontsize',12); box off
    title('Micro-saccades','fontsize',16);
    subplot(3,1,3);hold on; grid on %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h1 = shadedErrorBar(lags,squeeze(nanmean(mu_simsac_data(:,1,:))),squeeze(nanstd(mu_simsac_data(:,1,:)))/sqrt(length(cur_mu_set)),{'color','r'});
    h2 = shadedErrorBar(lags,squeeze(nanmean(mu_simsac_data(:,2,:))),squeeze(nanstd(mu_simsac_data(:,2,:)))/sqrt(length(cur_mu_set)),{'color','b'});
    h3 = shadedErrorBar(lags,squeeze(nanmean(mu_simsac_data(:,3,:))),squeeze(nanstd(mu_simsac_data(:,3,:)))/sqrt(length(cur_mu_set)),{'color','k'});
    legend([h1.mainLine h2.mainLine h3.mainLine],{'1 deg','2 deg','4 deg'});legend('boxoff');
    xlim([-0.2 0.5]);
    ylim([0.9 1.1]);
    xl = xlim(); yl = ylim();
    line(xl,[1 1],'color','k','linestyle','--');
    line([0 0],yl,'color','k','linestyle','--');
    xlabel('Time since saccade onset (s)','fontsize',14)
    ylabel('Relative firing rate','fontsize',14)
    set(gca,'fontsize',12); box off
    title('Simulated Saccades','fontsize',16);
    
    fname = 'mu_stripewidth_trigavgs_trialnorm';
    fillPage(gcf,'papersize',[6 15]);
    print(fname,'-dpdf','-painters');close

%     figure(f4);
%     subplot(2,1,1);hold on; grid on
%     plot(lags,squeeze(nanmean(mu_gsac_gray_data(:,1,:))),'r');
%     plot(lags,squeeze(nanmean(mu_gsac_gray_data(:,2,:))),'b');
%     plot(lags,squeeze(nanmean(mu_gsac_gray_data(:,3,:))),'k');
%     xlim([-0.2 0.5]);
%     ylim([0.65 1.4]);
%     xl = xlim(); yl = ylim();
%     line(xl,[1 1],'color','k','linestyle','--');
%     line([0 0],yl,'color','k','linestyle','--');
%     subplot(2,1,2);hold on; grid on;
%     plot(lags,squeeze(nanmean(mu_gsac_im_data(:,1,:))),'r--');
%      plot(lags,squeeze(nanmean(mu_gsac_im_data(:,2,:))),'b--');
%     plot(lags,squeeze(nanmean(mu_gsac_im_data(:,3,:))),'k--');
%    xlim([-0.2 0.5]);
%     ylim([0.65 1.4]);
%     xl = xlim(); yl = ylim();
%     line(xl,[1 1],'color','k','linestyle','--');
%     line([0 0],yl,'color','k','linestyle','--');
% 
%     figure(f5);
%     subplot(2,1,1);hold on; grid on
%     plot(lags,squeeze(nanmean(mu_msac_gray_data(:,1,:))),'r');
%     plot(lags,squeeze(nanmean(mu_msac_gray_data(:,2,:))),'b');
%     plot(lags,squeeze(nanmean(mu_msac_gray_data(:,3,:))),'k');
%    xlim([-0.2 0.5]);
%     ylim([0.65 1.4]);
%     xl = xlim(); yl = ylim();
%     line(xl,[1 1],'color','k','linestyle','--');
%     line([0 0],yl,'color','k','linestyle','--');
%     subplot(2,1,2);hold on; grid on
%     plot(lags,squeeze(nanmean(mu_msac_im_data(:,1,:))),'r--');
%      plot(lags,squeeze(nanmean(mu_msac_im_data(:,2,:))),'b--');
%     plot(lags,squeeze(nanmean(mu_msac_im_data(:,3,:))),'k--');
%    xlim([-0.2 0.5]);
%     ylim([0.65 1.4]);
%     xl = xlim(); yl = ylim();
%     line(xl,[1 1],'color','k','linestyle','--');
%     line([0 0],yl,'color','k','linestyle','--');
%     pause
%     figure(f1); clf;
%     figure(f2); clf;
%     figure(f3); clf;
% end

%% view all ind mu sacmod
bad_probes = [16; 92];
close all
n_un_wi = 3;
ex = 1;
fprintf('Expt %d\n',Expt_nums(ex));
cur_mu_set = find(all_mu_exnums == Expt_nums(ex));
cur_su_set = find(all_su_exnums == Expt_nums(ex));
if length(cur_mu_set) ~= 96
    error('Wrong mu Num');
end
exclude = [bad_probes; all_su_probes(cur_su_set)];
cur_mu_set(exclude) = [];

mu_gsac_data = get_struct_data(all_mua_data,cur_mu_set,'wi_gsac_avg');
mu_msac_data = get_struct_data(all_mua_data,cur_mu_set,'wi_msac_avg');
mu_simsac_data = get_struct_data(all_mua_data,cur_mu_set,'wi_simsac_avg');
mu_gsac_im_data = get_struct_data(all_mua_data,cur_mu_set,'wi_gsac_im_avg');
mu_gsac_gray_data = get_struct_data(all_mua_data,cur_mu_set,'wi_gsac_gray_avg');
mu_msac_im_data = get_struct_data(all_mua_data,cur_mu_set,'wi_msac_im_avg');
mu_msac_gray_data = get_struct_data(all_mua_data,cur_mu_set,'wi_msac_gray_avg');

mu_gsac_data = reshape(mu_gsac_data,[length(cur_mu_set) n_un_wi length(lags)]);
mu_msac_data = reshape(mu_msac_data,[length(cur_mu_set) n_un_wi length(lags)]);
mu_simsac_data = reshape(mu_simsac_data,[length(cur_mu_set) n_un_wi length(lags)]);
mu_gsac_im_data = reshape(mu_gsac_im_data,[length(cur_mu_set) n_un_wi length(lags)]);
mu_gsac_gray_data = reshape(mu_gsac_gray_data,[length(cur_mu_set) n_un_wi length(lags)]);
mu_msac_im_data = reshape(mu_msac_im_data,[length(cur_mu_set) n_un_wi length(lags)]);
mu_msac_gray_data = reshape(mu_msac_gray_data,[length(cur_mu_set) n_un_wi length(lags)]);

all_cond_avgs = get_struct_data(all_mua_data,cur_mu_set,'cond_avg_rate');
mu_gsac_data = bsxfun(@rdivide,mu_gsac_data,all_cond_avgs);
mu_msac_data = bsxfun(@rdivide,mu_msac_data,all_cond_avgs);
mu_simsac_data = bsxfun(@rdivide,mu_simsac_data,all_cond_avgs);
mu_gsac_im_data = bsxfun(@rdivide,mu_gsac_im_data,all_cond_avgs);
mu_gsac_gray_data = bsxfun(@rdivide,mu_gsac_gray_data,all_cond_avgs);
mu_msac_im_data = bsxfun(@rdivide,mu_msac_im_data,all_cond_avgs);
mu_msac_gray_data = bsxfun(@rdivide,mu_msac_gray_data,all_cond_avgs);

poss_lags = find(lags > 0 & lags < 0.2);

[mua_gsac_sup_rate,mua_gsac_sup_time] = min(mu_gsac_data(:,:,poss_lags),[],3);
mua_gsac_sup_time = lags(mua_gsac_sup_time + poss_lags(1) - 1);
[mua_gsac_enh_rate,mua_gsac_enh_time] = max(mu_gsac_data(:,:,poss_lags),[],3);
mua_gsac_enh_time = lags(mua_gsac_enh_time + poss_lags(1) - 1);
mua_gsac_sup_rate = 1-mua_gsac_sup_rate; mua_gsac_enh_rate = mua_gsac_enh_rate - 1;

[mua_gsac_im_sup_rate,mua_gsac_im_sup_time] = min(mu_gsac_im_data(:,:,poss_lags),[],3);
mua_gsac_im_sup_time = lags(mua_gsac_im_sup_time + poss_lags(1) - 1);
[mua_gsac_im_enh_rate,mua_gsac_im_enh_time] = max(mu_gsac_im_data(:,:,poss_lags),[],3);
mua_gsac_im_enh_time = lags(mua_gsac_im_enh_time + poss_lags(1) - 1);
[mua_gsac_gray_sup_rate,mua_gsac_gray_sup_time] = min(mu_gsac_gray_data(:,:,poss_lags),[],3);
mua_gsac_gray_sup_time = lags(mua_gsac_gray_sup_time + poss_lags(1) - 1);
[mua_gsac_gray_enh_rate,mua_gsac_gray_enh_time] = max(mu_gsac_gray_data(:,:,poss_lags),[],3);
mua_gsac_gray_enh_time = lags(mua_gsac_gray_enh_time + poss_lags(1) - 1);
mua_gsac_im_sup_rate = 1-mua_gsac_im_sup_rate; mua_gsac_im_enh_rate = mua_gsac_im_enh_rate - 1;
mua_gsac_gray_sup_rate = 1-mua_gsac_gray_sup_rate; mua_gsac_gray_enh_rate = mua_gsac_gray_enh_rate - 1;

[mua_msac_sup_rate,mua_msac_sup_time] = min(mu_msac_data(:,:,poss_lags),[],3);
mua_msac_sup_time = lags(mua_msac_sup_time + poss_lags(1) - 1);
[mua_msac_enh_rate,mua_msac_enh_time] = max(mu_msac_data(:,:,poss_lags),[],3);
mua_msac_enh_time = lags(mua_msac_enh_time + poss_lags(1) - 1);
mua_msac_sup_rate = 1-mua_msac_sup_rate; mua_msac_enh_rate = mua_msac_enh_rate - 1;

% poss_lags = find(lags > 0.0 & lags < 0.2);
[mua_simsac_sup_rate,mua_simsac_sup_time] = min(mu_simsac_data(:,:,poss_lags),[],3);
mua_simsac_sup_time = lags(mua_simsac_sup_time + poss_lags(1) - 1);
[mua_simsac_enh_rate,mua_simsac_enh_time] = max(mu_simsac_data(:,:,poss_lags),[],3);
mua_simsac_enh_time = lags(mua_simsac_enh_time + poss_lags(1) - 1);
mua_simsac_sup_rate = 1-mua_simsac_sup_rate; mua_simsac_enh_rate = mua_simsac_enh_rate -1;

strp_widths = [1 2 4];
subplot(3,1,1);hold on
errorbar(strp_widths,mean(mua_gsac_sup_rate),std(mua_gsac_sup_rate)/sqrt(length(cur_mu_set)),'o-');
errorbar(strp_widths,mean(mua_gsac_enh_rate),std(mua_gsac_enh_rate)/sqrt(length(cur_mu_set)),'ro-');
xlim([0.5 4.5]);
legend('Suppression','Enhancement'); legend('Boxoff');
xlabel('Stripe width (deg)','fontsize',14)
ylabel('Modulation index','fontsize',14)
set(gca,'fontsize',12); box off
title('Guided Saccades','fontsize',14);
subplot(3,1,2);hold on
errorbar(strp_widths,mean(mua_msac_sup_rate),std(mua_msac_sup_rate)/sqrt(length(cur_mu_set)),'o-');
errorbar(strp_widths,mean(mua_msac_enh_rate),std(mua_msac_enh_rate)/sqrt(length(cur_mu_set)),'ro-');
xlim([0.5 4.5]);
legend('Suppression','Enhancement'); legend('Boxoff');
xlabel('Stripe width (deg)','fontsize',14)
ylabel('Modulation index','fontsize',14)
set(gca,'fontsize',12); box off
title('Micro-Saccades','fontsize',14);
subplot(3,1,3);hold on
errorbar(strp_widths,mean(mua_simsac_sup_rate),std(mua_simsac_sup_rate)/sqrt(length(cur_mu_set)),'o-');
errorbar(strp_widths,mean(mua_simsac_enh_rate),std(mua_simsac_enh_rate)/sqrt(length(cur_mu_set)),'ro-');
xlim([0.5 4.5]);
legend('Suppression','Enhancement'); legend('Boxoff');
xlabel('Stripe width (deg)','fontsize',14)
ylabel('Modulation index','fontsize',14)
set(gca,'fontsize',12); box off
title('Simulated Saccades','fontsize',14);

fname = 'mu_stripewidth_modstrength';
fillPage(gcf,'papersize',[6 16]);
print(fname,'-dpdf','-painters');close


subplot(2,1,1);hold on
errorbar(strp_widths,mean(mua_gsac_gray_sup_rate),std(mua_gsac_gray_sup_rate)/sqrt(length(cur_mu_set)),'o-');
errorbar(strp_widths,mean(mua_gsac_im_sup_rate),std(mua_gsac_im_sup_rate)/sqrt(length(cur_mu_set)),'ro-');
xlim([0.5 4.5]);
legend('Gray back','Image back'); legend('Boxoff');
xlabel('Stripe width (deg)','fontsize',14)
ylabel('Modulation index','fontsize',14)
set(gca,'fontsize',12); box off
title('Suppression','fontsize',14);
subplot(2,1,2);hold on
errorbar(strp_widths,mean(mua_gsac_gray_enh_rate),std(mua_gsac_gray_enh_rate)/sqrt(length(cur_mu_set)),'o-');
errorbar(strp_widths,mean(mua_gsac_im_enh_rate),std(mua_gsac_im_enh_rate)/sqrt(length(cur_mu_set)),'ro-');
xlim([0.5 4.5]);
legend('Gray back','Image back'); legend('Boxoff');
xlabel('Stripe width (deg)','fontsize',14)
ylabel('Modulation index','fontsize',14)
set(gca,'fontsize',12); box off
title('Enhancement','fontsize',14);
fname = 'mu_stripewidth_gsac_backdep';
fillPage(gcf,'papersize',[6 8]);
print(fname,'-dpdf','-painters');close

%% view all expt-avg su sacmod
close all
n_un_wi = 3;
f1 = figure();
f2 = figure();
f3 = figure();
    cur_su_set = 1:length(all_sua_data);
    
    su_gsac_data = get_struct_data(all_sua_data,cur_su_set,'wi_gsac_avg');
    su_msac_data = get_struct_data(all_sua_data,cur_su_set,'wi_msac_avg');
    su_simsac_data = get_struct_data(all_sua_data,cur_su_set,'wi_simsac_avg');

    su_gsac_data = reshape(su_gsac_data,[length(cur_su_set) n_un_wi length(lags)]);
    su_msac_data = reshape(su_msac_data,[length(cur_su_set) n_un_wi length(lags)]);
    su_simsac_data = reshape(su_simsac_data,[length(cur_su_set) n_un_wi length(lags)]);
    
    all_cond_avgs = get_struct_data(all_sua_data,cur_su_set,'cond_avg_rate');
    su_gsac_data = bsxfun(@rdivide,su_gsac_data,all_cond_avgs);
    su_msac_data = bsxfun(@rdivide,su_msac_data,all_cond_avgs);
    su_simsac_data = bsxfun(@rdivide,su_simsac_data,all_cond_avgs);
    
    
    figure(f1);hold on; grid on
    shadedErrorBar(lags,squeeze(nanmean(su_gsac_data(:,1,:))),squeeze(nanstd(su_gsac_data(:,1,:)))/sqrt(length(cur_su_set)),{'color','r'});
    shadedErrorBar(lags,squeeze(nanmean(su_gsac_data(:,2,:))),squeeze(nanstd(su_gsac_data(:,2,:)))/sqrt(length(cur_su_set)),{'color','b'});
    shadedErrorBar(lags,squeeze(nanmean(su_gsac_data(:,3,:))),squeeze(nanstd(su_gsac_data(:,3,:)))/sqrt(length(cur_su_set)),{'color','k'});
    xlim([-0.2 0.5]);
    xl = xlim(); yl = ylim();
    line(xl,[1 1],'color','k','linestyle','--');
    line([0 0],yl,'color','k','linestyle','--');
    ylim([0.7 1.4]);

    figure(f2);hold on; grid on
    shadedErrorBar(lags,squeeze(nanmean(su_msac_data(:,1,:))),squeeze(nanstd(su_msac_data(:,1,:)))/sqrt(length(cur_su_set)),{'color','r'});
    shadedErrorBar(lags,squeeze(nanmean(su_msac_data(:,2,:))),squeeze(nanstd(su_msac_data(:,2,:)))/sqrt(length(cur_su_set)),{'color','b'});
    shadedErrorBar(lags,squeeze(nanmean(su_msac_data(:,3,:))),squeeze(nanstd(su_msac_data(:,3,:)))/sqrt(length(cur_su_set)),{'color','k'});
    xlim([-0.2 0.5]);
    xl = xlim(); yl = ylim();
    line(xl,[1 1],'color','k','linestyle','--');
    line([0 0],yl,'color','k','linestyle','--');
    ylim([0.7 1.3]);

    figure(f3);hold on; grid on
    shadedErrorBar(lags,squeeze(nanmean(su_simsac_data(:,1,:))),squeeze(nanstd(su_simsac_data(:,1,:)))/sqrt(length(cur_su_set)),{'color','r'});
    shadedErrorBar(lags,squeeze(nanmean(su_simsac_data(:,2,:))),squeeze(nanstd(su_simsac_data(:,2,:)))/sqrt(length(cur_su_set)),{'color','b'});
    shadedErrorBar(lags,squeeze(nanmean(su_simsac_data(:,3,:))),squeeze(nanstd(su_simsac_data(:,3,:)))/sqrt(length(cur_su_set)),{'color','k'});
    xlim([-0.2 0.5]);
    xl = xlim(); yl = ylim();
    line(xl,[1 1],'color','k','linestyle','--');
    line([0 0],yl,'color','k','linestyle','--');
    ylim([0.85 1.15]);

