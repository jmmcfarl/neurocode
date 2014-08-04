clear all
% close all

data_dir_base = '~/Analysis/bruce';
Expt_nums =  [86 88 89];
%% LOAD OVERALL SU DATA
load ~/Analysis/bruce/summary_analysis/su_data.mat
mahal_thresh = su_data.mah_thresh;

all_mua_data = [];
all_sua_data = [];
all_mu_probes = [];
all_mu_exnums = [];
all_su_probes = [];
all_su_exnums = [];
all_su_avgspks = [];
all_su_stdspks = [];
all_su_avgmahal = [];
for ex = 1:length(Expt_nums)
    fprintf('Loading data from expt %d\n',Expt_nums(ex));
    fname = [data_dir_base sprintf('/G0%d',Expt_nums(ex)) '/sac_mod/full_parorth_data.mat'];
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
bad_probes = [16 92];
close all
f1 = figure();
f2 = figure();
for ex = 1:length(Expt_nums)
    fprintf('Expt %d\n',Expt_nums(ex));
    cur_mu_set = find(all_mu_exnums == Expt_nums(ex));
    cur_su_set = find(all_su_exnums == Expt_nums(ex));
    if length(cur_mu_set) ~= 96
        error('Wrong mu Num');
    end
    exclude = [bad_probes; all_su_probes(cur_su_set)];
    cur_mu_set(exclude) = [];
    
    mu_p100_gsac_avgs = get_struct_data(all_mua_data,cur_mu_set,'p100_gsac_avg');
    mu_p30_gsac_avgs = get_struct_data(all_mua_data,cur_mu_set,'p30_gsac_avg');
    mu_o100_gsac_avgs = get_struct_data(all_mua_data,cur_mu_set,'o100_gsac_avg');
    mu_o30_gsac_avgs = get_struct_data(all_mua_data,cur_mu_set,'o30_gsac_avg');

    figure(f1);hold on; grid on
    shadedErrorBar(lags,nanmean(mu_p100_gsac_avgs),nanstd(mu_p100_gsac_avgs)/sqrt(length(cur_mu_set)),{'color','r'});
    shadedErrorBar(lags,nanmean(mu_p30_gsac_avgs),nanstd(mu_p30_gsac_avgs)/sqrt(length(cur_mu_set)),{'color','b'});
    shadedErrorBar(lags,nanmean(mu_o100_gsac_avgs),nanstd(mu_o100_gsac_avgs)/sqrt(length(cur_mu_set)),{'color','k'});
    shadedErrorBar(lags,nanmean(mu_o30_gsac_avgs),nanstd(mu_o30_gsac_avgs)/sqrt(length(cur_mu_set)),{'color','g'});
    xlim([-0.2 0.5]);
    xl = xlim(); yl = ylim();
    line(xl,[1 1],'color','k','linestyle','--');
    line([0 0],yl,'color','k','linestyle','--');
    ylim([0.7 1.4]);

    mu_p100_msac_avgs = get_struct_data(all_mua_data,cur_mu_set,'p100_msac_avg');
    mu_p30_msac_avgs = get_struct_data(all_mua_data,cur_mu_set,'p30_msac_avg');
    mu_o100_msac_avgs = get_struct_data(all_mua_data,cur_mu_set,'o100_msac_avg');
    mu_o30_msac_avgs = get_struct_data(all_mua_data,cur_mu_set,'o30_msac_avg');

    figure(f2);hold on; grid on
    shadedErrorBar(lags,nanmean(mu_p100_msac_avgs),nanstd(mu_p100_msac_avgs)/sqrt(length(cur_mu_set)),{'color','r'});
    shadedErrorBar(lags,nanmean(mu_p30_msac_avgs),nanstd(mu_p30_msac_avgs)/sqrt(length(cur_mu_set)),{'color','b'});
    shadedErrorBar(lags,nanmean(mu_o100_msac_avgs),nanstd(mu_o100_msac_avgs)/sqrt(length(cur_mu_set)),{'color','k'});
    shadedErrorBar(lags,nanmean(mu_o30_msac_avgs),nanstd(mu_o30_msac_avgs)/sqrt(length(cur_mu_set)),{'color','g'});
    xlim([-0.2 0.5]);
    xl = xlim(); yl = ylim();
    line(xl,[1 1],'color','k','linestyle','--');
    line([0 0],yl,'color','k','linestyle','--');
    ylim([0.7 1.2]);
   
    pause
    figure(f1); clf;
    figure(f2); clf;
end

%% view all expt-avg mu sacmod
poss_lags = find(lags > 0 & lags < 0.25);

bad_probes = [16; 92];
close all
f1 = figure();
f2 = figure();
for ex = 1:length(Expt_nums)
    fprintf('Expt %d\n',Expt_nums(ex));
    cur_mu_set = find(all_mu_exnums == Expt_nums(ex));
    cur_su_set = find(all_su_exnums == Expt_nums(ex));
    if length(cur_mu_set) ~= 96
        error('Wrong mu Num');
    end
    exclude = [bad_probes; all_su_probes(cur_su_set)];
    cur_mu_set(exclude) = [];
    
    mu_p100_gsac_avgs = get_struct_data(all_mua_data,cur_mu_set,'p100_gsac_avg');
    mu_p30_gsac_avgs = get_struct_data(all_mua_data,cur_mu_set,'p30_gsac_avg');
    mu_o100_gsac_avgs = get_struct_data(all_mua_data,cur_mu_set,'o100_gsac_avg');
    mu_o30_gsac_avgs = get_struct_data(all_mua_data,cur_mu_set,'o30_gsac_avg');
    
    all_cond_avgs = get_struct_data(all_mua_data,cur_mu_set,'cond_avg');
    mu_p100_gsac_avgs = bsxfun(@rdivide,mu_p100_gsac_avgs,all_cond_avgs(:,1));
    mu_p30_gsac_avgs = bsxfun(@rdivide,mu_p30_gsac_avgs,all_cond_avgs(:,2));
    mu_o100_gsac_avgs = bsxfun(@rdivide,mu_o100_gsac_avgs,all_cond_avgs(:,3));
    mu_o30_gsac_avgs = bsxfun(@rdivide,mu_o30_gsac_avgs,all_cond_avgs(:,4));
    
    
    mua_gsac_p100_sup_rate = 1-min(mu_p100_gsac_avgs(:,poss_lags),[],2);
    mua_gsac_p100_enh_rate = max(mu_p100_gsac_avgs(:,poss_lags),[],2)-1;
    mua_gsac_p30_sup_rate = 1-min(mu_p30_gsac_avgs(:,poss_lags),[],2);
    mua_gsac_p30_enh_rate = max(mu_p30_gsac_avgs(:,poss_lags),[],2)-1;
    mua_gsac_o100_sup_rate = 1-min(mu_o100_gsac_avgs(:,poss_lags),[],2);
    mua_gsac_o100_enh_rate = max(mu_o100_gsac_avgs(:,poss_lags),[],2)-1;
    mua_gsac_o30_sup_rate = 1-min(mu_o30_gsac_avgs(:,poss_lags),[],2);
    mua_gsac_o30_enh_rate = max(mu_o30_gsac_avgs(:,poss_lags),[],2)-1;
    
    mua_msac_p100_sup_rate = 1-min(mu_p100_msac_avgs(:,poss_lags),[],2);
    mua_msac_p100_enh_rate = max(mu_p100_msac_avgs(:,poss_lags),[],2)-1;
    mua_msac_p30_sup_rate = 1-min(mu_p30_msac_avgs(:,poss_lags),[],2);
    mua_msac_p30_enh_rate = max(mu_p30_msac_avgs(:,poss_lags),[],2)-1;
    mua_msac_o100_sup_rate = 1-min(mu_o100_msac_avgs(:,poss_lags),[],2);
    mua_msac_o100_enh_rate = max(mu_o100_msac_avgs(:,poss_lags),[],2)-1;
    mua_msac_o30_sup_rate = 1-min(mu_o30_msac_avgs(:,poss_lags),[],2);
    mua_msac_o30_enh_rate = max(mu_o30_msac_avgs(:,poss_lags),[],2)-1;

    figure(f1);hold on; grid on
    shadedErrorBar(lags,nanmean(mu_p100_gsac_avgs),nanstd(mu_p100_gsac_avgs)/sqrt(length(cur_mu_set)),{'color','r'});
    shadedErrorBar(lags,nanmean(mu_p30_gsac_avgs),nanstd(mu_p30_gsac_avgs)/sqrt(length(cur_mu_set)),{'color','b'});
    shadedErrorBar(lags,nanmean(mu_o100_gsac_avgs),nanstd(mu_o100_gsac_avgs)/sqrt(length(cur_mu_set)),{'color','k'});
    shadedErrorBar(lags,nanmean(mu_o30_gsac_avgs),nanstd(mu_o30_gsac_avgs)/sqrt(length(cur_mu_set)),{'color','g'});
    xlim([-0.2 0.5]);
    xl = xlim(); yl = ylim();
    line(xl,[1 1],'color','k','linestyle','--');
    line([0 0],yl,'color','k','linestyle','--');
    ylim([0.7 1.4]);

    mu_p100_msac_avgs = get_struct_data(all_mua_data,cur_mu_set,'p100_msac_avg');
    mu_p30_msac_avgs = get_struct_data(all_mua_data,cur_mu_set,'p30_msac_avg');
    mu_o100_msac_avgs = get_struct_data(all_mua_data,cur_mu_set,'o100_msac_avg');
    mu_o30_msac_avgs = get_struct_data(all_mua_data,cur_mu_set,'o30_msac_avg');
    mu_p100_msac_avgs = bsxfun(@rdivide,mu_p100_msac_avgs,all_cond_avgs(:,1));
    mu_p30_msac_avgs = bsxfun(@rdivide,mu_p30_msac_avgs,all_cond_avgs(:,2));
     mu_o100_msac_avgs = bsxfun(@rdivide,mu_o100_msac_avgs,all_cond_avgs(:,3));
    mu_o30_msac_avgs = bsxfun(@rdivide,mu_o30_msac_avgs,all_cond_avgs(:,4));

    figure(f2);hold on; grid on
    shadedErrorBar(lags,nanmean(mu_p100_msac_avgs),nanstd(mu_p100_msac_avgs)/sqrt(length(cur_mu_set)),{'color','r'});
    shadedErrorBar(lags,nanmean(mu_p30_msac_avgs),nanstd(mu_p30_msac_avgs)/sqrt(length(cur_mu_set)),{'color','b'});
    shadedErrorBar(lags,nanmean(mu_o100_msac_avgs),nanstd(mu_o100_msac_avgs)/sqrt(length(cur_mu_set)),{'color','k'});
    shadedErrorBar(lags,nanmean(mu_o30_msac_avgs),nanstd(mu_o30_msac_avgs)/sqrt(length(cur_mu_set)),{'color','g'});
    xlim([-0.2 0.5]);
    xl = xlim(); yl = ylim();
    line(xl,[1 1],'color','k','linestyle','--');
    line([0 0],yl,'color','k','linestyle','--');
    ylim([0.7 1.2]);
   
    pause
    figure(f1); clf;
    figure(f2); clf;
end

%% PLOT INDIVIDUAL EXPT
cd ~/Analysis/bruce/summary_analysis/

poss_lags = find(lags > 0 & lags < 0.25);
bad_probes = [16; 92];
close all

ex = 3;
fprintf('Expt %d\n',Expt_nums(ex));
cur_mu_set = find(all_mu_exnums == Expt_nums(ex));
cur_su_set = find(all_su_exnums == Expt_nums(ex));
if length(cur_mu_set) ~= 96
    error('Wrong mu Num');
end
exclude = [bad_probes; all_su_probes(cur_su_set)];
cur_mu_set(exclude) = [];

mu_p100_gsac_avgs = get_struct_data(all_mua_data,cur_mu_set,'p100_gsac_avg');
mu_p30_gsac_avgs = get_struct_data(all_mua_data,cur_mu_set,'p30_gsac_avg');
mu_o100_gsac_avgs = get_struct_data(all_mua_data,cur_mu_set,'o100_gsac_avg');
mu_o30_gsac_avgs = get_struct_data(all_mua_data,cur_mu_set,'o30_gsac_avg');
mu_p100_gsac_avgs = bsxfun(@rdivide,mu_p100_gsac_avgs,all_cond_avgs(:,1));
mu_p30_gsac_avgs = bsxfun(@rdivide,mu_p30_gsac_avgs,all_cond_avgs(:,2));
mu_o100_gsac_avgs = bsxfun(@rdivide,mu_o100_gsac_avgs,all_cond_avgs(:,3));
mu_o30_gsac_avgs = bsxfun(@rdivide,mu_o30_gsac_avgs,all_cond_avgs(:,4));

mu_p100_msac_avgs = get_struct_data(all_mua_data,cur_mu_set,'p100_msac_avg');
mu_p30_msac_avgs = get_struct_data(all_mua_data,cur_mu_set,'p30_msac_avg');
mu_o100_msac_avgs = get_struct_data(all_mua_data,cur_mu_set,'o100_msac_avg');
mu_o30_msac_avgs = get_struct_data(all_mua_data,cur_mu_set,'o30_msac_avg');
mu_p100_msac_avgs = bsxfun(@rdivide,mu_p100_msac_avgs,all_cond_avgs(:,1));
mu_p30_msac_avgs = bsxfun(@rdivide,mu_p30_msac_avgs,all_cond_avgs(:,2));
mu_o100_msac_avgs = bsxfun(@rdivide,mu_o100_msac_avgs,all_cond_avgs(:,3));
mu_o30_msac_avgs = bsxfun(@rdivide,mu_o30_msac_avgs,all_cond_avgs(:,4));



figure(f1);
subplot(2,1,1);hold on; grid on
h1 = shadedErrorBar(lags,nanmean(mu_p100_gsac_avgs),nanstd(mu_p100_gsac_avgs)/sqrt(length(cur_mu_set)),{'color','r'});
h2 = shadedErrorBar(lags,nanmean(mu_p30_gsac_avgs),nanstd(mu_p30_gsac_avgs)/sqrt(length(cur_mu_set)),{'color','b'});
h3 = shadedErrorBar(lags,nanmean(mu_o100_gsac_avgs),nanstd(mu_o100_gsac_avgs)/sqrt(length(cur_mu_set)),{'color','k'});
h4 = shadedErrorBar(lags,nanmean(mu_o30_gsac_avgs),nanstd(mu_o30_gsac_avgs)/sqrt(length(cur_mu_set)),{'color','g'});
legend([h1.mainLine h2.mainLine h3.mainLine h4.mainLine],{'Parallel 100','Parallel 30','Orthoganol 100','Orthoganol 30'});legend('boxoff');
xlim([-0.2 0.5]);
ylim([0.7 1.4]);
xl = xlim(); yl = ylim();
line(xl,[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
set(gca,'fontsize',12); box off
xlabel('Time since saccade onset (s)','fontsize',14)
ylabel('Relative firing rate','fontsize',14)

subplot(2,1,2);hold on; grid on
h1 = shadedErrorBar(lags,nanmean(mu_p100_msac_avgs),nanstd(mu_p100_msac_avgs)/sqrt(length(cur_mu_set)),{'color','r'});
h2 = shadedErrorBar(lags,nanmean(mu_p30_msac_avgs),nanstd(mu_p30_msac_avgs)/sqrt(length(cur_mu_set)),{'color','b'});
h3 = shadedErrorBar(lags,nanmean(mu_o100_msac_avgs),nanstd(mu_o100_msac_avgs)/sqrt(length(cur_mu_set)),{'color','k'});
h4 = shadedErrorBar(lags,nanmean(mu_o30_msac_avgs),nanstd(mu_o30_msac_avgs)/sqrt(length(cur_mu_set)),{'color','g'});
legend([h1.mainLine h2.mainLine h3.mainLine h4.mainLine],{'Parallel 100','Parallel 30','Orthoganol 100','Orthoganol 30'});legend('boxoff');
xlim([-0.2 0.5]);
ylim([0.7 1.4]);
xl = xlim(); yl = ylim();
line(xl,[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
set(gca,'fontsize',12); box off
xlabel('Time since saccade onset (s)','fontsize',14)
ylabel('Relative firing rate','fontsize',14)

fname = 'mu_parorth_avgs';
fillPage(gcf,'papersize',[6 10]);
print(fname,'-dpdf','-painters');close


%% view all individual su sacmod
close all
n_sus = length(all_sua_data);
f1 = figure();
f2 = figure();
for ss = 1:n_sus
    fprintf('Expt %d probe %d\n',all_su_exnums(ss),all_su_probes(ss));
%     fprintf('Avg rate:%.3f Nspikes: %d\n',all_sua_data(ss).avg_rate/dt,all_sua_data(ss).tot_nspikes);
%     fprintf('%d gsacs,  %d msacs,  %d simsacs \n',all_sua_data(ss).nused_gsac,all_sua_data(ss).nused_msac,all_sua_data(ss).nused_simsac);
    
    figure(f1);hold on; grid on
    shadedErrorBar(lags,all_sua_data(ss).gsac_avg,all_sua_data(ss).gsac_sem,{'color','r'});
    shadedErrorBar(lags,all_sua_data(ss).msac_avg,all_sua_data(ss).msac_sem,{'color','b'});
    shadedErrorBar(lags,all_sua_data(ss).simsac_avg,all_sua_data(ss).simsac_sem,{'color','k'});
    xlim([-0.2 0.5]);
    xl = xlim(); yl = ylim();
    line(xl,[1 1],'color','k','linestyle','--');
    line([0 0],yl,'color','k','linestyle','--');
    
    figure(f2)
    subplot(2,1,1);hold on; grid on
    shadedErrorBar(lags,all_sua_data(ss).gsac_gray_avg,all_sua_data(ss).gsac_gray_sem,{'color','r'});
    shadedErrorBar(lags,all_sua_data(ss).gsac_im_avg,all_sua_data(ss).gsac_im_sem,{'color','b'});
    xlim([-0.2 0.5]);
    xl = xlim(); yl = ylim();
    line(xl,[1 1],'color','k','linestyle','--');
    line([0 0],yl,'color','k','linestyle','--');
    subplot(2,1,2);hold on; grid on
    shadedErrorBar(lags,all_sua_data(ss).msac_gray_avg,all_sua_data(ss).msac_gray_sem,{'color','r'});
    shadedErrorBar(lags,all_sua_data(ss).msac_im_avg,all_sua_data(ss).msac_im_sem,{'color','b'});
    xlim([-0.2 0.5]);
    xl = xlim(); yl = ylim();
    line(xl,[1 1],'color','k','linestyle','--');
    line([0 0],yl,'color','k','linestyle','--');

    pause
    figure(f1);clf;
    figure(f2);clf;
end

%% average of SU data
n_sus = length(all_sua_data);
su_avgrate = get_struct_data(all_sua_data,1:n_sus,'avg_rate')/dt;
su_nspks = get_struct_data(all_sua_data,1:n_sus,'tot_nspikes');
su_ngsac = get_struct_data(all_sua_data,1:n_sus,'nused_gsac');
su_nmsac = get_struct_data(all_sua_data,1:n_sus,'nused_msac');
su_nsimsac = get_struct_data(all_sua_data,1:n_sus,'nused_simsac');

usable_sus = find(su_avgrate > 1 & su_nspks > 500 & all_su_avgmahal > 2);
used_sus = usable_sus(su_ngsac(usable_sus) > 50 & su_nmsac(usable_sus) > 50 & su_nsimsac(usable_sus) > 50);
su_gsac = get_struct_data(all_sua_data,used_sus,'gsac_avg');
su_msac = get_struct_data(all_sua_data,used_sus,'msac_avg');
su_simsac = get_struct_data(all_sua_data,used_sus,'simsac_avg');
figure
hold on
shadedErrorBar(lags,mean(su_gsac),std(su_gsac)/sqrt(length(used_sus)),{'color','r'});
shadedErrorBar(lags,mean(su_msac),std(su_msac)/sqrt(length(used_sus)),{'color','b'});
shadedErrorBar(lags,mean(su_simsac),std(su_simsac)/sqrt(length(used_sus)),{'color','k'});
xlim([-0.2 0.5]);
xl = xlim(); yl = ylim();
line(xl,[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');

%% compare individual SUS with averages
close all
used_sus = usable_sus(su_ngsac(usable_sus) > 100 & su_nmsac(usable_sus) > 100 & su_nsimsac(usable_sus) > 100);
su_gsac = get_struct_data(all_sua_data,used_sus,'gsac_avg');
su_msac = get_struct_data(all_sua_data,used_sus,'msac_avg');
su_simsac = get_struct_data(all_sua_data,used_sus,'simsac_avg');
figure
for ii = 1:length(used_sus)
    ii
subplot(1,3,1)
plot(lags,su_gsac,'k')
hold on
shadedErrorBar(lags,mean(su_gsac),std(su_gsac)/sqrt(length(used_sus)),{'color','b'});
plot(lags,su_gsac(ii,:),'r','linewidth',4)
xlim([-0.2 0.5]);
xl = xlim(); yl = ylim();
line(xl,[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
subplot(1,3,2)
plot(lags,su_msac,'k')
hold on
shadedErrorBar(lags,mean(su_msac),std(su_msac)/sqrt(length(used_sus)),{'color','b'});
plot(lags,su_msac(ii,:),'r','linewidth',4)
xlim([-0.2 0.5]);
xl = xlim(); yl = ylim();
line(xl,[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
subplot(1,3,3)
plot(lags,su_simsac,'k')
hold on
shadedErrorBar(lags,mean(su_simsac),std(su_simsac)/sqrt(length(used_sus)),{'color','b'});
plot(lags,su_simsac(ii,:),'r','linewidth',4)
xlim([-0.2 0.5]);
ylim([0.4 1.6])
xl = xlim(); yl = ylim();
line(xl,[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');

pause
clf
end


%% view all individual MUA
close all
for ex = 1:length(Expt_nums)
    fprintf('Expt %d\n',Expt_nums(ex));
    cur_mu_set = find(all_mu_exnums == Expt_nums(ex));
    cur_su_set = find(all_su_exnums == Expt_nums(ex));
     exclude = [bad_probes; all_su_probes(cur_su_set)];
    cur_mu_set(exclude) = [];

    mu_gsac_avgs = get_struct_data(all_mua_data,cur_mu_set,'gsac_avg');
    mu_msac_avgs = get_struct_data(all_mua_data,cur_mu_set,'msac_avg');
    mu_simsac_avgs = get_struct_data(all_mua_data,cur_mu_set,'simsac_avg');
    
    for cc = 1:length(cur_mu_set)
        fprintf('Expt %d, MU %d\n',Expt_nums(ex),cur_mu_set(cc) - 96*(ex-1));
        plot(lags,mu_gsac_avgs(cc,:),'r');
    hold on; grid on
        plot(lags,mu_msac_avgs(cc,:),'b');
        plot(lags,mu_simsac_avgs(cc,:),'k');
        xlim([-0.2 0.5]);
        xl = xlim(); yl = ylim();
        line(xl,[1 1],'color','k','linestyle','--');
        line([0 0],yl,'color','k','linestyle','--');
        ylim([0.7 1.4]);
        pause
        clf
    end
end