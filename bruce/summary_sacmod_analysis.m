clear all
% close all
clc

data_dir_base = '~/Analysis/bruce';
Expt_nums = [81 85 86 87 88 89 91 92 93 95];
% Expt_nums = [95];

bad_probes = [16; 92];
%% LOAD OVERALL SU DATA
load ~/Analysis/bruce/summary_analysis/su_data.mat
mahal_thresh = su_data.mah_thresh;

all_mua_data = [];
all_sua_data = [];
all_mu_probes = [];
all_mu_canuse = [];
all_mu_exnums = [];
all_su_probes = [];
all_su_exnums = [];
all_su_avgspks = [];
all_su_stdspks = [];
all_su_avgmahal = [];
all_gen_data = [];
for ex = 1:length(Expt_nums)
    fprintf('Loading data from expt %d\n',Expt_nums(ex));
    fname = [data_dir_base sprintf('/G0%d',Expt_nums(ex)) '/sac_mod/full_sacmod_data.mat'];
    load(fname);
    
    all_mua_data = cat(1,all_mua_data,mua_data');
    all_mu_probes = cat(1,all_mu_probes,(1:96)');
    all_mu_exnums = cat(1,all_mu_exnums,Expt_nums(ex)*ones(96,1));
    cur_canuse = true(96,1);
    cur_canuse(su_probes) = false;
    cur_canuse(bad_probes) = false;
    all_mu_canuse = cat(1,all_mu_canuse,cur_canuse);
    
    use_sua = find(get_struct_data(sua_data,1:length(sua_data),'used'));
    all_sua_data = cat(1,all_sua_data,sua_data(use_sua)');
    all_su_probes = cat(1,all_su_probes,su_probes(use_sua)');
    all_su_exnums = cat(1,all_su_exnums,Expt_nums(ex)*ones(length(use_sua),1));
    
    cur_su_avgmahal = su_data.avg_mahal(ex,su_probes(use_sua));
    all_su_avgmahal = [all_su_avgmahal; cur_su_avgmahal'];
    
    all_gen_data = cat(1,all_gen_data,gen_data');
    
%     fname = ['~/Data/bruce' sprintf('/G0%d', Expt_nums(ex)) '/fin_aclust_data.mat'];
%     load(fname);
%     cur_used = arrayfun(@(x) length(x.probenum),autoclust);
%     temp = find(cur_used(:,1) > 0,1,'first');
%     aclust_probenums = [autoclust(temp,:).probenum];
%     autoclust = autoclust(:,ismember(aclust_probenums,su_probes));
%     su_avg_wvfrms = zeros(length(su_probes),40);
%     su_std_wvfrms = zeros(length(su_probes),40);
%     for pp = 1:length(su_probes)
%         cur_avg_wvfrms = [];
%         for bb = 1:size(autoclust,1)
%             if autoclust(bb,pp).mahal_d > mahal_thresh | autoclust(bb,pp).man_code == 4
%                 cur_avg_wvfrms = [cur_avg_wvfrms autoclust(bb,pp).avg_wvfrm(:,1)];
%             end
%         end
%         su_avg_wvfrms(pp,:) = mean(cur_avg_wvfrms,2);
%         su_std_wvfrms(pp,:) = std(cur_avg_wvfrms,[],2);
%     end
%     all_su_avgspks = cat(1,all_su_avgspks,su_avg_wvfrms);
%     all_su_stdspks = cat(1,all_su_stdspks,su_std_wvfrms);
end
lags = anal_params.lags*anal_params.dt;
dt = anal_params.dt;
sac_bin_cents = anal_params.sac_bin_cents;
trial_lags = anal_params.trial_lags*dt;
acorr_lags = anal_params.acorr_lags*dt;
%% EXTRACT SAC-MOD STATS
poss_lags = find(lags > 0 & lags < 0.2);
poss_presac_lags = find(lags > -0.1 & lags < 0.03);

all_mua_gsac = get_struct_data(all_mua_data,find(all_mu_canuse),'gsac_avg');
all_mua_msac = get_struct_data(all_mua_data,find(all_mu_canuse),'msac_avg');
all_mua_simsac = get_struct_data(all_mua_data,find(all_mu_canuse),'simsac_avg');

%FOR MU
[mua_gsac_sup_rate,mua_gsac_sup_time] = min(all_mua_gsac(:,poss_lags),[],2); 
mua_gsac_sup_time = lags(mua_gsac_sup_time + poss_lags(1) - 1);
[mua_gsac_enh_rate,mua_gsac_enh_time] = max(all_mua_gsac(:,poss_lags),[],2); 
mua_gsac_enh_time = lags(mua_gsac_enh_time + poss_lags(1) - 1);
mua_gsac_sup_rate = 1-mua_gsac_sup_rate; mua_gsac_enh_rate = mua_gsac_enh_rate - 1;
[mua_gsac_presup_rate,mua_gsac_presup_ind] = min(all_mua_gsac(:,poss_presac_lags),[],2); 
mua_gsac_presup_ind = mua_gsac_presup_ind + poss_presac_lags(1)-1;
mua_gsac_presup_time = lags(mua_gsac_presup_ind);

[mua_msac_sup_rate,mua_msac_sup_time] = min(all_mua_msac(:,poss_lags),[],2); 
mua_msac_sup_time = lags(mua_msac_sup_time + poss_lags(1) - 1);
[mua_msac_enh_rate,mua_msac_enh_time] = max(all_mua_msac(:,poss_lags),[],2); 
mua_msac_enh_time = lags(mua_msac_enh_time + poss_lags(1) - 1);
mua_msac_sup_rate = 1-mua_msac_sup_rate; mua_msac_enh_rate = mua_msac_enh_rate - 1;
[mua_msac_presup_rate,mua_msac_presup_ind] = min(all_mua_msac(:,poss_presac_lags),[],2); 
mua_msac_presup_ind = mua_msac_presup_ind + poss_presac_lags(1)-1;
mua_msac_presup_time = lags(mua_msac_presup_ind);

[mua_simsac_sup_rate,mua_simsac_sup_time] = min(all_mua_simsac(:,poss_lags),[],2); 
mua_simsac_sup_time = lags(mua_simsac_sup_time + poss_lags(1) - 1);
[mua_simsac_enh_rate,mua_simsac_enh_time] = max(all_mua_simsac(:,poss_lags),[],2); 
mua_simsac_enh_time = lags(mua_simsac_enh_time + poss_lags(1) - 1);
mua_simsac_sup_rate = 1-mua_simsac_sup_rate; mua_simsac_enh_rate = mua_simsac_enh_rate - 1;

%FOR SU
all_sua_gsac = get_struct_data(all_sua_data,1:length(all_sua_data),'gsac_avg');
all_sua_msac = get_struct_data(all_sua_data,1:length(all_sua_data),'msac_avg');
all_sua_simsac = get_struct_data(all_sua_data,1:length(all_sua_data),'simsac_avg');

[sua_gsac_sup_rate,sua_gsac_sup_ind] = min(all_sua_gsac(:,poss_lags),[],2); 
sua_gsac_sup_ind = sua_gsac_sup_ind + poss_lags(1) - 1;
sua_gsac_sup_time = lags(sua_gsac_sup_ind);
[sua_gsac_enh_rate,sua_gsac_enh_ind] = max(all_sua_gsac(:,poss_lags),[],2);
sua_gsac_enh_ind = sua_gsac_enh_ind + poss_lags(1) - 1;
sua_gsac_enh_time = lags(sua_gsac_enh_ind);
sua_gsac_sup_rate = 1-sua_gsac_sup_rate; sua_gsac_enh_rate = sua_gsac_enh_rate - 1;
[sua_gsac_presup_rate,sua_gsac_presup_ind] = min(all_sua_gsac(:,poss_presac_lags),[],2); 
sua_gsac_presup_ind = sua_gsac_presup_ind + poss_presac_lags(1)-1;
sua_gsac_presup_time = lags(sua_gsac_presup_ind);

[sua_msac_sup_rate,sua_msac_sup_ind] = min(all_sua_msac(:,poss_lags),[],2); 
sua_msac_sup_ind = sua_msac_sup_ind + poss_lags(1) - 1;
sua_msac_sup_time = lags(sua_msac_sup_ind);
[sua_msac_enh_rate,sua_msac_enh_ind] = max(all_sua_msac(:,poss_lags),[],2); 
sua_msac_enh_ind = sua_msac_enh_ind + poss_lags(1) - 1;
sua_msac_enh_time = lags(sua_msac_enh_ind);
sua_msac_sup_rate = 1-sua_msac_sup_rate; sua_msac_enh_rate = sua_msac_enh_rate - 1;
[sua_msac_presup_rate,sua_msac_presup_ind] = min(all_sua_msac(:,poss_presac_lags),[],2); 
sua_msac_presup_ind = sua_msac_presup_ind + poss_presac_lags(1)-1;
sua_msac_presup_time = lags(sua_gsac_presup_ind);

[sua_simsac_sup_rate,sua_simsac_sup_ind] = min(all_sua_simsac(:,poss_lags),[],2); 
sua_simsac_sup_ind = sua_simsac_sup_ind + poss_lags(1) - 1;
sua_simsac_sup_time = lags(sua_simsac_sup_ind);
[sua_simsac_enh_rate,sua_simsac_enh_ind] = max(all_sua_simsac(:,poss_lags),[],2);
sua_simsac_enh_ind = sua_simsac_enh_ind + poss_lags(1) - 1;
sua_simsac_enh_time = lags(sua_simsac_enh_ind);
sua_simsac_sup_rate = 1-sua_simsac_sup_rate; sua_simsac_enh_rate = sua_simsac_enh_rate - 1;

n_sus = length(all_sua_data);
su_avgrate = get_struct_data(all_sua_data,1:n_sus,'avg_rate')/dt;
su_nspks = get_struct_data(all_sua_data,1:n_sus,'tot_nspikes');
su_ngsac = get_struct_data(all_sua_data,1:n_sus,'nused_gsac');
su_nmsac = get_struct_data(all_sua_data,1:n_sus,'nused_msac');
su_nsimsac = get_struct_data(all_sua_data,1:n_sus,'nused_simsac');

%which are non SU MU probes
all_use_mu = [];
for ex = 1:length(Expt_nums)
    cur_mu_set = find(all_mu_exnums == Expt_nums(ex));
    cur_su_set = find(all_su_exnums == Expt_nums(ex));
    exclude = [bad_probes; all_su_probes(cur_su_set)];
    cur_mu_set(exclude) = [];
    all_use_mu = [all_use_mu; cur_mu_set];
end

usable_sus = find(su_avgrate > 5 & su_nspks > 500 & all_su_avgmahal > 2.25);

%% SCATTERPLOTS
print_on = 0;
cd ~/Analysis/bruce/summary_analysis/

%SUP VS ENH MAG 
figure; hold on
plot(mua_gsac_enh_rate,mua_gsac_sup_rate,'.')
plot(sua_gsac_enh_rate(usable_sus),sua_gsac_sup_rate(usable_sus),'ro','markersize',4,'linewidth',1.5)
legend('MUs','SUs'); legend('boxoff');
xlabel('Enhancement index','fontsize',14)
ylabel('Suppression index','fontsize',14)
set(gca,'fontsize',12);
if print_on == 1
    fname = 'scatter_supvsenh';
    fillPage(gcf,'papersize',[5 5]);
    print(fname,'-dpdf','-painters');close
end

figure; hold on
plot(mua_gsac_enh_time*1e3,mua_gsac_sup_time*1e3,'.')
plot(sua_gsac_enh_time(usable_sus)*1e3,sua_gsac_sup_time(usable_sus)*1e3,'ro','markersize',4,'linewidth',1.5)
legend('MUs','SUs'); legend('boxoff');
xlabel('Enhancement time (ms)','fontsize',14)
ylabel('Suppression time (ms)','fontsize',14)
set(gca,'fontsize',12);
if print_on == 1
    fname = 'scatter_supvsenh_timing';
    fillPage(gcf,'papersize',[5 5]);
    print(fname,'-dpdf','-painters');close
end

figure;
subplot(2,2,1); hold on
plot(mua_gsac_enh_rate,mua_msac_enh_rate,'.')
plot(sua_gsac_enh_rate(usable_sus),sua_msac_enh_rate(usable_sus),'ro','markersize',4,'linewidth',1.5)
legend('MUs','SUs','Location','northwest'); legend('boxoff');
xlabel('Guided sac index','fontsize',14)
ylabel('Micro-sac index','fontsize',14)
xl = xlim(); yl = xlim();
line(xl,yl,'color','k');
set(gca,'fontsize',12);
title('Enhancement','fontsize',14);
subplot(2,2,2); hold on
plot(mua_gsac_sup_rate,mua_msac_sup_rate,'.')
plot(sua_gsac_sup_rate(usable_sus),sua_msac_sup_rate(usable_sus),'ro','markersize',4,'linewidth',1.5)
legend('MUs','SUs','Location','northwest'); legend('boxoff');
xlabel('Guided sac index','fontsize',14)
ylabel('Micro-sac index','fontsize',14)
xl = xlim(); yl = xlim();
line(xl,yl,'color','k');
set(gca,'fontsize',12);
title('Suppression','fontsize',14);
subplot(2,2,3); hold on
plot(mua_gsac_enh_time*1e3,mua_msac_enh_time*1e3,'.')
plot(sua_gsac_enh_time(usable_sus)*1e3,sua_msac_enh_time(usable_sus)*1e3,'ro','markersize',4,'linewidth',1.5)
legend('MUs','SUs','Location','northwest'); legend('boxoff');
xlabel('Guided sac time','fontsize',14)
ylabel('Micro-sac time','fontsize',14)
xl = xlim(); yl = xlim();
line(xl,yl,'color','k');
set(gca,'fontsize',12);
title('Enhancement','fontsize',14);
subplot(2,2,4); hold on
plot(mua_gsac_sup_time*1e3,mua_msac_sup_time*1e3,'.')
plot(sua_gsac_sup_time(usable_sus)*1e3,sua_msac_sup_time(usable_sus)*1e3,'ro','markersize',4,'linewidth',1.5)
legend('MUs','SUs','Location','northwest'); legend('boxoff');
xlabel('Guided sac time','fontsize',14)
ylabel('Micro-sac time','fontsize',14)
xl = xlim(); yl = xlim();
line(xl,yl,'color','k');
set(gca,'fontsize',12);
title('Suppression','fontsize',14);
if print_on == 1
    fname = 'scatter_gsac_vs_msac';
    fillPage(gcf,'papersize',[10 10]);
    print(fname,'-dpdf','-painters');close
end


figure;
subplot(2,2,1); hold on
plot(mua_gsac_enh_rate,mua_simsac_enh_rate,'.')
plot(sua_gsac_enh_rate(usable_sus),sua_simsac_enh_rate(usable_sus),'ro','markersize',4,'linewidth',1.5)
legend('MUs','SUs','Location','northwest'); legend('boxoff');
xlabel('Guided sac index','fontsize',14)
ylabel('Sim-sac index','fontsize',14)
xl = xlim(); yl = xlim();
line(xl,yl,'color','k');
set(gca,'fontsize',12);
title('Enhancement','fontsize',14);
subplot(2,2,2); hold on
plot(mua_gsac_sup_rate,mua_simsac_sup_rate,'.')
plot(sua_gsac_sup_rate(usable_sus),sua_simsac_sup_rate(usable_sus),'ro','markersize',4,'linewidth',1.5)
legend('MUs','SUs','Location','northwest'); legend('boxoff');
xlabel('Guided sac index','fontsize',14)
ylabel('Sim-sac index','fontsize',14)
xl = xlim(); yl = xlim();
line(xl,yl,'color','k');
set(gca,'fontsize',12);
title('Suppression','fontsize',14);
subplot(2,2,3); hold on
plot(mua_gsac_enh_time*1e3,mua_simsac_enh_time*1e3,'.')
plot(sua_gsac_enh_time(usable_sus)*1e3,sua_simsac_enh_time(usable_sus)*1e3,'ro','markersize',4,'linewidth',1.5)
legend('MUs','SUs','Location','northwest'); legend('boxoff');
xlabel('Guided sac time','fontsize',14)
ylabel('Sim-sac time','fontsize',14)
xl = xlim(); yl = xlim();
line(xl,yl,'color','k');
set(gca,'fontsize',12);
title('Enhancement','fontsize',14);
subplot(2,2,4); hold on
plot(mua_gsac_sup_time*1e3,mua_simsac_sup_time*1e3,'.')
plot(sua_gsac_sup_time(usable_sus)*1e3,sua_simsac_sup_time(usable_sus)*1e3,'ro','markersize',4,'linewidth',1.5)
legend('MUs','SUs','Location','northwest'); legend('boxoff');
xlabel('Guided sac time','fontsize',14)
ylabel('Sim-sac time','fontsize',14)
xl = xlim(); yl = xlim();
line(xl,yl,'color','k');
set(gca,'fontsize',12);
title('Suppression','fontsize',14);
if print_on == 1
    fname = 'scatter_gsac_vs_simsac';
    fillPage(gcf,'papersize',[10 10]);
    print(fname,'-dpdf','-painters');close
end

%% PLOT INDIVIDUAL SU EXAMPLE
su_num = 5; % [5 32 42 53]  [3 17 20 23 29 34 47 55]
print_on = 0;
cd ~/Analysis/bruce/summary_analysis/

for su_num = [60:69]
% for su_num = usable_sus'
    mu_gsac = get_struct_data(all_mua_data,all_use_mu,'gsac_avg');
    
    figure;hold on; grid on
    h1 = shadedErrorBar(lags,all_sua_data(su_num).gsac_avg,all_sua_data(su_num).gsac_sem,{'color','r'});
    h2 = shadedErrorBar(lags,all_sua_data(su_num).msac_avg,all_sua_data(su_num).msac_sem,{'color','b'});
    h3 = shadedErrorBar(lags,all_sua_data(su_num).simsac_avg,all_sua_data(su_num).simsac_sem,{'color','k'});
    h4 = plot(lags,mean(mu_gsac),'k--','linewidth',1.5);
    
    legend([h1.mainLine h2.mainLine h3.mainLine h4],{'Guided','Micro','Simulated','MU avg Guided'});legend('boxoff');
    xlim([-0.15 0.35]);
    % ylim([0.5 2.2]);
    xl = xlim(); yl = ylim();
    line(xl,[1 1],'color','k','linestyle','--');
    line([0 0],yl,'color','k','linestyle','--');
    set(gca,'fontsize',12); box off
    xlabel('Time since saccade onset (s)','fontsize',14)
    ylabel('Relative firing rate','fontsize',14)
    if print_on == 1
        fname = sprintf('examp_SU%d_sactrigavg',su_num);
        fillPage(gcf,'papersize',[6 5]);
        print(fname,'-dpdf','-painters');close
    else
        pause
        close
    end
end

%% TEST FOR STAT SIGNIFICANCE
z_thresh = 2;
for ii = usable_sus'
    fprintf('SU %d of %d\n',ii,length(usable_sus));
   sua_gsac_sup_zval(ii) = (all_sua_data(ii).gsac_avg(sua_gsac_sup_ind(ii)) - 1)/all_sua_data(ii).gsac_sem(sua_gsac_sup_ind(ii));
   sua_gsac_enh_zval(ii) = (all_sua_data(ii).gsac_avg(sua_gsac_enh_ind(ii)) - 1)/all_sua_data(ii).gsac_sem(sua_gsac_enh_ind(ii));
   sua_gsac_presup_zval(ii) = (all_sua_data(ii).gsac_avg(sua_gsac_presup_ind(ii)) - 1)/all_sua_data(ii).gsac_sem(sua_gsac_presup_ind(ii));
  
   sua_msac_sup_zval(ii) = (all_sua_data(ii).msac_avg(sua_msac_sup_ind(ii)) - 1)/all_sua_data(ii).msac_sem(sua_msac_sup_ind(ii));
   sua_msac_enh_zval(ii) = (all_sua_data(ii).msac_avg(sua_msac_enh_ind(ii)) - 1)/all_sua_data(ii).msac_sem(sua_msac_enh_ind(ii));
   sua_msac_presup_zval(ii) = (all_sua_data(ii).msac_avg(sua_msac_presup_ind(ii)) - 1)/all_sua_data(ii).msac_sem(sua_msac_presup_ind(ii));

   subplot(2,1,1)
   shadedErrorBar(lags,all_sua_data(ii).gsac_avg,z_thresh*all_sua_data(ii).gsac_sem,{'color','r'});
   hold on
   plot(lags(sua_gsac_enh_ind(ii)),all_sua_data(ii).gsac_avg(sua_gsac_enh_ind(ii)),'ko')
   plot(lags(sua_gsac_sup_ind(ii)),all_sua_data(ii).gsac_avg(sua_gsac_sup_ind(ii)),'ko')
   plot(lags(sua_gsac_presup_ind(ii)),all_sua_data(ii).gsac_avg(sua_gsac_presup_ind(ii)),'ko')
   xl = xlim();
   line(xl,[1 1],'color','k');
   subplot(2,1,2)
   shadedErrorBar(lags,all_sua_data(ii).msac_avg,z_thresh*all_sua_data(ii).msac_sem,{'color','r'});
   hold on
   plot(lags(sua_msac_enh_ind(ii)),all_sua_data(ii).msac_avg(sua_msac_enh_ind(ii)),'ko')
   plot(lags(sua_msac_sup_ind(ii)),all_sua_data(ii).msac_avg(sua_msac_sup_ind(ii)),'ko')
   plot(lags(sua_msac_presup_ind(ii)),all_sua_data(ii).msac_avg(sua_msac_presup_ind(ii)),'ko')
   xl = xlim();
   line(xl,[1 1],'color','k');
   
   fprintf('GSAC  Enh: %.3fz  Sup: %.3fz  Presup: %.3f\n',sua_gsac_enh_zval(ii),sua_gsac_sup_zval(ii),sua_gsac_presup_zval(ii));
   fprintf('MSAC  Enh: %.3fz  Sup: %.3fz  Presup: %.3f\n',sua_msac_enh_zval(ii),sua_msac_sup_zval(ii),sua_msac_presup_zval(ii));
   pause
   clf
   
end
n_sig_gsac_sup = sum(abs(sua_gsac_sup_zval(usable_sus)) > 3);
n_sig_gsac_enh = sum(abs(sua_gsac_enh_zval(usable_sus)) > 3);
n_sig_msac_sup = sum(abs(sua_msac_sup_zval(usable_sus)) > 3);
n_sig_msac_enh = sum(abs(sua_msac_enh_zval(usable_sus)) > 3);


print_on =1;
cd ~/Analysis/bruce/summary_analysis/

%SUP VS ENH Z 
figure; hold on
plot(abs(sua_gsac_enh_zval(usable_sus)),abs(sua_gsac_sup_zval(usable_sus)),'ro','markersize',4,'linewidth',1.5)
hold on
plot(abs(sua_msac_enh_zval(usable_sus)),abs(sua_msac_sup_zval(usable_sus)),'bo','markersize',4,'linewidth',1.5)
xl = xlim(); xlim([0 xl(2)]);xl = xlim();
yl = ylim(); ylim([0 yl(2)]);yl = ylim();
line(xl,[3 3],'color','k','linestyle','--')
line([3 3],yl,'color','k','linestyle','--')
legend('Guided sac','Micro sac'); %egend('boxoff');
xlabel('Enhancement (z)','fontsize',14)
ylabel('Suppression (z)','fontsize',14)
set(gca,'fontsize',12);
if print_on == 1
    fname = 'scatter_sua_magz';
    fillPage(gcf,'papersize',[5 5]);
    print(fname,'-dpdf','-painters');close
end


%% view all individual su sacmod
close all
n_sus = length(all_sua_data);
f1 = figure();
% f2 = figure();
for ss = 1:n_sus
    ss
    fprintf('Expt %d probe %d\n',all_su_exnums(ss),all_su_probes(ss));
    fprintf('Avg rate:%.3f Nspikes: %d\n',all_sua_data(ss).avg_rate/dt,all_sua_data(ss).tot_nspikes);
    fprintf('%d gsacs,  %d msacs,  %d simsacs \n',all_sua_data(ss).nused_gsac,all_sua_data(ss).nused_msac,all_sua_data(ss).nused_simsac);
    fprintf('gsac enh %.2f sup %.2f,  simsac enh %.2f sup %.2f\n',sua_gsac_enh_rate(ss),sua_gsac_sup_rate(ss),sua_simsac_enh_rate(ss),sua_simsac_sup_rate(ss));
    
    figure(f1);hold on; grid on
    shadedErrorBar(lags,all_sua_data(ss).gsac_avg,all_sua_data(ss).gsac_sem,{'color','r'});
    shadedErrorBar(lags,all_sua_data(ss).msac_avg,all_sua_data(ss).msac_sem,{'color','b'});
    shadedErrorBar(lags,all_sua_data(ss).simsac_avg,all_sua_data(ss).simsac_sem,{'color','k'});
    xlim([-0.2 0.5]);
    xl = xlim(); yl = ylim();
    line(xl,[1 1],'color','k','linestyle','--');
    line([0 0],yl,'color','k','linestyle','--');
    
%     figure(f2)
%     subplot(2,1,1);hold on; grid on
%     shadedErrorBar(lags,all_sua_data(ss).gsac_gray_avg,all_sua_data(ss).gsac_gray_sem,{'color','r'});
%     shadedErrorBar(lags,all_sua_data(ss).gsac_im_avg,all_sua_data(ss).gsac_im_sem,{'color','b'});
%     xlim([-0.2 0.5]);
%     xl = xlim(); yl = ylim();
%     line(xl,[1 1],'color','k','linestyle','--');
%     line([0 0],yl,'color','k','linestyle','--');
%     subplot(2,1,2);hold on; grid on
%     shadedErrorBar(lags,all_sua_data(ss).msac_gray_avg,all_sua_data(ss).msac_gray_sem,{'color','r'});
%     shadedErrorBar(lags,all_sua_data(ss).msac_im_avg,all_sua_data(ss).msac_im_sem,{'color','b'});
%     xlim([-0.2 0.5]);
%     xl = xlim(); yl = ylim();
%     line(xl,[1 1],'color','k','linestyle','--');
%     line([0 0],yl,'color','k','linestyle','--');
% 
%     figure(f2)
%     subplot(2,1,1);hold on; grid on
%     shadedErrorBar(lags,all_sua_data(ss).gsac_ver_avg,all_sua_data(ss).gsac_ver_sem,{'color','r'});
%     shadedErrorBar(lags,all_sua_data(ss).gsac_hor_avg,all_sua_data(ss).gsac_hor_sem,{'color','b'});
%     xlim([-0.2 0.5]);
%     xl = xlim(); yl = ylim();
%     line(xl,[1 1],'color','k','linestyle','--');
%     line([0 0],yl,'color','k','linestyle','--');
%     subplot(2,1,2);hold on; grid on
%     shadedErrorBar(lags,all_sua_data(ss).msac_ver_avg,all_sua_data(ss).msac_ver_sem,{'color','r'});
%     shadedErrorBar(lags,all_sua_data(ss).msac_hor_avg,all_sua_data(ss).msac_hor_sem,{'color','b'});
%     xlim([-0.2 0.5]);
%     xl = xlim(); yl = ylim();
%     line(xl,[1 1],'color','k','linestyle','--');
%     line([0 0],yl,'color','k','linestyle','--');

    pause
    figure(f1);clf;
%     figure(f2);clf;
end

%% OVERALL AVERAGES
cd ~/Analysis/bruce/summary_analysis/
print_on = 0;

close all

%For SU %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
usable_sus = find(su_avgrate > 1 & su_nspks > 500 & all_su_avgmahal > 2.25);
used_sus = usable_sus(su_ngsac(usable_sus) > 50 & su_nmsac(usable_sus) > 50 & su_nsimsac(usable_sus) > 50);
su_gsac = get_struct_data(all_sua_data,used_sus,'gsac_avg');
su_msac = get_struct_data(all_sua_data,used_sus,'msac_avg');
su_simsac = get_struct_data(all_sua_data,used_sus,'simsac_avg');
su_gsac_gray = get_struct_data(all_sua_data,used_sus,'gsac_gray_avg');
su_gsac_im = get_struct_data(all_sua_data,used_sus,'gsac_im_avg');
su_msac_gray = get_struct_data(all_sua_data,used_sus,'msac_gray_avg');
su_msac_im = get_struct_data(all_sua_data,used_sus,'msac_im_avg');

%AVG SACMOD
figure
hold on
h1 = shadedErrorBar(lags,mean(su_gsac),std(su_gsac)/sqrt(length(used_sus)),{'color','r'});
h2 = shadedErrorBar(lags,mean(su_msac),std(su_msac)/sqrt(length(used_sus)),{'color','b'});
h3 = shadedErrorBar(lags,mean(su_simsac),std(su_simsac)/sqrt(length(used_sus)),{'color','k'});
xlim([-0.2 0.5]);
ylim([0.65 1.4]);
xl = xlim(); yl = ylim();
line(xl,[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
legend([h1.mainLine h2.mainLine h3.mainLine],{'Guided','Micro','Simulated'});legend('boxoff');
set(gca,'fontsize',12); box off
xlabel('Time since saccade onset (s)','fontsize',14)
ylabel('Relative firing rate','fontsize',14)
if print_on == 1
    fname = 'su_sactrigavgs';
    fillPage(gcf,'papersize',[6 5]);
    print(fname,'-dpdf','-painters');close
end


%GRAY-BACK VS IMAGE-BACK
figure
subplot(2,1,1); hold on
h1 = shadedErrorBar(lags,nanmean(su_gsac_gray),nanstd(su_gsac_gray)/sqrt(length(used_sus)),{'color','r'});
h2 = shadedErrorBar(lags,nanmean(su_gsac_im),nanstd(su_gsac_im)/sqrt(length(used_sus)),{'color','k'});
legend([h1.mainLine h2.mainLine],{'Gray-back','Image-back'});legend('boxoff');
xlim([-0.2 0.5]);
ylim([0.65 1.35]);
xl = xlim(); yl = ylim();
line(xl,[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time since saccade onset (s)','fontsize',14)
ylabel('Relative firing rate','fontsize',14)
title('Guided sacs','fontsize',14);
set(gca,'fontsize',12); box off
subplot(2,1,2); hold on
h1 = shadedErrorBar(lags,nanmean(su_msac_gray),nanstd(su_msac_gray)/sqrt(length(used_sus)),{'color','r'});
h2 = shadedErrorBar(lags,nanmean(su_msac_im),nanstd(su_msac_im)/sqrt(length(used_sus)),{'color','k'});
legend([h1.mainLine h2.mainLine],{'Gray-back','Image-back'});legend('boxoff');
xlim([-0.2 0.5]);
ylim([0.65 1.35]);
xl = xlim(); yl = ylim();
line(xl,[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time since saccade onset (s)','fontsize',14)
ylabel('Relative firing rate','fontsize',14)
set(gca,'fontsize',12);
title('Micro-sacs','fontsize',14);
if print_on == 1
    fname = 'su_sta_backtype';
    fillPage(gcf,'papersize',[6 10]);
    print(fname,'-dpdf','-painters');close
end


% For MUA %%%%%%%%%%%%%%%%%%%%%%
mu_gsac = get_struct_data(all_mua_data,all_use_mu,'gsac_avg');
mu_msac = get_struct_data(all_mua_data,all_use_mu,'msac_avg');
mu_simsac = get_struct_data(all_mua_data,all_use_mu,'simsac_avg');
mu_gsac_gray = get_struct_data(all_mua_data,all_use_mu,'gsac_gray_avg');
mu_gsac_im = get_struct_data(all_mua_data,all_use_mu,'gsac_im_avg');
mu_msac_gray = get_struct_data(all_mua_data,all_use_mu,'msac_gray_avg');
mu_msac_im = get_struct_data(all_mua_data,all_use_mu,'msac_im_avg');

figure
hold on
h1 = shadedErrorBar(lags,mean(mu_gsac),std(mu_gsac)/sqrt(length(all_use_mu)),{'color','r'});
h2 = shadedErrorBar(lags,mean(mu_msac),std(mu_msac)/sqrt(length(all_use_mu)),{'color','b'});
h3 = shadedErrorBar(lags,mean(mu_simsac),std(mu_simsac)/sqrt(length(all_use_mu)),{'color','k'});
legend([h1.mainLine h2.mainLine h3.mainLine],{'Guided','Micro','Simulated'});legend('boxoff');
xlim([-0.2 0.5]);
ylim([0.7 1.4]);
xl = xlim(); yl = ylim();
line(xl,[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
set(gca,'fontsize',12); box off
xlabel('Time since saccade onset (s)','fontsize',14)
ylabel('Relative firing rate','fontsize',14)
if print_on == 1
    fname = 'mu_sactrigavgs';
    fillPage(gcf,'papersize',[6 5]);
    print(fname,'-dpdf','-painters');close
end


%GRAY-BACK VS IMAGE-BACK
figure
subplot(2,1,1); hold on
h1 = shadedErrorBar(lags,nanmean(mu_gsac_gray),nanstd(mu_gsac_gray)/sqrt(length(all_use_mu)),{'color','r'});
h2 = shadedErrorBar(lags,nanmean(mu_gsac_im),nanstd(mu_gsac_im)/sqrt(length(all_use_mu)),{'color','k'});
legend([h1.mainLine h2.mainLine],{'Gray-back','Image-back'});legend('boxoff');
xlim([-0.2 0.5]);
ylim([0.7 1.4]);
xl = xlim(); yl = ylim();
line(xl,[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time since saccade onset (s)','fontsize',14)
ylabel('Relative firing rate','fontsize',14)
set(gca,'fontsize',12); box off
title('Guided sacs','fontsize',14);
subplot(2,1,2); hold on
h1 = shadedErrorBar(lags,nanmean(mu_msac_gray),nanstd(mu_msac_gray)/sqrt(length(all_use_mu)),{'color','r'});
h2 = shadedErrorBar(lags,nanmean(mu_msac_im),nanstd(mu_msac_im)/sqrt(length(all_use_mu)),{'color','k'});
legend([h1.mainLine h2.mainLine],{'Gray-back','Image-back'});legend('boxoff');
xlim([-0.2 0.5]);
ylim([0.7 1.3]);
xl = xlim(); yl = ylim();
line(xl,[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time since saccade onset (s)','fontsize',14)
ylabel('Relative firing rate','fontsize',14)
set(gca,'fontsize',12); box off
title('Micro-sacs','fontsize',14);
if print_on == 1
    fname = 'mu_sta_backtype';
    fillPage(gcf,'papersize',[6 10]);
    print(fname,'-dpdf','-painters');close
end


figure
imagesc(lags,1:length(all_use_mu),mu_gsac)
xlim([-0.2 0.5]);
xl = xlim(); yl = ylim();
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time since saccade onset (s)','fontsize',14)
ylabel('MU Number','fontsize',14);
colorbar;
caxis([0.7 1.4])
set(gca,'fontsize',12); box off
title('Guided sacs','fontsize',14);
if print_on == 1
    fname = 'mu_all_gsta';
    fillPage(gcf,'papersize',[10 10]);
    print(fname,'-dpdf','-painters');close
end

figure
imagesc(lags,1:length(used_sus),su_gsac)
xlim([-0.2 0.5]);
xl = xlim(); yl = ylim();
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time since saccade onset (s)','fontsize',20)
ylabel('SU Number','fontsize',20);
colorbar;
caxis([0.5 1.7])
set(gca,'fontsize',20); box off
title('Guided sacs','fontsize',20);
if print_on == 1
    fname = 'su_all_gsta';
    fillPage(gcf,'papersize',[6 4]);
    print(fname,'-dpdf','-zbuffer');close
end

%% FLASHED VS SIM SACS
close all
cur_mu_use_flash = all_use_mu(all_mu_exnums(all_use_mu) == 81);
mu_simsac_flashed = get_struct_data(all_mua_data,cur_mu_use_flash,'simsac_avg');

cur_mu_use_sim = all_use_mu(all_mu_exnums(all_use_mu) > 81);
mu_simsac_sim = get_struct_data(all_mua_data,cur_mu_use_sim,'simsac_avg');

figure
hold on; grid on
h1 = shadedErrorBar(lags,mean(mu_simsac_flashed),std(mu_simsac_flashed)/sqrt(length(cur_mu_use_flash)),{'color','r'});
h2 = shadedErrorBar(lags,mean(mu_simsac_sim),std(mu_simsac_sim)/sqrt(length(cur_mu_use_sim)),{'color','b'});
legend([h1.mainLine h2.mainLine],{'Flashed simsac','Translated simsac'});legend('boxoff');
xlim([-0.05 0.25]);
ylim([0.9 1.1]);
xl = xlim(); yl = ylim();
line(xl,[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time since saccade onset (s)','fontsize',14)
ylabel('Relative firing rate','fontsize',14)
set(gca,'fontsize',12); box off
if print_on == 1
    fname = 'mu_flashed_vs_transl_simsac';
    fillPage(gcf,'papersize',[6 5]);
    print(fname,'-dpdf','-painters');close
end

%% SB VS SPARSE VS DENSE
close all
sb_mu_use = all_use_mu(all_mu_exnums(all_use_mu) == 81);
sb_mu_gsac = get_struct_data(all_mua_data,sb_mu_use,'gsac_avg');
sb_mu_msac = get_struct_data(all_mua_data,sb_mu_use,'msac_avg');
sb_mu_simsac = get_struct_data(all_mua_data,sb_mu_use,'simsac_avg');

sparse_mu_use = all_use_mu(all_mu_exnums(all_use_mu) > 81 & all_mu_exnums(all_use_mu) < 91);
sparse_mu_gsac = get_struct_data(all_mua_data,sparse_mu_use,'gsac_avg');
sparse_mu_msac = get_struct_data(all_mua_data,sparse_mu_use,'msac_avg');
sparse_mu_simsac = get_struct_data(all_mua_data,sparse_mu_use,'simsac_avg');

dense_mu_use = all_use_mu(all_mu_exnums(all_use_mu) >= 91);
dense_mu_gsac = get_struct_data(all_mua_data,dense_mu_use,'gsac_avg');
dense_mu_msac = get_struct_data(all_mua_data,dense_mu_use,'msac_avg');
dense_mu_simsac = get_struct_data(all_mua_data,dense_mu_use,'simsac_avg');

figure
subplot(3,1,1); %MU GSAC
hold on
h1 = shadedErrorBar(lags,mean(sb_mu_gsac),std(sb_mu_gsac)/sqrt(length(sb_mu_use)),{'color','r'});
h2 = shadedErrorBar(lags,mean(sparse_mu_gsac),std(sparse_mu_gsac)/sqrt(length(sb_mu_use)),{'color','b'});
h3 = shadedErrorBar(lags,mean(dense_mu_gsac),std(dense_mu_gsac)/sqrt(length(sb_mu_use)),{'color','k'});
legend([h1.mainLine h2.mainLine h3.mainLine],{'Single-bar','Sparse','Dense'});legend('boxoff');
xlim([-0.05 0.3]);
ylim([0.65 1.4]);
xl = xlim(); yl = ylim();
line(xl,[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time since saccade onset (s)','fontsize',14)
ylabel('Relative firing rate','fontsize',14)
title('Guided saccades','fontsize',16);
set(gca,'fontsize',12); box off

subplot(3,1,2);
hold on
h1 = shadedErrorBar(lags,mean(sb_mu_msac),std(sb_mu_msac)/sqrt(length(sb_mu_use)),{'color','r'});
h2 = shadedErrorBar(lags,mean(sparse_mu_msac),std(sparse_mu_msac)/sqrt(length(sb_mu_use)),{'color','b'});
h3 = shadedErrorBar(lags,mean(dense_mu_msac),std(dense_mu_msac)/sqrt(length(sb_mu_use)),{'color','k'});
legend([h1.mainLine h2.mainLine h3.mainLine],{'Single-bar','Sparse','Dense'});legend('boxoff');
xlim([-0.05 0.3]);
ylim([0.7 1.3]);
xl = xlim(); yl = ylim();
line(xl,[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time since saccade onset (s)','fontsize',14)
ylabel('Relative firing rate','fontsize',14)
title('Micro-saccades','fontsize',16);
set(gca,'fontsize',12); box off

subplot(3,1,3);
hold on
h1 = shadedErrorBar(lags,mean(sb_mu_simsac),std(sb_mu_simsac)/sqrt(length(sb_mu_use)),{'color','r'});
h2 = shadedErrorBar(lags,mean(sparse_mu_simsac),std(sparse_mu_simsac)/sqrt(length(sb_mu_use)),{'color','b'});
h3 = shadedErrorBar(lags,mean(dense_mu_simsac),std(dense_mu_simsac)/sqrt(length(sb_mu_use)),{'color','k'});
legend([h1.mainLine h2.mainLine h3.mainLine],{'Single-bar','Sparse','Dense'});legend('boxoff');
xlim([-0.05 0.3]);
ylim([0.85 1.15]);
xl = xlim(); yl = ylim();
line(xl,[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time since saccade onset (s)','fontsize',14)
ylabel('Relative firing rate','fontsize',14)
title('Simulated saccades','fontsize',16);
set(gca,'fontsize',12); box off
if print_on == 1
    fname = 'mu_sta_bardensity';
    fillPage(gcf,'papersize',[6 15]);
    print(fname,'-dpdf','-painters');close
end


%% SAC TIME CORRELATIONS
close all
gsac_grayback_tral_avg = get_struct_data(all_gen_data,2:length(Expt_nums),'gsac_grayback_trial_avg');
gsac_imback_tral_avg = get_struct_data(all_gen_data,2:length(Expt_nums),'gsac_imback_trial_avg');
msac_simsac_tral_avg = get_struct_data(all_gen_data,2:length(Expt_nums),'msac_simsac_trial_avg');
msac_grayback_tral_avg = get_struct_data(all_gen_data,2:length(Expt_nums),'msac_grayback_trial_avg');
msac_imback_tral_avg = get_struct_data(all_gen_data,2:length(Expt_nums),'msac_imback_trial_avg');
gsac_grayback_tral_avg = bsxfun(@rdivide,gsac_grayback_tral_avg,sum(gsac_grayback_tral_avg,2));
gsac_imback_tral_avg = bsxfun(@rdivide,gsac_imback_tral_avg,sum(gsac_imback_tral_avg,2));
msac_simsac_tral_avg = bsxfun(@rdivide,msac_simsac_tral_avg,sum(msac_simsac_tral_avg,2));
msac_grayback_tral_avg = bsxfun(@rdivide,msac_grayback_tral_avg,sum(msac_grayback_tral_avg,2));
msac_imback_tral_avg = bsxfun(@rdivide,msac_imback_tral_avg,sum(msac_imback_tral_avg,2));

sac_times = [0.7 1.4 2.1 2.8 3.5];

figure
subplot(2,1,1) %%%%%%%%%%%%%%%%%%%%
plot(trial_lags,jmm_smooth_1d_cor(nanmean(gsac_grayback_tral_avg),3),'b')
hold on
plot(trial_lags,jmm_smooth_1d_cor(nanmean(gsac_imback_tral_avg),3),'r')
xlabel('Time since trial onset (s)','fontsize',14)
ylabel('Relative frequency','fontsize',14)
legend('Gray back','Image back','Location','northwest'); legend('boxoff');
yl = ylim();
for ii = 1:length(sac_times); line(sac_times([ii ii]),yl,'color','k','linestyle','--'); end
set(gca,'fontsize',12); box off
title('Guided saccades','fontsize',14)
subplot(2,1,2) %%%%%%%%%%%%%%%%%%
plot(trial_lags,jmm_smooth_1d_cor(nanmean(msac_grayback_tral_avg),5),'b')
hold on
plot(trial_lags,jmm_smooth_1d_cor(nanmean(msac_imback_tral_avg),5),'r')
plot(trial_lags,jmm_smooth_1d_cor(nanmean(msac_simsac_tral_avg),5),'k')
legend('Gray back','Image back','Simulated saccade','Location','northwest'); legend('boxoff');
yl = ylim();
for ii = 1:length(sac_times); line(sac_times([ii ii]),yl,'color','k','linestyle','--'); end
xlabel('Time since trial onset (s)','fontsize',14)
ylabel('Relative frequency','fontsize',14)
set(gca,'fontsize',12); box off
title('Micro-saccades','fontsize',14)
if print_on == 1
    fname = 'relativeSaccadeTiming';
    fillPage(gcf,'papersize',[8 8]);
    print(fname,'-dpdf','-painters');close
end

%ACORR PLOT
msac_acorrs = get_struct_data(all_gen_data,2:length(Expt_nums),'msac_acorr');
msac_gsac_xcorr = get_struct_data(all_gen_data,2:length(Expt_nums),'msac_gsac_xcorr');
figure; hold on
h1 = shadedErrorBar(acorr_lags,mean(msac_acorrs),std(msac_acorrs)/sqrt(length(Expt_nums)),{'color','r'});
h2 = shadedErrorBar(acorr_lags,mean(msac_gsac_xcorr),std(msac_gsac_xcorr)/sqrt(length(Expt_nums)),{'color','k'});
legend([h1.mainLine h2.mainLine],{'Micro autocorr','Micro-guided xcorr'});legend('boxoff');
ylim([-0.05 0.1]);
xl = xlim();
line(xl,[0 0],'color','k','linestyle','--');
xlabel('Time lag (s)','fontsize',16);
ylabel('Autocorrelation','fontsize',16);
set(gca,'fontsize',12); box off
if print_on == 1
    fname = 'saccadeAcorrs';
    fillPage(gcf,'papersize',[6 5]);
    print(fname,'-dpdf','-painters');close
end

%% TRIG-AVG VS KERNEL ESTIMATES

genmu_gsac_avg = get_struct_data(all_gen_data,1:length(Expt_nums),'mua_gsac_avg');
genmu_msac_avg = get_struct_data(all_gen_data,1:length(Expt_nums),'mua_msac_avg');
genmu_msac_nb_avg = get_struct_data(all_gen_data,1:length(Expt_nums),'mua_nb_msac');
genmu_msac_b_avg = get_struct_data(all_gen_data,1:length(Expt_nums),'mua_b_msac');
genmu_simsac_avg = get_struct_data(all_gen_data,1:length(Expt_nums),'mua_simsac_avg');
genmu_gsac_kern = get_struct_data(all_gen_data,1:length(Expt_nums),'mua_gsac_kern');
genmu_msac_kern = get_struct_data(all_gen_data,1:length(Expt_nums),'mua_msac_kern');
genmu_simsac_kern = get_struct_data(all_gen_data,1:length(Expt_nums),'mua_simsac_kern');

figure
subplot(3,1,1);hold on %%%%%%%%%%%%%%%%%%%%
h1 = shadedErrorBar(lags,mean(genmu_gsac_avg),std(genmu_gsac_avg)/sqrt(length(Expt_nums)),{'color','r'});
h2 = shadedErrorBar(lags,mean(genmu_msac_avg),std(genmu_msac_avg)/sqrt(length(Expt_nums)),{'color','b'});
h3 = shadedErrorBar(lags,mean(genmu_simsac_avg),std(genmu_simsac_avg)/sqrt(length(Expt_nums)),{'color','k'});
legend([h1.mainLine h2.mainLine h3.mainLine],{'Guided','Micro','Simulated'});legend('boxoff');
xlim([-0.2 0.5]); yl = ylim();
line([-0.2 0.5],[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time since saccade onset (s)','fontsize',14);
ylabel('Relative firing rate','fontsize',14);
set(gca,'fontsize',12); box off
title('Triggered averages','fontsize',16);
subplot(3,1,2); hold on %%%%%%%%%%%%%%%%%%%%%%%%
h1 = shadedErrorBar(lags,mean(genmu_msac_nb_avg),std(genmu_msac_nb_avg)/sqrt(length(Expt_nums)),{'color','r'});
h2 = shadedErrorBar(lags,mean(genmu_msac_b_avg),std(genmu_msac_b_avg)/sqrt(length(Expt_nums)),{'color','b'});
legend([h1.mainLine h2.mainLine],{'No micro-bursts','Micro-burst-only'});legend('boxoff');
xlim([-0.2 0.5]); yl = ylim();
line([-0.2 0.5],[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time since saccade onset (s)','fontsize',14);
ylabel('Relative firing rate','fontsize',14);
set(gca,'fontsize',12); box off
title('10Hz micro-saccade bursting (trig-averages)','fontsize',16);
subplot(3,1,3); hold on %%%%%%%%%%%%%%%%%%%%%%%%%
h1 = shadedErrorBar(sac_bin_cents,mean(genmu_gsac_kern),std(genmu_gsac_kern)/sqrt(length(Expt_nums)),{'color','r'});
h2 = shadedErrorBar(sac_bin_cents,mean(genmu_msac_kern),std(genmu_msac_kern)/sqrt(length(Expt_nums)),{'color','b'});
h3 = shadedErrorBar(sac_bin_cents,mean(genmu_simsac_kern),std(genmu_simsac_kern)/sqrt(length(Expt_nums)),{'color','k'});
legend([h1.mainLine h2.mainLine h3.mainLine],{'Guided','Micro','Simulated'});legend('boxoff');
xlim([-0.2 0.5]); yl = ylim();
line([-0.2 0.5],[0 0],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time since saccade onset (s)','fontsize',14);
ylabel('Filter amplitude','fontsize',14);
set(gca,'fontsize',12); box off
title('Likelihood based filter estimation','fontsize',16);
if print_on == 1
    fname = 'sacta_vs_kernels';
    fillPage(gcf,'papersize',[6 15]);
    print(fname,'-dpdf','-painters');close
end


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

%% view all expt-avg mu sacmod
bad_probes = 16;
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
    mu_gsac_avgs = get_struct_data(all_mua_data,cur_mu_set,'gsac_avg');
    mu_msac_avgs = get_struct_data(all_mua_data,cur_mu_set,'msac_avg');
    mu_simsac_avgs = get_struct_data(all_mua_data,cur_mu_set,'simsac_avg');

    figure(f1);hold on; grid on
    shadedErrorBar(lags,nanmean(mu_gsac_avgs),nanstd(mu_gsac_avgs)/sqrt(length(cur_mu_set)),{'color','r'});
    shadedErrorBar(lags,nanmean(mu_msac_avgs),nanstd(mu_msac_avgs)/sqrt(length(cur_mu_set)),{'color','b'});
    shadedErrorBar(lags,nanmean(mu_simsac_avgs),nanstd(mu_simsac_avgs)/sqrt(length(cur_mu_set)),{'color','k'});
    xlim([-0.2 0.5]);
    xl = xlim(); yl = ylim();
    line(xl,[1 1],'color','k','linestyle','--');
    line([0 0],yl,'color','k','linestyle','--');
    ylim([0.7 1.4]);

    mu_gsac_grayback_avgs = get_struct_data(all_mua_data,cur_mu_set,'gsac_gray_avg');
    mu_gsac_imback_avgs = get_struct_data(all_mua_data,cur_mu_set,'gsac_im_avg');
    mu_msac_grayback_avgs = get_struct_data(all_mua_data,cur_mu_set,'msac_gray_avg');
    mu_msac_imback_avgs = get_struct_data(all_mua_data,cur_mu_set,'msac_im_avg');
    figure(f2);
    subplot(2,1,1);hold on; grid on
    shadedErrorBar(lags,nanmean(mu_gsac_grayback_avgs),nanstd(mu_gsac_grayback_avgs)/sqrt(length(cur_mu_set)),{'color','r'});
    shadedErrorBar(lags,nanmean(mu_gsac_imback_avgs),nanstd(mu_gsac_imback_avgs)/sqrt(length(cur_mu_set)),{'color','b'});
    xlim([-0.2 0.5]);
    xl = xlim(); yl = ylim();
    line(xl,[1 1],'color','k','linestyle','--');
    line([0 0],yl,'color','k','linestyle','--');
    ylim([0.7 1.4]);
    subplot(2,1,2);hold on; grid on
    shadedErrorBar(lags,nanmean(mu_msac_grayback_avgs),nanstd(mu_msac_grayback_avgs)/sqrt(length(cur_mu_set)),{'color','r'});
    shadedErrorBar(lags,nanmean(mu_msac_imback_avgs),nanstd(mu_msac_imback_avgs)/sqrt(length(cur_mu_set)),{'color','b'});
    xlim([-0.2 0.5]);
    xl = xlim(); yl = ylim();
    line(xl,[1 1],'color','k','linestyle','--');
    line([0 0],yl,'color','k','linestyle','--');
    ylim([0.7 1.4]);

    pause
    figure(f1); clf;
    figure(f2); clf;
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