%%
close all
clear all

dt = 0.005;
sm_sig = 0.01/dt;

fig_dir = '/home/james/Analysis/bruce/saccade_modulation/';
%%
Expt_list = {'M266','M270','M275','M277','M281','M287','M289','M294'};
Expt_RFPos = [1.09 3.57 0.96 3.17 3.85 2.69 3.44 4.17];
% Expt_list = {'M266','M270','M275','M277','M281','M287','M294'};
% Expt_list = {'M232','M235','M239','M266','M270','M275','M277','M281','M287','M294'};
% Expt_RFPos = [5.22 8.13 5.25 1.09 3.57 0.96 3.17 3.85 2.69 4.17];
% Expt_list = {'M266','M270','M275','M277','M294'};
% Expt_RFPos = [1.09 3.57 0.96 3.17 4.17];
% Expt_list = {'M266','M275','M281','M287','M289'};
% Expt_RFPos = [1.09  0.96 3.85 2.69 3.44];

sname = 'sac_trig_avg_data4';
lem_sua_exptnum = [];
lem_mua_exptnum = [];
for ee = 1:length(Expt_list)
    Expt_name = Expt_list{ee};
    Expt_num = str2num(Expt_name(2:end));
    save_dir = ['~/Analysis/bruce/' Expt_name '/sac_mod'];
    cd(save_dir)
    load(sname)
    all_sua_data(ee) = sua_data;
    lem_sua_exptnum = cat(1,lem_sua_exptnum,ones(size(sua_data.msac_avg,2),1)*Expt_num);
    all_mua_data(ee) = mua_data;
    lem_mua_exptnum = cat(1,lem_mua_exptnum,ones(size(mua_data.msac_avg,2),1)*Expt_num);
end
%%
lem.mua_gsac_avg = cat(2,all_mua_data(:).gsac_avg);
lem.mua_msac_avg = cat(2,all_mua_data(:).msac_avg);
lem.mua_simsac_avg = cat(2,all_mua_data(:).simsac_avg);

lem.mua_gsac_out_avg = cat(2,all_mua_data(:).gsac_out_avg);
lem.mua_gsac_in_avg = cat(2,all_mua_data(:).gsac_in_avg);
lem.mua_gsac_pos_avg = cat(2,all_mua_data(:).gsac_pos_avg);
lem.mua_gsac_neg_avg = cat(2,all_mua_data(:).gsac_neg_avg);

lem.mua_gsac_im_avg = cat(2,all_mua_data(:).gsac_im_avg);
lem.mua_gsac_gray_avg = cat(2,all_mua_data(:).gsac_gray_avg);
lem.mua_msac_im_avg = cat(2,all_mua_data(:).msac_im_avg);
lem.mua_msac_gray_avg = cat(2,all_mua_data(:).msac_gray_avg);

lem.mua_msac_end_avg = cat(2,all_mua_data(:).msac_end_avg);
lem.mua_gsac_end_avg = cat(2,all_mua_data(:).gsac_end_avg);

%%
lem.sua_gsac_avg = cat(2,all_sua_data(:).gsac_avg);
lem.sua_msac_avg = cat(2,all_sua_data(:).msac_avg);
lem.sua_simsac_avg = cat(2,all_sua_data(:).simsac_avg);
lem.sua_gsac_im_avg = cat(2,all_sua_data(:).gsac_im_avg);
lem.sua_gsac_gray_avg = cat(2,all_sua_data(:).gsac_gray_avg);
lem.sua_msac_im_avg = cat(2,all_sua_data(:).msac_im_avg);
lem.sua_msac_gray_avg = cat(2,all_sua_data(:).msac_gray_avg);

%% SMOOTH SUA SAC-TAs
lem.sua_gsac_sm = lem.sua_gsac_avg;
lem.sua_msac_sm = lem.sua_gsac_avg;
lem.sua_gsac_im_sm = lem.sua_gsac_avg;
lem.sua_msac_im_sm = lem.sua_gsac_avg;
lem.sua_simsac_sm = lem.sua_gsac_avg;
for ii = 1:size(lem.sua_gsac_avg,2)
    lem.sua_gsac_im_sm(:,ii) = jmm_smooth_1d_cor(lem.sua_gsac_im_avg(:,ii),sm_sig);
    lem.sua_msac_im_sm(:,ii) = jmm_smooth_1d_cor(lem.sua_msac_im_avg(:,ii),sm_sig);
    lem.sua_gsac_sm(:,ii) = jmm_smooth_1d_cor(lem.sua_gsac_gray_avg(:,ii),sm_sig);
    lem.sua_msac_sm(:,ii) = jmm_smooth_1d_cor(lem.sua_msac_gray_avg(:,ii),sm_sig);
    lem.sua_simsac_sm(:,ii) = jmm_smooth_1d_cor(lem.sua_simsac_avg(:,ii),sm_sig);
end

lem.mua_gsac_sm = lem.mua_gsac_avg;
lem.mua_msac_sm = lem.mua_gsac_avg;
lem.mua_gsac_im_sm = lem.mua_gsac_avg;
lem.mua_msac_im_sm = lem.mua_gsac_avg;
lem.mua_simsac_sm = lem.mua_gsac_avg;
for ii = 1:size(lem.mua_gsac_avg,2)
    lem.mua_gsac_im_sm(:,ii) = jmm_smooth_1d_cor(lem.mua_gsac_im_avg(:,ii),sm_sig);
    lem.mua_msac_im_sm(:,ii) = jmm_smooth_1d_cor(lem.mua_msac_im_avg(:,ii),sm_sig);
    lem.mua_gsac_sm(:,ii) = jmm_smooth_1d_cor(lem.mua_gsac_gray_avg(:,ii),sm_sig);
    lem.mua_msac_sm(:,ii) = jmm_smooth_1d_cor(lem.mua_msac_gray_avg(:,ii),sm_sig);
    lem.mua_simsac_sm(:,ii) = jmm_smooth_1d_cor(lem.mua_simsac_avg(:,ii),sm_sig);
end

%%
lagrange = find(lags*dt >= 0 & lags*dt <= 0.4);
lem.sua_avg_rates = cat(2,all_sua_data(:).avg_rates);
lem.sua_tot_nspikes = cat(2,all_sua_data(:).tot_nspikes);
[lem.sua_gsac_minvals,sua_gsac_mintime] = min(lem.sua_gsac_sm(lagrange,:));
[lem.sua_gsac_maxvals,sua_gsac_maxtime] = max(lem.sua_gsac_sm(lagrange,:));
[lem.sua_msac_minvals,sua_msac_mintime] = min(lem.sua_msac_sm(lagrange,:));
[lem.sua_msac_maxvals,sua_msac_maxtime] = max(lem.sua_msac_sm(lagrange,:));
lem.sua_gsac_mintime = lags(lagrange(sua_gsac_mintime))*dt;
lem.sua_msac_mintime = lags(lagrange(sua_msac_mintime))*dt;
lem.sua_gsac_maxtime = lags(lagrange(sua_gsac_maxtime))*dt;
lem.sua_msac_maxtime = lags(lagrange(sua_msac_maxtime))*dt;

lem.mua_avg_rates = cat(2,all_mua_data(:).avg_rates);
lem.mua_tot_nspikes = cat(2,all_mua_data(:).tot_nspikes);
[lem.mua_gsac_minvals,mua_gsac_mintime] = min(lem.mua_gsac_sm(lagrange,:));
[lem.mua_gsac_maxvals,mua_gsac_maxtime] = max(lem.mua_gsac_sm(lagrange,:));
[lem.mua_msac_minvals,mua_msac_mintime] = min(lem.mua_msac_sm(lagrange,:));
[lem.mua_msac_maxvals,mua_msac_maxtime] = max(lem.mua_msac_sm(lagrange,:));
lem.mua_gsac_mintime = lags(lagrange(mua_gsac_mintime))*dt;
lem.mua_msac_mintime = lags(lagrange(mua_msac_mintime))*dt;
lem.mua_gsac_maxtime = lags(lagrange(mua_gsac_maxtime))*dt;
lem.mua_msac_maxtime = lags(lagrange(mua_msac_maxtime))*dt;

%%
Expt_list = {'G085','G086','G087','G088','G089','G091','G093','G095'};
sname = 'sac_trig_avg_data4';
clear all_sua_data all_mua_data
jbe_sua_exptnum = [];
jbe_mua_exptnum = [];
for ee = 1:length(Expt_list)
    Expt_name = Expt_list{ee};
    Expt_num = str2num(Expt_name(2:end));
    save_dir = ['~/Analysis/bruce/' Expt_name '/sac_mod'];
    cd(save_dir)
    load(sname)
    all_sua_data(ee) = sua_data;
    jbe_sua_exptnum = cat(1,jbe_sua_exptnum,ones(size(sua_data.msac_avg,2),1)*Expt_num);
    all_mua_data(ee) = mua_data;
    jbe_mua_exptnum = cat(1,jbe_mua_exptnum,ones(size(mua_data.msac_avg,2),1)*Expt_num);
end

jbe.mua_gsac_avg = cat(2,all_mua_data(:).gsac_avg);
jbe.mua_msac_avg = cat(2,all_mua_data(:).msac_avg);
jbe.mua_simsac_avg = cat(2,all_mua_data(:).simsac_avg);

jbe.mua_gsac_out_avg = cat(2,all_mua_data(:).gsac_out_avg);
jbe.mua_gsac_in_avg = cat(2,all_mua_data(:).gsac_in_avg);
jbe.mua_gsac_pos_avg = cat(2,all_mua_data(:).gsac_pos_avg);
jbe.mua_gsac_neg_avg = cat(2,all_mua_data(:).gsac_neg_avg);

jbe.mua_gsac_im_avg = cat(2,all_mua_data(:).gsac_im_avg);
jbe.mua_gsac_gray_avg = cat(2,all_mua_data(:).gsac_gray_avg);

jbe.mua_msac_im_avg = cat(2,all_mua_data(:).msac_im_avg);
jbe.mua_msac_gray_avg = cat(2,all_mua_data(:).msac_gray_avg);

jbe.mua_msac_end_avg = cat(2,all_mua_data(:).msac_end_avg);
jbe.mua_gsac_end_avg = cat(2,all_mua_data(:).gsac_end_avg);

% jbe.mua_gsac_filts = cat(1,all_mua_data(:).gsac_filts)';
% jbe.mua_msac_filts = cat(1,all_mua_data(:).msac_filts)';
% jbe.sua_gsac_filts = cat(1,all_sua_data(:).gsac_filts)';
% jbe.sua_msac_filts = cat(1,all_sua_data(:).msac_filts)';

jbe.sua_gsac_avg = cat(2,all_sua_data(:).gsac_avg);
jbe.sua_msac_avg = cat(2,all_sua_data(:).msac_avg);
jbe.sua_simsac_avg = cat(2,all_sua_data(:).simsac_avg);
jbe.sua_gsac_im_avg = cat(2,all_sua_data(:).gsac_im_avg);
jbe.sua_gsac_gray_avg = cat(2,all_sua_data(:).gsac_gray_avg);
jbe.sua_msac_im_avg = cat(2,all_sua_data(:).msac_im_avg);
jbe.sua_msac_gray_avg = cat(2,all_sua_data(:).msac_gray_avg);

%%
lagrange = find(lags*dt >= 0 & lags*dt <= 0.4);
jbe.sua_avg_rates = cat(2,all_sua_data(:).avg_rates);
jbe.sua_tot_nspikes = cat(2,all_sua_data(:).tot_nspikes);
[jbe.sua_gsac_minvals,sua_gsac_mintime] = min(jbe.sua_gsac_avg(lagrange,:));
[jbe.sua_gsac_maxvals,sua_gsac_maxtime] = max(jbe.sua_gsac_avg(lagrange,:));
[jbe.sua_msac_minvals,sua_msac_mintime] = min(jbe.sua_msac_avg(lagrange,:));
[jbe.sua_msac_maxvals,sua_msac_maxtime] = max(jbe.sua_msac_avg(lagrange,:));
jbe.sua_gsac_mintime = lags(lagrange(sua_gsac_mintime))*dt;
jbe.sua_msac_mintime = lags(lagrange(sua_msac_mintime))*dt;
jbe.sua_gsac_maxtime = lags(lagrange(sua_gsac_maxtime))*dt;
jbe.sua_msac_maxtime = lags(lagrange(sua_msac_maxtime))*dt;

jbe.mua_avg_rates = cat(2,all_mua_data(:).avg_rates);
jbe.mua_tot_nspikes = cat(2,all_mua_data(:).tot_nspikes);
[jbe.mua_gsac_minvals,mua_gsac_mintime] = min(jbe.mua_gsac_avg(lagrange,:));
[jbe.mua_gsac_maxvals,mua_gsac_maxtime] = max(jbe.mua_gsac_avg(lagrange,:));
[jbe.mua_msac_minvals,mua_msac_mintime] = min(jbe.mua_msac_avg(lagrange,:));
[jbe.mua_msac_maxvals,mua_msac_maxtime] = max(jbe.mua_msac_avg(lagrange,:));
jbe.mua_gsac_mintime = lags(lagrange(mua_gsac_mintime))*dt;
jbe.mua_msac_mintime = lags(lagrange(mua_msac_mintime))*dt;
jbe.mua_gsac_maxtime = lags(lagrange(mua_gsac_maxtime))*dt;
jbe.mua_msac_maxtime = lags(lagrange(mua_msac_maxtime))*dt;

%% SMOOTH SUA SAC-TAs
jbe.sua_gsac_sm = jbe.sua_gsac_avg;
jbe.sua_msac_sm = jbe.sua_gsac_avg;
jbe.sua_gsac_im_sm = jbe.sua_gsac_avg;
jbe.sua_msac_im_sm = jbe.sua_gsac_avg;
jbe.sua_simsac_sm = jbe.sua_gsac_avg;
for ii = 1:size(jbe.sua_gsac_avg,2)
    jbe.sua_gsac_im_sm(:,ii) = jmm_smooth_1d_cor(jbe.sua_gsac_im_avg(:,ii),sm_sig);
    jbe.sua_msac_im_sm(:,ii) = jmm_smooth_1d_cor(jbe.sua_msac_im_avg(:,ii),sm_sig);
    jbe.sua_gsac_sm(:,ii) = jmm_smooth_1d_cor(jbe.sua_gsac_gray_avg(:,ii),sm_sig);
    jbe.sua_msac_sm(:,ii) = jmm_smooth_1d_cor(jbe.sua_msac_gray_avg(:,ii),sm_sig);
    jbe.sua_simsac_sm(:,ii) = jmm_smooth_1d_cor(jbe.sua_simsac_avg(:,ii),sm_sig);
end

jbe.mua_gsac_sm = jbe.mua_gsac_avg;
jbe.mua_msac_sm = jbe.mua_gsac_avg;
jbe.mua_gsac_im_sm = jbe.mua_gsac_avg;
jbe.mua_msac_im_sm = jbe.mua_gsac_avg;
jbe.mua_simsac_sm = jbe.mua_gsac_avg;
for ii = 1:size(jbe.mua_gsac_avg,2)
    jbe.mua_gsac_im_sm(:,ii) = jmm_smooth_1d_cor(jbe.mua_gsac_im_avg(:,ii),sm_sig);
    jbe.mua_msac_im_sm(:,ii) = jmm_smooth_1d_cor(jbe.mua_msac_im_avg(:,ii),sm_sig);
    jbe.mua_gsac_sm(:,ii) = jmm_smooth_1d_cor(jbe.mua_gsac_gray_avg(:,ii),sm_sig);
    jbe.mua_msac_sm(:,ii) = jmm_smooth_1d_cor(jbe.mua_msac_gray_avg(:,ii),sm_sig);
    jbe.mua_simsac_sm(:,ii) = jmm_smooth_1d_cor(jbe.mua_simsac_avg(:,ii),sm_sig);
%     jbe.mua_gsac_sm(:,ii) = jmm_smooth_1d_cor(jbe.mua_gsac_im_avg(:,ii),sm_sig);
%     jbe.mua_msac_sm(:,ii) = jmm_smooth_1d_cor(jbe.mua_msac_im_avg(:,ii),sm_sig);
%     jbe.mua_gsac_sm(:,ii) = jmm_smooth_1d_cor(jbe.mua_simsac_avg(:,ii),sm_sig);
end

%% FOR SUA
close all
min_nspks = 1e4;
min_rate = 5;
lem_use_sus = find(lem.sua_tot_nspikes >= min_nspks & lem.sua_avg_rates/dt >= min_rate);
jbe_use_sus = find(jbe.sua_tot_nspikes >= min_nspks & jbe.sua_avg_rates/dt >= min_rate);

jit_amp = 0.002;
f1 = figure; hold on
plot(lem.sua_gsac_mintime(lem_use_sus)+randn(length(lem_use_sus),1)*jit_amp,lem.sua_gsac_maxtime(lem_use_sus)+randn(length(lem_use_sus),1)*jit_amp,'.')
plot(jbe.sua_gsac_mintime(jbe_use_sus)+randn(length(jbe_use_sus),1)*jit_amp,jbe.sua_gsac_maxtime(jbe_use_sus)+randn(length(jbe_use_sus),1)*jit_amp,'r.')
line([0 0.3],[0 0.3],'color','k');
xlim([0 0.3]);
ylim([0 0.3]);
xlabel('Peak suppression time (s)');
ylabel('Peak enhancement time (s)');

% figure; hold on
% plot(1-lem.sua_gsac_minvals(lem_use_sus),lem.sua_gsac_maxvals(lem_use_sus)-1,'o')
% plot(1-jbe.sua_gsac_minvals(jbe_use_sus),jbe.sua_gsac_maxvals(jbe_use_sus)-1,'.')

% figure; hold on
% plot(lem.sua_gsac_mintime(lem_use_sus),lem.sua_msac_mintime(lem_use_sus),'o')
% plot(lem.sua_gsac_maxtime(lem_use_sus),lem.sua_msac_maxtime(lem_use_sus),'ro')
% plot(jbe.sua_gsac_mintime(jbe_use_sus),jbe.sua_msac_mintime(jbe_use_sus),'.')
% plot(jbe.sua_gsac_maxtime(jbe_use_sus),jbe.sua_msac_maxtime(jbe_use_sus),'r.')

% figure; hold on
% plot(1-lem.sua_gsac_minvals(lem_use_sus),lem.sua_gsac_mintime(lem_use_sus),'o')
% plot(1-jbe.sua_gsac_minvals(jbe_use_sus),jbe.sua_gsac_mintime(jbe_use_sus),'.')
% plot(lem.sua_gsac_maxvals(lem_use_sus)-1,lem.sua_gsac_maxtime(lem_use_sus),'ro')
% plot(jbe.sua_gsac_maxvals(jbe_use_sus)-1,jbe.sua_gsac_maxtime(jbe_use_sus),'r.')


f2 = figure;
ca = [0.4 1.6];
xl = [-0.1 0.35];
subplot(2,2,1)
% imagesc(lags*dt,1:length(jbe_use_sus),jbe.sua_gsac_gray_avg(:,jbe_use_sus)');
imagesc(lags*dt,1:length(jbe_use_sus),jbe.sua_gsac_sm(:,jbe_use_sus)');
caxis(ca);
xlim(xl);
xlabel('Time (s)');
title('JBE GSac');

subplot(2,2,2);
% imagesc(lags*dt,1:length(lem_use_sus),lem.sua_gsac_gray_avg(:,lem_use_sus)')
imagesc(lags*dt,1:length(lem_use_sus),lem.sua_gsac_sm(:,lem_use_sus)')
caxis(ca);
xlim(xl);
xlabel('Time (s)');
title('LEM Gsac');

ca = [0.5 1.5];
xl = [-0.1 0.35];
subplot(2,2,3)
% imagesc(lags*dt,1:length(jbe_use_sus),jbe.sua_msac_gray_avg(:,jbe_use_sus)');
imagesc(lags*dt,1:length(jbe_use_sus),jbe.sua_msac_sm(:,jbe_use_sus)');
caxis(ca);
xlim(xl);
xlabel('Time (s)');
title('JBE MSac');

subplot(2,2,4);
% imagesc(lags*dt,1:length(lem_use_sus),lem.sua_msac_gray_avg(:,lem_use_sus)')
imagesc(lags*dt,1:length(lem_use_sus),lem.sua_msac_sm(:,lem_use_sus)')
caxis(ca);
xlim(xl);
xlabel('Time (s)');
title('LEM MSac');


xl = [-0.15 0.4];
f3 = figure; hold on
h1=shadedErrorBar(lags*dt,mean(lem.sua_gsac_sm(:,lem_use_sus),2),std(lem.sua_gsac_sm(:,lem_use_sus),[],2)/sqrt(length(lem_use_sus)),{'color','r'});
h2=shadedErrorBar(lags*dt,mean(lem.sua_msac_sm(:,lem_use_sus),2),std(lem.sua_msac_sm(:,lem_use_sus),[],2)/sqrt(length(lem_use_sus)),{'color','b'});
h3=shadedErrorBar(lags*dt,mean(jbe.sua_gsac_sm(:,jbe_use_sus),2),std(jbe.sua_gsac_sm(:,jbe_use_sus),[],2)/sqrt(length(jbe_use_sus)),{'color','k'});
h4=shadedErrorBar(lags*dt,mean(jbe.sua_msac_sm(:,jbe_use_sus),2),std(jbe.sua_msac_sm(:,jbe_use_sus),[],2)/sqrt(length(jbe_use_sus)),{'color','m'});
legend([h1.mainLine h2.mainLine h3.mainLine h4.mainLine],{'LEM Gsac','LEM Msac','JBE Gsac','JBE Gsac'});
xlim(xl);
yl = ylim();
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time (s)');
ylabel('Relative rate');

xl = [-0.15 0.4];
f4 = figure; hold on
h1=shadedErrorBar(lags*dt,mean(lem.sua_gsac_sm(:,lem_use_sus),2),std(lem.sua_gsac_sm(:,lem_use_sus),[],2)/sqrt(length(lem_use_sus)),{'color','r'});
h2=shadedErrorBar(lags*dt,mean(lem.sua_simsac_sm(:,lem_use_sus),2),std(lem.sua_simsac_sm(:,lem_use_sus),[],2)/sqrt(length(lem_use_sus)),{'color','b'});
h3=shadedErrorBar(lags*dt,mean(jbe.sua_gsac_sm(:,jbe_use_sus),2),std(jbe.sua_gsac_sm(:,jbe_use_sus),[],2)/sqrt(length(jbe_use_sus)),{'color','k'});
h4=shadedErrorBar(lags*dt,nanmean(jbe.sua_simsac_sm(:,jbe_use_sus),2),nanstd(jbe.sua_simsac_sm(:,jbe_use_sus),[],2)/sqrt(length(jbe_use_sus)),{'color','m'});
legend([h1.mainLine h2.mainLine h3.mainLine h4.mainLine],{'LEM Gsac','LEM Simsac','JBE Gsac','JBE Simsac'});
xlim(xl);
yl = ylim();
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time (s)');
ylabel('Relative rate');


xl = [-0.15 0.4];
f5 = figure; hold on
h1=shadedErrorBar(lags*dt,mean(lem.sua_gsac_sm(:,lem_use_sus),2),std(lem.sua_gsac_sm(:,lem_use_sus),[],2)/sqrt(length(lem_use_sus)),{'color','r'});
h2=shadedErrorBar(lags*dt,mean(lem.sua_gsac_im_sm(:,lem_use_sus),2),std(lem.sua_gsac_im_sm(:,lem_use_sus),[],2)/sqrt(length(lem_use_sus)),{'color','b'});
h3=shadedErrorBar(lags*dt,mean(jbe.sua_gsac_sm(:,jbe_use_sus),2),std(jbe.sua_gsac_sm(:,jbe_use_sus),[],2)/sqrt(length(jbe_use_sus)),{'color','k'});
h4=shadedErrorBar(lags*dt,nanmean(jbe.sua_gsac_im_sm(:,jbe_use_sus),2),nanstd(jbe.sua_gsac_im_sm(:,jbe_use_sus),[],2)/sqrt(length(jbe_use_sus)),{'color','m'});
legend([h1.mainLine h2.mainLine h3.mainLine h4.mainLine],{'LEM Gsac','LEM Isac','JBE Gsac','JBE Isac'});
xlim(xl);
yl = ylim();
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time (s)');
ylabel('Relative rate');

%%
fig_width = 4; rel_height = 0.8;
figufy(f1);
fname = [fig_dir 'SUA_Enh_Sup_timing.pdf'];
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f1);

fig_width = 8; rel_height = 1;
figufy(f2);
fname = [fig_dir 'SUA_Allmod.pdf'];
exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f2);

fig_width = 4; rel_height = 0.8;
figufy(f3);
fname = [fig_dir 'SUA_avgs.pdf'];
exportfig(f3,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f3);

fig_width = 4; rel_height = 0.8;
figufy(f4);
fname = [fig_dir 'SUA_simsac_compare.pdf'];
exportfig(f4,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f4);

fig_width = 4; rel_height = 0.8;
figufy(f5);
fname = [fig_dir 'SUA_IMsac_compare.pdf'];
exportfig(f5,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f5);

%% FOR MUA
close all
min_nspks = 1e4;
min_rate = 10;
lem_use_mus = find(lem.mua_tot_nspikes >= min_nspks & lem.mua_avg_rates/dt >= min_rate);
jbe_use_mus = find(jbe.mua_tot_nspikes >= min_nspks & jbe.mua_avg_rates/dt >= min_rate);

jit_amp = 0.002;
f1 = figure; hold on
plot(lem.mua_gsac_mintime(lem_use_mus)+randn(length(lem_use_mus),1)*jit_amp,lem.mua_gsac_maxtime(lem_use_mus)+randn(length(lem_use_mus),1)*jit_amp,'.')
plot(jbe.mua_gsac_mintime(jbe_use_mus)+randn(length(jbe_use_mus),1)*jit_amp,jbe.mua_gsac_maxtime(jbe_use_mus)+randn(length(jbe_use_mus),1)*jit_amp,'r.')
line([0 0.3],[0 0.3],'color','k');
xlim([0 0.3]); ylim([0 0.3]);
xlabel('Peak suppression time (s)');
ylabel('Peak enhancement time (s)');

% figure; hold on
% plot(1-lem.mua_gsac_minvals(lem_use_mus),lem.mua_gsac_maxvals(lem_use_mus)-1,'o')
% plot(1-jbe.mua_gsac_minvals(jbe_use_mus),jbe.mua_gsac_maxvals(jbe_use_mus)-1,'.')
% 
% figure; hold on
% plot(lem.mua_gsac_mintime(lem_use_mus),lem.mua_msac_mintime(lem_use_mus),'o')
% plot(lem.mua_gsac_maxtime(lem_use_mus),lem.mua_msac_maxtime(lem_use_mus),'ro')
% plot(jbe.mua_gsac_mintime(jbe_use_mus),jbe.mua_msac_mintime(jbe_use_mus),'.')
% plot(jbe.mua_gsac_maxtime(jbe_use_mus),jbe.mua_msac_maxtime(jbe_use_mus),'r.')

% figure; hold on
% plot(1-lem.mua_gsac_minvals(lem_use_mus),lem.mua_gsac_mintime(lem_use_mus),'o')
% plot(1-jbe.mua_gsac_minvals(jbe_use_mus),jbe.mua_gsac_mintime(jbe_use_mus),'.')
% plot(lem.mua_gsac_maxvals(lem_use_mus)-1,lem.mua_gsac_maxtime(lem_use_mus),'ro')
% plot(jbe.mua_gsac_maxvals(jbe_use_mus)-1,jbe.mua_gsac_maxtime(jbe_use_mus),'r.')

% figure; hold on
% plot(lem.mua_gsac_mintime(lem_use_mus),lem.mua_gsac_maxtime(lem_use_mus),'o')
% plot(jbe.mua_gsac_mintime(jbe_use_mus),jbe.mua_gsac_maxtime(jbe_use_mus),'.')

f2 = figure;
ca = [0.4 1.6];
xl = [-0.1 0.35];
subplot(2,2,1)
% imagesc(lags*dt,1:length(jbe_use_mus),jbe.mua_gsac_gray_avg(:,jbe_use_mus)');
imagesc(lags*dt,1:length(jbe_use_mus),jbe.mua_gsac_sm(:,jbe_use_mus)');
caxis(ca);
xlim(xl);
xlabel('Time (s)');
title('JBE GSac');
subplot(2,2,2);
% imagesc(lags*dt,1:length(lem_use_mus),lem.mua_gsac_gray_avg(:,lem_use_mus)')
imagesc(lags*dt,1:length(lem_use_mus),lem.mua_gsac_sm(:,lem_use_mus)')
caxis(ca);
xlim(xl);
xlabel('Time (s)');
title('LEM GSac');

ca = [0.5 1.5];
xl = [-0.1 0.35];
subplot(2,2,3)
% imagesc(lags*dt,1:length(jbe_use_mus),jbe.mua_msac_gray_avg(:,jbe_use_mus)');
imagesc(lags*dt,1:length(jbe_use_mus),jbe.mua_msac_sm(:,jbe_use_mus)');
caxis(ca);
xlim(xl);
xlabel('Time (s)');
title('JBE MSac');
subplot(2,2,4);
% imagesc(lags*dt,1:length(lem_use_mus),lem.mua_msac_gray_avg(:,lem_use_mus)')
imagesc(lags*dt,1:length(lem_use_mus),lem.mua_msac_sm(:,lem_use_mus)')
caxis(ca);
xlim(xl);
xlabel('Time (s)');
title('LEM MSac');


xl = [-0.15 0.4];
f3 = figure; hold on
h1=shadedErrorBar(lags*dt,mean(lem.mua_gsac_sm(:,lem_use_mus),2),std(lem.mua_gsac_sm(:,lem_use_mus),[],2)/sqrt(length(lem_use_mus)),{'color','r'});
h2=shadedErrorBar(lags*dt,mean(lem.mua_msac_sm(:,lem_use_mus),2),std(lem.mua_msac_sm(:,lem_use_mus),[],2)/sqrt(length(lem_use_mus)),{'color','b'});
h3=shadedErrorBar(lags*dt,mean(jbe.mua_gsac_sm(:,jbe_use_mus),2),std(jbe.mua_gsac_sm(:,jbe_use_mus),[],2)/sqrt(length(jbe_use_mus)),{'color','k'});
h4=shadedErrorBar(lags*dt,mean(jbe.mua_msac_sm(:,jbe_use_mus),2),std(jbe.mua_msac_sm(:,jbe_use_mus),[],2)/sqrt(length(jbe_use_mus)),{'color','m'});
legend([h1.mainLine h2.mainLine h3.mainLine h4.mainLine],{'LEM Gsac','LEM Msac','JBE Gsac','JBE Msac'});
xlim(xl);
yl = ylim();
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time (s)');
ylabel('Relative rate');


xl = [-0.15 0.4];
f4 = figure; hold on
h1=shadedErrorBar(lags*dt,mean(lem.mua_gsac_sm(:,lem_use_mus),2),std(lem.mua_gsac_sm(:,lem_use_mus),[],2)/sqrt(length(lem_use_mus)),{'color','r'});
h2=shadedErrorBar(lags*dt,mean(lem.mua_simsac_sm(:,lem_use_mus),2),std(lem.mua_simsac_sm(:,lem_use_mus),[],2)/sqrt(length(lem_use_mus)),{'color','b'});
h3=shadedErrorBar(lags*dt,mean(jbe.mua_gsac_sm(:,jbe_use_mus),2),std(jbe.mua_gsac_sm(:,jbe_use_mus),[],2)/sqrt(length(jbe_use_mus)),{'color','k'});
h4=shadedErrorBar(lags*dt,nanmean(jbe.mua_simsac_sm(:,jbe_use_mus),2),nanstd(jbe.mua_simsac_sm(:,jbe_use_mus),[],2)/sqrt(length(jbe_use_mus)),{'color','m'});
legend([h1.mainLine h2.mainLine h3.mainLine h4.mainLine],{'LEM Gsac','LEM Simsac','JBE Gsac','JBE Simsac'});
xlim(xl);
yl = ylim();
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time (s)');
ylabel('Relative rate');


xl = [-0.15 0.4];
f5 = figure; hold on
h1=shadedErrorBar(lags*dt,mean(lem.mua_gsac_sm(:,lem_use_mus),2),std(lem.mua_gsac_sm(:,lem_use_mus),[],2)/sqrt(length(lem_use_mus)),{'color','r'});
h2=shadedErrorBar(lags*dt,mean(lem.mua_gsac_im_sm(:,lem_use_mus),2),std(lem.mua_gsac_im_sm(:,lem_use_mus),[],2)/sqrt(length(lem_use_mus)),{'color','b'});
h3=shadedErrorBar(lags*dt,mean(jbe.mua_gsac_sm(:,jbe_use_mus),2),std(jbe.mua_gsac_sm(:,jbe_use_mus),[],2)/sqrt(length(jbe_use_mus)),{'color','k'});
h4=shadedErrorBar(lags*dt,nanmean(jbe.mua_gsac_im_sm(:,jbe_use_mus),2),nanstd(jbe.mua_gsac_im_sm(:,jbe_use_mus),[],2)/sqrt(length(jbe_use_mus)),{'color','m'});
legend([h1.mainLine h2.mainLine h3.mainLine h4.mainLine],{'LEM Gsac','LEM IMsac','JBE Gsac','JBE IMsac'});
xlim(xl);
yl = ylim();
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time (s)');
ylabel('Relative rate');

%%
fig_width = 4; rel_height = 0.8;
figufy(f1);
fname = [fig_dir 'MUA_Enh_Sup_timing.pdf'];
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f1);

fig_width = 8; rel_height = 1;
figufy(f2);
fname = [fig_dir 'MUA_Allmod.pdf'];
exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f2);

fig_width = 4; rel_height = 0.8;
figufy(f3);
fname = [fig_dir 'MUA_avgs.pdf'];
exportfig(f3,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f3);

fig_width = 4; rel_height = 0.8;
figufy(f4);
fname = [fig_dir 'MUA_simsac_compare.pdf'];
exportfig(f4,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f4);


fig_width = 4; rel_height = 0.8;
figufy(f5);
fname = [fig_dir 'MUA_IMsac_compare.pdf'];
exportfig(f5,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f5);

%%
un_lem_exptnums = unique(lem_mua_exptnum);
for ii = 1:length(un_lem_exptnums)
    cur_set = find(lem_mua_exptnum == un_lem_exptnums(ii));
    expt_avg_gsac_mua(ii,:) = mean(lem.mua_gsac_sm(:,cur_set),2);
    expt_avg_msac_mua(ii,:) = mean(lem.mua_msac_sm(:,cur_set),2);
end

[~,expt_Ecc_ord] = sort(Expt_RFPos);
expt_isfov = Expt_RFPos < 2;
expt_isperi = Expt_RFPos > 4;
f1=figure; hold on
cmap = jet(length(Expt_RFPos));
for ii = 1:length(Expt_RFPos)
    %     plot(lags*dt,expt_avg_gsac_mua(expt_Ecc_ord(ii),:),'color',cmap(ii,:),'linewidth',2);
    if expt_isfov(ii)
        plot(lags*dt,expt_avg_gsac_mua(ii,:),'b','linewidth',2);
    else
        plot(lags*dt,expt_avg_gsac_mua(ii,:),'r','linewidth',1);
    end
    if expt_isperi(ii)
        plot(lags*dt,expt_avg_gsac_mua(ii,:),'g','linewidth',1);
        
    end
end
plot(lags*dt,nanmean(jbe.mua_gsac_sm,2),'k--','linewidth',4);
xl = [-0.15 0.5];
xlim(xl);
xlabel('Time (s)');
ylabel('Relative rate');
yl = ylim();
line([0 0],yl,'color','k','linestyle','--');

% fig_width = 4; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'Expt_avg_fov_compare.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
% 



% figure; hold on
% cmap = jet(length(Expt_RFPos));
% for ii = 1:length(Expt_RFPos)
%     plot(lags*dt,expt_avg_msac_mua(expt_Ecc_ord(ii),:),'color',cmap(ii,:),'linewidth',2);
% end
% plot(lags*dt,mean(jbe.mua_msac_sm,2),'k--','linewidth',4);

%%
% close all
% ca = [0.6 1.4];
% 
% figure;
% imagesc(lags*dt,1:size(mua_gsac_avg,2),mua_gsac_gray_avg');
% xl = xlim();
% n_probes = 24;
% for ii = 1:length(all_mua_data)
%     line(xl,[0.5 0.5]+ii*n_probes,'color','w');
% end
% caxis(ca);
% 
% figure;
% imagesc(lags*dt,1:size(mua_gsac_avg,2),mua_gsac_im_avg');
% xl = xlim();
% n_probes = 24;
% for ii = 1:length(all_mua_data)
%     line(xl,[0.5 0.5]+ii*n_probes,'color','w');
% end
% caxis(ca);
% 
% figure;
% imagesc(lags*dt,1:size(mua_gsac_avg,2),mua_simsac_avg');
% xl = xlim();
% n_probes = 24;
% for ii = 1:length(all_mua_data)
%     line(xl,[0.5 0.5]+ii*n_probes,'color','w');
% end
% caxis(ca);
% 
% figure;
% imagesc(lags*dt,1:size(mua_gsac_avg,2),mua_msac_gray_avg');
% xl = xlim();
% n_probes = 24;
% for ii = 1:length(all_mua_data)
%     line(xl,[0.5 0.5]+ii*n_probes,'color','w');
% end
% caxis(ca);
% 
% %%
% mua_gsac_means = cell2mat(arrayfun(@(x) mean(x.gsac_avg,2),all_mua_data,'uniformoutput',0));
% mua_gsac_end_means = cell2mat(arrayfun(@(x) mean(x.gsac_end_avg,2),all_mua_data,'uniformoutput',0));
