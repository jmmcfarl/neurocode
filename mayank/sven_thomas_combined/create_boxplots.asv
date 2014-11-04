clear all
cd C:\WC_Germany\sven_thomas_combined
load ./combined_dir_nd.mat
uset = sort([l3mec l3lec l3mec_np l3lec_np]);
all_cells = 1:length(combined_dir);
l3mec = find(ismember(all_cells(uset),l3mec));
l3lec = find(ismember(all_cells(uset),l3lec));
l3mec_np = find(ismember(all_cells(uset),l3mec_np));
l3lec_np = find(ismember(all_cells(uset),l3lec_np));
combined_dir = combined_dir(uset);

load ./combined_core_analysis_fin_nd_np
l3mec_np(l3mec_np == 62) = []; %interneuron
l3mec_np(l3mec_np == 59) = [];
l3mec_np(l3mec_np == 60) = [];
int = 62;

%%
ec_rt2_ups = fract_rt2_ups;
% ec_downlag = median_downlag;
% ec_uplag = median_uplag;
ec_downlag = mean_downlag;
ec_uplag = mean_uplag;
load ./combined_core_analysis_dg.mat
dg = (length(uset)+1):(length(uset)+10);
% fract_rt2_ups = [ec_rt2_ups fract_rt2_ups];
% median_downlag = [ec_downlag median_downlag];
% median_uplag = [ec_uplag median_uplag];
% [use,ord] = sort([l3mec l3mec_np l3lec l3lec_np int dg]);

fract_rt2_ups = [ec_rt2_ups];
% median_downlag = [ec_downlag];
mean_downlag = [ec_downlag];
% median_uplag = [ec_uplag];
mean_uplag = [ec_uplag];
[use,ord] = sort([l3mec l3mec_np l3lec l3lec_np int]);


% [use,ord] = sort([l3mec l3mec_np l3lec l3lec_np int]);
Y = fract_rt2_ups(use);

% G = [ones(length(l3mec),1); 2*ones(length(l3mec_np),1); 4*ones(length(l3lec),1); 5*ones(length(l3lec_np),1); 3; 6*ones(length(dg),1)];
G = [ones(length(l3mec),1); 2*ones(length(l3mec_np),1); 4*ones(length(l3lec),1); 5*ones(length(l3lec_np),1); 3];
G = G(ord);
% Y = Y(ord);
boxplot(Y*100,G);
set(gca,'xtick',1:6, 'xticklabel',{'MEC pyr','MEC non-pyr','MEC int','LEC pyr','LEC non-pyr','DG'})
ylabel('Persistence (%)','fontsize',14)

% Y = median_downlag(use)/Fsd;
Y = mean_downlag(use)/Fsd;
% Y = Y(ord);
figure
boxplot(Y,G);
set(gca,'xtick',1:6, 'xticklabel',{'MEC pyr','MEC non-pyr','MEC int','LEC pyr','LEC non-pyr','DG'})
ylabel('Down-transition lag (s)','fontsize',14)

% Y = median_uplag(use)/Fsd;
Y = mean_uplag(use)/Fsd;
% Y = Y(ord);
figure
boxplot(Y,G);
set(gca,'xtick',1:6, 'xticklabel',{'MEC pyr','MEC non-pyr','MEC int','LEC pyr','LEC non-pyr','DG'})
ylabel('Up-transition lag (s)','fontsize',14)

%%
clear all
load ./combined_dir_nd.mat
uset = sort([l3mec l3lec l3mec_np l3lec_np]);
all_cells = 1:length(combined_dir);
l3mec = find(ismember(all_cells(uset),l3mec));
l3lec = find(ismember(all_cells(uset),l3lec));
l3mec_np = find(ismember(all_cells(uset),l3mec_np));
l3lec_np = find(ismember(all_cells(uset),l3lec_np));
combined_dir = combined_dir(uset);

l3mec_np(l3mec_np == 62) = []; %interneuron
l3mec_np(l3mec_np == 59) = [];
l3mec_np(l3mec_np == 60) = [];
int = 62;
load ./spike_rate_data_fin_nd_np

[use,ord] = sort([l3mec l3mec_np l3lec l3lec_np int]);
Y = mp_up_rate(use);
G = [ones(length(l3mec),1); 2*ones(length(l3mec_np),1); 4*ones(length(l3lec),1); 5*ones(length(l3lec_np),1); 3];
G = G(ord);
% Y = Y(ord)

boxplot(Y,G);
set(gca,'xtick',1:5, 'xticklabel',{'MEC pyr','MEC non-pyr','MEC int','LEC pyr','LEC non-pyr'})
ylabel('UP state firing rate (Hz)','fontsize',14)


%%
clear all
close all
load combined_heka_UDS_data_allcells_fin_nd_np
load ./combined_dir_nd.mat
l3mec_np(ismember(l3mec_np,[64])) = [];
uset = sort([l3mec l3lec l3mec_np l3lec_np]);
all_cells = 1:length(combined_dir);
l3mec = find(ismember(all_cells(uset),l3mec));
l3lec = find(ismember(all_cells(uset),l3lec));
l3mec_np = find(ismember(all_cells(uset),l3mec_np));
l3lec_np = find(ismember(all_cells(uset),l3lec_np));
combined_dir = combined_dir(uset);

l3lec(5) = []; %get rid of this outlier in terms of MP amplitude

l3mec_np(l3mec_np==61) = [];
l3mec_np(l3mec_np == 58) = [];
l3mec_np(l3mec_np == 59) = [];
int = 61;
[use,ord] = sort([l3mec l3mec_np l3lec l3lec_np int]);
Y = upstate_mean_spksub_beg(use)-7;
% Y = upstate_mode_spksub_beg(:);
G = [ones(length(l3mec),1); 2*ones(length(l3mec_np),1); 4*ones(length(l3lec),1); 5*ones(length(l3lec_np),1); 3];
G = G(ord);
boxplot(Y,G);
set(gca,'xtick',1:5, 'xticklabel',{'MEC pyr','MEC non-pyr','MEC int','LEC pyr','LEC non-pyr'})
ylabel('UP state MP (mV)','fontsize',14)

figure
% Y = downstate_mean_spksub_beg(:);
Y = downstate_mode_spksub_beg(use)-7;
boxplot(Y,G);
set(gca,'xtick',1:5, 'xticklabel',{'MEC pyr','MEC non-pyr','MEC int','LEC pyr','LEC non-pyr'})
ylabel('DOWN state MP (mV)','fontsize',14)

figure
% Y = downstate_mean_spksub_beg(:);
Y = upstate_mean_spksub_beg(use) - downstate_mode_spksub_beg(use);
boxplot(Y,G);
set(gca,'xtick',1:5, 'xticklabel',{'MEC pyr','MEC non-pyr','MEC int','LEC pyr','LEC non-pyr'})
ylabel('UDS Amp (mV)','fontsize',14)

%%
cd C:\WC_Germany\sven_thomas_combined
load ./combined_core_analysis_fin_nd_np
load combined_heka_UDS_data_allcells_fin_nd_np
load ./combined_dir_nd.mat
l3mec_np(ismember(l3mec_np,[64])) = [];
uset = sort([l3mec l3lec l3mec_np l3lec_np]);
all_cells = 1:length(combined_dir);
l3mec = find(ismember(all_cells(uset),l3mec));
l3lec = find(ismember(all_cells(uset),l3lec));
l3mec_np = find(ismember(all_cells(uset),l3mec_np));
l3lec_np = find(ismember(all_cells(uset),l3lec_np));
l3mec_np(l3mec_np == 62) = [];
l3mec_np(l3mec_np == 59) = [];
l3mec_np(l3mec_np == 60) = [];
combined_dir = combined_dir(uset);

l3lec(5) = []; %get rid of this outlier in terms of MP

uds_amps = upstate_mean_spksub_beg-downstate_mode_spksub_beg;

figure
plot(upstate_mean_spksub_beg(l3mec)-7,100*fract_rt2_ups(l3mec),'ro')
hold on
plot(upstate_mean_spksub_beg(l3lec)-7,100*fract_rt2_ups(l3lec),'o')
% plot(upstate_mean_spksub_beg(l3mec_np)-7,100*fract_rt2_ups(l3mec_np),'ko')
% plot(upstate_mean_spksub_beg(l3lec_np)-7,100*fract_rt2_ups(l3lec_np),'go')
xlabel('Up state amp (mV)','fontsize',16)
ylabel('Persistence (%)','fontsize',16)
set(gca,'fontsize',14,'fontname','arial')
legend('MEC pyr','LEC pyr','MEC non-pyr','LEC non-pyr')

figure
plot(upstate_mean_spksub_beg(l3mec)-7,100*fract_rt2_ups(l3mec),'ro')
hold on
plot(upstate_mean_spksub_beg(l3lec)-7,100*fract_rt2_ups(l3lec),'o')
% plot(upstate_mean_spksub_beg(l3mec_np)-7,100*fract_rt2_ups(l3mec_np),'ko')
% plot(upstate_mean_spksub_beg(l3lec_np)-7,100*fract_rt2_ups(l3lec_np),'go')
xlabel('Up state amp (mV)','fontsize',16)
ylabel('Persistence (%)','fontsize',16)
set(gca,'fontsize',14,'fontname','arial')
legend('MEC pyr','LEC pyr','MEC non-pyr','LEC non-pyr')

figure
plot(downstate_mode_spksub_beg(l3mec)-7,uds_amps(l3mec),'ro')
hold on
plot(downstate_mode_spksub_beg(l3lec)-7,uds_amps(l3lec),'o')
xlabel('Down state amp (mV)','fontsize',16)
ylabel('UDS amp (mV)','fontsize',16)
set(gca,'fontsize',14,'fontname','arial')

figure
plot(upstate_mean_spksub_beg(l3mec)-7,uds_amps(l3mec),'ro')
hold on
plot(upstate_mean_spksub_beg(l3lec)-7,uds_amps(l3lec),'o')
xlabel('Up state amp (mV)','fontsize',16)
ylabel('UDS amp (mV)','fontsize',16)
set(gca,'fontsize',14,'fontname','arial')

figure
plot(uds_amps(l3mec),100*fract_rt2_ups(l3mec),'ro')
hold on
plot(uds_amps(l3lec),100*fract_rt2_ups(l3lec),'o')
% plot(uds_amps(l3mec_np),100*fract_rt2_ups(l3mec_np),'ko')
% plot(uds_amps(l3lec_np),100*fract_rt2_ups(l3lec_np),'go')
xlabel('UDS amp amp (mV)','fontsize',16)
ylabel('Persistence (%)','fontsize',16)
set(gca,'fontsize',14,'fontname','arial')

load ./spike_rate_data_fin_nd_np mp_up_rate
figure
plot(mp_up_rate(l3mec),100*fract_rt2_ups(l3mec),'ro')
hold on
plot(mp_up_rate(l3lec),100*fract_rt2_ups(l3lec),'o')
% plot(mp_up_rate(l3mec_np),100*fract_rt2_ups(l3mec_np),'ko')
% plot(mp_up_rate(l3lec_np),100*fract_rt2_ups(l3lec_np),'go')
xlabel('Up state firing rate (Hz)','fontsize',16)
ylabel('Persistence (%)','fontsize',16)
set(gca,'fontsize',14,'fontname','arial')
legend('MEC pyr','LEC pyr','MEC non-pyr','LEC non-pyr')

figure
plot(mp_up_rate(l3mec),mean_downlag(l3mec)/Fsd,'ro')
hold on
plot(mp_up_rate(l3lec),mean_downlag(l3lec)/Fsd,'o')
% plot(mp_up_rate(l3mec_np),mean_downlag(l3mec_np)/Fsd,'ko')
% plot(mp_up_rate(l3lec_np),mean_downlag(l3lec_np)/Fsd,'go')
xlabel('Up state firing rate (Hz)','fontsize',16)
ylabel('Down transition lag (s)','fontsize',16)
set(gca,'fontsize',14,'fontname','arial')
legend('MEC pyr','LEC pyr','MEC non-pyr','LEC non-pyr')

%%
figure
plot(upstate_mean_spksub_beg(l3mec)-7,mean_downlag(l3mec)/Fsd,'ro')
hold on
plot(upstate_mean_spksub_beg(l3lec)-7,mean_downlag(l3lec)/Fsd,'o')
% plot(upstate_mean_spksub_beg(l3mec_np)-7,mean_downlag(l3mec_np)/Fsd,'ko')
% plot(upstate_mean_spksub_beg(l3lec_np)-7,mean_downlag(l3lec_np)/Fsd,'go')
xlabel('Up state amp (mV)','fontsize',16)
ylabel('Down-transition lag (s)','fontsize',16)
set(gca,'fontsize',14,'fontname','arial')
legend('MEC pyr','LEC pyr','MEC non-pyr','LEC non-pyr')

figure
plot(downstate_mode_spksub_beg(l3mec)-7,mean_downlag(l3mec)/Fsd,'ro')
hold on
plot(downstate_mode_spksub_beg(l3lec)-7,mean_downlag(l3lec)/Fsd,'o')
% plot(downstate_mode_spksub_beg(l3mec_np)-7,mean_downlag(l3mec_np)/Fsd,'ko')
% plot(downstate_mode_spksub_beg(l3lec_np)-7,mean_downlag(l3lec_np)/Fsd,'go')
xlabel('Down state amp (mV)','fontsize',16)
ylabel('Down-transition lag (s)','fontsize',16)
set(gca,'fontsize',14,'fontname','arial')
legend('MEC pyr','LEC pyr','MEC non-pyr','LEC non-pyr')

figure
plot(uds_amps(l3mec),mean_downlag(l3mec)/Fsd,'ro')
hold on
plot(uds_amps(l3lec),mean_downlag(l3lec)/Fsd,'o')
xlabel('UDS amp amp (mV)','fontsize',16)
ylabel('Down-transition lag (s)','fontsize',16)
set(gca,'fontsize',14,'fontname','arial')

%%
figure
plot(upstate_mean_spksub_beg(l3mec)-7,median_uplag(l3mec)/Fsd,'o')
hold on
plot(upstate_mean_spksub_beg(l3lec)-7,median_uplag(l3lec)/Fsd,'ro')
plot(upstate_mean_spksub_beg(l3mec_np)-7,median_uplag(l3mec_np)/Fsd,'ko')
plot(upstate_mean_spksub_beg(l3lec_np)-7,median_uplag(l3lec_np)/Fsd,'go')
xlabel('Up state amp (mV)','fontsize',16)
ylabel('Down-transition lag (s)','fontsize',16)
set(gca,'fontsize',14,'fontname','arial')
legend('MEC pyr','LEC pyr','MEC non-pyr','LEC non-pyr')

figure
plot(downstate_mode_spksub_beg(l3mec)-7,median_uplag(l3mec)/Fsd,'o')
hold on
plot(downstate_mode_spksub_beg(l3lec)-7,median_uplag(l3lec)/Fsd,'ro')
plot(downstate_mode_spksub_beg(l3mec_np)-7,median_uplag(l3mec_np)/Fsd,'ko')
plot(downstate_mode_spksub_beg(l3lec_np)-7,median_uplag(l3lec_np)/Fsd,'go')
xlabel('Down state amp (mV)','fontsize',16)
ylabel('Down-transition lag (s)','fontsize',16)
set(gca,'fontsize',14,'fontname','arial')
legend('MEC pyr','LEC pyr','MEC non-pyr','LEC non-pyr')
