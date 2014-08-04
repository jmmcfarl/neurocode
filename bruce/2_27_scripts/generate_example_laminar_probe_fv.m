cd ~/Data/bruce/2_27_12/
load Blocks.mat

load ./examp_data

suprobes = Blocks{1}.suprobes;
muprobes = Blocks{1}.muprobes;
sutimes = Blocks{1}.spktimes;
mutimes = Blocks{1}.mutimes;
all_probes = [suprobes muprobes];
all_spktimes = [sutimes mutimes];
[~,probe_ord] = sort(all_probes);
all_spktimes = all_spktimes(probe_ord);

%%
close all
sup_ch = 4;
gran_ch = 14;
inf_ch = 23;

sup_unit = 4;
gran_unit = 14;
inf_unit = 23;

offset = 4;
tick_size = 0.6;
xr = [41 44]; %42 44

figure

subplot(4,1,[1 2])
plot(lfp_timed,lfp_sampsd(:,sup_ch)+offset)
hold on
plot(lfp_timed,lfp_sampsd(:,gran_ch),'r')
plot(lfp_timed,lfp_sampsd(:,inf_ch)-offset,'k')
xlim(xr)
set(gca,'ytick',[])
xlabel('Time (s)','fontsize',16)

subplot(4,1,3)
cur_spktimes = all_spktimes{sup_unit};
cur_spktimes = cur_spktimes(cur_spktimes >= xr(1) & cur_spktimes <= xr(2));
for ii = 1:length(cur_spktimes)
    line(cur_spktimes([ii ii]),1 + [0 tick_size])
end
cur_spktimes = all_spktimes{gran_unit};
cur_spktimes = cur_spktimes(cur_spktimes >= xr(1) & cur_spktimes <= xr(2));
for ii = 1:length(cur_spktimes)
    line(cur_spktimes([ii ii]),0 + [0 tick_size],'color','r')
end
cur_spktimes = all_spktimes{inf_unit};
cur_spktimes = cur_spktimes(cur_spktimes >= xr(1) & cur_spktimes <= xr(2));
for ii = 1:length(cur_spktimes)
    line(cur_spktimes([ii ii]),-1 + [0 tick_size],'color','k')
end
xlim(xr)
ylim([-1.2 1+tick_size + 0.2])
set(gca,'ytick',[])
xlabel('Time (s)','fontsize',16)

subplot(4,1,4)
plot(eyets,eye_speed,'r')
xlim(xr)
xlabel('Time (s)','fontsize',16)
ylabel('Eye speed (deg/s)','fontsize',16)


%%
%%
close all
sup_ch = 4;
gran_ch = 14;
inf_ch = 23;

use_units = 1:24;

offset = 4;
tick_size = 0.8;
xr = [41 44]; %42 44

figure

subplot(4,1,[1 2])
plot(lfp_timed,lfp_sampsd(:,sup_ch)+offset)
hold on
plot(lfp_timed,lfp_sampsd(:,gran_ch),'r')
plot(lfp_timed,lfp_sampsd(:,inf_ch)-offset,'k')
xlim(xr)
set(gca,'ytick',[])
xlabel('Time (s)','fontsize',16)

subplot(4,1,3)
cnt = 0;
for cc = 1:length(use_units)
cur_spktimes = all_spktimes{use_units(cc)};
cur_spktimes = cur_spktimes(cur_spktimes >= xr(1) & cur_spktimes <= xr(2));
for ii = 1:length(cur_spktimes)
    line(cur_spktimes([ii ii]),cnt + [0 tick_size],'color','k')
end
cnt = cnt + 1;
end
xlim(xr)
ylim([-0.2 cnt + 0.2])
set(gca,'ytick',[])
xlabel('Time (s)','fontsize',16)

subplot(4,1,4)
plot(eyets,eye_speed,'r')
xlim(xr)
xlabel('Time (s)','fontsize',16)G
ylabel('Eye speed (deg/s)','fontsize',16)

%%
t_range = [80 83];



target_eye_inds = find(eyets >= t_range(1) & eyets < t_range(2));
target_lfp_inds = find(lfp_timed >= t_range(1) & lfp_timed < t_range(2));

use_lfp_ch = 2;

f1 = figure;
subplot(2,1,1)
plot(eyets(target_eye_inds),eye_speed(target_eye_inds))
xlabel('Time (s)','fontsize',16)
ylabel('Eye speed (deg/s)','fontsize',16)
title('Eye speed','fontsize',18)
subplot(2,1,2)
plot(lfp_timed(target_lfp_inds),lfp_sampsd(target_lfp_inds,use_lfp_ch))
xlabel('Time (s)','fontsize',16)
ylabel('Amplitude','fontsize',16)
title('LFP signal','fontsize',18)

f2 = figure;
subplot(2,1,1)
pcolor(lfp_timed(target_lfp_inds),wfreqs,squeeze(all_filt_bank(use_lfp_ch,target_lfp_inds,:))');shading interp
set(gca,'yscale','log');
xlabel('Time (s)','fontsize',16)
ylabel('Frequency (Hz)','fontsize',16)
title('Wavelet transformed LFP (real component)','fontsize',18)
subplot(2,1,2)
pcolor(lfp_timed(target_lfp_inds),wfreqs,squeeze(all_cwt_amp(use_lfp_ch,target_lfp_inds,:))');shading interp
set(gca,'yscale','log')
xlabel('Time (s)','fontsize',16)
ylabel('Frequency (Hz)','fontsize',16)
title('Wavelet amplitude','fontsize',18)
