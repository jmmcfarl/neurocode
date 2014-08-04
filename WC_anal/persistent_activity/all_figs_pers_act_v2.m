%% generate figures

%% Amplitude Histogram Example
clear all
close all
load C:\WC_Germany\Persistent_activity\heka_amp\heka_amp_data

example_cell = 5;

figure
bar(x{example_cell},y{example_cell},'k')
xlabel('Membrane Potential (mV)','FontSize',14)
ylabel('Probability Density','FontSize',14)
xlim([-70 0])
cd C:\WC_Germany\Persistent_activity\potential_figures

%% Example neuralynx mp trace
clear all
% close all
load C:\WC_Germany\Persistent_activity\dir_tree_update

%specify session number
d = 14;
niqf = 2016/2;
lcf = 0.05/niqf;
hcf = 4/niqf;
[b,a] = butter(2,[lcf hcf]);

dsf = 8;
Fsd = 2016/dsf;

cd(dir_array{d})
pwd

load used_data lf8 wcv wcv_minus_spike
load spike_time_jmm

spike_ids = round(spkid/dsf);

lf8_f = filtfilt(b,a,lf8);
wcv_f = filtfilt(b,a,wcv_minus_spike);

lf8_d = downsample(lf8_f,dsf);
wcv_d = downsample(wcv_f,dsf);

lf8_d = zscore(lf8_d);
wcv_d = zscore(wcv_d);

wcv = zscore(wcv);
lf8 = zscore(lf8);
t = (1:length(wcv_d))/2016*dsf;

plot(t-66,wcv_d,'linewidth',3)
xlabel('Time (s)','FontSize',14)
ylabel('Amplitude (zscore)','FontSize',14)
title('MEC Membrane Potential','FontSize',16)
hold on
plot(t-66,lf8_d,'r','linewidth',3)
xlabel('Time (s)','FontSize',14)
ylabel('Amplitude (zscore)','FontSize',14)
title('Cortical LFP','FontSize',16)
legend('MP','LFP')
xlim([66-66 82-66])
ylim([-4 4])

used_spikes = find(spike_ids/Fsd > 66 & spike_ids/Fsd < 82);
% plot(t(spike_ids)-66,ones(size(spike_ids))*3,'ro')
line_bottom = 2.8;
line_top = 3.2;
for i = 1:length(used_spikes)
   line([t(spike_ids(used_spikes(i)))-66 t(spike_ids(used_spikes(i)))-66],...
       [line_bottom line_top],'Color','k','linewidth',2) 
end
shg

% figure
% plot(t,wcv_d,'linewidth',3)
% xlabel('Time (s)','FontSize',14)
% ylabel('Amplitude (zscore)','FontSize',14)
% title('MEC Membrane Potential','FontSize',16)
% hold on
% plot(t,lf8_d,'r','linewidth',3)
% xlabel('Time (s)','FontSize',14)
% ylabel('Amplitude (zscore)','FontSize',14)
% title('Cortical LFP','FontSize',16)
% legend('MP','LFP')
% xlim([0 70])
% ylim([-4 4])
% shg


%% Example neuralynx mp trace
clear all

load C:\WC_Germany\Persistent_activity\dir_tree_update

%specify session number
d = 14;

dsf = 1;
Fsd = 2016/dsf;

cd(dir_array{d})
pwd

load used_data lf8 wcv wcv_minus_spike

wcv_d = zscore(wcv);

t = (1:length(wcv_d))/2016*dsf;

% plot(t-66,wcv_d,'linewidth',3)
% xlabel('Time (s)','FontSize',14)
% ylabel('Amplitude (zscore)','FontSize',14)
% title('MEC Membrane Potential','FontSize',16)
% hold on
% plot(t-66,lf8_d,'r','linewidth',3)
% xlabel('Time (s)','FontSize',14)
% ylabel('Amplitude (zscore)','FontSize',14)
% title('Cortical LFP','FontSize',16)
% legend('MP','LFP')
% xlim([66-66 84-66])
% ylim([-4 4])
% shg

figure
plot(t-94,wcv_d,'linewidth',3)
xlabel('Time (s)','FontSize',14)
ylabel('Amplitude (zscore)','FontSize',14)
title('MEC Membrane Potential','FontSize',16)
hold on
xlim([94-94 115-94])
ylim([-3 8])
shg


%% Example MP Trace
close all
clear all
load C:\WC_Germany\Persistent_activity\pyr_heka_dir

%specify session number
d = 14;

load(f_loc{d})
    dat_name = [f_loc{d} '_MP'];
    dat_name(1:24) = [];
i=7    

zscore_data = zscore(A2007_07_04_CWC_LFP_C_MP);

timing_offset = 3;
    max_time = max(A2007_07_04_CWC_LFP_C_sampletimes);
 sweep_dur=find(A2007_07_04_CWC_LFP_C_sampletimes ==max_time,1,'first');
 plot(A2007_07_04_CWC_LFP_C_sampletimes(1+sweep_dur*(i-1):i*sweep_dur)-timing_offset,...
     100*A2007_07_04_CWC_LFP_C_MP(1+sweep_dur*(i-1):i*sweep_dur),'linewidth',2)
 xlim([0 5])
 ylim([-80 10])
 
 figure
  plot(A2007_07_04_CWC_LFP_C_sampletimes(1+sweep_dur*(i-1):i*sweep_dur)-timing_offset,...
     zscore_data(1+sweep_dur*(i-1):i*sweep_dur),'linewidth',2)
 xlim([0 5])

 
 xlabel('Time (s)','FontSize',14)
ylabel('Amplitude (mV)','FontSize',14)
 cd C:\WC_Germany\Persistent_activity\potential_figures
 
 %% Example MP Trace
close all
clear all
load C:\WC_Germany\Persistent_activity\pyr_heka_dir

%specify session number
d = 14;

load(f_loc{d})
    dat_name = [f_loc{d} '_MP'];
    dat_name(1:24) = [];
i=7    

zscore_data = zscore(A2007_07_04_CWC_LFP_C_MP);

timing_offset = 3;
    max_time = max(A2007_07_04_CWC_LFP_C_sampletimes);
 sweep_dur=find(A2007_07_04_CWC_LFP_C_sampletimes ==max_time,1,'first');
 plot(A2007_07_04_CWC_LFP_C_sampletimes(1+sweep_dur*(i-1):i*sweep_dur)-timing_offset,...
     100*A2007_07_04_CWC_LFP_C_MP(1+sweep_dur*(i-1):i*sweep_dur),'linewidth',2)
 xlim([0 5])
 ylim([-80 10])
 
 figure
  plot(A2007_07_04_CWC_LFP_C_sampletimes(1+sweep_dur*(i-1):i*sweep_dur)-timing_offset,...
     zscore_data(1+sweep_dur*(i-1):i*sweep_dur),'linewidth',2)
 xlim([0 5])

 
 xlabel('Time (s)','FontSize',14)
ylabel('Amplitude (mV)','FontSize',14)
 cd C:\WC_Germany\Persistent_activity\potential_figures

%% Example MP Trace2
close all
clear all
load C:\WC_Germany\Persistent_activity\pyr_heka_dir

%specify session number
d = 5;

load(f_loc{d})
    dat_name = [f_loc{d} '_MP'];
    dat_name(1:24) = [];
for i=1:200   

zscore_data = zscore(A2007_05_30_CWC_LFP_C_MP);

timing_offset = 3;
    max_time = max(A2007_05_30_CWC_LFP_C_sampletimes);
 sweep_dur=find(A2007_05_30_CWC_LFP_C_sampletimes ==max_time,1,'first');
 plot(A2007_05_30_CWC_LFP_C_sampletimes(1+sweep_dur*(i-1):i*sweep_dur)-timing_offset,...
     100*A2007_05_30_CWC_LFP_C_MP(1+sweep_dur*(i-1):i*sweep_dur),'linewidth',2)
 ylim([-65 10])
 pause
 clf
end
 figure
  plot(A2007_05_30_CWC_LFP_C_sampletimes(1+sweep_dur*(i-1):i*sweep_dur)-timing_offset,...
     zscore_data(1+sweep_dur*(i-1):i*sweep_dur),'linewidth',2)

 
 xlabel('Time (s)','FontSize',14)
ylabel('Amplitude (mV)','FontSize',14)
 cd C:\WC_Germany\Persistent_activity\potential_figures

%% All Distributions LFP 

clear all
close all
load C:\WC_Germany\Persistent_activity\neural_amp\low_pass_data

% %for individuals and mean
% plot(gridv,lf8_dist')
% hold on
% plot(gridv,mean(lf8_dist),'r','linewidth',3)
% xlim([-3 3])
% xlabel('Amplitude (z)','FontSize',14)
% ylabel('Probability Density','FontSize',14)
% cd C:\WC_Germany\Persistent_activity\potential_figures

%for mean and std error
m_dist = mean(lf8_dist);
u_dist = m_dist+2.1*std(lf8_dist)/sqrt(17);
l_dist = m_dist-2.1*std(lf8_dist)/sqrt(17);
plot(gridv,m_dist,'r','linewidth',2)
hold on
plot(gridv,u_dist,'r--')
plot(gridv,l_dist,'r--')
xlabel('Amplitude (z)','FontSize',14)
ylabel('Probability Density','FontSize',14)
xlim([-3 3])
ylim([0 0.7])

%% All Distributions MP

clear all
close all
load C:\WC_Germany\Persistent_activity\neural_amp\low_pass_data_3

% %for individuals and mean
% plot(gridv,wcv_dist')
% hold on
% plot(gridv,mean(wcv_dist),'linewidth',3)
% xlim([-3 3])
% xlabel('Amplitude (z)','FontSize',14)
% ylabel('Probability Density','FontSize',14)
% ylim([0 1.1])
% xlim([-3 3])
% cd C:\WC_Germany\Persistent_activity\potential_figures

%for mean and std error
m_dist = mean(wcv_dist);
u_dist = m_dist+1*std(wcv_dist)/sqrt(17);
l_dist = m_dist-1*std(wcv_dist)/sqrt(17);
plot(gridv,m_dist,'linewidth',2)
hold on
plot(gridv,u_dist,'--')
plot(gridv,l_dist,'--')
xlabel('Amplitude (z)','FontSize',14)
ylabel('Probability Density','FontSize',14)
xlim([-3 3])
ylim([0 0.5])

begp = find(gridv < -3,1,'last');
endp = find(gridv > 3,1,'first');

X = [gridv(begp:endp)' fliplr(gridv(begp:endp)')];
Y = [u_dist(begp:endp) fliplr(l_dist(begp:endp))];
fill(X,Y,'b')

%% Overall Power Spectra
clear all 
close all
load C:\WC_Germany\Persistent_activity\power_spectra\power_spec_data

max_f = 1;

%convert to log power
Sw = 10*log10(Sw);
S8 = 10*log10(S8);

mean_pow = mean(Sw);
mean_pow8 = mean(S8);
u_pow = mean_pow+2*std(Sw)/sqrt(17);
l_pow = mean_pow-2*std(Sw)/sqrt(17);
u_pow8 = mean_pow8+2*std(S8)/sqrt(17);
l_pow8 = mean_pow8-2*std(S8)/sqrt(17);

plot(f,mean_pow,'linewidth',2)
hold on
plot(f,mean_pow8,'r','linewidth',2)
legend('MP','LFP')
plot(f,u_pow,'--')
plot(f,l_pow,'--')
plot(f,u_pow8,'r--')
plot(f,l_pow8,'r--')
xlim([0 max_f])
ylabel('Log Power (au)','FontSize',14)
xlabel('Frequency (Hz)','FontSize',14)
 
cd C:\WC_Germany\Persistent_activity\potential_figures

%% Overall Coherence
clear all
close all
load C:\WC_Germany\Persistent_activity\coherence\coherence_data
max_f = 1; %maximum frequency of interest in Hz
max_f_loc = find(f < max_f,1,'last');
mean_coh = mean(Cmn);
u_coh = mean_coh+1*std(Cmn)/sqrt(17);
l_coh = mean_coh-1*std(Cmn)/sqrt(17);

for i = 1:17
   sig_check(i,:) = Cerr(i,1,:) > ConfC(i);   
end
frac_sig = sum(sig_check)/17;

plot(f,mean_coh,'k','linewidth',2)
hold on
stairs(f,frac_sig,'r')
% plot(f,jmm_smooth_1d(frac_sig,1),'r')
% stairs(f,frac_sig,'m')
m_confc = mean(ConfC);
s_confc = std(ConfC);
u_confc = m_confc+1*s_confc/sqrt(17);
l_confc = m_confc-1*s_confc/sqrt(17);
X = [0 max_f max_f 0];
Y = [l_confc l_confc u_confc u_confc];
fill(X,Y,'b')
% legend('Coherence','Fraction of Cells with Significant Coherence','p=0.05 Sig Level')
% plot(f,u_coh,'k--')
% plot(f,l_coh,'k--')
X = [f(1:max_f_loc) fliplr(f(1:max_f_loc))];
Y = [u_coh(1:max_f_loc) fliplr(l_coh(1:max_f_loc))];
fill(X,Y,'k')
xlim([0 max_f])
xlabel('Frequency (Hz)','FontSize',14)

cd C:\WC_Germany\Persistent_activity\potential_figures


%% NEWEST (3-5-09) Up State Distribution Compare
clear all
close all
load C:\WC_Germany\Persistent_activity\lec_state_dur_data_new_method 

u_pyr_up = avg_pyr_up + se_pyr_up;
d_pyr_up = avg_pyr_up - se_pyr_up;

u_pyr_up8 = avg_pyr_up8+se_pyr_up8;
d_pyr_up8 = avg_pyr_up8-se_pyr_up8;

u_lec_up = avg_lec_up+se_lec_up;
d_lec_up = avg_lec_up-se_lec_up;

u_lec_up8 = avg_lec_up8+se_lec_up8;
d_lec_up8 = avg_lec_up8-se_lec_up8;

u_adj_up = avg_adj_up + se_adj_up;
d_adj_up = avg_adj_up - se_adj_up;

% upsamp_grid = linspace(log(0.3),log(maxdur),1000);
% upsamp_grid = exp(upsamp_grid);

% upsamp_grid = linspace(0.3,maxdur,1000);
% 
% us_avg_pyr_up = spline(up_grid,avg_pyr_up,upsamp_grid);
% us_avg_lec_up = spline(up_grid,avg_lec_up,upsamp_grid);
% us_avg_pyr8_up = spline(up_grid,avg_pyr_up8,upsamp_grid);
% us_avg_lec8_up = spline(up_grid,avg_lec_up8,upsamp_grid);
% 
% us_u_pyr_up = spline(up_grid,u_pyr_up,upsamp_grid);
% us_d_pyr_up = spline(up_grid,d_pyr_up,upsamp_grid);
% 
% us_u_lec_up = spline(up_grid,u_lec_up,upsamp_grid);
% us_d_lec_up = spline(up_grid,d_lec_up,upsamp_grid);
% 
% us_u_pyr8_up = spline(up_grid,u_pyr_up8,upsamp_grid);
% us_d_pyr8_up = spline(up_grid,d_pyr_up8,upsamp_grid);
% 
% us_u_lec8_up = spline(up_grid,u_lec_up8,upsamp_grid);
% us_d_lec8_up = spline(up_grid,d_lec_up8,upsamp_grid);

upsamp_grid = up_grid;

us_avg_pyr_up = avg_pyr_up;
us_avg_lec_up = avg_lec_up;
us_avg_pyr8_up = avg_pyr_up8;
us_avg_lec8_up = avg_lec_up8;
us_avg_adj_up = avg_adj_up;

us_u_pyr_up = u_pyr_up;
us_d_pyr_up = d_pyr_up;

us_u_lec_up = u_lec_up;
us_d_lec_up = d_lec_up;

us_u_pyr8_up = u_pyr_up8;
us_d_pyr8_up = d_pyr_up8;

us_u_lec8_up = u_lec_up8;
us_d_lec8_up = d_lec_up8;

us_u_adj_up = u_adj_up;
us_d_adj_up = d_adj_up;

epsilon = 1e-7;
us_avg_pyr_up(us_avg_pyr_up <= 0) = epsilon;
us_avg_lec_up(us_avg_lec_up <= 0) = epsilon;
us_avg_pyr8_up(us_avg_pyr8_up <= 0) = epsilon;
us_avg_lec8_up(us_avg_lec8_up <= 0) = epsilon;
us_avg_adj_up(us_avg_adj_up <= 0) = epsilon;
us_u_pyr_up(us_u_pyr_up <= 0) = epsilon;
us_u_lec_up(us_u_lec_up <= 0) = epsilon;
us_u_pyr8_up(us_u_pyr8_up <= 0) = epsilon;
us_u_lec8_up(us_u_lec8_up <= 0) = epsilon;
us_u_adj_up(us_u_adj_up <= 0) = epsilon;
us_d_pyr_up(us_d_pyr_up <= 0) = epsilon;
us_d_lec_up(us_d_lec_up <= 0) = epsilon;
us_d_pyr8_up(us_d_pyr8_up <= 0) = epsilon;
us_d_lec8_up(us_d_lec8_up <= 0) = epsilon;
us_d_adj_up(us_d_adj_up <= 0) = epsilon;

% 
% plot(upsamp_grid,us_avg_pyr_up,'linewidth',2)
% hold on
% plot(upsamp_grid,us_avg_pyr8_up,'r','linewidth',2)
% plot(upsamp_grid,us_avg_lec_up,'g','linewidth',2)
% % plot(upsamp_grid,us_avg_lec8_up,'k','linewidth',2)

% stairs(upsamp_grid,us_avg_pyr_up,'linewidth',2)
% hold on
% stairs(upsamp_grid,us_avg_pyr8_up,'r','linewidth',2)
% stairs(upsamp_grid,us_avg_lec_up,'g','linewidth',2)
% stairs(upsamp_grid,us_avg_adj_up,'k','linewidth',2)
% plot(upsamp_grid,us_avg_lec8_up,'k','linewidth',2)



legend('MEC','LFP','LEC')


% for i = 1:length(up_grid)-1
%     clear X Y
%     bin_width = (upsamp_grid(i+1)-upsamp_grid(i));
%     X = [upsamp_grid(i) upsamp_grid(i)+bin_width];
%     Y = [us_u_adj_up(i) us_u_adj_up(i)];
%     X = [X fliplr(X)];
%     Y = [Y us_d_adj_up(i) us_d_adj_up(i)];
%     fill(X,Y,'k','EdgeColor','none');
%     hold on
% end

for i = 1:length(up_grid)-1
    clear X Y
    bin_width = (upsamp_grid(i+1)-upsamp_grid(i));
    X = [upsamp_grid(i) upsamp_grid(i)+bin_width];
    Y = [us_u_lec_up(i) us_u_lec_up(i)];
    X = [X fliplr(X)];
    Y = [Y us_d_lec_up(i) us_d_lec_up(i)];
    fill(X,Y,'g','EdgeColor','none');
    hold on
end

for i = 1:length(up_grid)-1
    clear X Y
    bin_width = (upsamp_grid(i+1)-upsamp_grid(i));
    X = [upsamp_grid(i) upsamp_grid(i)+bin_width];
    Y = [us_u_pyr8_up(i) us_u_pyr8_up(i)];
    X = [X fliplr(X)];
    Y = [Y us_d_pyr8_up(i) us_d_pyr8_up(i)];
    fill(X,Y,'r','EdgeColor','none');
    hold on
end

for i = 1:length(up_grid)-1
    clear X Y
    bin_width = (upsamp_grid(i+1)-upsamp_grid(i));
    X = [upsamp_grid(i) upsamp_grid(i)+bin_width];
    Y = [us_u_pyr_up(i) us_u_pyr_up(i)];
    X = [X fliplr(X)];
    Y = [Y us_d_pyr_up(i) us_d_pyr_up(i)];
    fill(X,Y,'b','EdgeColor','none');
    hold on
end

set(gca,'yscale','log')
set(gca,'YMinorGrid','off')
set(gca,'YGrid','off')
% set(gca,'XMinorGrid','off')
ylim([1e-4 0.2])
xlim([0 10])
shg
xlabel('Up State Duration (s)','FontSize',14)
ylabel('Percent of Data','FontSize',14)
cd C:\WC_Germany\Persistent_activity\potential_figures

%% NEWEST (3-5-09) Down State Distribution Compare
clear all
close all
load C:\WC_Germany\Persistent_activity\lec_state_dur_data_new_method 

u_pyr_down = avg_pyr_down + se_pyr_down;
d_pyr_down = avg_pyr_down - se_pyr_down;

u_pyr_down8 = avg_pyr_down8+se_pyr_down8;
d_pyr_down8 = avg_pyr_down8-se_pyr_down8;

u_lec_down = avg_lec_down+se_lec_down;
d_lec_down = avg_lec_down-se_lec_down;

u_lec_down8 = avg_lec_down8+se_lec_down8;
d_lec_down8 = avg_lec_down8-se_lec_down8;

u_adj_down = avg_adj_down + se_adj_down;
d_adj_down = avg_adj_down - se_adj_down;

% upsamp_grid = linspace(log(0.3),log(maxdur),1000);
% upsamp_grid = exp(upsamp_grid);

% upsamp_grid = linspace(0.3,maxdur,1000);
% 
% us_avg_pyr_up = spline(up_grid,avg_pyr_up,upsamp_grid);
% us_avg_lec_up = spline(up_grid,avg_lec_up,upsamp_grid);
% us_avg_pyr8_up = spline(up_grid,avg_pyr_up8,upsamp_grid);
% us_avg_lec8_up = spline(up_grid,avg_lec_up8,upsamp_grid);
% 
% us_u_pyr_up = spline(up_grid,u_pyr_up,upsamp_grid);
% us_d_pyr_up = spline(up_grid,d_pyr_up,upsamp_grid);
% 
% us_u_lec_up = spline(up_grid,u_lec_up,upsamp_grid);
% us_d_lec_up = spline(up_grid,d_lec_up,upsamp_grid);
% 
% us_u_pyr8_up = spline(up_grid,u_pyr_up8,upsamp_grid);
% us_d_pyr8_up = spline(up_grid,d_pyr_up8,upsamp_grid);
% 
% us_u_lec8_up = spline(up_grid,u_lec_up8,upsamp_grid);
% us_d_lec8_up = spline(up_grid,d_lec_up8,upsamp_grid);

upsamp_grid = down_grid;

us_avg_pyr_down = avg_pyr_down;
us_avg_lec_down = avg_lec_down;
us_avg_pyr8_down = avg_pyr_down8;
us_avg_lec8_down = avg_lec_down8;
us_avg_adj_down = avg_adj_down;

us_u_pyr_down = u_pyr_down;
us_d_pyr_down = d_pyr_down;

us_u_lec_down = u_lec_down;
us_d_lec_down = d_lec_down;

us_u_pyr8_down = u_pyr_down8;
us_d_pyr8_down = d_pyr_down8;

us_u_lec8_down = u_lec_down8;
us_d_lec8_down = d_lec_down8;

us_u_adj_down = u_adj_down;
us_d_adj_down = d_adj_down;

epsilon = 1e-7;
us_avg_pyr_down(us_avg_pyr_down <= 0) = epsilon;
us_avg_lec_down(us_avg_lec_down <= 0) = epsilon;
us_avg_pyr8_down(us_avg_pyr8_down <= 0) = epsilon;
us_avg_lec8_down(us_avg_lec8_down <= 0) = epsilon;
us_avg_adj_down(us_avg_adj_down <= 0) = epsilon;
us_u_pyr_down(us_u_pyr_down <= 0) = epsilon;
us_u_lec_down(us_u_lec_down <= 0) = epsilon;
us_u_pyr8_down(us_u_pyr8_down <= 0) = epsilon;
us_u_lec8_down(us_u_lec8_down <= 0) = epsilon;
us_u_adj_down(us_u_adj_down <= 0) = epsilon;
us_d_pyr_down(us_d_pyr_down <= 0) = epsilon;
us_d_lec_down(us_d_lec_down <= 0) = epsilon;
us_d_pyr8_down(us_d_pyr8_down <= 0) = epsilon;
us_d_lec8_down(us_d_lec8_down <= 0) = epsilon;
us_d_adj_down(us_d_adj_down <= 0) = epsilon;

% 
% plot(upsamp_grid,us_avg_pyr_up,'linewidth',2)
% hold on
% plot(upsamp_grid,us_avg_pyr8_up,'r','linewidth',2)
% plot(upsamp_grid,us_avg_lec_up,'g','linewidth',2)
% % plot(upsamp_grid,us_avg_lec8_up,'k','linewidth',2)

% stairs(upsamp_grid,us_avg_pyr_down,'linewidth',2)
% hold on
% stairs(upsamp_grid,us_avg_pyr8_down,'r','linewidth',2)
% stairs(upsamp_grid,us_avg_lec_down,'g','linewidth',2)
% stairs(upsamp_grid,us_avg_adj_down,'k','linewidth',2)
% plot(upsamp_grid,us_avg_lec8_down,'k','linewidth',2)



legend('MEC','LFP','LEC')

% X = [upsamp_grid fliplr(upsamp_grid)];
% 
% % Y = [us_u_lec8_up fliplr(us_d_lec8_up)];
% % cur_fill = fill(X,Y,'k')
% 
% Y = [us_u_adj_up fliplr(us_d_adj_up)];
% cur_fill = fill(X,Y,'k')
% 
% Y = [us_u_lec_up fliplr(us_d_lec_up)];
% cur_fill = fill(X,Y,'g')
% 
% Y = [us_u_pyr8_up fliplr(us_d_pyr8_up)];
% cur_fill = fill(X,Y,'r')
% 
% Y = [us_u_pyr_up fliplr(us_d_pyr_up)];
% cur_fill = fill(X,Y,'b')


% for i = 1:length(down_grid)-1
%     clear X Y
%     bin_width = (upsamp_grid(i+1)-upsamp_grid(i));
%     X = [upsamp_grid(i) upsamp_grid(i)+bin_width];
%     Y = [us_u_adj_down(i) us_u_adj_down(i)];
%     X = [X fliplr(X)];
%     Y = [Y us_d_adj_down(i) us_d_adj_down(i)];
%     fill(X,Y,'k','EdgeColor','none');
%     hold on
% end

for i = 1:length(down_grid)-1
    clear X Y
    bin_width = (upsamp_grid(i+1)-upsamp_grid(i));
    X = [upsamp_grid(i) upsamp_grid(i)+bin_width];
    Y = [us_u_lec_down(i) us_u_lec_down(i)];
    X = [X fliplr(X)];
    Y = [Y us_d_lec_down(i) us_d_lec_down(i)];
    fill(X,Y,'g','EdgeColor','none');
    hold on
end

for i = 1:length(down_grid)-1
    clear X Y
    bin_width = (upsamp_grid(i+1)-upsamp_grid(i));
    X = [upsamp_grid(i) upsamp_grid(i)+bin_width];
    Y = [us_u_pyr8_down(i) us_u_pyr8_down(i)];
    X = [X fliplr(X)];
    Y = [Y us_d_pyr8_down(i) us_d_pyr8_down(i)];
    fill(X,Y,'r','EdgeColor','none');
    hold on
end

for i = 1:length(down_grid)-1
    clear X Y
    bin_width = (upsamp_grid(i+1)-upsamp_grid(i));
    X = [upsamp_grid(i) upsamp_grid(i)+bin_width];
    Y = [us_u_pyr_down(i) us_u_pyr_down(i)];
    X = [X fliplr(X)];
    Y = [Y us_d_pyr_down(i) us_d_pyr_down(i)];
    fill(X,Y,'b','EdgeColor','none');
    hold on
end

set(gca,'yscale','log')
set(gca,'YMinorGrid','off')
set(gca,'YGrid','off')
% set(gca,'XMinorGrid','off')
ylim([1e-4 0.2])
xlim([0 10])
shg
xlabel('Down State Duration (s)','FontSize',14)
ylabel('Percent of Data','FontSize',14)
cd C:\WC_Germany\Persistent_activity\potential_figures


%% MP Triggered UP Trans Overall
clear all
close all
load C:\WC_Germany\Persistent_activity\persistent_analysis\wcv_up_ctrig_data_new_method

[dummy,peaks_mp] = max(m_wcv_up_ctrig_mat,[],2);
[dummy,peaks_lfp] = max(m_lf8_up_ctrig_mat,[],2);
peak_lags = (peaks_mp-peaks_lfp)/2016*8;
half_pt = round(size(m_lf8_up_ctrig_mat,2)/2);
minlfp = min(m_lf8_up_ctrig_mat(:,1:half_pt),[],2);
maxlfp = max(m_lf8_up_ctrig_mat,[],2);
midpt = (minlfp+maxlfp)/2;
for i = 1:17
lfp_half_loc(i) = find(m_lf8_up_ctrig_mat(i,1:half_pt) < midpt(i),1,'last');
end
lfp_half_loc = lags(lfp_half_loc)/Fsd;
% subplot(2,1,1)
% plot(lags/Fsd,m_wcv_up_ctrig_mat)
% hold on
% plot(lags/Fsd,mean(m_wcv_up_ctrig_mat),'linewidth',3)
% xlabel('Time (s)','FontSize',14)
% ylabel('Amplitude (zscore)','FontSize',14)
% title('Average MP','FontSize',16)
% ylim([-2 3])
% hold on
% xlim([-2 2])
% line([0 0],[-2 3],'Color','k','linewidth',2)
% subplot(2,1,2)

m_lup = mean(m_lf8_up_ctrig_mat);
u_lup = m_lup+1*std(m_lf8_up_ctrig_mat)/sqrt(17);
l_lup = m_lup-1*std(m_lf8_up_ctrig_mat)/sqrt(17);

[min_avg_lfp,min_loc] = min(m_lup(1:half_pt));
[max_avg_lfp,max_loc] = max(m_lup);
midpt_avg = (min_avg_lfp+max_avg_lfp)/2;
half_loc_avg = find(m_lup(1:half_pt) < midpt_avg,1,'last');
half_loc_avg = lags(half_loc_avg)/Fsd;
u_half_avg = find(u_lup(1:half_pt) < midpt_avg,1,'last');
u_half_avg = lags(u_half_avg)/Fsd;
l_half_avg = find(l_lup(1:half_pt) < midpt_avg,1,'last');
l_half_avg = lags(l_half_avg)/Fsd;

min_u_lfp = min(u_lup(1:half_pt));
max_u_lfp = max(u_lup);
midpt_u = (min_u_lfp + max_u_lfp)/2;
half_loc_u = find(u_lup(1:half_pt) < midpt_u,1,'last');
half_loc_u = lags(half_loc_u)/Fsd;

min_l_lfp = min(l_lup(1:half_pt));
max_l_lfp = max(l_lup);
midpt_l = (min_l_lfp+max_l_lfp)/2;
half_loc_l = find(l_lup(1:half_pt) < midpt_l,1,'last');
half_loc_l = lags(half_loc_l)/Fsd;

m_wup = mean(m_wcv_up_ctrig_mat);
u_wup = m_wup+1*std(m_wcv_up_ctrig_mat)/sqrt(17);
l_wup = m_wup-1*std(m_wcv_up_ctrig_mat)/sqrt(17);
startp = find(lags/Fsd<-2,1,'last');
endp = find(lags/Fsd>2,1,'first');

plot(lags(startp:endp)/Fsd,m_lup(startp:endp),'r','linewidth',3)
hold on
plot(lags(min_loc)/Fsd,m_lup(min_loc),'ro')
plot(lags(max_loc)/Fsd,m_lup(max_loc),'ro')
plot(half_loc_avg,midpt_avg,'ro')
plot(lags(startp:endp)/Fsd,m_wup(startp:endp),'linewidth',3)
% plot(lags(startp:endp)/Fsd,u_lup(startp:endp),'r--')
% plot(lags(startp:endp)/Fsd,l_lup(startp:endp),'r--')
% plot(lags(startp:endp)/Fsd,u_wup(startp:endp),'--')
% plot(lags(startp:endp)/Fsd,l_wup(startp:endp),'--')
xlabel('Time (s)','FontSize',14)
ylabel('Amplitude (zscore)','FontSize',14)
title('Average LFP','FontSize',16)
xlim([-2 2])
ylim([-1.5 1.5])
line([0 0],[-1.5 2],'Color','k','linewidth',2)
% line([-0.2897 -0.2897],[-1.5 2],'Color','k')
% line([-2 2],[0 0],'Color','k','linewidth',2)
X = [lags(startp:endp)/Fsd fliplr(lags(startp:endp))/Fsd];
Y = [u_lup(startp:endp) fliplr(l_lup(startp:endp))];
fill(X,Y,'r')
Y = [u_wup(startp:endp) fliplr(l_wup(startp:endp))];
line([-2 2],[midpt_avg midpt_avg],'Color','k')
fill(X,Y,'b')
% plot(lags/Fsd,m_lf8_up_ctrig_mat)
% hold on
% plot(lags/Fsd,mean(m_lf8_up_ctrig_mat),'r','linewidth',4)
% xlabel('Time (s)','FontSize',14)
% ylabel('Amplitude (zscore)','FontSize',14)
% title('Average LFP','FontSize',16)
% xlim([-2 2])
% ylim([-1 2])
% line([0 0],[-1 2],'Color','k','linewidth',2)
% line([-0.2897 -0.2897],[-1 2],'Color','k')

cd C:\WC_Germany\Persistent_activity\potential_figures

%% MP Triggered Down Trans Overall
clear all
close all
load C:\WC_Germany\Persistent_activity\persistent_analysis\wcv_up_ctrig_data_new_method
half_pt = round(size(m_lf8_down_ctrig_mat,2)/2);
minlfp = min(m_lf8_down_ctrig_mat,[],2);
maxlfp = max(m_lf8_down_ctrig_mat(:,1:half_pt),[],2);
midpt = (minlfp+maxlfp)/2;
for i = 1:17
lfp_half_loc(i) = find(m_lf8_down_ctrig_mat(i,1:half_pt) > midpt(i),1,'last');
end
lfp_half_loc = lags(lfp_half_loc)/Fsd;

% subplot(2,1,1)
% plot(lags/Fsd,m_wcv_down_ctrig_mat)
% hold on
% plot(lags/Fsd,mean(m_wcv_down_ctrig_mat),'linewidth',3)
% xlabel('Time (s)','FontSize',14)
% ylabel('Amplitude (zscore)','FontSize',14)
% title('Average MP','FontSize',16)
% xlim([-2 2])
% ylim([-2.5 2.5])
% line([0 0],[-2.5 2.5],'Color','k','linewidth',2)
% subplot(2,1,2)
% plot(lags/Fsd,m_lf8_down_ctrig_mat)
% hold on
% plot(lags/Fsd,mean(m_lf8_down_ctrig_mat),'r','linewidth',3)
% xlabel('Time (s)','FontSize',14)
% ylabel('Amplitude (zscore)','FontSize',14)
% title('Average LFP','FontSize',16)
% line([0 0],[-1.5 1.5],'Color','k','linewidth',2)
% xlim([-2 2])
% ylim([-1.5 1.5])

m_ldown = mean(m_lf8_down_ctrig_mat);
u_ldown = m_ldown+1*std(m_lf8_down_ctrig_mat)/sqrt(17);
l_ldown = m_ldown-1*std(m_lf8_down_ctrig_mat)/sqrt(17);

[min_avg_lfp,min_loc] = min(m_ldown(1:half_pt));
[max_avg_lfp,max_loc] = max(m_ldown(1:half_pt));
midpt_avg = (min_avg_lfp+max_avg_lfp)/2;
half_loc_avg = find(m_ldown(1:half_pt) > midpt_avg,1,'last');
half_loc_avg = lags(half_loc_avg)/Fsd;
u_half_loc_avg = find(u_ldown(1:half_pt) > midpt_avg,1,'last');
u_half_loc_avg = lags(u_half_loc_avg)/Fsd;
l_half_loc_avg = find(l_ldown(1:half_pt) > midpt_avg,1,'last');
l_half_loc_avg = lags(l_half_loc_avg)/Fsd;

min_u_lfp = min(u_ldown(1:half_pt));
max_u_lfp = max(u_ldown);
midpt_u = (min_u_lfp + max_u_lfp)/2;
half_loc_u = find(u_ldown(1:half_pt) > midpt_u,1,'last');
half_loc_u = lags(half_loc_u)/Fsd;

min_l_lfp = min(l_ldown(1:half_pt));
max_l_lfp = max(l_ldown);
midpt_l = (min_l_lfp+max_l_lfp)/2;
half_loc_l = find(l_ldown(1:half_pt) > midpt_l,1,'last');
half_loc_l = lags(half_loc_l)/Fsd;


m_wdown = mean(m_wcv_down_ctrig_mat);
u_wdown = m_wdown+1*std(m_wcv_down_ctrig_mat)/sqrt(17);
l_wdown = m_wdown-1*std(m_wcv_down_ctrig_mat)/sqrt(17);
startp = find(lags/Fsd<-2,1,'last');
endp = find(lags/Fsd>2,1,'first');

plot(lags(startp:endp)/Fsd,m_ldown(startp:endp),'r','linewidth',3)
hold on
plot(lags(min_loc)/Fsd,m_ldown(min_loc),'ro')
plot(lags(max_loc)/Fsd,m_ldown(max_loc),'ro')
plot(lags(half_loc_avg)/Fsd,m_ldown(half_loc_avg),'ro')

plot(lags(startp:endp)/Fsd,m_wdown(startp:endp),'linewidth',3)
% plot(lags(startp:endp)/Fsd,u_ldown(startp:endp),'r--')
% plot(lags(startp:endp)/Fsd,l_ldown(startp:endp),'r--')
% plot(lags(startp:endp)/Fsd,u_wdown(startp:endp),'--')
% plot(lags(startp:endp)/Fsd,l_wdown(startp:endp),'--')
xlabel('Time (s)','FontSize',14)
ylabel('Amplitude (zscore)','FontSize',14)
title('Average LFP','FontSize',16)
xlim([-2 2])
ylim([-1.5 1.5])
line([0 0],[-1.5 2],'Color','k','linewidth',2)
thresh_cross = startp+find(m_ldown(startp:end) < -0.1400,1,'first');
line([lags(thresh_cross)/Fsd lags(thresh_cross)/Fsd],[-1.5 2],'Color','k','linestyle','--')

X = [lags(startp:endp)/Fsd fliplr(lags(startp:endp))/Fsd];
Y = [u_ldown(startp:endp) fliplr(l_ldown(startp:endp))];
fill(X,Y,'r')
Y = [u_wdown(startp:endp) fliplr(l_wdown(startp:endp))];
fill(X,Y,'b')


cd C:\WC_Germany\Persistent_activity\potential_figures

%% Quantization Figure Overall
clear all
% close all
load C:\WC_Germany\Persistent_activity\persistent_analysis\wcv_up_ctrig_data_new_method

% for i = 1:17
%     cor_mp_lperiod_dur(i,:) = cor_mp_lperiod_dur(i,:)/sum(cor_mp_lperiod_dur(i,:));
%     cor_mp_lperiod_dur(i,:) = jmm_smooth_1d(cor_mp_lperiod_dur(i,:),2)/mean(diff(hist_range));
% end

hist_range = linspace(0,7,100);
binsize = hist_range(2)-hist_range(1);
for i = 1:17
    
    cell_hist(i,:) = hist(lfp_period_dur{i},hist_range);
    cell_hist(i,:) = cell_hist(i,:)/sum(cell_hist(i,:))/binsize;
    cell_hist(i,:) = jmm_smooth_1d_cor(cell_hist(i,:),2);
    
end

% sm_hist = jmm_smooth_1d(mean(cell_hist),2);
sm_hist = mean(cell_hist);
u = mean(cell_hist)+1*std(cell_hist)/sqrt(17);
l = mean(cell_hist)-1*std(cell_hist)/sqrt(17);


epsilon = 1e-10;
upsamp_grid = linspace(0,7,1000);
sm_hist = spline(hist_range,sm_hist,upsamp_grid);
u = spline(hist_range,u,upsamp_grid);
l = spline(hist_range,l,upsamp_grid);
l(l < epsilon) = epsilon;
u(u < epsilon) = epsilon;
sm_hist(sm_hist < epsilon) = epsilon;
plot(upsamp_grid,sm_hist)
hold on
X = [upsamp_grid fliplr(upsamp_grid)];
Y = [u fliplr(l)];
fill(X,Y,'b')
% plot(hist_range,u,'r')
% plot(hist_range,l,'r')

set(gca,'yscale','log')
xlim([0 6.2])
ylim([5e-4 2])

[a,b] = findpeaks(sm_hist);
peak_locs = upsamp_grid(b);

for i = 1:length(peak_locs)
    line([peak_locs(i) peak_locs(i)],[5e-4 2],'Color','k')
end

cd C:\WC_Germany\Persistent_activity\potential_figures
shg
%% Quantization Figure Single Example
clear all
close all
load C:\WC_Germany\Persistent_activity\persistent_analysis\wcv_up_ctrig_data_new_method

examp_ind = 14;
x = linspace(0,7,100);
binsize = x(2)-x(1);

y = hist(cor_lfp_period_dur{examp_ind},x);
% [yg,xg] = gpkde(cor_lfp_period_dur{14}',[0 6 numBins],-3);
% y = y/sum(y)/binsize;

bar(x,y,1,'k')
hold on
temp = jmm_smooth_1d_cor(y,2);
fine_axis = linspace(0,7,1000);
smooth_temp = spline(x,temp,fine_axis);
% plot(x,jmm_smooth_1d_cor(y,2),'r','linewidth',2)
plot(fine_axis,smooth_temp,'r','linewidth',2)
xlim([0 4.2])
% ylim([0 1.7])
shg

%% Quantization Figure Second Example
clear all
close all
load C:\WC_Germany\Persistent_activity\persistent_analysis\wcv_up_ctrig_data_new_method

%examp 1 = 14, examp2 = 8, examp3 = 1, examp4 = 5, examp5 = 15, examp6 = 2
examp_ind = 10;
x = linspace(0,7,100);
binsize = x(2)-x(1);

y = hist(lfp_period_dur{examp_ind},x);
% [yg,xg] = gpkde(cor_lfp_period_dur{14}',[0 6 numBins],-3);
% y = y/sum(y)/binsize;

bar(x,y,1,'k')
hold on
temp = jmm_smooth_1d_cor(y,2);
fine_axis = linspace(0,7,1000);
smooth_temp = spline(x,temp,fine_axis);
% plot(x,jmm_smooth_1d_cor(y,2),'r','linewidth',2)
plot(fine_axis,smooth_temp,'r','linewidth',2)
xlim([0 4.2])
% ylim([0 1.7])
shg

%% Histogram of Persistent Activity
clear all
close all
load C:\WC_Germany\Persistent_activity\persistent_analysis\wcv_up_ctrig_data_new_method pers_fract*

numBins = 10;
% pers_g = linspace(0,0.6,numBins)
[pers_h,pers_g] = hist(pers_fract_noshift,numBins);
xlim([0 60])
ylim([0 4])
xlabel('Fraction of Up States Persisting','FontSize',14)
ylabel('Number of Cells','FontSize',14)

load C:\WC_Germany\overall_calcs\persistence_calcs\wcv_up_ctrig_data_new_method pers_fract 
load C:\WC_Germany\overall_calcs\overall_dir 

bad_cells = find(cell_type == 8 | cell_type == 8);
pers_fract(bad_cells) = [];
cell_type(bad_cells) = [];
cell_type(cell_type == 6) = 4;
cell_type(cell_type == 5) = 4;

lec_cells = find(cell_type==4);
cell_type(lec_cells([2 3 5 6 7 8 9])) = [];
pers_fract(lec_cells([2 3 5 6 7 8 9])) = [];
unused_cells = find(cell_type ~= 4);
pers_fract(unused_cells) = [];
pers_lec = hist(pers_fract,pers_g);
% stairs(100*pers_g,pers_h,'linewidth',2)
% hold on
% stairs(100*pers_g,pers_lec,'k','linewidth',2)
% bar(100*pers_g,[pers_h' pers_lec'])
figure
bar(100*pers_g,pers_h,'k')
xlim([0 60])
ylim([0 4])
xlabel('Fraction of Up States Persisting','FontSize',14)
ylabel('Number of Cells','FontSize',14)


%% Depth from pia vs up state duration or freq
clear all
close all
load C:\WC_Germany\Persistent_activity\mean_state_dur_data
load C:\WC_Germany\overall_data\meta_data_use_2
load C:\WC_Germany\Persistent_activity\power_spectra\com_freq_data

depths = depth_from_pia(17:end);
figure
scatter_with_cor(depths,log(mean_wcv_up))
figure
scatter_with_cor(depths,com_freq)

%% WCV Up Trig Mat example
clear all
close all

%%for new method
load C:\WC_Germany\Persistent_activity\persistent_analysis\sing_examp_new_method
load C:\WC_Germany\Persistent_activity\UDS_dur_raw_no_cut\UDS_raw_data_no_cut
load C:\WC_Germany\Persistent_activity\UDS_synch_state_dur\UDS_synch_state_dur_data_new_method

d = 14;
dsf = 8;
Fsd = 2016/dsf;

cur_mp_mat = wcv_up_ctrig_mat;
cur_lfp_mat = lf8_up_ctrig_mat;

cur_mp_mat(bad_rows,:) = [];
cur_lfp_mat(bad_rows,:) = [];
synch_up_dur{d}(bad_rows) = [];
[dummy,up_order] = sort(synch_up_dur{d});

%     Fig = figure(1)
%     clf
%     set(Fig,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
%     imagesc(lags/Fsd,(1:length(synch_ups{d})) ...
%         ,(cur_mp_mat(up_order,:)));shading flat;
%     caxis([-3 3]);colorbar
%         xlim([-2 10])
% 
%     hold on
%     xlim([-2 10])
%     xlabel('Time (s)','FontSize',14)
%     ylabel('MP Up Duration Percentile','FontSize',14)
%     Fig2 = figure(2)
%     clf
%     set(Fig2,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(Fig2,'PaperPosition',[0,0,(get(Fig2,'PaperSize'))])
%     imagesc(lags/Fsd,(1:length(synch_ups{d}))...
%         ,(cur_lfp_mat(up_order,:)));shading flat;
%     caxis([-2 4]);colorbar
%         xlim([-2 10])
% 
% figure
%     plot(synch_up_dur{d}(up_order),...
%         1-(1:length(synch_up_dur{d}))/length(synch_up_dur{d}),'r','linewidth',2)
%     line([0 0],[0 1],'Color','k')
%     xlim([-2 10])
%     xlabel('Time (s)','FontSize',14)
%     ylabel('MP Up Duration Percentile','FontSize',14)
% 
    
load C:\WC_Germany\Persistent_activity\dir_tree_update

%specify session number
d = 14;
niqf = 2016/2;
lcf = 0.1/niqf;
hcf = 6/niqf;
[b,a] = butter(2,[lcf hcf]);

dsf = 8;
Fsd = 2016/dsf;

cd(dir_array{d})
pwd

load used_data lf8 wcv_minus_spike
load spike_time_jmm
spike_ids = round(spkid/dsf);

lf8_f = filtfilt(b,a,lf8);
wcv_f = filtfilt(b,a,wcv_minus_spike);

lf8_d = downsample(lf8_f,dsf);
wcv_d = downsample(wcv_f,dsf);

lf8_d = zscore(lf8_d);
wcv_d = zscore(wcv_d);
t = (1:length(lf8_d))/Fsd;


    %at 17%
    trans_num = 24;
    prev_lfp_tran = find(up_trans8{d}<up_trans{d}(trans_num),1,'last');
    cur_lag = up_trans{d}(trans_num)-up_trans8{d}(prev_lfp_tran);
    st = up_trans{d}(trans_num)/Fsd-4;
    et = st+9;
    figure
    plot(t,wcv_d,'linewidth',3)
    hold on
    
    %if you want to lag lfp
%     plot(t+cur_lag/Fsd,lf8_d-1,'r','linewidth',3)
%else
    plot(t,lf8_d-1,'r','linewidth',3)
    hold on
    xlim([100.8 103.5])
    ylim([-3.5 3])
    used_spikes = find(spike_ids/Fsd > 100.8 & spike_ids/Fsd < 103.5);
line_bottom = 2.2;
line_top = 2.6;
for i = 1:length(used_spikes)
   line([t(spike_ids(used_spikes(i))) t(spike_ids(used_spikes(i)))],...
       [line_bottom line_top],'Color','k','linewidth',2) 
end
    
    shg
    
    %at 60%
    trans_num = 94;
    prev_lfp_tran = find(up_trans8{d}<up_trans{d}(trans_num),1,'last');
    cur_lag = up_trans{d}(trans_num)-up_trans8{d}(prev_lfp_tran);
    st = up_trans{d}(trans_num)/Fsd-3;
    et = st+10;
    figure
    plot(t,wcv_d,'linewidth',3)
    hold on
%     plot(t+cur_lag/Fsd,lf8_d-1,'r','linewidth',3)
    plot(t,lf8_d-1,'r','linewidth',3)
    xlim([495 499.3])
    ylim([-3.5 3])
        used_spikes = find(spike_ids/Fsd > 495 & spike_ids/Fsd < 499.3);

    line_bottom = 2.2;
line_top = 2.6;
for i = 1:length(used_spikes)
   line([t(spike_ids(used_spikes(i))) t(spike_ids(used_spikes(i)))],...
       [line_bottom line_top],'Color','k','linewidth',2) 
end

    shg
    
       %at 80%
       trans_num = 31;
      prev_lfp_tran = find(up_trans8{d}<up_trans{d}(trans_num),1,'last');
    cur_lag = up_trans{d}(trans_num)-up_trans8{d}(prev_lfp_tran);     
    st = up_trans{d}(trans_num)/Fsd-6;
    et = st+15;
    figure
    plot(t,wcv_d,'linewidth',3)
    hold on
%     plot(t+cur_lag/Fsd,lf8_d-1,'r','linewidth',3)
    plot(t,lf8_d-1,'r','linewidth',3)
    xlim([125.3 131.1])
    ylim([-3.5 3])
        used_spikes = find(spike_ids/Fsd > 125.3 & spike_ids/Fsd < 131.1);
line_bottom = 2.2;
line_top = 2.6;
for i = 1:length(used_spikes)
   line([t(spike_ids(used_spikes(i))) t(spike_ids(used_spikes(i)))],...
       [line_bottom line_top],'Color','k','linewidth',2) 
end

    shg

       %at 95%
       trans_num = 116;
      prev_lfp_tran = find(up_trans8{d}<up_trans{d}(trans_num),1,'last');
    cur_lag = up_trans{d}(trans_num)-up_trans8{d}(prev_lfp_tran);     
    st = up_trans{d}(trans_num)/Fsd-6;
    et = st+15;
    figure
    plot(t,wcv_d,'linewidth',3)
    hold on
%     plot(t+cur_lag/Fsd,lf8_d-1,'r','linewidth',3)
plot(t,lf8_d-1,'r','linewidth',3)
%     xlim([609 616.4])
    xlim([609.1 616.4])

    ylim([-3.5 3])
        used_spikes = find(spike_ids/Fsd > 609.1 & spike_ids/Fsd < 616.4);
line_bottom = 2.2;
line_top = 2.6;
for i = 1:length(used_spikes)
   line([t(spike_ids(used_spikes(i))) t(spike_ids(used_spikes(i)))],...
       [line_bottom line_top],'Color','k','linewidth',2) 
end

    shg
   
    cd C:\WC_Germany\Persistent_activity\potential_figures

%% persistent activity examples
clear all
load C:\WC_Germany\Persistent_activity\dir_tree_update
load C:\WC_Germany\Persistent_activity\UDS_dist\synch_UDS_inds
load C:\WC_Germany\Persistent_activity\UDS_dur_raw\UDS_raw_data

%specify session number
d = 14;
niqf = 2016/2;
lcf = 0.05/niqf;
hcf = 4/niqf;
[b,a] = butter(2,[lcf hcf]);

dsf = 1;
Fsd = 2016/dsf;

cd(dir_array{d})
pwd

load used_data lf8 wcv_minus_spike wcv

lf8_f = filtfilt(b,a,lf8);
% wcv_f = filtfilt(b,a,wcv_minus_spike);
wcv_f = wcv;

lf8_d = downsample(lf8_f,dsf);
wcv_d = downsample(wcv_f,dsf);

lf8_d = zscore(lf8_d);
wcv_d = zscore(wcv_d);

t = (1:length(lf8_d))/Fsd;

plot(t,wcv_d,'linewidth',3)
xlabel('Time (s)','FontSize',14)
ylabel('Amplitude (zscore)','FontSize',14)
title('MEC Membrane Potential','FontSize',16)
hold on
plot(t,lf8_d,'r','linewidth',3)
xlabel('Time (s)','FontSize',14)
ylabel('Amplitude (zscore)','FontSize',14)
title('Cortical LFP','FontSize',16)

% %example of no PA
% x1s = 314;
% x1e = 324;
% xlim([x1s x1e])
% shg
% 
%example1 PA
x2s = 170;
x2e = 180;
xlim([x2s x2e])
shg

% %example2 PA
x3s = 740;
x3e = 749;
xlim([x3s x3e])
% shg


%% persistent activity examples
clear all
load C:\WC_Germany\Persistent_activity\dir_tree_update

%specify session number
d =15;
niqf = 2016/2;
lcf = 0.05/niqf;
hcf = 4/niqf;
[b,a] = butter(2,[lcf hcf]);

dsf = 1;
Fsd = 2016/dsf;

cd(dir_array{d})
pwd

load used_data lf8 wcv_minus_spike wcv

lf8_f = filtfilt(b,a,lf8);
% wcv_f = filtfilt(b,a,wcv_minus_spike);
wcv_f = wcv;

lf8_d = downsample(lf8_f,dsf);
wcv_d = downsample(wcv_f,dsf);

lf8_d = zscore(lf8_d);
wcv_d = zscore(wcv_d);

t = (1:length(lf8_d))/Fsd;

plot(t,wcv_d,'linewidth',3)
xlabel('Time (s)','FontSize',14)
ylabel('Amplitude (zscore)','FontSize',14)
title('MEC Membrane Potential','FontSize',16)
hold on
plot(t,lf8_d,'r','linewidth',3)
xlabel('Time (s)','FontSize',14)
ylabel('Amplitude (zscore)','FontSize',14)
title('Cortical LFP','FontSize',16)

% %example1 noPA
x2s = 14;
x2e = 22;
xlim([x2s x2e])
shg

%example1 PA
x2s = 726;
x2e = 743;
xlim([x2s x2e])
shg

%% SUPER persistent activity examples
clear all, close all
load C:\WC_Germany\Persistent_activity\dir_tree_update

%specify session number
d =5;
niqf = 2016/2;
lcf = 0.05/niqf;
hcf = 4/niqf;
[b,a] = butter(2,hcf,'low');

dsf = 1;
Fsd = 2016/dsf;

cd(dir_array{d})
pwd

load used_data lf8 wcv_minus_spike wcv

lf8_f = filtfilt(b,a,lf8);
% wcv_f = filtfilt(b,a,wcv_minus_spike);
wcv_f = wcv;

lf8_d = downsample(lf8_f,dsf);
wcv_d = downsample(wcv_f,dsf);
% wcv_d = downsample(wcv,dsf);

lf8_d = zscore(lf8_d);
wcv_d = zscore(wcv_d);

t = (1:length(lf8_d))/Fsd;


% %example1 noPA
% x2s = 448;
% x2e = 465;
% xlim([x2s x2e])
% shg

%example1 PA
x2s = 220;
x2e = 275;


% subplot(2,1,1)
plot(t-x2s,wcv_d,'linewidth',1)
xlabel('Time (s)','FontSize',14)
ylabel('Amplitude (zscore)','FontSize',14)
title('MEC Membrane Potential','FontSize',16)
% subplot(2,1,2)
figure
% hold on
plot(t-x2s,lf8_d,'r','linewidth',1)
xlabel('Time (s)','FontSize',14)
ylabel('Amplitude (zscore)','FontSize',14)
title('Cortical LFP','FontSize',16)


% subplot(2,1,1)
figure(1)
xlim([x2s-x2s x2e-x2s])
figure(2)
xlim([x2s-x2s x2e-x2s])

% ylim([-3 2])
% subplot(2,1,2)

shg
% 
% %example2 PA
% x2s = 225;
% x2e = 273;
% subplot(2,1,1)
% xlim([x2s-225 x2e-225])
% subplot(2,1,2)
% xlim([x2s-225 x2e-225])
% shg
%% SUPER persistent activity examples
clear all
load C:\WC_Germany\Persistent_activity\dir_tree_update

%specify session number
d =15;
niqf = 2016/2;
lcf = 0.05/niqf;
hcf = 4/niqf;
[b,a] = butter(2,hcf,'low');

dsf = 8;
Fsd = 2016/dsf;

cd(dir_array{d})
pwd

load used_data lf8 wcv_minus_spike wcv

lf8_f = filtfilt(b,a,lf8);
wcv_f = filtfilt(b,a,wcv_minus_spike);

lf8_d = downsample(lf8_f,dsf);
wcv_d = downsample(wcv_f,dsf);
% wcv_d = downsample(wcv,dsf);

lf8_d = zscore(lf8_d);
wcv_d = zscore(wcv_d);

t = (1:length(lf8_d))/Fsd;


% %example1 noPA
% x2s = 11;
% x2e = 19.5;
% xlim([x2s x2e])
% shg

%example1 PA
x2s = 724.5;
x2e = 746;

subplot(2,1,1)
plot(t-x2s,wcv_d,'linewidth',1)
xlabel('Time (s)','FontSize',14)
ylabel('Amplitude (zscore)','FontSize',14)
title('MEC Membrane Potential','FontSize',16)
subplot(2,1,2)
hold on
plot(t-x2s,lf8_d,'r','linewidth',1)
xlabel('Time (s)','FontSize',14)
ylabel('Amplitude (zscore)','FontSize',14)
title('Cortical LFP','FontSize',16)


subplot(2,1,1)
xlim([x2s-x2s x2e-x2s])
ylim([-2 1.5])
subplot(2,1,2)
ylim([-2 2])
xlim([x2s-x2s x2e-x2s])
shg
% 
% %example2 PA
% x2s = 225;
% x2e = 273;
% subplot(2,1,1)
% xlim([x2s-225 x2e-225])
% subplot(2,1,2)
% xlim([x2s-225 x2e-225])
% shg

%% SUPER persistent activity examples
clear all
load C:\WC_Germany\Persistent_activity\dir_tree_update

%specify session number
d =14;
niqf = 2016/2;
lcf = 0.05/niqf;
hcf = 4/niqf;
[b,a] = butter(2,hcf,'low');

dsf = 1;
Fsd = 2016/dsf;

cd(dir_array{d})
pwd

load used_data lf8 wcv_minus_spike wcv
load spike_time_jmm

spike_ids = round(spkid/dsf);

lf8_f = filtfilt(b,a,lf8);
% wcv_f = filtfilt(b,a,wcv_minus_spike);
wcv_f = wcv;

lf8_d = downsample(lf8_f,dsf);
wcv_d = downsample(wcv_f,dsf);
% wcv_d = downsample(wcv,dsf);

lf8_d = zscore(lf8_d);
wcv_d = zscore(wcv_d);

t = (1:length(lf8_d))/Fsd;


% %example1 noPA
x2s = 94.5+760;
x2e = 94.5+768.5;
% xlim([x2s x2e])
% shg

%example1 PA
x2s = 424.5;
x2e = 449.8;
% 

subplot(2,1,1)
plot(t-x2s,wcv_d,'linewidth',1)
xlabel('Time (s)','FontSize',14)
ylabel('Amplitude (zscore)','FontSize',14)
title('MEC Membrane Potential','FontSize',16)
used_spikes = find(spike_ids/Fsd > x2s & spike_ids/Fsd < x2e);
line_bottom = 1.2;
line_top = 1.6;
for i = 1:length(used_spikes)
   line([t(spike_ids(used_spikes(i))) - x2s t(spike_ids(used_spikes(i))) - x2s],...
       [line_bottom line_top],'Color','k','linewidth',2) 
end


subplot(2,1,2)
hold on
plot(t-x2s,lf8_d,'r','linewidth',1)
xlabel('Time (s)','FontSize',14)
ylabel('Amplitude (zscore)','FontSize',14)
title('Cortical LFP','FontSize',16)


subplot(2,1,1)
xlim([x2s-x2s x2e-x2s])
% ylim([-2 2])
subplot(2,1,2)
ylim([-1.5 2.5])
xlim([x2s-x2s x2e-x2s])
shg

%% Supplementary figure all MP autocorr
clear all
close all
load C:\WC_Germany\Persistent_activity\corr_data
% load corr_data
m_corr = mean(w_acorr);
u_corr = m_corr+1*std(w_acorr)/sqrt(17);
l_corr = m_corr-1*std(w_acorr)/sqrt(17);

m_lcorr = mean(l_acorr);
u_lcorr = m_lcorr+1*std(l_acorr)/sqrt(17);
l_lcorr = m_lcorr-1*std(l_acorr)/sqrt(17);

% plot(lags/Fsd,w_acorr)
% hold on
% plot(lags/Fsd,mean(w_acorr),'linewidth',4)
plot(lags/Fsd,m_corr,'linewidth',2)
hold on
plot(lags/Fsd,m_lcorr,'r','linewidth',2)
X = [lags/Fsd fliplr(lags)/Fsd];
Y = [l_corr fliplr(u_corr)];
fill(X,Y,'b')
Y = [l_lcorr fliplr(u_lcorr)];
fill(X,Y,'r')
line([0 0],[-0.4 1],'Color','k')
line([-8 8],[0 0],'Color','k')
xlim([-8 8])
ylim([-0.4 1])
ylabel('Correlation Coefficient')
xlabel('Time (s)')
%% Supplementary figure for MP/LFP xcorr
clear all
close all
load C:\WC_Germany\Persistent_activity\corr_data
[peak_cor,peak_loc] = max(x_cor(:,2521:end),[],2);
[ov_peak_cor,ov_peak_loc] = max(x_cor,[],2);
ov_peak_loc = lags(ov_peak_loc)/Fsd;
peak_loc = lags(peak_loc+2520)/Fsd;

mean_xcorr = mean(x_cor);
u_xcorr = mean_xcorr+std(x_cor)/sqrt(17);
l_xcorr = mean_xcorr-std(x_cor)/sqrt(17);
plot(lags/Fsd,mean_xcorr,'linewidth',2)
hold on
X = [lags/Fsd fliplr(lags)/Fsd];
Y = [u_xcorr fliplr(l_xcorr)];
fill(X,Y,'b')
% plot(lags/Fsd,x_cor)
% hold on
% plot(lags/Fsd,mean(x_cor),'linewidth',4)
ylabel('Correlation Coefficient')
xlabel('Time (s)')
ylim([-0.4 0.4])
xlim([-4 4])
line([0 0],[-0.5 0.6],'Color','k')
%% Supplementary figure for LFP trig avg
clear all
close all
load C:\WC_Germany\Persistent_activity\lfp_trig_avg_data
plot(lags/Fsd,mean(mean_lfp_trig_lfp),'r','linewidth',4)
hold on
plot(lags/Fsd,mean(mean_lfp_trig_mp),'linewidth',4)
plot(lags/Fsd,mean_lfp_trig_mp)
ylabel('Amplitude (z)')
xlabel('Time (s)')
ylim([-1 1.7])
xlim([-4 4])
line([0 0],[-1 1.7],'Color','k')

%% heka amplitude supp
clear all
close all
load('C:\WC_Germany\Persistent_activity\heka_amp\heka_amp_data.mat')

for i = 1:17
    cell_amp_dist(i,:) = y{i};
end
grid = x{1}';

mean_amp = mean(cell_amp_dist);
u_amp = mean_amp+std(cell_amp_dist)/sqrt(17);
l_amp = mean_amp-std(cell_amp_dist)/sqrt(17);

% plot(grid,mean_amp,'linewidth',2)
% hold on
X = [grid fliplr(grid)];
Y = [u_amp fliplr(l_amp)];
% fill(X,Y,'b')

grid2 = grid-7;
X = [grid2 fliplr(grid2)];
plot(grid2,mean_amp,'r','linewidth',2)
hold on
fill(X,Y,'r')
xlim([-100 -20])

