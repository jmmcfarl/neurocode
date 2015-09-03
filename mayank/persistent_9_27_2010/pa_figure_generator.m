%% Example MEC MP trace
clear all
close all
load G:\WC_Germany\overall_EC\overall_EC_dir
used_data = [l3mec_p l3lec_p];
sess_data = sess_data(used_data);
mec_cells = 1:length(l3mec_p);
lec_cells = (mec_cells(end)+1):(mec_cells(end)+length(l3lec_p));

%specify session number
d = 14;
niqf = 2016/2;
lcf = 0.05/niqf;
hcf = 4/niqf;
[b,a] = butter(2,[lcf hcf]);

dsf = 1;
Fsd = 2016/dsf;

cd(sess_data(d).directory)
pwd

load ./used_data lf8 lf3 wcv wcv_minus_spike

lf8_f = filtfilt(b,a,lf8);
lf3_f = filtfilt(b,a,lf3);
% wcv_f = filtfilt(b,a,wcv_minus_spike);

lf8_d = downsample(lf8_f,dsf);
lf3_d = downsample(lf3_f,dsf);
wcv_d = downsample(wcv,dsf);

lf8_d = zscore(lf8_d);
lf3_d = zscore(lf3_d);
wcv_d = zscore(wcv_d);
wcv_std = std(wcv)*1e3; %in mv
twenty = 20/wcv_std;

% wcv = zscore(wcv);
% lf8 = zscore(lf8);
t = (1:length(wcv_d))/2016*dsf;

plot(t-66,wcv_d,'linewidth',2)
xlabel('Time (s)','FontSize',14)
ylabel('Amplitude (zscore)','FontSize',14)
title('MEC Membrane Potential','FontSize',16)
xlim([66-66 82-66])
ylim([-3 8])
line([5 5],[0 0+twenty],'color','k')

figure
plot(t-66,lf8_d,'r','linewidth',3)
hold on
plot(t-66,lf3_d,'k','linewidth',3)
xlabel('Time (s)','FontSize',14)
ylabel('Amplitude (zscore)','FontSize',14)
title('Cortical LFP','FontSize',16)
legend('MP','LFP')
xlim([66-66 82-66])
ylim([-2 3])

% used_spikes = find(spike_ids/Fsd > 66 & spike_ids/Fsd < 82);
% % plot(t(spike_ids)-66,ones(size(spike_ids))*3,'ro')
% line_bottom = 2.8;
% line_top = 3.2;
% for i = 1:length(used_spikes)
%    line([t(spike_ids(used_spikes(i)))-66 t(spike_ids(used_spikes(i)))-66],...
%        [line_bottom line_top],'Color','k','linewidth',2)
% end
shg

% load ./aligned_heka.mat
% figure
% plot(dc_time,dc_data*10)
% xlim([66 82])

figure
plot(t,wcv_d), hold on
plot(t,lf8_d,'r','linewidth',2)
plot(t,lf3_d,'k','linewidth',2)
xlim([299 325])
xlim([130 170])
xlim([400 430])
%% Example LEC MP trace
clear all
close all
load F:\WC_Germany\overall_EC\overall_EC_dir
used_data = [l3mec_p l3lec_p];
sess_data = sess_data(used_data);
mec_cells = 1:length(l3mec_p);
lec_cells = (mec_cells(end)+1):(mec_cells(end)+length(l3lec_p));

%specify session number
d = 34;
niqf = 2016/2;
lcf = 0.05/niqf;
hcf = 4/niqf;
[b,a] = butter(2,[lcf hcf]);

dsf = 1;
Fsd = 2016/dsf;

cd(sess_data(d).directory)
pwd

load ./used_data lf8 wcv wcv_minus_spike

lf8_f = filtfilt(b,a,lf8);

lf8_d = downsample(lf8_f,dsf);
wcv_d = downsample(wcv,dsf);
wcv_std = std(wcv)*1e3; %in mv
twenty = 20/wcv_std;

lf8_d = zscore(lf8_d);
wcv_d = zscore(wcv_d);

wcv = zscore(wcv);
lf8 = zscore(lf8);
t = (1:length(wcv_d))/2016*dsf;

plot(t,wcv_d,'linewidth',2), hold on
xlabel('Time (s)','FontSize',14)
ylabel('Amplitude (zscore)','FontSize',14)
line([5 5]+430,[0 0+twenty],'color','k')
% title('MEC Membrane Potential','FontSize',16)
% xlim([214 230])
% ylim([-3 8])
% xlim([970 986])
% xlim([430 446])
% xlim([633 649])
% xlim([1332 1348]) %for d = 28
% xlim([616.5 632.5]) %for d = 34
xlim([390 406])
xlim([430 446])

% figure
plot(t,lf8_d,'r','linewidth',2)
xlabel('Time (s)','FontSize',14)
ylabel('Amplitude (zscore)','FontSize',14)
% title('Cortical LFP','FontSize',16)
% legend('MP','LFP')
% xlim([214 230])
% ylim([-2 3])
% xlim([430 446])
% xlim([633 649])
% xlim([1332 1348]) %for d = 28
% xlim([616.5 632.5]) %for d = 34
xlim([390 406])
xlim([430 446])

load ./aligned_heka.mat
figure
plot(dc_time,dc_data)
xlim([430 446])

%% Up State Distribution Compare
clear all
close all
load G:\WC_Germany\overall_EC\overall_EC_dir

used_data = [l3mec_p l3lec_p];
sess_data = sess_data(used_data);

mec_cells = 1:length(l3mec_p);
lec_cells = (mec_cells(end)+1):(mec_cells(end)+length(l3lec_p));

cd G:\WC_Germany\persistent_9_27_2010
% load ./pa_state_dur_stats2
load ./pa_state_dur_stats_new2.mat

%use logarithmic hists
up_hist = up_loghist;
up_hist8 = up_loghist8;
% up_hist = up_hist;
% up_hist8 = up_hist8;

%normalize state duration distributions 
up_hist = up_hist./repmat(sum(up_hist,2),1,size(up_hist,2));
up_hist8 = up_hist8./repmat(sum(up_hist8,2),1,size(up_hist8,2));

up_hist = cumsum(up_hist,2);
up_hist8 = cumsum(up_hist8,2);

mean_up_mec = mean(up_hist(mec_cells,:));
u_up_mec = mean_up_mec + 1*std(up_hist(mec_cells,:))/sqrt(length(mec_cells));
l_up_mec = mean_up_mec - 1*std(up_hist(mec_cells,:))/sqrt(length(mec_cells));

mean_up_lfp = mean(up_hist8([mec_cells lec_cells],:));
u_up_lfp = mean_up_lfp + 1*std(up_hist8([mec_cells lec_cells],:))/sqrt(length(mec_cells)+length(lec_cells));
l_up_lfp = mean_up_lfp - 1*std(up_hist8([mec_cells lec_cells],:))/sqrt(length(mec_cells)+length(lec_cells));

mean_up_lec = mean(up_hist(lec_cells,:));
u_up_lec = mean_up_lec + 1*std(up_hist(lec_cells,:))/sqrt(length(lec_cells));
l_up_lec = mean_up_lec - 1*std(up_hist(lec_cells,:))/sqrt(length(lec_cells));

figure
stairs(log_grid,mean_up_lfp,'r','linewidth',2)
hold on
stairs(log_grid,mean_up_mec,'linewidth',2)
stairs(log_grid,mean_up_lec,'g','linewidth',2)

% legend('LFP','MEC','LEC')

epsilon = 1e-7;
u_up_lfp(u_up_lfp < epsilon) = epsilon;
l_up_lfp(l_up_lfp < epsilon) = epsilon;
mean_up_lfp(mean_up_lfp < epsilon) = epsilon;
u_up_mec(u_up_mec < epsilon) = epsilon;
l_up_mec(l_up_mec < epsilon) = epsilon;
mean_up_mec(mean_up_mec < epsilon) = epsilon;
u_up_lec(u_up_lec < epsilon) = epsilon;
l_up_lec(l_up_lec < epsilon) = epsilon;
mean_up_lec(mean_up_lec < epsilon) = epsilon;

for i = 1:length(log_grid)-1
    clear X Y
    bin_width = (log_grid(i+1)-log_grid(i));
    X = [log_grid(i) log_grid(i)+bin_width];
    Y = [u_up_lec(i) u_up_lec(i)];
    X = [X fliplr(X)];
    Y = [Y l_up_lec(i) l_up_lec(i)];
    fill(X,Y,'g','EdgeColor','none');
    hold on
end

for i = 1:length(log_grid)-1
    clear X Y
    bin_width = (log_grid(i+1)-log_grid(i));
    X = [log_grid(i) log_grid(i)+bin_width];
    Y = [u_up_lfp(i) u_up_lfp(i)];
    X = [X fliplr(X)];
    Y = [Y l_up_lfp(i) l_up_lfp(i)];
    fill(X,Y,'r','EdgeColor','none');
    hold on
end

for i = 1:length(log_grid)-1
    clear X Y
    bin_width = (log_grid(i+1)-log_grid(i));
    X = [log_grid(i) log_grid(i)+bin_width];
    Y = [u_up_mec(i) u_up_mec(i)];
    X = [X fliplr(X)];
    Y = [Y l_up_mec(i) l_up_mec(i)];
    fill(X,Y,'b','EdgeColor','none');
    hold on
end


% set(gca,'yscale','log')
% set(gca,'YMinorGrid','off')
% set(gca,'YGrid','off')
% set(gca,'XMinorGrid','off')
% ylim([1e-4 0.2])
% xlim([0.05 10])
ylim([0 1])
xlim([0 2.5])
shg
% set(gca,'xscale','log')
xlabel('Up State Duration (s)','FontSize',14)
ylabel('Percent of Data','FontSize',14)

% figure
% stairs(lin_grid,mean_up_lfp,'r','linewidth',2)
% % set(gca,'XMinorGrid','off')
% ylim([0 1])
% xlim([0 2.5])
% shg
% % set(gca,'xscale','log')
% xlabel('Down State Duration (s)','FontSize',14)
% ylabel('Percent of Data','FontSize',14)
% line([0.5 0.5],[0 0.1],'color','k')
% line([0 0.5],[0.1 0.1],'color','k')

%% Down State Distribution Compare
clear all
close all
load G:\WC_Germany\overall_EC\overall_EC_dir

used_data = [l3mec_p l3lec_p];
sess_data = sess_data(used_data);

mec_cells = 1:length(l3mec_p);
lec_cells = (mec_cells(end)+1):(mec_cells(end)+length(l3lec_p));

cd G:\WC_Germany\persistent_9_27_2010
load ./pa_state_dur_stats_new2.mat

%use logarithmic hists
down_hist = down_loghist;
down_hist8 = down_loghist8;

%normalize state duration distributions 
down_hist = down_hist./repmat(sum(down_hist,2),1,size(down_hist,2));
down_hist8 = down_hist8./repmat(sum(down_hist8,2),1,size(down_hist8,2));

% down_hist = cumsum(down_hist,2);
% down_hist8 = cumsum(down_hist8,2);

mean_down_mec = mean(down_hist(mec_cells,:));
u_down_mec = mean_down_mec + 1*std(down_hist(mec_cells,:))/sqrt(length(mec_cells));
l_down_mec = mean_down_mec - 1*std(down_hist(mec_cells,:))/sqrt(length(mec_cells));

mean_down_lfp = mean(down_hist8([mec_cells lec_cells],:));
u_down_lfp = mean_down_lfp + 1*std(down_hist8([mec_cells lec_cells],:))/sqrt(length(mec_cells)+length(lec_cells));
l_down_lfp = mean_down_lfp - 1*std(down_hist8([mec_cells lec_cells],:))/sqrt(length(mec_cells)+length(lec_cells));

mean_down_lec = mean(down_hist(lec_cells,:));
u_down_lec = mean_down_lec + 1*std(down_hist(lec_cells,:))/sqrt(length(lec_cells));
l_down_lec = mean_down_lec - 1*std(down_hist(lec_cells,:))/sqrt(length(lec_cells));

figure
stairs(log_grid,mean_down_lfp,'r','linewidth',2)
hold on
stairs(log_grid,mean_down_mec,'linewidth',2)
stairs(log_grid,mean_down_lec,'g','linewidth',2)

legend('LFP','MEC','LEC')

epsilon = 1e-7;
u_down_lfp(u_down_lfp < epsilon) = epsilon;
l_down_lfp(l_down_lfp < epsilon) = epsilon;
mean_down_lfp(mean_down_lfp < epsilon) = epsilon;
u_down_mec(u_down_mec < epsilon) = epsilon;
l_down_mec(l_down_mec < epsilon) = epsilon;
mean_down_mec(mean_down_mec < epsilon) = epsilon;
u_down_lec(u_down_lec < epsilon) = epsilon;
l_down_lec(l_down_lec < epsilon) = epsilon;
mean_down_lec(mean_down_lec < epsilon) = epsilon;

for i = 1:length(log_grid)-1
    clear X Y
    bin_width = (log_grid(i+1)-log_grid(i));
    X = [log_grid(i) log_grid(i)+bin_width];
    Y = [u_down_lec(i) u_down_lec(i)];
    X = [X fliplr(X)];
    Y = [Y l_down_lec(i) l_down_lec(i)];
    fill(X,Y,'g','EdgeColor','none');
    hold on
end

for i = 1:length(log_grid)-1
    clear X Y
    bin_width = (log_grid(i+1)-log_grid(i));
    X = [log_grid(i) log_grid(i)+bin_width];
    Y = [u_down_lfp(i) u_down_lfp(i)];
    X = [X fliplr(X)];
    Y = [Y l_down_lfp(i) l_down_lfp(i)];
    fill(X,Y,'r','EdgeColor','none');
    hold on
end

for i = 1:length(log_grid)-1
    clear X Y
    bin_width = (log_grid(i+1)-log_grid(i));
    X = [log_grid(i) log_grid(i)+bin_width];
    Y = [u_down_mec(i) u_down_mec(i)];
    X = [X fliplr(X)];
    Y = [Y l_down_mec(i) l_down_mec(i)];
    fill(X,Y,'b','EdgeColor','none');
    hold on
end


set(gca,'yscale','log')
set(gca,'YMinorGrid','off')
set(gca,'YGrid','off')
% set(gca,'XMinorGrid','off')
ylim([1e-4 0.1])
xlim([0.05 10])
% ylim([0 1])
% xlim([0 2.5])
shg
% set(gca,'xscale','log')
xlabel('Down State Duration (s)','FontSize',14)
ylabel('Percent of Data','FontSize',14)


%% Quantization Figure Overall
clear all
% close all
load G:\WC_Germany\overall_EC\overall_EC_dir

used_data = [l3mec_p l3lec_p];
sess_data = sess_data(used_data);

mec_cells = 1:length(l3mec_p);
lec_cells = (mec_cells(end)+1):(mec_cells(end)+length(l3lec_p));
n_mec = length(mec_cells);
n_lec = length(lec_cells);

cd G:\WC_Germany\persistent_9_27_2010
% load ./pa_corresponding_lfp_state_data_rtest
load ./pa_corresponding_lfp_revised_simp_new2.mat

hist_range = linspace(0,7,1000);
binsize = hist_range(2)-hist_range(1);

for i = 1:n_mec+n_lec
    useable_ups = find(corresp_lfp_dmindowndur{i} > 0.5);
    cell_hist(i,:) = hist(mp_updur_lfpc_delay{i}(useable_ups),hist_range);
%     cell_hist(i,:) = hist(mp_updur_lfpc{i}(useable_ups),hist_range);
    cell_hist(i,:) = cell_hist(i,:)/sum(cell_hist(i,:))/binsize;
    cell_hist(i,:) = jmm_smooth_1d_cor(cell_hist(i,:),25);
end

% load F:\WC_Germany\persistent_9_27_2010\pa_state_dur_stats
% for i = 1:n_mec+n_lec
%     p25d = prctile(lfp_state_durations{i}{1},25); %25th percentile of LFP down state durations
%     used_corr_lfps = find(corresp_lfp_downdur{i} > p25d); %only count LFP down states which are in the upper 3/4
%     cell_hist(i,:) = hist(mp_updur_lfpc_delay{i}(used_corr_lfps),hist_range);
%     cell_hist(i,:) = cell_hist(i,:)/sum(cell_hist(i,:))/binsize;
%     cell_hist(i,:) = jmm_smooth_1d_cor(cell_hist(i,:),25);
% end
% 

mean_hist_mec = nanmean(cell_hist(mec_cells,:));
u_mec = mean_hist_mec+1*nanstd(cell_hist(mec_cells,:))./sqrt(sum(~isnan(cell_hist(mec_cells,:))));
l_mec = mean_hist_mec-1*nanstd(cell_hist(mec_cells,:))./sqrt(sum(~isnan(cell_hist(mec_cells,:))));

mean_hist_lec = nanmean(cell_hist(lec_cells,:));
u_lec = mean_hist_lec+1*nanstd(cell_hist(lec_cells,:))./sqrt(sum(~isnan(cell_hist(lec_cells,:))));
l_lec = mean_hist_lec-1*nanstd(cell_hist(lec_cells,:))./sqrt(sum(~isnan(cell_hist(lec_cells,:))));

figure
plot(hist_range,mean_hist_mec,'linewidth',2)
hold on
plot(hist_range,mean_hist_lec,'g','linewidth',2)
legend('MEC','LEC')

epsilon = 1e-10;
cell_hist(cell_hist < epsilon) = epsilon;
u_mec(u_mec < epsilon) = epsilon;
l_mec(l_mec < epsilon) = epsilon;
u_lec(u_lec < epsilon) = epsilon;
l_lec(l_lec < epsilon) = epsilon;

X = [hist_range fliplr(hist_range)];
Y = [u_mec fliplr(l_mec)];
fill(X,Y,'b')
Y = [u_lec fliplr(l_lec)];
fill(X,Y,'g')
% 
set(gca,'yscale','log')
xlim([0 6.2])
ylim([5e-4 3])

[peak_amps,b] = findpeaks(mean_hist_mec);
peak_locs = hist_range(b);
plot(peak_locs,peak_amps,'ro')
% for i = 1:length(peak_locs)
%     line([peak_locs(i) peak_locs(i)],[5e-4 2],'Color','k')
% end
temp_p = polyfit(peak_locs,log10(peak_amps),1);
pred_p = 10.^polyval(temp_p,hist_range);
plot(hist_range,pred_p,'k','linewidth',2)
xlabel('MP Up state duration (LFP cycles)')
ylabel('Probability')

k = 1:7;
pooled_mec_ups = [];
for i = 1:length(mec_cells)
    useable_ups = find(corresp_lfp_mindowndur{i} > 0.5);
    cur_ups = mp_updur_lfpc_delay{i}(useable_ups);
    cur_ups = ceil(cur_ups);
    cur_ups(isnan(cur_ups)) = [];
    cur_ups(cur_ups > 12) = [];
    emp(i,:) = hist(cur_ups,k);
%     emp(i,:) = emp(i,:)/sum(emp(i,:));
    p(i) = 1/mean(cur_ups);
    geom(i,:) = (1-p(i)).^(k-1).*p(i);
    pred_geom(i,:) = geom(i,:)*length(cur_ups);
    chisqr(i) = sum((emp(i,:)-pred_geom(i,:)).^2./pred_geom(i,:));
    pooled_mec_ups = [pooled_mec_ups; cur_ups];
end

p_pooled = 1/mean(pooled_mec_ups);
geom_pooled = (1-p_pooled).^(k-1).*p_pooled;
pred_geom_pooled = geom_pooled*length(pooled_mec_ups);
emp_pooled = hist(pooled_mec_ups,k);
chisqr_pooled = sum((emp_pooled(1:5)-pred_geom_pooled(1:5)).^2./pred_geom_pooled(1:5));
emp = emp/sum(emp);

% figure
% plot(k,emp,'o')
% hold on
% plot(k,geom,'r.-')
%% Quantization Figure Single Example
clear all
% close all
load G:\WC_Germany\overall_EC\overall_EC_dir

used_data = [l3mec_p l3lec_p];
sess_data = sess_data(used_data);

mec_cells = 1:length(l3mec_p);
lec_cells = (mec_cells(end)+1):(mec_cells(end)+length(l3lec_p));
n_mec = length(mec_cells);
n_lec = length(lec_cells);
load G:\WC_Germany\persistent_9_27_2010\pa_corresponding_lfp_revised_simp_new2

%good examples are 1,2,5,8,14,15,19, 21
hist_range = linspace(0,5,1000);
binsize = hist_range(2)-hist_range(1);
for i = 1:n_mec+n_lec
    cell_hist(i,:) = hist(rmp_updurs_lfpc{i},hist_range);
%     cell_hist(i,:) = hist(mp_updur_lfpc{i}(useable_ups),hist_range);
    cell_hist(i,:) = cell_hist(i,:)/sum(cell_hist(i,:))/binsize;
    cell_hist(i,:) = jmm_smooth_1d_cor(cell_hist(i,:),20);
end

examp_ind = 21;
x = linspace(0,7,100);
binsize = x(2)-x(1);

y = hist(rmp_updurs_lfpc{examp_ind},x);
sy = sum(y);

figure
set(gca,'fontsize',14,'fontname','arial')
bar(x,y/sy/binsize,1,'k')
hold on
temp = jmm_smooth_1d_cor(y,2);
temp = temp/sy;
% fine_axis = linspace(0,7,1000);
% smooth_temp = spline(x,temp,fine_axis);
% plot(x,jmm_smooth_1d_cor(y,2),'r','linewidth',2)
plot(hist_range,cell_hist(examp_ind,:),'r','linewidth',2)
xlim([0 4.2])
% ylim([0 0.15])
xlabel('Duration (Ncx UDS Cycles)','fontsize',16,'fontname','arial')
ylabel('Probability Density','fontsize',16,'fontname','arial')
shg

%% Histogram of Type 1 Persistent Activity
clear all
close all
load G:\WC_Germany\overall_EC\overall_EC_dir

used_data = [l3mec_p l3lec_p];
sess_data = sess_data(used_data);

mec_cells = 1:length(l3mec_p);
lec_cells = (mec_cells(end)+1):(mec_cells(end)+length(l3lec_p));
n_mec = length(mec_cells);
n_lec = length(lec_cells);

cd G:\WC_Germany\persistent_9_27_2010
% load ./pa_corresponding_lfp_state_data_rtest
load ./pa_corresponding_lfp_revised_simp_new2

numBins = 20;
pers_range = linspace(0,1.,numBins+1);

persd_within = fract_rt1_ups;
persd_within_np = fract_rt1_ups_nt2;

pers_mec = histc(persd_within(mec_cells),pers_range);
pers_lec = histc(persd_within(lec_cells),pers_range);
pers_mec_np = histc(persd_within_np(mec_cells),pers_range);


stairs(pers_range*100,pers_mec,'linewidth',2)
hold on
stairs(pers_range*100,pers_lec,'g','linewidth',2)
stairs(pers_range*100,pers_mec_np,'y','linewidth',2)
% xlim([0 65])
ylim([0 8])
xlabel('Fraction of Up States Persisting','FontSize',14)
ylabel('Number of Cells','FontSize',14)
legend('MEC','LEC','MEC Nt2')

%% Histogram of Type 1 UP-skipping
clear all
close all
load G:\WC_Germany\overall_EC\overall_EC_dir

used_data = [l3mec_p l3lec_p];
sess_data = sess_data(used_data);

mec_cells = 1:length(l3mec_p);
lec_cells = (mec_cells(end)+1):(mec_cells(end)+length(l3lec_p));
n_mec = length(mec_cells);
n_lec = length(lec_cells);

cd G:\WC_Germany\persistent_9_27_2010
load ./pa_corresponding_lfp_state_data_downs

numBins = 20;
pers_range = linspace(0,1.,numBins+1);

pers_mec = histc(down_pers_within(mec_cells),pers_range);
pers_lec = histc(down_pers_within(lec_cells),pers_range);


stairs(pers_range*100,pers_mec,'linewidth',2)
hold on
stairs(pers_range*100,pers_lec,'g','linewidth',2)
% xlim([0 65])
ylim([0 8])
xlabel('Fraction of Up States Persisting','FontSize',14)
ylabel('Number of Cells','FontSize',14)
legend('MEC','LEC','MEC Nt2')

%% Histogram of Type 2 Persistent Activity
clear all
% close all
load G:\WC_Germany\overall_EC\overall_EC_dir

used_data = [l3mec_p l3lec_p];
sess_data = sess_data(used_data);

mec_cells = 1:length(l3mec_p);
lec_cells = (mec_cells(end)+1):(mec_cells(end)+length(l3lec_p));
n_mec = length(mec_cells);
n_lec = length(lec_cells);

cd G:\WC_Germany\persistent_9_27_2010
load ./pa_corresponding_lfp_revised_simp_new2.mat

numBins = 15;
pers_range = linspace(0,0.4,numBins+1);

% pers_mec = histc(pers_across(mec_cells),pers_range);
% pers_lec = histc(pers_across(lec_cells),pers_range);
% pers_mec = histc(rpersd_across(mec_cells),pers_range);
% pers_lec = histc(rpersd_across(lec_cells),pers_range);
pers_mec = histc(fract_rt2_ups(mec_cells),pers_range);
pers_lec = histc(fract_rt2_ups(lec_cells),pers_range);

figure
stairs(pers_range*100,pers_mec,'linewidth',2)
hold on
stairs(pers_range*100,pers_lec,'g','linewidth',2)
xlim([0 40])
% ylim([0 5])
xlabel('Fraction of Up States Persisting','FontSize',14)
ylabel('Number of Cells','FontSize',14)
legend('MEC','LEC')

% figure
% plot(mean_down_lag(mec_cells),100*pers_across(mec_cells),'o')
% hold on
% plot(mean_down_lag(lec_cells),100*pers_across(lec_cells),'go')
% xlabel('Down-transition lag (s)')
% ylabel('Percent Type 2 persistence')

% figure
% plot(mean_down_lag(mec_cells),100*rpersd_across(mec_cells),'o')
% hold on
% plot(mean_down_lag(lec_cells),100*rpersd_across(lec_cells),'go')
% xlabel('Down-transition lag (s)')
% ylabel('Percent Type 2 persistence')
% lag_ax1 = linspace(0.3,1.4,400);
% temp_p = polyfit(mean_down_lag(mec_cells),100*rpersd_across(mec_cells),1);
% plot(lag_ax1,polyval(temp_p,lag_ax1))
% lag_ax2 = linspace(-0.3,0.3,400);
% temp_p = polyfit(mean_down_lag(lec_cells),100*rpersd_across(lec_cells),1);
% plot(lag_ax2,polyval(temp_p,lag_ax2),'g')
% ylim([0 50])

%% Histogram of Type 2 up-skipping
clear all
% close all
load G:\WC_Germany\overall_EC\overall_EC_dir

used_data = [l3mec_p l3lec_p];
sess_data = sess_data(used_data);

mec_cells = 1:length(l3mec_p);
lec_cells = (mec_cells(end)+1):(mec_cells(end)+length(l3lec_p));
n_mec = length(mec_cells);
n_lec = length(lec_cells);

cd G:\WC_Germany\persistent_9_27_2010
load ./pa_corresponding_lfp_state_data_downs

numBins = 15;
pers_range = linspace(0,0.70,numBins+1);

% pers_mec = histc(pers_across(mec_cells),pers_range);
% pers_lec = histc(pers_across(lec_cells),pers_range);
pers_mec = histc(rdown_pers_across(mec_cells),pers_range);
pers_lec = histc(rdown_pers_across(lec_cells),pers_range);

figure
stairs(pers_range*100,pers_mec,'linewidth',2)
hold on
stairs(pers_range*100,pers_lec,'g','linewidth',2)
xlim([0 70])
% ylim([0 5])
xlabel('Fraction of Up States Persisting','FontSize',14)
ylabel('Number of Cells','FontSize',14)
legend('MEC','LEC')


%% UP TRAN LAG
clear all
close all
load F:\WC_Germany\overall_EC\overall_EC_dir

used_data = [l3mec_p l3lec_p];
sess_data = sess_data(used_data);

mec_cells = 1:length(l3mec_p);
lec_cells = (mec_cells(end)+1):(mec_cells(end)+length(l3lec_p));
n_mec = length(mec_cells);
n_lec = length(lec_cells);

cd F:\WC_Germany\persistent_9_27_2010
load ./pa_corresponding_lfp_state_data

lag_range = linspace(-3,3,1000);
dlag = lag_range(2)-lag_range(1);
for d = 1:length(sess_data)
    lag_hist(d,:) = histc(lfp_up_lag{d},lag_range);
    lag_hist(d,:) = lag_hist(d,:)/sum(lag_hist(d,:))/dlag;
        non_pers_states = find(mp_updur_lfpc{d} < 1);
    lag_hist_np(d,:) = histc(lfp_up_lag{d}(non_pers_states),lag_range);
    lag_hist_np(d,:) = lag_hist_np(d,:)/sum(lag_hist_np(d,:))/dlag;
    sm_lag_hist(d,:) = jmm_smooth_1d_cor(lag_hist(d,:),10);
    sm_lag_hist_np(d,:) = jmm_smooth_1d_cor(lag_hist_np(d,:),10);
%     phase_hist(d,:) = fgsmooth(phase_hist(d,:),4);
end

mean_lag_mec = mean(sm_lag_hist(mec_cells,:));
u_lag_mec = mean_lag_mec + 1*std(sm_lag_hist(mec_cells,:))/sqrt(length(mec_cells));
l_lag_mec = mean_lag_mec - 1*std(sm_lag_hist(mec_cells,:))/sqrt(length(mec_cells));
mean_lag_lec = mean(sm_lag_hist(lec_cells,:));
u_lag_lec = mean_lag_lec + 1*std(sm_lag_hist(lec_cells,:))/sqrt(length(lec_cells));
l_lag_lec = mean_lag_lec - 1*std(sm_lag_hist(lec_cells,:))/sqrt(length(lec_cells));
mean_lag_np_mec = mean(sm_lag_hist_np(mec_cells,:));
u_lag_np_mec = mean_lag_np_mec + 1*std(sm_lag_hist_np(mec_cells,:))/sqrt(length(mec_cells));
l_lag_np_mec = mean_lag_np_mec - 1*std(sm_lag_hist_np(mec_cells,:))/sqrt(length(mec_cells));

plot(lag_range,mean_lag_mec,'linewidth',2)
hold on
plot(lag_range,mean_lag_lec,'g','linewidth',2)
plot(lag_range,mean_lag_np_mec,'y','linewidth',2)
% plot(lag_range,mean_lag_pre,'k','linewidth',2)
% plot(lag_range,mean_lag_fro,'c','linewidth',2)
% legend('MEC','LEC')

X = [lag_range fliplr(lag_range)];
Y = [u_lag_lec fliplr(l_lag_lec)];
fill(X,Y,'g')
X = [lag_range fliplr(lag_range)];
Y = [u_lag_mec fliplr(l_lag_mec)];
fill(X,Y,'b')
X = [lag_range fliplr(lag_range)];
Y = [u_lag_np_mec fliplr(l_lag_np_mec)];
fill(X,Y,'c')
xlim([-0.5 1.])
% ylim([0 3])
line([0 0],[0 3],'Color','k')
xlabel('Lag (s)')
ylabel('Probability')

%% DOWN TRAN LAG
clear all
close all
load F:\WC_Germany\overall_EC\overall_EC_dir

used_data = [l3mec_p l3lec_p];
sess_data = sess_data(used_data);

mec_cells = 1:length(l3mec_p);
lec_cells = (mec_cells(end)+1):(mec_cells(end)+length(l3lec_p));
n_mec = length(mec_cells);
n_lec = length(lec_cells);

cd F:\WC_Germany\persistent_9_27_2010
load ./pa_corresponding_lfp_state_data

lag_range = linspace(-1,5,1000);
dlag = lag_range(2)-lag_range(1);
for d = 1:length(sess_data)
    lag_hist(d,:) = histc(lfp_down_lag{d},lag_range);
    lag_hist(d,:) = lag_hist(d,:)/sum(lag_hist(d,:))/dlag;
    useable_ups = find(corresp_lfp_downdur{d} > 0.5);
    non_pers_states = find(mp_updur_lfpc_delay{d} < 1);
    pers_states = useable_ups(mp_updur_lfpc_delay{d}(useable_ups) > 1);
    mean_down_pers_lag_np(d) = mean(lfp_down_lag{d}(non_pers_states));
    mean_down_pers_lag_p(d) = mean(lfp_down_lag{d}(pers_states));
    lag_hist_np(d,:) = histc(lfp_down_lag{d}(non_pers_states),lag_range);
    lag_hist_np(d,:) = lag_hist_np(d,:)/sum(lag_hist_np(d,:))/dlag;
    up_lag_hist(d,:) = histc(lfp_up_lag{d},lag_range);
    up_lag_hist(d,:) = up_lag_hist(d,:)/sum(up_lag_hist(d,:))/dlag;
    sm_lag_hist(d,:) = jmm_smooth_1d_cor(lag_hist(d,:),20);
    sm_lag_hist_np(d,:) = jmm_smooth_1d_cor(lag_hist_np(d,:),20);
    sm_up_lag_hist(d,:) = jmm_smooth_1d_cor(up_lag_hist(d,:),20);
%     phase_hist(d,:) = fgsmooth(phase_hist(d,:),4);
end

mean_lag_mec = mean(sm_lag_hist(mec_cells,:));
u_lag_mec = mean_lag_mec + 1*std(sm_lag_hist(mec_cells,:))/sqrt(length(mec_cells));
l_lag_mec = mean_lag_mec - 1*std(sm_lag_hist(mec_cells,:))/sqrt(length(mec_cells));
mean_lag_lec = mean(sm_lag_hist(lec_cells,:));
u_lag_lec = mean_lag_lec + 1*std(sm_lag_hist(lec_cells,:))/sqrt(length(lec_cells));
l_lag_lec = mean_lag_lec - 1*std(sm_lag_hist(lec_cells,:))/sqrt(length(lec_cells));
mean_uplag_mec = mean(sm_up_lag_hist(mec_cells,:));
u_uplag_mec = mean_uplag_mec + 1*std(sm_up_lag_hist(mec_cells,:))/sqrt(length(mec_cells));
l_uplag_mec = mean_uplag_mec - 1*std(sm_up_lag_hist(mec_cells,:))/sqrt(length(mec_cells));
mean_uplag_lec = mean(sm_up_lag_hist(lec_cells,:));
u_uplag_lec = mean_uplag_lec + 1*std(sm_up_lag_hist(lec_cells,:))/sqrt(length(lec_cells));
l_uplag_lec = mean_uplag_lec - 1*std(sm_up_lag_hist(lec_cells,:))/sqrt(length(lec_cells));
mean_lag_np_mec = mean(sm_lag_hist_np(mec_cells,:));
u_lag_np_mec = mean_lag_np_mec + 1*std(sm_lag_hist_np(mec_cells,:))/sqrt(length(mec_cells));
l_lag_np_mec = mean_lag_np_mec - 1*std(sm_lag_hist_np(mec_cells,:))/sqrt(length(mec_cells));

plot(lag_range,mean_lag_mec,'linewidth',2)
hold on
plot(lag_range,mean_lag_lec,'g','linewidth',2)
plot(lag_range,mean_uplag_mec,'r','linewidth',2)
% plot(lag_range,mean_lag_np_mec,'y','linewidth',2)
% plot(lag_range,mean_lag_pre,'k','linewidth',2)
% plot(lag_range,mean_lag_fro,'c','linewidth',2)

X = [lag_range fliplr(lag_range)];
Y = [u_lag_lec fliplr(l_lag_lec)];
fill(X,Y,'g')
X = [lag_range fliplr(lag_range)];
Y = [u_lag_mec fliplr(l_lag_mec)];
fill(X,Y,'b')
X = [lag_range fliplr(lag_range)];
Y = [u_uplag_mec fliplr(l_uplag_mec)];
fill(X,Y,'r')
% X = [lag_range fliplr(lag_range)];
% Y = [u_lag_np_mec fliplr(l_lag_np_mec)];
% fill(X,Y,'c')
xlim([-0.5 2.])
line([0 0],[0 2],'Color','k')
xlabel('Lag (s)')
ylabel('Probability')

%% Supplementary figure all MP autocorr
clear all
close all
load G:\WC_Germany\overall_EC\overall_EC_dir

used_data = [l3mec_p l3lec_p];
sess_data = sess_data(used_data);

mec_cells = 1:length(l3mec_p);
lec_cells = (mec_cells(end)+1):(mec_cells(end)+length(l3lec_p));
n_mec = length(mec_cells);
n_lec = length(lec_cells);

cd G:\WC_Germany\persistent_9_27_2010
load ./time_domain_data.mat

zero_lag = find(lags/Fsd >= 0,1,'first');

m_corr_mec = mean(tot_wcv_acorr(mec_cells,:));
u_corr_mec = m_corr_mec + 1*std(tot_wcv_acorr(mec_cells,:))/sqrt(length(mec_cells));
l_corr_mec = m_corr_mec - 1*std(tot_wcv_acorr(mec_cells,:))/sqrt(length(mec_cells));

% [min_mec,min_loc] = min(m_corr_mec(zero_lag:end));
% [mec_peaks,mec_peak_locs] = findpeaks(m_corr_mec(zero_lag+min_loc:end));
% for i = 1:length(mec_peaks)
%     peak_test_mec(i) = ttest(tot_wcv_acorr(mec_cells,zero_lag+min_loc+mec_peak_locs(i)));
% end
[min_mec,min_mec_loc] = min(m_corr_mec(zero_lag:end));
[peak_mec,peak_loc_mec] = max(m_corr_mec(zero_lag+min_mec_loc:end));
peak_loc_mec = peak_loc_mec + zero_lag+min_mec_loc;
mec_peaks = tot_wcv_acorr(mec_cells,peak_loc_mec);

m_corr_lec = mean(tot_wcv_acorr(lec_cells,:));
u_corr_lec = m_corr_lec + 1*std(tot_wcv_acorr(lec_cells,:))/sqrt(length(lec_cells));
l_corr_lec = m_corr_lec - 1*std(tot_wcv_acorr(lec_cells,:))/sqrt(length(lec_cells));

[min_lec,min_lec_loc] = min(m_corr_lec(zero_lag:end));
[peak_lec,peak_loc_lec] = max(m_corr_lec(zero_lag+min_lec_loc:end));
peak_loc_lec = peak_loc_lec + zero_lag+min_lec_loc;
lec_peaks = tot_wcv_acorr(lec_cells,peak_loc_lec);

m_corr_lfp = mean(tot_lf8_acorr([mec_cells lec_cells],:));
u_corr_lfp = m_corr_lfp + 1*std(tot_lf8_acorr([mec_cells lec_cells],:))/sqrt(length(mec_cells)+length(lec_cells));
l_corr_lfp = m_corr_lfp - 1*std(tot_lf8_acorr([mec_cells lec_cells],:))/sqrt(length(mec_cells)+length(lec_cells));

[min_lfp,min_lfp_loc] = min(m_corr_lfp(zero_lag:end));
[peak_lfp,peak_loc_lfp] = max(m_corr_lfp(zero_lag+min_lfp_loc:end));
peak_loc_lfp = peak_loc_lfp + zero_lag+min_lfp_loc;
lfp_peaks = tot_lf8_acorr(:,peak_loc_lfp);

[a_lec,b_lec] = ttest(tot_wcv_acorr(lec_cells,peak_loc_lfp),tot_lf8_acorr(lec_cells,peak_loc_lfp))
[a_mec,b_mec] = ttest(tot_wcv_acorr(mec_cells,peak_loc_lfp),tot_lf8_acorr(mec_cells,peak_loc_lfp))


plot(lags/Fsd,m_corr_mec,'linewidth',2)
hold on
plot(lags/Fsd,m_corr_lec,'g','linewidth',2)
plot(lags/Fsd,m_corr_lfp,'r','linewidth',2)

X = [lags/Fsd fliplr(lags)/Fsd];
Y = [l_corr_mec fliplr(u_corr_mec)];
fill(X,Y,'b')
Y = [l_corr_lec fliplr(u_corr_lec)];
fill(X,Y,'g')
Y = [l_corr_lfp fliplr(u_corr_lfp)];
fill(X,Y,'r')

line([0 0],[-0.4 1],'Color','k')
line([-8 8],[0 0],'Color','k')
xlim([-6 6])
ylim([-0.4 0.4])
ylabel('Correlation Coefficient')
xlabel('Time (s)')

%% Supplementary figure for MP/LFP xcorr
clear all
close all
load G:\WC_Germany\overall_EC\overall_EC_dir

used_data = [l3mec_p l3lec_p];
sess_data = sess_data(used_data);

mec_cells = 1:length(l3mec_p);
lec_cells = (mec_cells(end)+1):(mec_cells(end)+length(l3lec_p));
n_mec = length(mec_cells);
n_lec = length(lec_cells);

cd G:\WC_Germany\persistent_9_27_2010
load ./time_domain_data.mat

first_point = find(lags/Fsd > -1,1,'first');

[all_peaks,all_peak_locs] = max(tot_w8_x(:,first_point:end),[],2);
peak_lags = lags(all_peak_locs+first_point)/Fsd;

m_corr_mec = mean(tot_w8_x(mec_cells,:));
u_corr_mec = m_corr_mec + 1*std(tot_w8_x(mec_cells,:))/sqrt(length(mec_cells));
l_corr_mec = m_corr_mec - 1*std(tot_w8_x(mec_cells,:))/sqrt(length(mec_cells));

[pk_mec,pk_mec_loc] = max(m_corr_mec);
corr_at_peak_mec = tot_w8_x(mec_cells,pk_mec_loc);

pk_mec_u = u_corr_mec(pk_mec_loc);
pk_mec_l = l_corr_mec(pk_mec_loc);
pk_shift_mec = lags(pk_mec_loc)/Fsd;

m_corr_lec = mean(tot_w8_x(lec_cells,:));
u_corr_lec = m_corr_lec + 1*std(tot_w8_x(lec_cells,:))/sqrt(length(lec_cells));
l_corr_lec = m_corr_lec - 1*std(tot_w8_x(lec_cells,:))/sqrt(length(lec_cells));

[pk_lec,pk_lec_loc] = max(m_corr_lec);
pk_lec_u = u_corr_lec(pk_lec_loc);
pk_lec_l = l_corr_lec(pk_lec_loc);
pk_shift_lec = lags(pk_lec_loc)/Fsd;
corr_at_peak_lec = tot_w8_x(lec_cells,pk_lec_loc);

plot(lags/Fsd,m_corr_mec,'linewidth',2)
hold on
plot(lags/Fsd,m_corr_lec,'g','linewidth',2)

X = [lags/Fsd fliplr(lags)/Fsd];
Y = [l_corr_mec fliplr(u_corr_mec)];
fill(X,Y,'b')
Y = [l_corr_lec fliplr(u_corr_lec)];
fill(X,Y,'g')

line([0 0],[-0.5 1],'Color','k')
line([-8 8],[0 0],'Color','k')
xlim([-2.5 2.5])
ylim([-0.4 0.8])
ylabel('Correlation Coefficient')
xlabel('Time (s)')

figure
set(gca,'fontsize',14,'fontname','arial')
plot(peak_lags(mec_cells),all_peaks(mec_cells),'o','markersize',8)
hold on
plot(peak_lags(lec_cells),all_peaks(lec_cells),'go','markersize',8)
xlabel('Peak lags (s)')
ylabel('Peak correlation')
temp_p = polyfit(peak_lags,all_peaks',1);
lag_x = linspace(-0.1,1,100);
plot(lag_x,polyval(temp_p,lag_x),'k')
xlim([-0.1 1])
ylim([0.1 0.8])

load ./simcort_timedomain_data
load G:\WC_Germany\persistent_revised\pers_sim_cortlfp_dir
pre_cells = 1:16;
fro_cells = 17:24;

m_corr_pre = mean(tot_w8_x(pre_cells,:));
u_corr_pre = m_corr_pre + 1*std(tot_w8_x(pre_cells,:))/sqrt(length(pre_cells));
l_corr_pre = m_corr_pre - 1*std(tot_w8_x(pre_cells,:))/sqrt(length(pre_cells));

m_corr_fro = mean(tot_w8_x(fro_cells,:));
u_corr_fro = m_corr_fro + 1*std(tot_w8_x(fro_cells,:))/sqrt(length(fro_cells));
l_corr_fro = m_corr_fro - 1*std(tot_w8_x(fro_cells,:))/sqrt(length(fro_cells));

[all_peaks,all_peak_locs] = max(tot_w8_x(:,first_point:end),[],2);
peak_lags = lags(all_peak_locs+first_point)/Fsd;

figure(1)
plot(lags/Fsd,m_corr_pre,'c','linewidth',2)
hold on
plot(lags/Fsd,m_corr_fro,'k','linewidth',2)

X = [lags/Fsd fliplr(lags)/Fsd];
Y = [l_corr_pre fliplr(u_corr_pre)];
fill(X,Y,'c')
Y = [l_corr_fro fliplr(u_corr_fro)];
fill(X,Y,'k')
ylim([-0.5 1])
%% Supplementary figure for MP/LFP xcorr LF3
clear all
close all
load F:\WC_Germany\overall_EC\overall_EC_dir

used_data = [l3mec_p l3lec_p];
sess_data = sess_data(used_data);

mec_cells = 1:length(l3mec_p);
lec_cells = (mec_cells(end)+1):(mec_cells(end)+length(l3lec_p));
n_mec = length(mec_cells);
n_lec = length(lec_cells);

cd F:\WC_Germany\persistent_9_27_2010
load ./time_domain_data.mat

first_point = find(lags/Fsd > -1,1,'first');

[all_peaks,all_peak_locs] = max(tot_w3_x(:,first_point:end),[],2);
peak_lags = lags(all_peak_locs+first_point)/Fsd;

m_corr_mec = mean(tot_w3_x(mec_cells,:));
u_corr_mec = m_corr_mec + 1*std(tot_w3_x(mec_cells,:))/sqrt(length(mec_cells));
l_corr_mec = m_corr_mec - 1*std(tot_w3_x(mec_cells,:))/sqrt(length(mec_cells));

[pk_mec,pk_mec_loc] = max(m_corr_mec);
corr_at_peak_mec = tot_w3_x(mec_cells,pk_mec_loc);

pk_mec_u = u_corr_mec(pk_mec_loc);
pk_mec_l = l_corr_mec(pk_mec_loc);
pk_shift_mec = lags(pk_mec_loc)/Fsd;

m_corr_lec = mean(tot_w3_x(lec_cells,:));
u_corr_lec = m_corr_lec + 1*std(tot_w3_x(lec_cells,:))/sqrt(length(lec_cells));
l_corr_lec = m_corr_lec - 1*std(tot_w3_x(lec_cells,:))/sqrt(length(lec_cells));

[pk_lec,pk_lec_loc] = max(m_corr_lec);
pk_lec_u = u_corr_lec(pk_lec_loc);
pk_lec_l = l_corr_lec(pk_lec_loc);
pk_shift_lec = lags(pk_lec_loc)/Fsd;
corr_at_peak_lec = tot_w8_x(lec_cells,pk_lec_loc);

plot(lags/Fsd,m_corr_mec,'linewidth',2)
hold on
plot(lags/Fsd,m_corr_lec,'g','linewidth',2)

X = [lags/Fsd fliplr(lags)/Fsd];
Y = [l_corr_mec fliplr(u_corr_mec)];
fill(X,Y,'b')
Y = [l_corr_lec fliplr(u_corr_lec)];
fill(X,Y,'g')

line([0 0],[-0.4 1],'Color','k')
line([-8 8],[0 0],'Color','k')
xlim([-4 4])
ylim([-0.4 0.8])
ylabel('Correlation Coefficient')
xlabel('Time (s)')

figure
plot(peak_lags(mec_cells),all_peaks(mec_cells),'o')
hold on
plot(peak_lags(lec_cells),all_peaks(lec_cells),'go')
xlabel('Peak lags (s)')
ylabel('Peak correlation')

load ./simcort_timedomain_data

%% Heka amplitude distributions
clear all
% close all
load F:\WC_Germany\overall_EC\overall_EC_dir

used_data = [l3mec_p l3lec_p];
sess_data = sess_data(used_data);

mec_cells = 1:length(l3mec_p);
lec_cells = (mec_cells(end)+1):(mec_cells(end)+length(l3lec_p));
n_mec = length(mec_cells);
n_lec = length(lec_cells);

cd F:\WC_Germany\persistent_9_27_2010
load ./pa_heka_UDS_data.mat

figure
h = errorbar(heka_amp_range,nanmean(upstate_heka_dist_spksub_sig(mec_cells,:)),nanstd(upstate_heka_dist_spksub_sig(mec_cells,:))/sqrt(length(mec_cells)));
errorbar_tick(h,.001,'units');
hold on
h = errorbar(heka_amp_range,nanmean(downstate_heka_dist_spksub_sig(mec_cells,:)),nanstd(downstate_heka_dist_spksub_sig(mec_cells,:))/sqrt(length(mec_cells)),'r');
errorbar_tick(h,.001,'units');
xlim([-90 -20])

figure
h = errorbar(heka_amp_range,nanmean(upstate_heka_dist_spksub_sig(lec_cells,:)),nanstd(upstate_heka_dist_spksub_sig(lec_cells,:))/sqrt(length(lec_cells)));
errorbar_tick(h,.001,'units');
hold on
h = errorbar(heka_amp_range,nanmean(downstate_heka_dist_spksub_sig(lec_cells,:)),nanstd(downstate_heka_dist_spksub_sig(lec_cells,:))/sqrt(length(lec_cells)),'r');
errorbar_tick(h,.001,'units');
xlim([-90 -20])

figure
h = errorbar(heka_amp_range,nanmean(ov_heka_dist(mec_cells,:)),nanstd(ov_heka_dist(mec_cells,:))/sqrt(length(mec_cells)));
errorbar_tick(h,.001,'units');
xlim([-90 -20])
ylim([0 0.12])

figure
h = errorbar(heka_amp_range,nanmean(ov_heka_dist(lec_cells,:)),nanstd(ov_heka_dist(lec_cells,:))/sqrt(length(lec_cells)),'g');
errorbar_tick(h,.001,'units');
xlim([-90 -20])
ylim([0 0.12])

%% Nlx amplitude distributions
clear all
close all
load G:\WC_Germany\overall_EC\overall_EC_dir

used_data = [l3mec_p l3lec_p];
sess_data = sess_data(used_data);

mec_cells = 1:length(l3mec_p);
lec_cells = (mec_cells(end)+1):(mec_cells(end)+length(l3lec_p));
n_mec = length(mec_cells);
n_lec = length(lec_cells);

cd G:\WC_Germany\persistent_9_27_2010
load ./pa_heka_UDS_data.mat

m_dist_mec = nanmean(ov_nlx_dist(mec_cells,:));
u_dist_mec = m_dist_mec+nanstd(ov_nlx_dist(mec_cells,:))/sqrt(length(mec_cells));
l_dist_mec = m_dist_mec-nanstd(ov_nlx_dist(mec_cells,:))/sqrt(length(mec_cells));

m_dist_lec = nanmean(ov_nlx_dist(lec_cells,:));
u_dist_lec = m_dist_lec+nanstd(ov_nlx_dist(lec_cells,:))/sqrt(length(lec_cells));
l_dist_lec = m_dist_lec-nanstd(ov_nlx_dist(lec_cells,:))/sqrt(length(lec_cells));

m_dist_lf8 = nanmean(ov_nlx_dist8);
u_dist_lf8 = m_dist_lf8+nanstd(ov_nlx_dist)/sqrt(length(used_data));
l_dist_lf8 = m_dist_lf8-nanstd(ov_nlx_dist)/sqrt(length(used_data));

figure
plot(nlx_amp_range,m_dist_mec,'linewidth',2)
hold on
plot(nlx_amp_range,m_dist_lec,'g','linewidth',2)
plot(nlx_amp_range,m_dist_lf8,'r','linewidth',2)
xlabel('Amplitude (z)','FontSize',14)
ylabel('Probability Density','FontSize',14)
xlim([-3 3])
ylim([0 1.0])
legend('MEC','LEC','LF8')

begp = find(nlx_amp_range < -3,1,'last');
endp = find(nlx_amp_range > 3,1,'first');

X = [nlx_amp_range(begp:endp) fliplr(nlx_amp_range(begp:endp))];
Y = [u_dist_mec(begp:endp) fliplr(l_dist_mec(begp:endp))];
fill(X,Y,'b')
X = [nlx_amp_range(begp:endp) fliplr(nlx_amp_range(begp:endp))];
Y = [u_dist_lec(begp:endp) fliplr(l_dist_lec(begp:endp))];
fill(X,Y,'g')
X = [nlx_amp_range(begp:endp) fliplr(nlx_amp_range(begp:endp))];
Y = [u_dist_lf8(begp:endp) fliplr(l_dist_lf8(begp:endp))];
fill(X,Y,'r')
xlabel('Amplitude (z)')
ylabel('Probability')

% figure
% h = errorbar(nlx_amp_range,nanmean(upstate_nlx_dist(mec_cells,:)),nanstd(upstate_nlx_dist(mec_cells,:))/sqrt(length(mec_cells)));
% errorbar_tick(h,.001,'units');
% hold on
% h = errorbar(nlx_amp_range,nanmean(downstate_nlx_dist(mec_cells,:)),nanstd(downstate_nlx_dist(mec_cells,:))/sqrt(length(mec_cells)),'r');
% errorbar_tick(h,.001,'units');
% xlim([-4 4])
% ylim([0 1.4])
% 
% figure
% h = errorbar(nlx_amp_range,nanmean(upstate_nlx_dist(lec_cells,:)),nanstd(upstate_nlx_dist(lec_cells,:))/sqrt(length(lec_cells)));
% errorbar_tick(h,.001,'units');
% hold on
% h = errorbar(nlx_amp_range,nanmean(downstate_nlx_dist(lec_cells,:)),nanstd(downstate_nlx_dist(lec_cells,:))/sqrt(length(lec_cells)),'r');
% errorbar_tick(h,.001,'units');
% xlim([-4 4])
% ylim([0 1.4])
% 
% figure
% h = errorbar(nlx_amp_range,nanmean(upstate_nlx_dist8),nanstd(upstate_nlx_dist)/sqrt(length(sess_data)));
% errorbar_tick(h,.001,'units');
% hold on
% h = errorbar(nlx_amp_range,nanmean(downstate_nlx_dist8),nanstd(downstate_nlx_dist)/sqrt(length(sess_data)),'r');
% errorbar_tick(h,.001,'units');
% xlim([-4 4])
% ylim([0 1.4])
% 
% figure
% h = errorbar(nlx_amp_range,nanmean(ov_nlx_dist(mec_cells,:)),nanstd(ov_nlx_dist(mec_cells,:))/sqrt(length(mec_cells)));
% errorbar_tick(h,.001,'units');
% xlim([-3 3])
% ylim([0 0.9])
% 
% figure
% h = errorbar(nlx_amp_range,nanmean(ov_nlx_dist(lec_cells,:)),nanstd(ov_nlx_dist(lec_cells,:))/sqrt(length(lec_cells)));
% errorbar_tick(h,.001,'units');
% xlim([-3 3])
% ylim([0 0.9])
% 
% figure
% h = errorbar(nlx_amp_range,nanmean(ov_nlx_dist8),nanstd(ov_nlx_dist8)/sqrt(length(sess_data)));
% errorbar_tick(h,.001,'units');
% xlim([-3 3])
% ylim([0 0.9])

%% Firing Rate Analysis
clear all
close all
load G:\WC_Germany\overall_EC\overall_EC_dir

used_data = [l3mec_p l3lec_p];
sess_data = sess_data(used_data);

mec_cells = 1:length(l3mec_p);
lec_cells = (mec_cells(end)+1):(mec_cells(end)+length(l3lec_p));
n_mec = length(mec_cells);
n_lec = length(lec_cells);

cd G:\WC_Germany\persistent_9_27_2010\
load ./spike_rate_data

plot(mp_up_rate_lfup(mec_cells),mp_up_rate_lfdown(mec_cells),'.','markersize',14)
% hold on
% plot(mp_up_rate_lfup(lec_cells),mp_up_rate_lfup(lec_cells),'g.')
legend('MEC')
line([0 12],[0 12],'Color','k')
xlabel('Rate during LFP up state (Hz)','Fontsize',14)
ylabel('Rate during LFP down state (Hz)','Fontsize',14)
P = polyfit(mp_up_rate_lfup(mec_cells),mp_up_rate_lfdown(mec_cells),1);
x_ax = linspace(0,12,1000);
predict = polyval(P,x_ax);
hold on
plot(x_ax,predict,'r')

%% Overall Power Spectra
clear all
close all
load G:\WC_Germany\overall_EC\overall_EC_dir

used_data = [l3mec_p l3lec_p];
sess_data = sess_data(used_data);

mec_cells = 1:length(l3mec_p);
lec_cells = (mec_cells(end)+1):(mec_cells(end)+length(l3lec_p));
n_mec = length(mec_cells);
n_lec = length(lec_cells);

cd G:\WC_Germany\persistent_9_27_2010\
load ./spectral_data_3

max_fr = 1.;
max_f = find(f > max_fr,1,'first');

% convert to log power
% Sw = 10*log10(zPww);
% S8 = 10*log10(zP88);
% S3 = 10*log10(zP33);
Sw = 10*(zPww);
S8 = 10*(zP88);
S3 = 10*(zP33);
S3h = 10*(zP33h);
% Sw = Pww;
% S8 = P88;

mean_S8 = mean(S8([mec_cells lec_cells],:));
u_S8 = mean_S8 + 1*std(S8([mec_cells lec_cells],:))/sqrt(length(mec_cells)+length(lec_cells));
l_S8 = mean_S8 - 1*std(S8([mec_cells lec_cells],:))/sqrt(length(mec_cells)+length(lec_cells));

mean_S3 = mean(S3([mec_cells lec_cells],:));
u_S3 = mean_S3 + 1*std(S3([mec_cells lec_cells],:))/sqrt(length(mec_cells)+length(lec_cells));
l_S3 = mean_S3 - 1*std(S3([mec_cells lec_cells],:))/sqrt(length(mec_cells)+length(lec_cells));

mean_S3h = mean(S3h([mec_cells lec_cells],:));
u_S3h = mean_S3h + 1*std(S3h([mec_cells lec_cells],:))/sqrt(length(mec_cells)+length(lec_cells));
l_S3h = mean_S3h - 1*std(S3h([mec_cells lec_cells],:))/sqrt(length(mec_cells)+length(lec_cells));

mean_Sw_mec = mean(Sw(mec_cells,:));
u_Sw_mec = mean_Sw_mec + 1*std(Sw(mec_cells,:))/sqrt(length(mec_cells));
l_Sw_mec = mean_Sw_mec - 1*std(Sw(mec_cells,:))/sqrt(length(mec_cells));

mean_Sw_lec = mean(Sw(lec_cells,:));
u_Sw_lec = mean_Sw_lec + 1*std(Sw(lec_cells,:))/sqrt(length(lec_cells));
l_Sw_lec = mean_Sw_lec - 1*std(Sw(lec_cells,:))/sqrt(length(lec_cells));

figure
set(gca,'fontsize',14,'fontname','arial')
plot(f(1:max_f),mean_Sw_mec(1:max_f),'linewidth',2)
hold on
plot(f(1:max_f),mean_Sw_lec(1:max_f),'g','linewidth',2)
plot(f(1:max_f),mean_S8(1:max_f),'r','linewidth',2)
plot(f(1:max_f),mean_S3(1:max_f),'k','linewidth',2)
% plot(f(1:max_f),mean_S3h(1:max_f),'c','linewidth',2)

X = [f(1:max_f) fliplr(f(1:max_f))];
Y = [u_Sw_mec(1:max_f) fliplr(l_Sw_mec(1:max_f))];
fill(X,Y,'b')
X = [f(1:max_f) fliplr(f(1:max_f))];
Y = [u_Sw_lec(1:max_f) fliplr(l_Sw_lec(1:max_f))];
fill(X,Y,'g')
X = [f(1:max_f) fliplr(f(1:max_f))];
Y = [u_S8(1:max_f) fliplr(l_S8(1:max_f))];
fill(X,Y,'r')
X = [f(1:max_f) fliplr(f(1:max_f))];
Y = [u_S3(1:max_f) fliplr(l_S3(1:max_f))];
fill(X,Y,'k')
% X = [f(1:max_f) fliplr(f(1:max_f))];
% Y = [u_S3h(1:max_f) fliplr(l_S3h(1:max_f))];
% fill(X,Y,'c')

xlabel('Frequency (Hz)')
ylabel('Power (V^2/Hz)')
xlim([0.0 max_fr])
ylim([-12 2])
% ax1 = gca;
% xlim([0 max_fr])
% % set(ax1,'yscale','log')
% ax2 = axes('Position',get(ax1,'Position'),'YAxisLocation','right','Color','none','YColor','r');
% line(f(1:max_f),mean_S8(1:max_f),'color','r','linewidth',2,'parent',ax2)
% % X = [f(1:max_f) fliplr(f(1:max_f))];
% % Y = [u_S8(1:max_f) fliplr(l_S8(1:max_f))];
% % fill(X,Y,'r')
% % set(ax2,'yscale','log')
% xlim([0 max_fr])

uds_freqs = find(f > 0 & f < 1);
[peak_Sw,peak_Sw_loc] = max(Sw(:,uds_freqs),[],2);
[peak_S8,peak_S8_loc] = max(S8(:,uds_freqs),[],2);
[peak_S3,peak_S3_loc] = max(S3(:,uds_freqs),[],2);
peak_Sw_loc = f(uds_freqs(peak_Sw_loc));
peak_S8_loc = f(uds_freqs(peak_S8_loc));
peak_S3_loc = f(uds_freqs(peak_S3_loc));

%% Overall Coherence
clear all
close all
load G:\WC_Germany\overall_EC\overall_EC_dir

used_data = [l3mec_p l3lec_p];
sess_data = sess_data(used_data);

mec_cells = 1:length(l3mec_p);
lec_cells = (mec_cells(end)+1):(mec_cells(end)+length(l3lec_p));
n_mec = length(mec_cells);
n_lec = length(lec_cells);

cd G:\WC_Germany\persistent_9_27_2010\
% load ./spectral_data
% load ./pa_coherence_data_2.mat
load ./pa_coherence_data_lf2.mat

f = f_i;
max_fr = 1; %maximum frequency of interest in Hz
max_f = find(f > max_fr,1,'first');

%convert to atanh to stabilize the variance
% Cmn = atanh(Cmn);

mean_coh_mec = nanmean(aCw8_cor(mec_cells,:));
u_coh_mec = mean_coh_mec+1*nanstd(aCw8_cor(mec_cells,:))/sqrt(length(mec_cells));
l_coh_mec = mean_coh_mec-1*nanstd(aCw8_cor(mec_cells,:))/sqrt(length(mec_cells));
mean_coh_lec = nanmean(aCw8_cor(lec_cells,:));
u_coh_lec = mean_coh_lec+1*nanstd(aCw8_cor(lec_cells,:))/sqrt(length(lec_cells));
l_coh_lec = mean_coh_lec-1*nanstd(aCw8_cor(lec_cells,:))/sqrt(length(lec_cells));
% mean_coh_mec = nanmean(partial_w8(mec_cells,:));
% u_coh_mec = mean_coh_mec+1*nanstd(partial_w8(mec_cells,:))/sqrt(length(mec_cells));
% l_coh_mec = mean_coh_mec-1*nanstd(partial_w8(mec_cells,:))/sqrt(length(mec_cells));
% mean_coh_lec = nanmean(partial_w8(lec_cells,:));
% u_coh_lec = mean_coh_lec+1*nanstd(partial_w8(lec_cells,:))/sqrt(length(lec_cells));
% l_coh_lec = mean_coh_lec-1*nanstd(partial_w8(lec_cells,:))/sqrt(length(lec_cells));

partial_w3 = aCw3_cor;

mean_coh3_mec = nanmean(partial_w3(mec_cells,:));
u_coh3_mec = mean_coh3_mec+1*nanstd(partial_w3(mec_cells,:))/sqrt(length(mec_cells));
l_coh3_mec = mean_coh3_mec-1*nanstd(partial_w3(mec_cells,:))/sqrt(length(mec_cells));

mean_coh3_lec = nanmean(partial_w3(lec_cells,:));
u_coh3_lec = mean_coh3_lec+1*nanstd(partial_w3(lec_cells,:))/sqrt(length(lec_cells));
l_coh3_lec = mean_coh3_lec-1*nanstd(partial_w3(lec_cells,:))/sqrt(length(lec_cells));

mean_coh3h_mec = nanmean(partial_w3h(mec_cells,:));
u_coh3h_mec = mean_coh3_mec+1*nanstd(partial_w3h(mec_cells,:))/sqrt(length(mec_cells));
l_coh3h_mec = mean_coh3_mec-1*nanstd(partial_w3h(mec_cells,:))/sqrt(length(mec_cells));

mean_coh3h_lec = nanmean(partial_w3h(lec_cells,:));
u_coh3h_lec = mean_coh3_lec+1*nanstd(partial_w3h(lec_cells,:))/sqrt(length(lec_cells));
l_coh3h_lec = mean_coh3_lec-1*nanstd(partial_w3h(lec_cells,:))/sqrt(length(lec_cells));

%convert back to 0-1 interval
mean_coh_mec = tanh(mean_coh_mec);
mean_coh_lec = tanh(mean_coh_lec);
u_coh_mec = tanh(u_coh_mec);
u_coh_lec = tanh(u_coh_lec);
l_coh_mec = tanh(l_coh_mec);
l_coh_lec = tanh(l_coh_lec);

mean_coh3_mec = tanh(mean_coh3_mec);
mean_coh3_lec = tanh(mean_coh3_lec);
u_coh3_mec = tanh(u_coh3_mec);
u_coh3_lec = tanh(u_coh3_lec);
l_coh3_mec = tanh(l_coh3_mec);
l_coh3_lec = tanh(l_coh3_lec);

mean_coh3h_mec = tanh(mean_coh3h_mec);
mean_coh3h_lec = tanh(mean_coh3h_lec);
u_coh3h_mec = tanh(u_coh3h_mec);
u_coh3h_lec = tanh(u_coh3h_lec);
l_coh3h_mec = tanh(l_coh3h_mec);
l_coh3h_lec = tanh(l_coh3h_lec);

plot(f(1:max_f),mean_coh_mec(1:max_f),'linewidth',2)
hold on
plot(f(1:max_f),mean_coh_lec(1:max_f),'g','linewidth',2)
plot(f(1:max_f),mean_coh3_mec(1:max_f),'c','linewidth',2)
plot(f(1:max_f),mean_coh3_lec(1:max_f),'r','linewidth',2)
% plot(f(1:max_f),mean_coh3h_mec(1:max_f),'k','linewidth',2)
% plot(f(1:max_f),mean_coh3h_lec(1:max_f),'y','linewidth',2)
legend('MEC','LEC','Pre','Fro')

X = [f(1:max_f) fliplr(f(1:max_f))];
Y = [u_coh_mec(1:max_f) fliplr(l_coh_mec(1:max_f))];
fill(X,Y,'b')
X = [f(1:max_f) fliplr(f(1:max_f))];
Y = [u_coh3_lec(1:max_f) fliplr(l_coh3_lec(1:max_f))];
fill(X,Y,'g')
X = [f(1:max_f) fliplr(f(1:max_f))];
Y = [u_coh3_mec(1:max_f) fliplr(l_coh3_mec(1:max_f))];
fill(X,Y,'b')
X = [f(1:max_f) fliplr(f(1:max_f))];
Y = [u_coh_lec(1:max_f) fliplr(l_coh_lec(1:max_f))];
fill(X,Y,'g')
xlim([0 max_fr])
xlabel('Frequency (Hz)','FontSize',14)
ylabel('Coherence','Fontsize',14)
ylim([0. 1])

% load ./simcort_spectral_data.mat
% load ./pa_simcort_dir.mat
% Cmn = atanh(Cmn);
% mean_coh_pre = mean(Cmn(pre,:));
% u_coh_pre = mean_coh_pre+1*std(Cmn(pre,:))/sqrt(length(pre));
% l_coh_pre = mean_coh_pre-1*std(Cmn(pre,:))/sqrt(length(pre));
% mean_coh_fro = mean(Cmn(fro,:));
% u_coh_fro = mean_coh_fro+1*std(Cmn(fro,:))/sqrt(length(fro));
% l_coh_fro = mean_coh_fro-1*std(Cmn(fro,:))/sqrt(length(fro));
% %convert back to 0-1 interval
% mean_coh_pre = tanh(mean_coh_pre);
% mean_coh_fro = tanh(mean_coh_fro);
% u_coh_pre = tanh(u_coh_pre);
% u_coh_fro = tanh(u_coh_fro);
% l_coh_pre = tanh(l_coh_pre);
% l_coh_fro = tanh(l_coh_fro);
% 
% plot(f(1:max_f),mean_coh_pre(1:max_f),'c','linewidth',2)
% plot(f(1:max_f),mean_coh_fro(1:max_f),'k','linewidth',2)
% 
% % X = [f(1:max_f) fliplr(f(1:max_f))];
% % Y = [u_coh_pre(1:max_f) fliplr(l_coh_pre(1:max_f))];
% % fill(X,Y,'c')
% % X = [f(1:max_f) fliplr(f(1:max_f))];
% % Y = [u_coh_fro(1:max_f) fliplr(l_coh_fro(1:max_f))];
% % fill(X,Y,'k')
% 
% 
% figure
% plot(f(1:max_f),mean_coh3_mec(1:max_f),'linewidth',2)
% hold on
% plot(f(1:max_f),mean_coh3_lec(1:max_f),'g','linewidth',2)
% legend('MEC','LEC')
% X = [f(1:max_f) fliplr(f(1:max_f))];
% Y = [u_coh3_mec(1:max_f) fliplr(l_coh3_mec(1:max_f))];
% fill(X,Y,'b')
% X = [f(1:max_f) fliplr(f(1:max_f))];
% Y = [u_coh3_lec(1:max_f) fliplr(l_coh3_lec(1:max_f))];
% fill(X,Y,'g')
% xlim([0 max_fr])
% xlabel('Frequency (Hz)','FontSize',14)
% ylabel('Coherence','Fontsize',14)
% ylim([0 1])
% 
% % figure
% % h = errorbar(f,nanmean(lPww(mec_cells,:)),nanstd(lPww(mec_cells,:))/sqrt(22));
% % errorbar_tick(h,.001,'units')
% % hold on
% % h = errorbar(f,nanmean(lPww(lec_cells,:)),nanstd(lPww(lec_cells,:))/sqrt(14),'g');
% % errorbar_tick(h,.001,'units')
% % h = errorbar(f,nanmean(lP88),nanstd(lP88)/sqrt(14),'r');
% % errorbar_tick(h,.001,'units')
% % h = errorbar(f,nanmean(lP33),nanstd(lP33)/sqrt(36),'k');
% % errorbar_tick(h,.001,'units')
% % xlim([0 max_fr])
% 
uds_freqs = find(f > 0 & f < 1);
[peak_Cw8,peak_Cw8_loc] = nanmax(tanh(aCw8_cor(:,uds_freqs)),[],2);
[peak_Cw3,peak_Cw3_loc] = nanmax(tanh(aCw3_cor(:,uds_freqs)),[],2);
[peak_partial_Cw3,peak_partial_Cw3_loc] = nanmax(partial_w3(:,uds_freqs),[],2);
peak_Cw8_loc = f(uds_freqs(peak_Cw8_loc));
peak_Cw3_loc = f(uds_freqs(peak_Cw3_loc));
peak_partial_Cw3_loc = f(uds_freqs(peak_partial_Cw3_loc));

%% LF8 Triggered UP Trans Overall
clear all
close all
load G:\WC_Germany\overall_EC\overall_EC_dir

used_data = [l3mec_p l3lec_p];
sess_data = sess_data(used_data);

mec_cells = 1:length(l3mec_p);
lec_cells = (mec_cells(end)+1):(mec_cells(end)+length(l3lec_p));
n_mec = length(mec_cells);
n_lec = length(lec_cells);

cd G:\WC_Germany\persistent_9_27_2010\
load ./lf8_trig_avg_data_new2

first_look_point = find(lags > 0,1,'first');

lfp_avg_utrig = nanmean(lf8_utrig_lf8([mec_cells lec_cells],:));
lfp_uc_utrig = lfp_avg_utrig+nanstd(lf8_utrig_lf8([mec_cells lec_cells],:))/sqrt(length(sess_data));
lfp_lc_utrig = lfp_avg_utrig-nanstd(lf8_utrig_lf8([mec_cells lec_cells],:))/sqrt(length(sess_data));

mec_avg_utrig = nanmean(lf8_utrig_mp(mec_cells,:));
mec_uc_utrig = mec_avg_utrig+nanstd(lf8_utrig_mp(mec_cells,:))/sqrt(length(mec_cells));
mec_lc_utrig = mec_avg_utrig-nanstd(lf8_utrig_mp(mec_cells,:))/sqrt(length(mec_cells));

[mval,mloc] = max(mec_avg_utrig);
[minval,minloc] = min(mec_avg_utrig(1:mloc));
mid_pt = (mval+minval)/2;
mid_loc_mec = lags(first_look_point+find(mec_avg_utrig(first_look_point:end) > mid_pt,1,'first'));
mid_loc_mec_u = lags(first_look_point+find(mec_uc_utrig(first_look_point:end) > mid_pt,1,'first'));
mid_loc_mec_l = lags(first_look_point+find(mec_lc_utrig(first_look_point:end) > mid_pt,1,'first'));

lec_avg_utrig = nanmean(lf8_utrig_mp(lec_cells,:));
lec_uc_utrig = lec_avg_utrig+nanstd(lf8_utrig_mp(lec_cells,:))/sqrt(length(lec_cells));
lec_lc_utrig = lec_avg_utrig-nanstd(lf8_utrig_mp(lec_cells,:))/sqrt(length(lec_cells));

[mval,mloc] = max(lec_avg_utrig);
[minval,minloc] = min(lec_avg_utrig(1:mloc));
mid_pt = (mval+minval)/2;
mid_loc_lec = lags(first_look_point+find(lec_avg_utrig(first_look_point:end) > mid_pt,1,'first'));
mid_loc_lec_u = lags(first_look_point+find(lec_uc_utrig(first_look_point:end) > mid_pt,1,'first'));
mid_loc_lec_l = lags(first_look_point+find(lec_lc_utrig(first_look_point:end) > mid_pt,1,'first'));

figure
plot(lags,mec_avg_utrig,'linewidth',2)
hold on
plot(lags,lec_avg_utrig,'g','linewidth',2)
plot(lags,lfp_avg_utrig,'r','linewidth',2)
legend('MEC','LEC','LFP')

X = [lags fliplr(lags)];
Y = [mec_uc_utrig fliplr(mec_lc_utrig)];
fill(X,Y,'b')

Y = [lec_uc_utrig fliplr(lec_lc_utrig)];
fill(X,Y,'g')
Y = [lfp_uc_utrig fliplr(lfp_lc_utrig)];
fill(X,Y,'r')

xlim([-2 2])
xlabel('Time Lag (s)')
ylabel('Amplitude (z)')
line([0 0],[-1 1.5],'Color','k')
line([-2 2],[0 0],'Color','k')
ylim([-1 1.6])

figure
plot(lags,lf8_utrig_mp(mec_cells,:))
hold on
plot(lags,lfp_avg_utrig,'r','linewidth',2)
xlim([-2 2])
xlabel('Time Lag (s)')
ylabel('Amplitude (z)')
line([0 0],[-1 1.5],'Color','k')
line([-2 2],[0 0],'Color','k')
ylim([-1 1.6])

figure
plot(lags,lf8_utrig_mp(lec_cells,:))
hold on
plot(lags,lfp_avg_utrig,'r','linewidth',2)
xlim([-2 2])
xlabel('Time Lag (s)')
ylabel('Amplitude (z)')
line([0 0],[-1 1.5],'Color','k')
line([-2 2],[0 0],'Color','k')
ylim([-1 1.6])


%% LF8 Triggered DOWN Trans Overall
clear all
close all
load G:\WC_Germany\overall_EC\overall_EC_dir

used_data = [l3mec_p l3lec_p];
sess_data = sess_data(used_data);

mec_cells = 1:length(l3mec_p);
lec_cells = (mec_cells(end)+1):(mec_cells(end)+length(l3lec_p));
n_mec = length(mec_cells);
n_lec = length(lec_cells);

cd G:\WC_Germany\persistent_9_27_2010\
load ./lf8_trig_avg_data_new2

first_look_point = find(lags > 0,1,'first');

lfp_avg_dtrig = nanmean(lf8_dtrig_lf8([lec_cells mec_cells],:));
lfp_uc_dtrig = lfp_avg_dtrig+nanstd(lf8_dtrig_lf8([lec_cells mec_cells],:))/sqrt(length(sess_data));
lfp_lc_dtrig = lfp_avg_dtrig-nanstd(lf8_dtrig_lf8([lec_cells mec_cells],:))/sqrt(length(sess_data));

mec_avg_dtrig = nanmean(lf8_dtrig_mp(mec_cells,:));
mec_uc_dtrig = mec_avg_dtrig+nanstd(lf8_dtrig_mp(mec_cells,:))/sqrt(length(mec_cells));
mec_lc_dtrig = mec_avg_dtrig-nanstd(lf8_dtrig_mp(mec_cells,:))/sqrt(length(mec_cells));

[mval,mloc] = min(mec_avg_dtrig);
[maxval,maxloc] = max(mec_avg_dtrig(1:mloc));
mid_pt = (mval+maxval)/2;
mid_loc_mec = lags(first_look_point+find(mec_avg_dtrig(first_look_point:end) < mid_pt,1,'first'));
mid_loc_mec_u = lags(first_look_point+find(mec_uc_dtrig(first_look_point:end) < mid_pt,1,'first'));
mid_loc_mec_l = lags(first_look_point+find(mec_lc_dtrig(first_look_point:end) < mid_pt,1,'first'));

lec_avg_dtrig = nanmean(lf8_dtrig_mp(lec_cells,:));
lec_uc_dtrig = lec_avg_dtrig+nanstd(lf8_dtrig_mp(lec_cells,:))/sqrt(length(lec_cells));
lec_lc_dtrig = lec_avg_dtrig-nanstd(lf8_dtrig_mp(lec_cells,:))/sqrt(length(lec_cells));

[mval,mloc] = min(lec_avg_dtrig);
[maxval,maxloc] = max(lec_avg_dtrig(1:mloc));
mid_pt = (mval+maxval)/2;
mid_loc_lec = lags(first_look_point+find(lec_avg_dtrig(first_look_point:end) < mid_pt,1,'first'));
mid_loc_lec_u = lags(first_look_point+find(lec_uc_dtrig(first_look_point:end) < mid_pt,1,'first'));
mid_loc_lec_l = lags(first_look_point+find(lec_lc_dtrig(first_look_point:end) < mid_pt,1,'first'));

figure
plot(lags,mec_avg_dtrig,'linewidth',2)
hold on
plot(lags,lec_avg_dtrig,'g','linewidth',2)
plot(lags,lfp_avg_dtrig,'r','linewidth',2)
legend('MEC','LEC')

X = [lags fliplr(lags)];
Y = [mec_uc_dtrig fliplr(mec_lc_dtrig)];
fill(X,Y,'b')

Y = [lec_uc_dtrig fliplr(lec_lc_dtrig)];
fill(X,Y,'g')
Y = [lfp_uc_dtrig fliplr(lfp_lc_dtrig)];
fill(X,Y,'r')

ylim([-1 1.5])
xlim([-2 2])
xlabel('Time Lag (s)')
ylabel('Amplitude (z)')
line([0 0],[-1 2],'Color','k')
line([-2 2],[0 0],'Color','k')
ylim([-1 1.6])

figure
plot(lags,lf8_dtrig_mp(mec_cells,:))
hold on
plot(lags,lfp_avg_dtrig,'r','linewidth',2)
xlim([-2 2])
xlabel('Time Lag (s)')
ylabel('Amplitude (z)')
line([0 0],[-1 1.5],'Color','k')
line([-2 2],[0 0],'Color','k')
ylim([-1 1.6])

figure
plot(lags,lf8_dtrig_mp(lec_cells,:))
hold on
plot(lags,lfp_avg_dtrig,'r','linewidth',2)
xlim([-2 2])
xlabel('Time Lag (s)')
ylabel('Amplitude (z)')
line([0 0],[-1 1.5],'Color','k')
line([-2 2],[0 0],'Color','k')
ylim([-1 1.6])

%% WCV Up Trig Mat example
clear all
close all
load G:\WC_Germany\overall_EC\overall_EC_dir
addpath('G:\WC_Germany\hsmm_state_detection\')

used_data = [l3mec_p l3lec_p];
sess_data = sess_data(used_data);

mec_cells = 1:length(l3mec_p);
lec_cells = (mec_cells(end)+1):(mec_cells(end)+length(l3lec_p));
n_mec = length(mec_cells);
n_lec = length(lec_cells);

%
load G:\WC_Germany\persistent_9_27_2010\pa_trigmat_cell35_new2
d = 35;
trans_num_1 = 24;
begt_1 = 100.8;
endt_1 = 103.5;
trans_num_2 = 115;
begt_2 = 495.1;
endt_2 = 499.3;
trans_num_3 = 30;
begt_3 = 125.3;
endt_3 = 131.1;
trans_num_4 = 35;
begt_4 = 145.2;
endt_4 = 153.2;
% load G:\WC_Germany\persistent_9_27_2010\pa_trigmat_cell1
% d = 1;
% trans_num_1 = 5;
% begt_1 = 20.5;
% endt_1 = 23.5;
% % trans_num_2 = 44;
% % begt_2 = 151;
% % endt_2 = 157;
% trans_num_2 = 3;
% begt_2 = 12;
% endt_2 = 17;
% trans_num_3 = 60;
% begt_3 = 224;
% endt_3 = 230.6;
% % trans_num_3 = 45;
% % begt_3 = 156;
% % endt_3 = 165;
% trans_num_4 = 128;
% begt_4 = 505;
% endt_4 = 514;

dsf = 8;
Fsd = 2016/dsf;

cur_mp_mat = mp_utrig_mp_mat;
cur_lfp_mat = mp_utrig_lf8_mat;
cur_lf3_mat = mp_utrig_lf3_mat;

%get rid of the very first row of the matrix 
% bad_rows = 1;
bad_rows = [1 size(cur_mp_mat,1) size(cur_mp_mat,1)-1];
cur_mp_mat(bad_rows,:) = [];
cur_lfp_mat(bad_rows,:) = [];
cur_lf3_mat(bad_rows,:) = [];
mp_updur(bad_rows) = [];

[dummy,up_order] = sort(mp_updur);

Fig = figure(1);
clf
set(Fig,'PaperUnits','centimeters');
set(Fig, 'PaperSize', [30 20]);% paper size is in [width height] format
set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
imagesc(lags,(1:length(mp_updur)),(cur_mp_mat(up_order,:)));shading flat;
hold on
plot(mp_updur(up_order),(1:length(mp_updur)),'r','linewidth',2)
line([0 0],[0 1],'Color','k')
caxis([-3 3]);colorbar
xlim([-2 10])
xlabel('Time (s)','FontSize',14)
ylabel('MP Up Duration Percentile','FontSize',14)

Fig2 = figure(2);
clf
set(Fig2,'PaperUnits','centimeters');
set(Fig2, 'PaperSize', [30 20]);% paper size is in [width height] format
set(Fig2,'PaperPosition',[0,0,(get(Fig2,'PaperSize'))])
imagesc(lags,(1:length(mp_updur)),(cur_lfp_mat(up_order,:)));shading flat;
hold on
plot(mp_updur(up_order),(1:length(mp_updur)),'r','linewidth',2)
line([0 0],[0 1],'Color','k')
caxis([-1.5 3]);colorbar
xlim([-2 10])

Fig3 = figure(3);
clf
set(Fig3,'PaperUnits','centimeters');
set(Fig3, 'PaperSize', [30 20]);% paper size is in [width height] format
set(Fig3,'PaperPosition',[0,0,(get(Fig3,'PaperSize'))])
imagesc(lags,(1:length(mp_updur)),(cur_lf3_mat(up_order,:)));shading flat;
hold on
plot(mp_updur(up_order),(1:length(mp_updur)),'r','linewidth',2)
line([0 0],[0 1],'Color','k')
caxis([-1.5 2]);colorbar
xlim([-2 10])

% figure
% plot(mp_updur(up_order),temp,'r','linewidth',2)
% line([0 0],temp([1 end]),'color','k')
% xlim([-2 10])

% temp = log(5+(length(mp_updur):-1:1));
% figure
% pcolor(lags,temp,(cur_lfp_mat(up_order,:)));shading flat;
% hold on
% plot(mp_updur(up_order),temp,'r','linewidth',2)
% line([0 0],temp([1 end]),'Color','k')
% caxis([-1.8 3]);colorbar
% xlim([-2 10])
% 
% temp = log(5+(length(mp_updur):-1:1));
% figure
% pcolor(lags,temp,(cur_mp_mat(up_order,:)));shading flat;
% hold on
% plot(mp_updur(up_order),temp,'r','linewidth',2)
% line([0 0],temp([1 end]),'Color','k')
% caxis([-2.5 2.5]);colorbar
% xlim([-2 10])

figure(4)
plot(mp_updur(up_order),1-(1:length(mp_updur))/length(mp_updur),'r','linewidth',2)
line([0 0],[0 1],'Color','k')
xlim([-2 10])
xlabel('Time (s)','FontSize',14)
ylabel('MP Up Duration Percentile','FontSize',14)

niqf = 2016/2;
lcf = 0.05/niqf;
hcf = 4/niqf;
[b,a] = butter(2,[lcf hcf]);

dsf = 8;
Fsd = 2016/dsf;

cdir = sess_data(d).directory;
cdir(1) = 'G';
cd(cdir)
pwd

load ./used_data lf8 wcv_minus_spike lf3
load ./spike_time_jmm
spike_ids = round(spkid/dsf);
load ./pa_hsmm_state_seq_new2
load ./pa_hsmm_state_seq8_new2

mp_state_seq_c =  hsmm_bbstate_seq;
lf8_state_seq_c = hsmm_bbstate_seq8;

lf3_f = filtfilt(b,a,lf3);
lf8_f = filtfilt(b,a,lf8);
wcv_f = filtfilt(b,a,wcv_minus_spike);

lf3_d = downsample(lf3_f,dsf);
lf8_d = downsample(lf8_f,dsf);
wcv_d = downsample(wcv_f,dsf);

lf3_d = zscore(lf3_d);
lf8_d = zscore(lf8_d);
wcv_d = zscore(wcv_d);
t = (1:length(lf8_d))/Fsd;
[new_mp_seg_inds] = round(resample_uds_seg_inds(hmm.UDS_segs,50.4,Fsd,length(t)));
[new_lf8_seg_inds] = round(resample_uds_seg_inds(hmm8.UDS_segs,50.4,Fsd,length(t)));
[up_trans_inds,down_trans_inds] = compute_state_transitions_seg(new_mp_seg_inds,mp_state_seq_c);
[up_trans_inds8,down_trans_inds8] = compute_state_transitions_seg(new_lf8_seg_inds,lf8_state_seq_c);

up_trans_inds(bad_rows) = []; down_trans_inds(bad_rows) = [];

% at 16 percentile, 43rd up trans in rank order
prev_lfp_tran = find(t(up_trans_inds8) < t(up_trans_inds(trans_num_1)),1,'last');
begpt = find(t > begt_1,1,'first');
endpt = find(t > endt_1,1,'first');
figure
plot(t(begpt:endpt),wcv_d(begpt:endpt),'linewidth',3)
hold on
plot(t(begpt:endpt),lf8_d(begpt:endpt)-1,'r','linewidth',3)
plot(t(begpt:endpt),lf3_d(begpt:endpt)-1,'k','linewidth',3)
xlim([begt_1 endt_1])
ylim([-3 2.2])
% ylim([-2.5 2.2])
used_spikes = find(spike_ids/Fsd > begt_1 & spike_ids/Fsd < endt_1);
% line_bottom = 2.2;
% line_top = 2.6;
line_bottom = 1.4;
line_top = 1.8;
for i = 1:length(used_spikes)
   line([t(spike_ids(used_spikes(i))) t(spike_ids(used_spikes(i)))],...
       [line_bottom line_top],'Color','k','linewidth',2)
end
line([t(up_trans_inds(trans_num_1)) t(up_trans_inds(trans_num_1))],[0 1],'color','k')
line([t(down_trans_inds(trans_num_1)) t(down_trans_inds(trans_num_1))],[0 1],'color','k')
line([t(up_trans_inds8(prev_lfp_tran)) t(up_trans_inds8(prev_lfp_tran))],[0 1],'color','g')
line([t(down_trans_inds8(prev_lfp_tran)) t(down_trans_inds8(prev_lfp_tran))],[0 1],'color','g')


%at 68%, rank 147
prev_lfp_tran = find(up_trans_inds8 < up_trans_inds(trans_num_2),1,'last');
cur_lag = up_trans_inds(trans_num_2) - up_trans_inds8(prev_lfp_tran);
begpt = find(t > begt_2,1,'first');
endpt = find(t > endt_2,1,'first');
figure
plot(t(begpt:endpt),wcv_d(begpt:endpt),'linewidth',2)
hold on
plot(t(begpt:endpt),lf8_d(begpt:endpt)-1,'r','linewidth',3)
plot(t(begpt:endpt),lf3_d(begpt:endpt)-1,'k','linewidth',3)
xlim([begt_2 endt_2])
ylim([-3 2])
% ylim([-2.5 3])
used_spikes = find(spike_ids/Fsd > begt_2 & spike_ids/Fsd < endt_2);
% line_bottom = 2.2;
% line_top = 2.6;
line_bottom = 1.4;
line_top = 1.8;
for i = 1:length(used_spikes)
    line([t(spike_ids(used_spikes(i))) t(spike_ids(used_spikes(i)))],...
        [line_bottom line_top],'Color','k','linewidth',2)
end
line([t(up_trans_inds(trans_num_2)) t(up_trans_inds(trans_num_2))],[0 1],'color','k')
line([t(down_trans_inds(trans_num_2)) t(down_trans_inds(trans_num_2))],[0 1],'color','k')
line([t(up_trans_inds8(prev_lfp_tran)) t(up_trans_inds8(prev_lfp_tran))],[0 1],'color','g')
line([t(down_trans_inds8(prev_lfp_tran)) t(down_trans_inds8(prev_lfp_tran))],[0 1],'color','g')
line([t(up_trans_inds8(prev_lfp_tran+1)) t(up_trans_inds8(prev_lfp_tran+1))],[0 1],'color','g')
line([t(down_trans_inds8(prev_lfp_tran+1)) t(down_trans_inds8(prev_lfp_tran+1))],[0 1],'color','g')


%at 84%, rank 183
prev_lfp_tran = find(up_trans_inds8 < up_trans_inds(trans_num_3),1,'last');
cur_lag = up_trans_inds(trans_num_3) - up_trans_inds8(prev_lfp_tran);
figure
begpt = find(t > begt_3,1,'first');
endpt = find(t > endt_3,1,'first');
plot(t(begpt:endpt),wcv_d(begpt:endpt),'linewidth',2)
hold on
plot(t(begpt:endpt),lf8_d(begpt:endpt)-1,'r','linewidth',3)
plot(t(begpt:endpt),lf3_d(begpt:endpt)-1,'k','linewidth',3)
xlim([begt_3 endt_3])
% ylim([-3 3])
ylim([-3 2.2])
used_spikes = find(spike_ids/Fsd > begt_3 & spike_ids/Fsd < endt_3);
% line_bottom = 2.2;
% line_top = 2.6;
line_bottom = 1.5;
line_top = 1.9;
for i = 1:length(used_spikes)
    line([t(spike_ids(used_spikes(i))) t(spike_ids(used_spikes(i)))],...
        [line_bottom line_top],'Color','k','linewidth',2)
end
line([t(up_trans_inds(trans_num_3)) t(up_trans_inds(trans_num_3))],[0 1],'color','k')
line([t(down_trans_inds(trans_num_3)) t(down_trans_inds(trans_num_3))],[0 1],'color','k')
line([t(up_trans_inds8(prev_lfp_tran)) t(up_trans_inds8(prev_lfp_tran))],[0 1],'color','g')
line([t(down_trans_inds8(prev_lfp_tran)) t(down_trans_inds8(prev_lfp_tran))],[0 1],'color','g')
line([t(up_trans_inds8(prev_lfp_tran+1)) t(up_trans_inds8(prev_lfp_tran+1))],[0 1],'color','g')
line([t(down_trans_inds8(prev_lfp_tran+1)) t(down_trans_inds8(prev_lfp_tran+1))],[0 1],'color','g')
line([t(up_trans_inds8(prev_lfp_tran+2)) t(up_trans_inds8(prev_lfp_tran+2))],[0 1],'color','g')
line([t(down_trans_inds8(prev_lfp_tran+2)) t(down_trans_inds8(prev_lfp_tran+2))],[0 1],'color','g')


%at 95%
prev_lfp_tran = find(up_trans_inds8 < up_trans_inds(trans_num_4),1,'last');
cur_lag = up_trans_inds(trans_num_4) - up_trans_inds8(prev_lfp_tran);
begpt = find(t > begt_4,1,'first');
endpt = find(t > endt_4,1,'first');
figure
plot(t(begpt:endpt),wcv_d(begpt:endpt),'linewidth',2)
hold on
plot(t(begpt:endpt),lf8_d(begpt:endpt)-1,'r','linewidth',3)
plot(t(begpt:endpt),lf3_d(begpt:endpt)-1,'k','linewidth',3)
xlim([begt_4 endt_4])
ylim([-3.3 2.5])
% ylim([-2.6 1.8])
used_spikes = find(spike_ids/Fsd > begt_4 & spike_ids/Fsd < endt_4);
line_bottom = 1.8;
line_top = 2.2;
% line_bottom = 1.1;
% line_top = 1.5;
for i = 1:length(used_spikes)
   line([t(spike_ids(used_spikes(i))) t(spike_ids(used_spikes(i)))],...
       [line_bottom line_top],'Color','k','linewidth',2)
end
line([t(up_trans_inds(trans_num_4)) t(up_trans_inds(trans_num_4))],[0 1],'color','k')
line([t(down_trans_inds(trans_num_4)) t(down_trans_inds(trans_num_4))],[0 1],'color','k')
line([t(up_trans_inds8(prev_lfp_tran)) t(up_trans_inds8(prev_lfp_tran))],[0 1],'color','g')
line([t(down_trans_inds8(prev_lfp_tran)) t(down_trans_inds8(prev_lfp_tran))],[0 1],'color','g')
line([t(up_trans_inds8(prev_lfp_tran+1)) t(up_trans_inds8(prev_lfp_tran+1))],[0 1],'color','g')
line([t(down_trans_inds8(prev_lfp_tran+1)) t(down_trans_inds8(prev_lfp_tran+1))],[0 1],'color','g')
line([t(up_trans_inds8(prev_lfp_tran+2)) t(up_trans_inds8(prev_lfp_tran+2))],[0 1],'color','g')
line([t(down_trans_inds8(prev_lfp_tran+2)) t(down_trans_inds8(prev_lfp_tran+2))],[0 1],'color','g')
line([t(up_trans_inds8(prev_lfp_tran+3)) t(up_trans_inds8(prev_lfp_tran+3))],[0 1],'color','g')
line([t(down_trans_inds8(prev_lfp_tran+3)) t(down_trans_inds8(prev_lfp_tran+3))],[0 1],'color','g')

%% Super persistent example 1
clear all
close all
load G:\WC_Germany\overall_EC\overall_EC_dir
addpath('G:\WC_Germany\hsmm_state_detection\')

used_data = [l3mec_p l3lec_p];
sess_data = sess_data(used_data);

mec_cells = 1:length(l3mec_p);
lec_cells = (mec_cells(end)+1):(mec_cells(end)+length(l3lec_p));
n_mec = length(mec_cells);
n_lec = length(lec_cells);

d = 1;
niqf = 2016/2;
lcf = 0.05/niqf;
hcf = 10/niqf;
[b,a] = butter(2,[lcf hcf]);
dsf = 8;
Fsd = 2016/dsf;

cdir = sess_data(d).directory;
cdir(1) = 'G';
cd(cdir)
pwd

load ./used_data lf8 wcv_minus_spike
load ./spike_time_jmm
spike_ids = round(spkid/dsf);


lf8_f = filtfilt(b,a,lf8);
wcv_f = filtfilt(b,a,wcv_minus_spike);
lf8_d = downsample(lf8_f,dsf);
wcv_d = downsample(wcv_f,dsf);
lf8_d = zscore(lf8_d);
wcv_d = zscore(wcv_d);

load ./pa_hsmm_state_seq_sm400.mat
load ./pa_hsmm_state_seq8_sm400.mat
mp_state_seq_c =  hsmm_bbstate_seq;
lf8_state_seq_c = hsmm_bbstate_seq8;
% mp_smoothed_seq = thresh_state_smooth_seg(mp_state_seq_c,Fsd,500,500); %smooth out excessively short state durations to improve model fits
% lf8_smoothed_seq = thresh_state_smooth_seg(lf8_state_seq_c,Fsd,500,500); %smooth out excessively short state durations to improve model fits
t = (1:length(lf8_d))/Fsd;
[new_mp_seg_inds] = round(resample_uds_seg_inds(hmm.UDS_segs,50.4,Fsd,length(t)));
[new_lf8_seg_inds] = round(resample_uds_seg_inds(hmm8.UDS_segs,50.4,Fsd,length(t)));
mp_state_vec = nan(size(t));
lf8_state_vec = nan(size(t));
% mp_state_sm_vec = nan(size(t));
% lf8_state_sm_vec = nan(size(t));
for i = 1:hmm.Nsegs
   mp_state_vec(new_mp_seg_inds(i,1):new_mp_seg_inds(i,2)) = mp_state_seq_c{i}; 
   lf8_state_vec(new_lf8_seg_inds(i,1):new_lf8_seg_inds(i,2)) = lf8_state_seq_c{i}; 
%    mp_state_sm_vec(new_mp_seg_inds(i,1):new_mp_seg_inds(i,2)) = mp_smoothed_seq{i}; 
%    lf8_state_sm_vec(new_lf8_seg_inds(i,1):new_lf8_seg_inds(i,2)) = lf8_smoothed_seq{i}; 
end

% xl1 = 463;
% xl2 = 520;

figure
plot(t,wcv_d)
hold on
plot(t,lf8_d+3,'r')
plot(t,mp_state_vec-1.1,'k')
plot(t,lf8_state_vec+2.1,'k')
ot = (1:hmm.UDS_segs(1,2))/50.4;
% plot(ot,hsmm.state(1).meanfun{1},'k')
% plot(ot,hsmm.state(2).meanfun{1},'k')

load ./pa_hsmm_state_seq8_new2.mat
load ./pa_hsmm_state_seq_new2.mat

mp_state_seq_c =  hsmm_bbstate_seq;
lf8_state_seq_c = hsmm_bbstate_seq8;
t = (1:length(lf8_d))/Fsd;
[new_mp_seg_inds] = round(resample_uds_seg_inds(hmm.UDS_segs,50.4,Fsd,length(t)));
[new_lf8_seg_inds] = round(resample_uds_seg_inds(hmm8.UDS_segs,50.4,Fsd,length(t)));
mp_state_vec = nan(size(t));
lf8_state_vec = nan(size(t));
for i = 1:hmm.Nsegs
   mp_state_vec(new_mp_seg_inds(i,1):new_mp_seg_inds(i,2)) = mp_state_seq_c{i}; 
   lf8_state_vec(new_lf8_seg_inds(i,1):new_lf8_seg_inds(i,2)) = lf8_state_seq_c{i}; 
end
plot(t,mp_state_vec-1.2,'g')
plot(t,lf8_state_vec+2.2,'g')

% load ./pa_hsmm_state_seq_2hz_80.mat
% % % load ./pa_hsmm_state_seq8_4hz_30.mat
% mp_state_seq_c =  hsmm_bbstate_seq;
% lf8_state_seq_c = hsmm_bbstate_seq8;
% 
% t = (1:length(lf8_d))/Fsd;
% [new_mp_seg_inds] = round(resample_uds_seg_inds(hmm.UDS_segs,50.4,Fsd,length(t)));
% [new_lf8_seg_inds] = round(resample_uds_seg_inds(hmm8.UDS_segs,50.4,Fsd,length(t)));
% mp_state_vec2 = nan(size(t));
% lf8_state_vec2 = nan(size(t));
% for i = 1:hmm.Nsegs
%    mp_state_vec2(new_mp_seg_inds(i,1):new_mp_seg_inds(i,2)) = mp_state_seq_c{i}; 
%    lf8_state_vec2(new_lf8_seg_inds(i,1):new_lf8_seg_inds(i,2)) = lf8_state_seq_c{i}; 
% end
% % 
% plot(t,mp_state_vec2-2.2,'g')
% plot(ot,hsmm.state(1).meanfun{1},'r')
% plot(ot,hsmm.state(2).meanfun{1},'r')
% 
% % plot(t,lf8_state_vec2+3.2,'b')
% % xlim([xl1 xl2])