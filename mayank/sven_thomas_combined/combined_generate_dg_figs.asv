%% UP STATE DURATION DIST
clear all
close all
cd C:\WC_Germany\sven_thomas_combined\
load ./combined_core_analysis_dg

%use logarithmic hists
up_hist = up_dur_loghist;
up_hist8 = up_dur_loghist8;

%normalize state duration distributions 
up_hist = up_hist./repmat(sum(up_hist,2),1,size(up_hist,2));
up_hist8 = up_hist8./repmat(sum(up_hist8,2),1,size(up_hist8,2));
 
mean_up_dg = mean(up_hist);
u_up_dg = mean_up_dg + 1*std(up_hist)/sqrt(10);
l_up_dg = mean_up_dg - 1*std(up_hist)/sqrt(10);

mean_up_lfp = mean(up_hist8);
u_up_lfp = mean_up_lfp + 1*std(up_hist8)/sqrt(10);
l_up_lfp = mean_up_lfp - 1*std(up_hist8)/sqrt(10);

load ./combined_core_analysis_fin_nd_np.mat
load ./combined_dir_nd.mat
uset = sort([l3mec l3lec]);
all_cells = 1:61;
l3mec = find(ismember(all_cells(uset),l3mec));
combined_dir = combined_dir(uset);
up_hist = up_dur_loghist;
up_hist8 = up_dur_loghist8;
%normalize state duration distributions 
up_hist = up_hist./repmat(sum(up_hist,2),1,size(up_hist,2));
up_hist8 = up_hist8./repmat(sum(up_hist8,2),1,size(up_hist8,2));
mean_up_mec = mean(up_hist(l3mec,:));
u_up_mec = mean_up_mec + 1*std(up_hist(l3mec,:))/sqrt(length(l3mec));
l_up_mec = mean_up_mec - 1*std(up_hist(l3mec,:))/sqrt(length(l3mec));

figure
stairs(log_dur_grid,mean_up_lfp,'color',[0.2 0.2 0.2],'linewidth',2)
hold on
stairs(log_dur_grid,mean_up_dg,'b','linewidth',2)
stairs(log_dur_grid,mean_up_mec,'r','linewidth',2)

epsilon = 1e-7;
u_up_lfp(u_up_lfp < epsilon) = epsilon;
l_up_lfp(l_up_lfp < epsilon) = epsilon;
mean_up_lfp(mean_up_lfp < epsilon) = epsilon;
u_up_dg(u_up_dg < epsilon) = epsilon;
l_up_dg(l_up_dg < epsilon) = epsilon;
mean_up_dg(mean_up_dg < epsilon) = epsilon;
u_up_mec(u_up_mec < epsilon) = epsilon;
l_up_mec(l_up_mec < epsilon) = epsilon;
mean_up_mec(mean_up_mec < epsilon) = epsilon;

for i = 1:length(log_dur_grid)-1
    clear X Y
    bin_width = (log_dur_grid(i+1)-log_dur_grid(i));
    X = [log_dur_grid(i) log_dur_grid(i)+bin_width];
    Y = [u_up_lfp(i) u_up_lfp(i)];
    X = [X fliplr(X)];
    Y = [Y l_up_lfp(i) l_up_lfp(i)];
    fill(X,Y,[0.2 0.2 0.2],'EdgeColor','none');
    hold on
end

for i = 1:length(log_dur_grid)-1
    clear X Y
    bin_width = (log_dur_grid(i+1)-log_dur_grid(i));
    X = [log_dur_grid(i) log_dur_grid(i)+bin_width];
    Y = [u_up_dg(i) u_up_dg(i)];
    X = [X fliplr(X)];
    Y = [Y l_up_dg(i) l_up_dg(i)];
    fill(X,Y,'b','EdgeColor','none');
    hold on
end

for i = 1:length(log_dur_grid)-1
    clear X Y
    bin_width = (log_dur_grid(i+1)-log_dur_grid(i));
    X = [log_dur_grid(i) log_dur_grid(i)+bin_width];
    Y = [u_up_mec(i) u_up_mec(i)];
    X = [X fliplr(X)];
    Y = [Y l_up_mec(i) l_up_mec(i)];
    fill(X,Y,'r','EdgeColor','none');
    hold on
end

set(gca,'fontname','arial')
set(gca,'fontsize',14)
set(gca,'yscale','log')
set(gca,'YMinorGrid','off')
set(gca,'YGrid','off')
% set(gca,'XMinorGrid','off')
ylim([1e-4 0.15])
xlim([0.05 10])
shg
% set(gca,'xscale','log')
xlabel('Up State Duration (s)','FontSize',16)
ylabel('Percent of Data','FontSize',16)

%% DOWN STATE DURATION DIST
clear all
close all
cd C:\WC_Germany\sven_thomas_combined\
load ./combined_core_analysis_dg

%use logarithmic hists
down_hist = down_dur_loghist;
down_hist8 = down_dur_loghist8;

%normalize state duration distributions 
down_hist = down_hist./repmat(sum(down_hist,2),1,size(down_hist,2));
down_hist8 = down_hist8./repmat(sum(down_hist8,2),1,size(down_hist8,2));
 
mean_down_dg = mean(down_hist);
u_down_dg = mean_down_dg + 1*std(down_hist)/sqrt(10);
l_down_dg = mean_down_dg - 1*std(down_hist)/sqrt(10);

mean_down_lfp = mean(down_hist8);
u_down_lfp = mean_down_lfp + 1*std(down_hist8)/sqrt(10);
l_down_lfp = mean_down_lfp - 1*std(down_hist8)/sqrt(10);

load ./combined_core_analysis_fin_nd_np.mat
load ./combined_dir_nd.mat
uset = sort([l3mec l3lec]);
all_cells = 1:61;
l3mec = find(ismember(all_cells(uset),l3mec));
combined_dir = combined_dir(uset);
down_hist = down_dur_loghist;
down_hist8 = down_dur_loghist8;
%normalize state duration distributions 
down_hist = down_hist./repmat(sum(down_hist,2),1,size(down_hist,2));
down_hist8 = down_hist8./repmat(sum(down_hist8,2),1,size(down_hist8,2));
mean_down_mec = mean(down_hist(l3mec,:));
u_down_mec = mean_down_mec + 1*std(down_hist(l3mec,:))/sqrt(length(l3mec));
l_down_mec = mean_down_mec - 1*std(down_hist(l3mec,:))/sqrt(length(l3mec));

figure
stairs(log_dur_grid,mean_down_lfp,'color',[0.2 0.2 0.2],'linewidth',2)
hold on
stairs(log_dur_grid,mean_down_dg,'b','linewidth',2)
stairs(log_dur_grid,mean_down_mec,'r','linewidth',2)

epsilon = 1e-7;
u_down_lfp(u_down_lfp < epsilon) = epsilon;
l_down_lfp(l_down_lfp < epsilon) = epsilon;
mean_down_lfp(mean_down_lfp < epsilon) = epsilon;
u_down_dg(u_down_dg < epsilon) = epsilon;
l_down_dg(l_down_dg < epsilon) = epsilon;
mean_down_dg(mean_down_dg < epsilon) = epsilon;
u_down_mec(u_down_mec < epsilon) = epsilon;
l_down_mec(l_down_mec < epsilon) = epsilon;
mean_down_mec(mean_down_mec < epsilon) = epsilon;

for i = 1:length(log_dur_grid)-1
    clear X Y
    bin_width = (log_dur_grid(i+1)-log_dur_grid(i));
    X = [log_dur_grid(i) log_dur_grid(i)+bin_width];
    Y = [u_down_lfp(i) u_down_lfp(i)];
    X = [X fliplr(X)];
    Y = [Y l_down_lfp(i) l_down_lfp(i)];
    fill(X,Y,[0.2 0.2 0.2],'EdgeColor','none');
    hold on
end

for i = 1:length(log_dur_grid)-1
    clear X Y
    bin_width = (log_dur_grid(i+1)-log_dur_grid(i));
    X = [log_dur_grid(i) log_dur_grid(i)+bin_width];
    Y = [u_down_dg(i) u_down_dg(i)];
    X = [X fliplr(X)];
    Y = [Y l_down_dg(i) l_down_dg(i)];
    fill(X,Y,'b','EdgeColor','none');
    hold on
end

for i = 1:length(log_dur_grid)-1
    clear X Y
    bin_width = (log_dur_grid(i+1)-log_dur_grid(i));
    X = [log_dur_grid(i) log_dur_grid(i)+bin_width];
    Y = [u_down_mec(i) u_down_mec(i)];
    X = [X fliplr(X)];
    Y = [Y l_down_mec(i) l_down_mec(i)];
    fill(X,Y,'r','EdgeColor','none');
    hold on
end

set(gca,'fontname','arial')
set(gca,'fontsize',14)
set(gca,'yscale','log')
set(gca,'YMinorGrid','off')
set(gca,'YGrid','off')
% set(gca,'XMinorGrid','off')
ylim([1e-4 0.15])
xlim([0.05 10])
shg
% set(gca,'xscale','log')
xlabel('down State Duration (s)','FontSize',16)
ylabel('Percent of Data','FontSize',16)

%% UP transition lag
clear all
close all
cd C:\WC_Germany\sven_thomas_combined\
load ./combined_core_analysis_dg

lag_range = linspace(-0.5,1.5,1000);
dlag = lag_range(2)-lag_range(1);
for d = 1:10
    lag_hist(d,:) = histc(mp_uplags{d}/Fsd,lag_range);
    lag_hist(d,:) = lag_hist(d,:)/sum(lag_hist(d,:))/dlag;
    dg_lag_hist(d,:) = jmm_smooth_1d_cor(lag_hist(d,:),20);
end
load ./combined_core_analysis_fin_nd_np.mat
load ./combined_dir_nd.mat
uset = sort([l3mec l3lec]);
all_cells = 1:61;
l3mec = find(ismember(all_cells(uset),l3mec));
combined_dir = combined_dir(uset);
for d = 1:length(combined_dir)
    lag_hist(d,:) = histc(mp_uplags{d}/Fsd,lag_range);
    lag_hist(d,:) = lag_hist(d,:)/sum(lag_hist(d,:))/dlag;
    sm_lag_hist(d,:) = jmm_smooth_1d_cor(lag_hist(d,:),20);
end

figure
set(gca,'fontsize',14,'fontname','arial')
plot(lag_range,mean(dg_lag_hist),'b','linewidth',2)
hold on
plot(lag_range,mean(sm_lag_hist(l3mec,:)),'r','linewidth',2)
shadedErrorBar(lag_range,mean(dg_lag_hist),std(dg_lag_hist)./sqrt(10),{'b'});
shadedErrorBar(lag_range,mean(sm_lag_hist(l3mec,:)),std(sm_lag_hist(l3mec,:))./sqrt(length(l3mec)),{'r'});

xlim([-0.4 1.2])
ylim([0 3.2])
line([0 0],[0 3.2],'Color','k')
xlabel('Lag (s)','fontsize',16,'fontname','arial')
ylabel('Probability','fontsize',16,'fontname','arial')

%% DOWN transition lag
clear all
close all
cd C:\WC_Germany\sven_thomas_combined\
load ./combined_core_analysis_dg

lag_range = linspace(-1,2.3,1000);
dlag = lag_range(2)-lag_range(1);
for d = 1:10
    lag_hist(d,:) = histc(mp_downlags{d}/Fsd,lag_range);
    lag_hist(d,:) = lag_hist(d,:)/sum(lag_hist(d,:))/dlag;
    dg_lag_hist(d,:) = jmm_smooth_1d_cor(lag_hist(d,:),20);
end
load ./combined_core_analysis_fin_nd_np.mat
load ./combined_dir_nd.mat
uset = sort([l3mec l3lec]);
all_cells = 1:61;
l3mec = find(ismember(all_cells(uset),l3mec));
combined_dir = combined_dir(uset);
for d = 1:length(combined_dir)
    lag_hist(d,:) = histc(mp_downlags{d}/Fsd,lag_range);
    lag_hist(d,:) = lag_hist(d,:)/sum(lag_hist(d,:))/dlag;
    sm_lag_hist(d,:) = jmm_smooth_1d_cor(lag_hist(d,:),20);
end

figure
set(gca,'fontsize',14,'fontname','arial')
plot(lag_range,mean(dg_lag_hist),'b','linewidth',2)
hold on
plot(lag_range,mean(sm_lag_hist(l3mec,:)),'r','linewidth',2)
shadedErrorBar(lag_range,mean(dg_lag_hist),std(dg_lag_hist)./sqrt(10),{'b'});
shadedErrorBar(lag_range,mean(sm_lag_hist(l3mec,:)),std(sm_lag_hist(l3mec,:))./sqrt(length(l3mec)),{'r'});

xlim([-0.6 1.5])
ylim([0 2.2])
line([0 0],[0 3],'Color','k')
xlabel('Lag (s)','fontsize',16,'fontname','arial')
ylabel('Probability','fontsize',16,'fontname','arial')

%% Histogram of Type 2 Persistent Activity
cd C:\WC_Germany\sven_thomas_combined\
load ./combined_dir_nd.mat
load ./combined_core_analysis_fin_nd_np.mat

numBins = 10;
pers_range = linspace(0,0.40,numBins+1);

pers_mec = histc(fract_rt2_ups(l3mec),pers_range);
pers_lec = histc(fract_rt2_ups(l3lec),pers_range);

load ./combined_core_analysis_dg
pers_dg = histc(fract_rt2_ups,pers_range);

figure
set(gca,'fontname','arial','fontsize',14)
stairs(pers_range*100,pers_mec,'r','linewidth',2)
hold on
stairs(pers_range*100,pers_dg,'b','linewidth',2)
xlim([0 38])
% ylim([0 5])
xlabel('Fraction of Up States Persisting','FontSize',16,'fontname','arial')
ylabel('Number of Cells','FontSize',16,'fontname','arial')
legend('MEC','DG')
ylim([0 13])


