clear all
close all
cd C:\WC_Germany\sven_thomas_combined\
load ./combined_dir_nd.mat
load ./combined_core_analysis_fin_nd.mat
uset = sort([l3mec l3lec]);
all_cells = 1:61;
l3mec = find(ismember(all_cells(uset),l3mec));
l3lec = find(ismember(all_cells(uset),l3lec));
combined_dir = combined_dir(uset);
hpc_mua = hpc_mua(uset);
hpc_lfp = hpc_lfp(uset);
ctx_lfp = ctx_lfp(uset);


%% sorted LFP up-trig avgs
load combined_trig_mua
% [~,mec_ord] = sort(mean_updur(l3mec));
% mec_ord = l3mec(mec_ord);
all_cells = [l3mec l3lec];
l3lec(l3lec==33) = [];
cent = find(lags == 0);
for i = 1:length(all_cells)
     temp = find(lf8_utrig_wcv(i,cent:end) > 0.,1,'first')+cent-1;
     ep(i) = find(lf8_utrig_wcv(i,temp:end) < 0,1,'first') + temp-1;
end
[~,mec_ord] = sort(ep(l3mec));
mec_ord = l3mec(mec_ord);
 [~,lec_ord] = sort(ep(l3lec));
lec_ord = l3lec(lec_ord);
   
full_ord = [mec_ord lec_ord];
figure; set(gca,'fontname','arial','fontsize',14)
imagesc(lags,1:length(full_ord),lf8_utrig_wcv(full_ord,:))
xl = xlim(); yl = ylim();
line(xl,[length(l3mec) length(l3mec)]+0.5,'color','w')
line([0 0],yl,'color','k')
xlabel('Time lag (s)','fontsize',14)
ylabel('Cell Number','fontsize',14)

for i = 1:length(all_cells)
     temp = find(lf8_dtrig_wcv(i,1:cent) < 0.,1,'last');
     ep(i) = find(lf8_dtrig_wcv(i,temp+1:end) < 0,1,'first') + temp;
end
[~,mec_ord] = sort(ep(l3mec));
mec_ord = l3mec(mec_ord);
 [~,lec_ord] = sort(ep(l3lec));
lec_ord = l3lec(lec_ord);

full_ord = [mec_ord lec_ord];
figure; set(gca,'fontname','arial','fontsize',14)
imagesc(lags,1:length(full_ord),lf8_dtrig_wcv(full_ord,:))
xl = xlim(); yl = ylim();
line(xl,[length(l3mec) length(l3mec)]+0.5,'color','w')
line([0 0],yl,'color','k')
xlabel('Time lag (s)','fontsize',14)
ylabel('Cell Number','fontsize',14)

% figure; set(gca,'fontname','arial','fontsize',14)
% imagesc(lags,1:length(mec_ord),lf8_utrig_wcv(mec_ord,:))
% 
% figure; set(gca,'fontname','arial','fontsize',14)
% imagesc(lags,1:length(lec_ord),lf8_utrig_wcv(lec_ord,:))

% [~,lf8_ord] = sort(median_updur8);
% imagesc(lags,1:length(lf8_ord),lf8_utrig_lf8(lf8_ord,:))
%% Up State Distribution Compare
% l3mec(l3mec > 22) = [];

%use logarithmic hists
up_hist = up_dur_loghist;
up_hist8 = up_dur_loghist8;
% up_hist = up_dur_hist;
% up_hist8 = up_dur_hist8;

%normalize state duration distributions 
up_hist = up_hist./repmat(sum(up_hist,2),1,size(up_hist,2));
up_hist8 = up_hist8./repmat(sum(up_hist8,2),1,size(up_hist8,2));

% up_hist = cumsum(up_hist,2);
% up_hist8 = cumsum(up_hist8,2);
 
mean_up_mec = mean(up_hist(l3mec,:));
u_up_mec = mean_up_mec + 1*std(up_hist(l3mec,:))/sqrt(length(l3mec));
l_up_mec = mean_up_mec - 1*std(up_hist(l3mec,:))/sqrt(length(l3mec));

mean_up_lfp = mean(up_hist8([l3mec l3lec],:));
u_up_lfp = mean_up_lfp + 1*std(up_hist8([l3mec l3lec],:))/sqrt(length(l3mec)+length(l3lec));
l_up_lfp = mean_up_lfp - 1*std(up_hist8([l3mec l3lec],:))/sqrt(length(l3mec)+length(l3lec));

mean_up_lec = mean(up_hist(l3lec,:));
u_up_lec = mean_up_lec + 1*std(up_hist(l3lec,:))/sqrt(length(l3lec));
l_up_lec = mean_up_lec - 1*std(up_hist(l3lec,:))/sqrt(length(l3lec));

figure
stairs(log_dur_grid,mean_up_lfp,'color',[0.2 0.2 0.2],'linewidth',2)
hold on
stairs(log_dur_grid,mean_up_mec,'r','linewidth',2)
stairs(log_dur_grid,mean_up_lec,'b','linewidth',2)

legend('Ctx LFP','MEC','LEC')

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

for i = 1:length(log_dur_grid)-1
    clear X Y
    bin_width = (log_dur_grid(i+1)-log_dur_grid(i));
    X = [log_dur_grid(i) log_dur_grid(i)+bin_width];
    Y = [u_up_lec(i) u_up_lec(i)];
    X = [X fliplr(X)];
    Y = [Y l_up_lec(i) l_up_lec(i)];
    fill(X,Y,'b','EdgeColor','none');
    hold on
end

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


%% DOWN State Distribution Compare
% l3mec(l3mec > 22) = [];

%use logarithmic hists
down_hist = down_dur_loghist;
down_hist8 = down_dur_loghist8;

%normalize state duration distributions 
down_hist = down_hist./repmat(sum(down_hist,2),1,size(down_hist,2));
down_hist8 = down_hist8./repmat(sum(down_hist8,2),1,size(down_hist8,2));

down_hist = cumsum(down_hist,2);
down_hist8 = cumsum(down_hist8,2);
%  
% down_hist8 = down_hist8*100;

mean_down_mec = mean(down_hist(l3mec,:));
u_down_mec = mean_down_mec + 1*std(down_hist(l3mec,:))/sqrt(length(l3mec));
l_down_mec = mean_down_mec - 1*std(down_hist(l3mec,:))/sqrt(length(l3mec));

mean_down_lfp = mean(down_hist8([l3mec l3lec],:));
u_down_lfp = mean_down_lfp + 1*std(down_hist8([l3mec l3lec],:))/sqrt(length(l3mec)+length(l3lec));
l_down_lfp = mean_down_lfp - 1*std(down_hist8([l3mec l3lec],:))/sqrt(length(l3mec)+length(l3lec));

mean_down_lec = mean(down_hist(l3lec,:));
u_down_lec = mean_down_lec + 1*std(down_hist(l3lec,:))/sqrt(length(l3lec));
l_down_lec = mean_down_lec - 1*std(down_hist(l3lec,:))/sqrt(length(l3lec));

figure
stairs(log_dur_grid,mean_down_lfp,'color',[0.2 0.2 0.2],'linewidth',2)
hold on
% stairs(log_dur_grid,mean_down_mec,'r','linewidth',2)
% stairs(log_dur_grid,mean_down_lec,'b','linewidth',2)

legend('Ctx LFP','MEC','LEC')

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

% for i = 1:length(log_dur_grid)-1
%     clear X Y
%     bin_width = (log_dur_grid(i+1)-log_dur_grid(i));
%     X = [log_dur_grid(i) log_dur_grid(i)+bin_width];
%     Y = [u_down_lec(i) u_down_lec(i)];
%     X = [X fliplr(X)];
%     Y = [Y l_down_lec(i) l_down_lec(i)];
%     fill(X,Y,'b','EdgeColor','none');
%     hold on
% end
% 
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

% for i = 1:length(log_dur_grid)-1
%     clear X Y
%     bin_width = (log_dur_grid(i+1)-log_dur_grid(i));
%     X = [log_dur_grid(i) log_dur_grid(i)+bin_width];
%     Y = [u_down_mec(i) u_down_mec(i)];
%     X = [X fliplr(X)];
%     Y = [Y l_down_mec(i) l_down_mec(i)];
%     fill(X,Y,'r','EdgeColor','none');
%     hold on
% end


% set(gca,'yscale','log')
set(gca,'YMinorGrid','off')
set(gca,'YGrid','off')
set(gca,'fontname','arial','fontsize',14)
% set(gca,'XMinorGrid','off')
% ylim([1e-4 0.1])
% ylim([0 100])
% xlim([0.05 10])
xlim([0.05 2.5])
shg
% set(gca,'xscale','log')
xlabel('Down State Duration (s)','FontSize',16,'fontname','arial')
ylabel('Percent of Data','FontSize',16,'fontname','arial')

%% UP transition lag
close all
clear lag_* sm_*
lag_range = linspace(-0.5,1.5,1000);
dlag = lag_range(2)-lag_range(1);
for d = 1:length(combined_dir)
    lag_hist(d,:) = histc(mp_uplags{d}/Fsd,lag_range);
    lag_hist(d,:) = lag_hist(d,:)/sum(lag_hist(d,:))/dlag;
    lag_hist_nt2(d,:) = histc(mp_uplags{d}(nrt2_ups{d})/Fsd,lag_range);
    lag_hist_nt2(d,:) = lag_hist_nt2(d,:)/sum(lag_hist_nt2(d,:))/dlag;
    lag_hist_t2(d,:) = histc(mp_uplags{d}(rt2_ups{d})/Fsd,lag_range);
    lag_hist_t2(d,:) = lag_hist_t2(d,:)/sum(lag_hist_t2(d,:))/dlag;
    sm_lag_hist(d,:) = jmm_smooth_1d_cor(lag_hist(d,:),20);
    sm_lag_hist_nt2(d,:) = jmm_smooth_1d_cor(lag_hist_nt2(d,:),20);
    sm_lag_hist_t2(d,:) = jmm_smooth_1d_cor(lag_hist_t2(d,:),20);
end

% mean_lag_mec = mean(sm_lag_hist(l3mec,:));
% u_lag_mec = mean_lag_mec + 1*std(sm_lag_hist(l3mec,:))/sqrt(length(l3mec));
% l_lag_mec = mean_lag_mec - 1*std(sm_lag_hist(l3mec,:))/sqrt(length(l3mec));
% mean_lag_lec = mean(sm_lag_hist(l3lec,:));
% u_lag_lec = mean_lag_lec + 1*std(sm_lag_hist(l3lec,:))/sqrt(length(l3lec));
% l_lag_lec = mean_lag_lec - 1*std(sm_lag_hist(l3lec,:))/sqrt(length(l3lec));
% mean_lag_np_mec = nanmean(sm_lag_hist_nt2(l3mec,:));
% u_lag_np_mec = mean_lag_np_mec + 1*nanstd(sm_lag_hist_nt2(l3mec,:))/sqrt(length(l3mec));
% l_lag_np_mec = mean_lag_np_mec - 1*nanstd(sm_lag_hist_nt2(l3mec,:))/sqrt(length(l3mec));
% mean_lag_p_mec = nanmean(sm_lag_hist_t2(l3mec,:));
% u_lag_p_mec = mean_lag_p_mec + 1*nanstd(sm_lag_hist_t2(l3mec,:))/sqrt(length(l3mec));
% l_lag_p_mec = mean_lag_p_mec - 1*nanstd(sm_lag_hist_t2(l3mec,:))/sqrt(length(l3mec));

figure
set(gca,'fontsize',14,'fontname','arial')
plot(lag_range,mean(sm_lag_hist(l3mec,:)),'r','linewidth',2)
hold on
plot(lag_range,mean(sm_lag_hist(l3lec,:)),'b','linewidth',2)
% plot(lag_range,mean_lag_np_mec,'c','linewidth',2)
% plot(lag_range,mean_lag_np_mec,'g','linewidth',2)
legend('MEC','LEC')
shadedErrorBar(lag_range,mean(sm_lag_hist(l3lec,:)),std(sm_lag_hist(l3lec,:))./sqrt(length(l3lec)),{'b'});
shadedErrorBar(lag_range,mean(sm_lag_hist(l3mec,:)),std(sm_lag_hist(l3mec,:))./sqrt(length(l3mec)),{'r'});

% X = [lag_range fliplr(lag_range)];
% Y = [u_lag_lec fliplr(l_lag_lec)];
% fill(X,Y,'b')
% X = [lag_range fliplr(lag_range)];
% Y = [u_lag_mec fliplr(l_lag_mec)];
% fill(X,Y,'r')
% X = [lag_range fliplr(lag_range)];
% Y = [u_lag_np_mec fliplr(l_lag_np_mec)];
% fill(X,Y,'c')
% X = [lag_range fliplr(lag_range)];
% Y = [u_lag_p_mec fliplr(l_lag_p_mec)];
% fill(X,Y,'g')
% % xlim([-0.4 1.2])
xlim([-0.4 1.2])
ylim([0 3.2])
line([0 0],[0 3.2],'Color','k')
xlabel('Lag (s)','fontsize',16,'fontname','arial')
ylabel('Probability','fontsize',16,'fontname','arial')

%% DOWN transition lag
close all
clear lag_* sm_*
lag_range = linspace(-1,2.3,1000);
dlag = lag_range(2)-lag_range(1);
for d = 1:length(combined_dir)
    lag_hist(d,:) = histc(mp_downlags{d}/Fsd,lag_range);
    lag_hist(d,:) = lag_hist(d,:)/sum(lag_hist(d,:))/dlag;
    lag_hist_nt2(d,:) = histc(mp_downlags{d}(nrt2_ups{d})/Fsd,lag_range);
    lag_hist_nt2(d,:) = lag_hist_nt2(d,:)/sum(lag_hist_nt2(d,:))/dlag;
    sm_lag_hist(d,:) = jmm_smooth_1d_cor(lag_hist(d,:),20);
    sm_lag_hist_nt2(d,:) = jmm_smooth_1d_cor(lag_hist_nt2(d,:),20);
    %     phase_hist(d,:) = fgsmooth(phase_hist(d,:),4);
end

% mean_lag_mec = mean(sm_lag_hist(l3mec,:));
% u_lag_mec = mean_lag_mec + 1*std(sm_lag_hist(l3mec,:))/sqrt(length(l3mec));
% l_lag_mec = mean_lag_mec - 1*std(sm_lag_hist(l3mec,:))/sqrt(length(l3mec));
% mean_lag_lec = mean(sm_lag_hist(l3lec,:));
% u_lag_lec = mean_lag_lec + 1*std(sm_lag_hist(l3lec,:))/sqrt(length(l3lec));
% l_lag_lec = mean_lag_lec - 1*std(sm_lag_hist(l3lec,:))/sqrt(length(l3lec));
% mean_lag_np_mec = nanmean(sm_lag_hist_nt2(l3mec,:));
% u_lag_np_mec = mean_lag_np_mec + 1*nanstd(sm_lag_hist_nt2(l3mec,:))/sqrt(length(l3mec));
% l_lag_np_mec = mean_lag_np_mec - 1*nanstd(sm_lag_hist_nt2(l3mec,:))/sqrt(length(l3mec));

figure
set(gca,'fontsize',14,'fontname','arial')
plot(lag_range,mean(sm_lag_hist(l3mec,:)),'r','linewidth',2)
hold on
plot(lag_range,mean(sm_lag_hist(l3lec,:)),'b','linewidth',2)
% plot(lag_range,mean_lag_np_mec,'c','linewidth',2)
% plot(lag_range,mean_lag_np_mec,'g','linewidth',2)
legend('MEC','LEC')
shadedErrorBar(lag_range,mean(sm_lag_hist(l3lec,:)),std(sm_lag_hist(l3lec,:))./sqrt(length(l3lec)),{'b'});
shadedErrorBar(lag_range,mean(sm_lag_hist(l3mec,:)),std(sm_lag_hist(l3mec,:))./sqrt(length(l3mec)),{'r'});

% X = [lag_range fliplr(lag_range)];
% Y = [u_lag_lec fliplr(l_lag_lec)];
% fill(X,Y,'b')
% X = [lag_range fliplr(lag_range)];
% Y = [u_lag_mec fliplr(l_lag_mec)];
% fill(X,Y,'r')
% % X = [lag_range fliplr(lag_range)];
% % Y = [u_lag_np_mec fliplr(l_lag_np_mec)];
% % fill(X,Y,'c')
xlim([-0.6 1.5])
ylim([0 2.2])
line([0 0],[0 3],'Color','k')
xlabel('Lag (s)','fontsize',16,'fontname','arial')
ylabel('Probability','fontsize',16,'fontname','arial')

%% RELATIVE UP STATE DURATION
close all
clear dur_* sm_*
dur_range = linspace(-3,6,1000);
ddur = dur_range(2)-dur_range(1);
for d = 1:length(combined_dir)
    dur_hist(d,:) = histc(mp_rel_updurs{d},dur_range);
    dur_hist_nt2(d,:) = histc(mp_rel_updurs{d}(nrt2_ups{d}),dur_range);
    dur_hist_nt2(d,:) = dur_hist_nt2(d,:)/sum(dur_hist(d,:))/ddur;
    dur_hist(d,:) = dur_hist(d,:)/sum(dur_hist(d,:))/ddur;
    sm_dur_hist(d,:) = jmm_smooth_1d_cor(dur_hist(d,:),20);
    sm_dur_hist_nt2(d,:) = jmm_smooth_1d_cor(dur_hist_nt2(d,:),20);
end
eps = 1e-5;
mean_dur_mec = mean(sm_dur_hist(l3mec,:));
u_dur_mec = mean_dur_mec + 1*std(sm_dur_hist(l3mec,:))/sqrt(length(l3mec));
l_dur_mec = mean_dur_mec - 1*std(sm_dur_hist(l3mec,:))/sqrt(length(l3mec));
mean_dur_mec(mean_dur_mec <= eps) = eps;
u_dur_mec(u_dur_mec <= eps) = eps;
l_dur_mec(l_dur_mec <= eps) = eps;
mean_dur_lec = mean(sm_dur_hist(l3lec,:));
u_dur_lec = mean_dur_lec + 1*std(sm_dur_hist(l3lec,:))/sqrt(length(l3lec));
l_dur_lec = mean_dur_lec - 1*std(sm_dur_hist(l3lec,:))/sqrt(length(l3lec));
mean_dur_lec(mean_dur_lec <= eps) = eps;
u_dur_lec(u_dur_lec <= eps) = eps;
l_dur_lec(l_dur_lec <= eps) = eps;
mean_dur_np_mec = mean(sm_dur_hist_nt2(l3mec,:));
u_dur_np_mec = mean_dur_np_mec + 1*std(sm_dur_hist_nt2(l3mec,:))/sqrt(length(l3mec));
l_dur_np_mec = mean_dur_np_mec - 1*std(sm_dur_hist_nt2(l3mec,:))/sqrt(length(l3mec));
mean_dur_np_mec(mean_dur_np_mec <= eps) = eps;
u_dur_np_mec(u_dur_np_mec <= eps) = eps;
l_dur_np_mec(l_dur_np_mec <= eps) = eps;
mean_dur_np_lec = mean(sm_dur_hist_nt2(l3lec,:));
u_dur_np_lec = mean_dur_np_lec + 1*std(sm_dur_hist_nt2(l3lec,:))/sqrt(length(l3lec));
l_dur_np_lec = mean_dur_np_lec - 1*std(sm_dur_hist_nt2(l3lec,:))/sqrt(length(l3lec));
mean_dur_np_lec(mean_dur_np_lec <= eps) = eps;
u_dur_np_lec(u_dur_np_lec <= eps) = eps;
l_dur_np_lec(l_dur_np_lec <= eps) = eps;

figure
set(gca,'fontsize',14,'fontname','arial')
plot(dur_range,mean_dur_mec,'r','linewidth',2)
hold on
plot(dur_range,mean_dur_lec,'b','linewidth',2)
plot(dur_range,mean_dur_np_mec,'color',[1 0.3 0.4],'linewidth',2)
plot(dur_range,mean_dur_np_lec,'c','linewidth',2)

X = [dur_range fliplr(dur_range)];
Y = [u_dur_lec fliplr(l_dur_lec)];
fill(X,Y,'b')
X = [dur_range fliplr(dur_range)];
Y = [u_dur_mec fliplr(l_dur_mec)];
fill(X,Y,'r')
X = [dur_range fliplr(dur_range)];
Y = [u_dur_np_mec fliplr(l_dur_np_mec)];
fill(X,Y,[1 0.3 0.4])
X = [dur_range fliplr(dur_range)];
Y = [u_dur_np_lec fliplr(l_dur_np_lec)];
fill(X,Y,'c')
xlim([-3 5])
xlabel('Relative Duration (s)','fontsize',16,'fontname','arial')
ylabel('Probability Density','fontsize',16,'fontname','arial')
% ylim([0 1.3])
% line([0 0],[0 1.3],'color','k')
set(gca,'yscale','log')
line([0 0],[1e-3 1.3],'color','k')
ylim([3e-3 1.3])

%% UP quantization
% hist_range = linspace(0,6,1000);
hist_range = linspace(0,8,1500);
% hist_range = linspace(0,8,100);
binsize = hist_range(2)-hist_range(1);
% l3mec(l3mec > 22) = [];
for i = 1:length(combined_dir)
    cell_hist(i,:) = hist(rmp_updurs_lfpc{i},hist_range);
    cell_hist(i,:) = cell_hist(i,:)/sum(cell_hist(i,:))/binsize;
%     cell_hist(i,:) = cell_hist(i,:)/sum(cell_hist(i,:));
    cell_hist(i,:) = jmm_smooth_1d_cor(cell_hist(i,:),20);
%     cell_hist(i,:) = jmm_smooth_1d_cor(cell_hist(i,:),1);
end

mean_hist_mec = nanmean(cell_hist(l3mec,:));
u_mec = mean_hist_mec+1*nanstd(cell_hist(l3mec,:))./sqrt(sum(~isnan(cell_hist(l3mec,:))));
l_mec = mean_hist_mec-1*nanstd(cell_hist(l3mec,:))./sqrt(sum(~isnan(cell_hist(l3mec,:))));

mean_hist_lec = nanmean(cell_hist(l3lec,:));
u_lec = mean_hist_lec+1*nanstd(cell_hist(l3lec,:))./sqrt(sum(~isnan(cell_hist(l3lec,:))));
l_lec = mean_hist_lec-1*nanstd(cell_hist(l3lec,:))./sqrt(sum(~isnan(cell_hist(l3lec,:))));

figure
set(gca,'fontsize',14,'fontname','arial')
plot(hist_range,mean_hist_mec,'r','linewidth',2)
hold on
% plot(hist_range,mean_hist_lec,'b','linewidth',2)
legend('MEC','LEC')

epsilon = 1e-10;
cell_hist(cell_hist < epsilon) = epsilon;
u_mec(u_mec < epsilon) = epsilon;
l_mec(l_mec < epsilon) = epsilon;
u_lec(u_lec < epsilon) = epsilon;
l_lec(l_lec < epsilon) = epsilon;

X = [hist_range fliplr(hist_range)];
Y = [u_mec fliplr(l_mec)];
fill(X,Y,'r')
Y = [u_lec fliplr(l_lec)];
% fill(X,Y,'b')
%
set(gca,'yscale','log')
xlim([0 5.3])
% xlim([0 7.3])
ylim([1e-4 4])
xlabel('Duration (Ncx UDS Cycles)','fontsize',16,'fontname','arial')
ylabel('Probability Density','fontsize',16,'fontname','arial')
[peak_amps,b] = findpeaks(mean_hist_mec);
peak_locs = hist_range(b);

% load ./combined_core_analysis_lfpmprev.mat
% for i = 1:length(combined_dir)
%     cell_hist(i,:) = hist(rmp_updurs_lfpc{i},hist_range);
%     cell_hist(i,:) = cell_hist(i,:)/sum(cell_hist(i,:))/binsize;
%     cell_hist(i,:) = jmm_smooth_1d_cor(cell_hist(i,:),20);
% end
% 
% mean_hist_mec = nanmean(cell_hist(l3mec,:));
% u_mec = mean_hist_mec+1*nanstd(cell_hist(l3mec,:))./sqrt(sum(~isnan(cell_hist(l3mec,:))));
% l_mec = mean_hist_mec-1*nanstd(cell_hist(l3mec,:))./sqrt(sum(~isnan(cell_hist(l3mec,:))));
% 
% mean_hist_lec = nanmean(cell_hist(l3lec,:));
% u_lec = mean_hist_lec+1*nanstd(cell_hist(l3lec,:))./sqrt(sum(~isnan(cell_hist(l3lec,:))));
% l_lec = mean_hist_lec-1*nanstd(cell_hist(l3lec,:))./sqrt(sum(~isnan(cell_hist(l3lec,:))));
% 
% plot(hist_range,mean_hist_mec,'k','linewidth',2)
% hold on
% plot(hist_range,mean_hist_lec,'g','linewidth',2)

%%
%good examples are 1,2,8,14,15,19,21,46,53,54,57
%1, 2, 8, 14, 21, 43, 
hist_range = linspace(0,7,1000);
binsize = hist_range(2)-hist_range(1);
clear cell_hist
for i = 1:length(combined_dir)
    cell_hist(i,:) = hist(rmp_updurs_lfpc{i},hist_range);
    cell_hist(i,:) = cell_hist(i,:)/sum(cell_hist(i,:))/binsize;
    cell_hist(i,:) = jmm_smooth_1d_cor(cell_hist(i,:),20);
end
% for i = 36:58
examp_ind = 43
x = linspace(0,7,100);
binsize = x(2)-x(1);

y = hist(rmp_updurs_lfpc{examp_ind},x);
sy = sum(y);

figure
set(gca,'fontname','arial','fontsize',14)
bar(x,y/sy/binsize,1,'k')
hold on
xlim([0 4.2])
xlabel('Neocortical UDS Cycles','fontsize',16,'fontname','arial')
ylabel('Relative frequency','fontsize',16,'fontname','arial')
shg

% pause
% close
% end

%%
figure
plot(median_uplag_nt2(l3mec)/Fsd,median_uplag_t2(l3mec)/Fsd,'ro')
% hold on
% plot(median_uplag_nt2(l3lec),median_uplag_t2(l3lec),'bo')
xlabel('Up-transition lag non-pers (s)','fontsize',16)
ylabel('Up-transition lag pers (s)','fontsize',16)
xl = xlim(); yl = ylim();
line([0 .5],[0 .5],'color','k')

figure
plot(median_reluplag_nt2(l3mec),median_reluplag_t2(l3mec),'ro')
% hold on
% plot(median_uplag_nt2(l3lec),median_uplag_t2(l3lec),'bo')
xlabel('Relative Up-transition lag non-pers','fontsize',16)
ylabel('Relative Up-transition lag pers','fontsize',16)
xl = xlim(); yl = ylim();
line([0 .5],[0 .5],'color','k')

figure
plot(median_downlag_nt2(l3mec)/Fsd,median_downlag_t2(l3mec)/Fsd,'ro')
% hold on
% plot(median_uplag_nt2(l3lec),median_uplag_t2(l3lec),'bo')
xlabel('Down-transition lag non-pers (s)','fontsize',16)
ylabel('Down-transition lag pers (s)','fontsize',16)
xl = xlim(); yl = ylim();
line([0 1.5],[0 1.5],'color','k')

figure
plot(median_reldownlag_nt2(l3mec),median_reldownlag_t2(l3mec),'ro')
% hold on
% plot(median_uplag_nt2(l3lec),median_uplag_t2(l3lec),'bo')
xlabel('Relative Down-transition lag non-pers','fontsize',16)
ylabel('Relative Down-transition lag pers','fontsize',16)
xl = xlim(); yl = ylim();
line([0 1.],[0 1.],'color','k')

%%
[dlag_mec_c,dlag_mec_p] = corrcoef(median_downlag(l3mec),fract_rt2_ups(l3mec));
[dlag_lec_c,dlag_lec_p] = corrcoef(median_downlag(l3lec),fract_rt2_ups(l3lec));
uset = find(fract_rt2_ups(l3mec) > 0);
[dlag_mec_cl,dlag_mec_pl] = corrcoef(median_downlag(l3mec(uset)),log(fract_rt2_ups(l3mec(uset))));
uset = find(fract_rt2_ups(l3lec) > 0);
[dlag_lec_cl,dlag_lec_pl] = corrcoef(median_downlag(l3lec(uset)),log(fract_rt2_ups(l3lec(uset))));
[ulag_mec_c,ulag_mec_p] = corrcoef(median_uplag(l3mec),fract_rt2_downs(l3mec));
[ulag_lec_c,ulag_lec_p] = corrcoef(median_uplag(l3lec),fract_rt2_downs(l3lec));
uset = find(fract_rt2_downs(l3mec) > 0);
[ulag_mec_cl,ulag_mec_pl] = corrcoef(median_uplag(l3mec(uset)),log(fract_rt2_downs(l3mec(uset))));
uset = find(fract_rt2_downs(l3lec) > 0);
[ulag_lec_cl,ulag_lec_pl] = corrcoef(median_uplag(l3lec(uset)),log(fract_rt2_downs(l3lec(uset))));
   
eps = 0;
figure
plot(median_downlag(l3mec)/Fsd,fract_rt2_ups(l3mec)+eps,'ro','markersize',8)
hold on
plot(median_downlag(l3lec)/Fsd,fract_rt2_ups(l3lec)+eps,'o','markersize',8)
used = find(fract_rt2_ups + eps > 0);
temp_all = polyfit(median_downlag(used)/Fsd,log(fract_rt2_ups(used)+eps),1);
temp_mec = polyfit(median_downlag(l3mec)/Fsd,log(fract_rt2_ups(l3mec)+eps),1);
temp_lec = polyfit(median_downlag(l3lec)/Fsd,log(fract_rt2_ups(l3lec)+eps),1);
x_ax = linspace(-0.2,1.25,100);
set(gca,'fontsize',14,'fontname','arial')
plot(x_ax,exp(polyval(temp_all,x_ax)),'k')
xlabel('Median down-transition lag (s)','fontsize',14,'fontname','arial')
ylabel('Log probability persistent up state','fontsize',14,'fontname','arial')
set(gca,'yscale','log')
xlim([-0.2 1.25])
ylim([4e-3 1])

figure
plot(median_reldownlag(l3mec),fract_rt2_ups(l3mec)+eps,'ro','markersize',8,'linewidth',2)
hold on
plot(median_reldownlag(l3lec),fract_rt2_ups(l3lec)+eps,'o','markersize',8,'linewidth',2)
used = find(fract_rt2_ups + eps > 0);
temp_all = polyfit(median_reldownlag(used),log(fract_rt2_ups(used)+eps),1);
temp_mec = polyfit(median_reldownlag(l3mec),log(fract_rt2_ups(l3mec)+eps),1);
temp_lec = polyfit(median_reldownlag(l3lec),log(fract_rt2_ups(l3lec)+eps),1);
x_ax = linspace(-0.2,1.25,100);
set(gca,'fontsize',14,'fontname','arial')
plot(x_ax,exp(polyval(temp_all,x_ax)),'k')
xlabel('Median relative down-transition lag','fontsize',14,'fontname','arial')
ylabel('Log probability persistent up state','fontsize',14,'fontname','arial')
set(gca,'yscale','log')
xlim([-0.1 0.8])
ylim([4e-3 1])

%%
%for state duration distribution calculations
dur_range = [0.02 15];
numBins = 60;
dur_grid = logspace(log10(dur_range(1)),log10(dur_range(2)),numBins+1);
% dur_grid = linspace(dur_range(1),dur_range(2),numBins+1);
% dur_grid = [-fliplr(dur_grid) dur_grid];
dur_grid_cents = (dur_grid(1:end-1)+dur_grid(2:end))/2;
clear temp temp2
for i = 1:61
%     temp(i,:) = histc(mp_rel_updurs{i},dur_grid);
    temp(i,:) = histc(mp_state_durations{i}{2},dur_grid);
    temp(i,:) = temp(i,:)/length(mp_rel_updurs{i});
    temp2(i,:) = histc(lfp_state_durations{i}{2},dur_grid);
    temp2(i,:) = temp2(i,:)/length(lfp_state_durations{i}{2});
    mm(i) = max(lfp_state_durations{i}{2});
end
temp(:,end) = [];
temp = bsxfun(@rdivide,temp,diff(dur_grid));
temp2(:,end) = [];
temp2 = bsxfun(@rdivide,temp2,diff(dur_grid));
shadedErrorBar(dur_grid_cents,mean(temp(l3mec,:)),std(temp(l3mec,:))/sqrt(length(l3mec)),'r')
hold on
shadedErrorBar(dur_grid_cents,mean(temp2(l3mec,:)),std(temp2(l3mec,:))/sqrt(length(l3mec)),'b')


%% Up State Distribution Compare
boundaries = prctile(median_downdur8(l3mec),[33 67]);
low_chunk = l3mec(median_downdur8(l3mec) < boundaries(1));
high_chunk = l3mec(median_downdur8(l3mec) > boundaries(2));

%use logarithmic hists
up_hist = up_dur_loghist;
up_hist8 = up_dur_loghist8;

%normalize state duration distributions 
up_hist = up_hist./repmat(sum(up_hist,2),1,size(up_hist,2));
up_hist8 = up_hist8./repmat(sum(up_hist8,2),1,size(up_hist8,2));
 
mean_up_mec_l = mean(up_hist(low_chunk,:));
s_up_mec_l = 1*std(up_hist(low_chunk,:))/sqrt(length(low_chunk));
mean_up_mec_h = mean(up_hist(high_chunk,:));
s_up_mec_h = 1*std(up_hist(high_chunk,:))/sqrt(length(high_chunk));

mean_up_lfp_l = mean(up_hist8(low_chunk,:));
s_up_lfp_l = 1*std(up_hist8(low_chunk,:))/sqrt(length(low_chunk));
mean_up_lfp_h = mean(up_hist8(high_chunk,:));
s_up_lfp_h = 1*std(up_hist8(high_chunk,:))/sqrt(length(high_chunk));

figure
h = errorbar(log_dur_grid,mean_up_lfp_l,s_up_lfp_l,'color',[0.2 0.2 0.2],'linewidth',2)
hold on
h = errorbar(log_dur_grid,mean_up_mec_l,s_up_mec_l,'r','linewidth',2)
h = errorbar(log_dur_grid,mean_up_lfp_h,s_up_lfp_h,'g','linewidth',2)
h = errorbar(log_dur_grid,mean_up_mec_h,s_up_mec_h,'c','linewidth',2)

set(gca,'yscale','log')
set(gca,'YMinorGrid','off')
set(gca,'YGrid','off')
ylim([1e-4 0.2])
xlim([0.05 10])
shg
xlabel('Up State Duration (s)','FontSize',14)
ylabel('Percent of Data','FontSize',14)

%% Histogram of Type 2 Persistent Activity
cd C:\WC_Germany\sven_thomas_combined\
load ./combined_dir_nd.mat
load ./combined_core_analysis_fin_nd.mat

numBins = 10;
pers_range = linspace(0,0.40,numBins+1);

pers_mec = histc(fract_rt2_ups(l3mec),pers_range);
pers_lec = histc(fract_rt2_ups(l3lec),pers_range);

figure
set(gca,'fontname','arial','fontsize',14)
stairs(pers_range*100,pers_mec,'r','linewidth',2)
hold on
stairs(pers_range*100,pers_lec,'b','linewidth',2)
xlim([0 38])
% ylim([0 5])
xlabel('Fraction of Up States Persisting','FontSize',16,'fontname','arial')
ylabel('Number of Cells','FontSize',16,'fontname','arial')
legend('MEC','LEC')

%%
% figure
% plot(median_downlag(l3mec)/Fsd,100*fract_rt2_ups(l3mec),'o')
% hold on
% plot(median_downlag(l3lec)/Fsd,100*fract_rt2_ups(l3lec),'go')
% xlabel('Down-transition lag (s)')
% ylabel('Percent Type 2 persistence')
% lag_ax1 = linspace(0.3,1.4,400);
% temp_p = polyfit(median_downlag(l3mec)/Fsd,100*fract_rt2_ups(l3mec),1);
% plot(lag_ax1,polyval(temp_p,lag_ax1))
% lag_ax2 = linspace(-0.3,0.3,400);
% temp_p = polyfit(median_downlag(l3lec)/Fsd,100*fract_rt2_ups(l3lec),1);
% plot(lag_ax2,polyval(temp_p,lag_ax2),'g')
% ylim([0 50])

figure
plot(median_reldownlag(l3mec),100*fract_rt2_ups(l3mec),'o')
hold on
plot(median_reldownlag(l3lec),100*fract_rt2_ups(l3lec),'go')
xlabel('Down-transition lag (s)')
ylabel('Percent Type 2 persistence')
ylim([0 50])



%% Depth of anesthesia figures (cortical down dur)
figure
plot(median_downdur8(l3mec),100*fract_rt2_ups(l3mec),'o')
xlabel('Ctx down state duration(s)','fontsize',16)
ylabel('Percent Type 2 persistence','fontsize',16)
ylim([0 50])

figure
plot(median_downdur8(l3mec),median_downlag(l3mec)/Fsd,'o')
x_ax = linspace(0.5,2,100);
temp_p = polyfit(median_downdur8(l3mec),median_downlag(l3mec)/Fsd,1);
hold on
plot(x_ax,polyval(temp_p,x_ax),'k')
xlabel('Ctx down state duration(s)','fontsize',16)
ylabel('Down-transition lag (s)','fontsize',16)

figure
plot(median_downdur8(l3mec),median_reldownlag(l3mec),'o')
x_ax = linspace(0.5,2,100);
temp_p = polyfit(median_downdur8(l3mec),median_reldownlag(l3mec),1);
hold on
plot(x_ax,polyval(temp_p,x_ax),'k')
xlabel('Ctx down state duration(s)','fontsize',16)
ylabel('Relative down-transition lag','fontsize',16)

%% Depth of anesthesia figures (cortical duty cycle)
figure
plot(median_dc8(l3mec),100*fract_rt2_ups(l3mec),'o')
xlabel('Ctx duty cycle','fontsize',16)
ylabel('Percent Type 2 persistence','fontsize',16)
ylim([0 50])

figure
plot(median_dc8(l3mec),median_downlag(l3mec)/Fsd,'o')
x_ax = linspace(0.2,0.7,100);
temp_p = polyfit(median_dc8(l3mec),median_downlag(l3mec)/Fsd,1);
hold on
plot(x_ax,polyval(temp_p,x_ax),'k')
xlabel('Ctx duty cycle','fontsize',16)
ylabel('Down-transition lag (s)','fontsize',16)

figure
plot(median_dc8(l3mec),median_reldownlag(l3mec),'o')
x_ax = linspace(0.2,0.7,100);
temp_p = polyfit(median_dc8(l3mec),median_reldownlag(l3mec),1);
hold on
plot(x_ax,polyval(temp_p,x_ax),'k')
xlabel('Ctx duty cycle','fontsize',16)
ylabel('Relative down-transition lag','fontsize',16)

%% Depth of anesthesia figures (cortical uds freq)
load ./combined_spectra
uds_freqs = find(f > 0.2 & f < 1);
lS = real(10*log10(S));
[peak_S,peakf] = max(lS(:,:,uds_freqs),[],3);
peak_uds_freq = f(uds_freqs(peakf));
uds8_freq = peak_uds_freq(:,8)';
figure
plot(uds8_freq(l3mec),100*fract_rt2_ups(l3mec),'o')
xlabel('Ctx UDS freq (Hz)','fontsize',16)
ylabel('Percent Type 2 persistence','fontsize',16)
ylim([0 50])

figure
plot(uds8_freq(l3mec),median_downlag(l3mec)/Fsd,'o')
x_ax = linspace(0.15,0.7,100);
temp_p = polyfit(uds8_freq(l3mec),median_downlag(l3mec)/Fsd,1);
hold on
plot(x_ax,polyval(temp_p,x_ax),'k')
xlabel('Ctx UDS freq (Hz)','fontsize',16)
ylabel('Down-transition lag (s)','fontsize',16)

figure
plot(uds8_freq(l3mec),median_reldownlag(l3mec),'o')
x_ax = linspace(0.15,0.7,100);
temp_p = polyfit(uds8_freq(l3mec),median_reldownlag(l3mec),1);
hold on
plot(x_ax,polyval(temp_p,x_ax),'k')
xlabel('Ctx UDS freq (Hz)','fontsize',16)
ylabel('Relative down-transition lag','fontsize',16)