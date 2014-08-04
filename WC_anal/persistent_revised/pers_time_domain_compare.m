load C:\WC_Germany\persistent_revised\time_domain\time_domain_data

pyr_cells = 1:17;
lec_cells = 18:21;
adj_cells = 22:27;

%% for MP acorr triggered average
pyr_avg_acor = nanmean(tot_wcv_acorr(pyr_cells,:));
pyr_uc_acor = pyr_avg_acor+nanstd(tot_wcv_acorr(pyr_cells,:))/sqrt(length(pyr_cells));
pyr_lc_acor = pyr_avg_acor-nanstd(tot_wcv_acorr(pyr_cells,:))/sqrt(length(pyr_cells));

pyr8_avg_acor = nanmean(tot_lf8_acorr(pyr_cells,:));
pyr8_uc_acor = pyr8_avg_acor+nanstd(tot_lf8_acorr(pyr_cells,:))/sqrt(length(pyr_cells));
pyr8_lc_acor = pyr8_avg_acor-nanstd(tot_lf8_acorr(pyr_cells,:))/sqrt(length(pyr_cells));

lec_avg_acor = nanmean(tot_wcv_acorr(lec_cells,:));
lec_uc_acor = lec_avg_acor+nanstd(tot_wcv_acorr(lec_cells,:))/sqrt(length(lec_cells));
lec_lc_acor = lec_avg_acor-nanstd(tot_wcv_acorr(lec_cells,:))/sqrt(length(lec_cells));

adj_avg_acor = nanmean(tot_wcv_acorr(adj_cells,:));
adj_uc_acor = adj_avg_acor+nanstd(tot_wcv_acorr(adj_cells,:))/sqrt(length(adj_cells));
adj_lc_acor = adj_avg_acor-nanstd(tot_wcv_acorr(adj_cells,:))/sqrt(length(adj_cells));



figure
plot(lags/Fsd,pyr_avg_acor,'linewidth',2)
hold on
plot(lags/Fsd,pyr8_avg_acor,'r','linewidth',2)
plot(lags/Fsd,lec_avg_acor,'g','linewidth',2)
plot(lags/Fsd,adj_avg_acor,'k','linewidth',2)

legend('MEC','LF8','LEC','ADJ')

X = [lags/Fsd fliplr(lags/Fsd)];
Y = [pyr_uc_acor fliplr(pyr_lc_acor)];

fill(X,Y,'b')

Y = [pyr8_uc_acor fliplr(pyr8_lc_acor)];
fill(X,Y,'r')

Y = [lec_uc_acor fliplr(lec_lc_acor)];
fill(X,Y,'g')

Y = [adj_uc_acor fliplr(adj_lc_acor)];
fill(X,Y,'k')
ylim([-0.5 1])
xlim([-5 5])
xlabel('Time Lag (s)')
ylabel('Correlation (z)')
line([0 0],[-1 2],'Color','k')
line([-5 5],[0 0],'Color','k')


%% for MP xcor triggered average
pyr_avg_xcor = nanmean(tot_w8_x(pyr_cells,:));
pyr_uc_xcor = pyr_avg_xcor+nanstd(tot_w8_x(pyr_cells,:))/sqrt(length(pyr_cells));
pyr_lc_xcor = pyr_avg_xcor-nanstd(tot_w8_x(pyr_cells,:))/sqrt(length(pyr_cells));

lec_avg_xcor = nanmean(tot_w8_x(lec_cells,:));
lec_uc_xcor = lec_avg_xcor+nanstd(tot_w8_x(lec_cells,:))/sqrt(length(lec_cells));
lec_lc_xcor = lec_avg_xcor-nanstd(tot_w8_x(lec_cells,:))/sqrt(length(lec_cells));

adj_avg_xcor = nanmean(tot_w8_x(adj_cells,:));
adj_uc_xcor = adj_avg_xcor+nanstd(tot_w8_x(adj_cells,:))/sqrt(length(adj_cells));
adj_lc_xcor = adj_avg_xcor-nanstd(tot_w8_x(adj_cells,:))/sqrt(length(adj_cells));



figure
plot(lags/Fsd,pyr_avg_xcor,'linewidth',2)
hold on
plot(lags/Fsd,lec_avg_xcor,'g','linewidth',2)
plot(lags/Fsd,adj_avg_xcor,'k','linewidth',2)

legend('MEC','LEC','ADJ')

X = [lags/Fsd fliplr(lags/Fsd)];
Y = [pyr_uc_xcor fliplr(pyr_lc_xcor)];

fill(X,Y,'b')

Y = [lec_uc_xcor fliplr(lec_lc_xcor)];
fill(X,Y,'g')

Y = [adj_uc_xcor fliplr(adj_lc_xcor)];
fill(X,Y,'k')
ylim([-0.5 1])
xlim([-5 5])
xlabel('Time Lag (s)')
ylabel('Correlation (z)')
line([0 0],[-1 2],'Color','k')
line([-5 5],[0 0],'Color','k')

