load C:\WC_Germany\persistent_revised\trig_avgs\trig_avg_data_wideband
load C:\WC_Germany\persistent_revised\pers_revised_dir

pyr_cells = 1:n_mec;
lec_cells = n_mec+1:n_mec+n_lec;
adj_cells = [];

%% for up triggered average
pyr_avg_utrig = nanmean(mp_utrig_lf8(pyr_cells,:));
pyr_uc_utrig = pyr_avg_utrig+nanstd(mp_utrig_lf8(pyr_cells,:))/sqrt(length(pyr_cells));
pyr_lc_utrig = pyr_avg_utrig-nanstd(mp_utrig_lf8(pyr_cells,:))/sqrt(length(pyr_cells));

lec_avg_utrig = nanmean(mp_utrig_lf8(lec_cells,:));
lec_uc_utrig = lec_avg_utrig+nanstd(mp_utrig_lf8(lec_cells,:))/sqrt(length(lec_cells));
lec_lc_utrig = lec_avg_utrig-nanstd(mp_utrig_lf8(lec_cells,:))/sqrt(length(lec_cells));

adj_avg_utrig = nanmean(mp_utrig_lf8(adj_cells,:));
adj_uc_utrig = adj_avg_utrig+nanstd(mp_utrig_lf8(adj_cells,:))/sqrt(length(adj_cells));
adj_lc_utrig = adj_avg_utrig-nanstd(mp_utrig_lf8(adj_cells,:))/sqrt(length(adj_cells));



figure
plot(lags,pyr_avg_utrig,'linewidth',2)
hold on
plot(lags,lec_avg_utrig,'g','linewidth',2)
plot(lags,adj_avg_utrig,'k','linewidth',2)

legend('MEC','LEC','ADJ')

X = [lags fliplr(lags)];
Y = [pyr_uc_utrig fliplr(pyr_lc_utrig)];

fill(X,Y,'b')

Y = [lec_uc_utrig fliplr(lec_lc_utrig)];
fill(X,Y,'g')

Y = [adj_uc_utrig fliplr(adj_lc_utrig)];
fill(X,Y,'k')

xlim([-2 2])
xlabel('Time Lag (s)')
ylabel('Amplitude (z)')
line([0 0],[-1 2],'Color','k')
line([-2 2],[0 0],'Color','k')



%% for down triggered average
pyr_avg_dtrig = nanmean(mp_dtrig_lf8(pyr_cells,:));
pyr_uc_dtrig = pyr_avg_dtrig+nanstd(mp_dtrig_lf8(pyr_cells,:))/sqrt(length(pyr_cells));
pyr_lc_dtrig = pyr_avg_dtrig-nanstd(mp_dtrig_lf8(pyr_cells,:))/sqrt(length(pyr_cells));

lec_avg_dtrig = nanmean(mp_dtrig_lf8(lec_cells,:));
lec_uc_dtrig = lec_avg_dtrig+nanstd(mp_dtrig_lf8(lec_cells,:))/sqrt(length(lec_cells));
lec_lc_dtrig = lec_avg_dtrig-nanstd(mp_dtrig_lf8(lec_cells,:))/sqrt(length(lec_cells));

adj_avg_dtrig = nanmean(mp_dtrig_lf8(adj_cells,:));
adj_uc_dtrig = adj_avg_dtrig+nanstd(mp_dtrig_lf8(adj_cells,:))/sqrt(length(adj_cells));
adj_lc_dtrig = adj_avg_dtrig-nanstd(mp_dtrig_lf8(adj_cells,:))/sqrt(length(adj_cells));



figure
plot(lags,pyr_avg_dtrig,'linewidth',2)
hold on
plot(lags,lec_avg_dtrig,'g','linewidth',2)
plot(lags,adj_avg_dtrig,'k','linewidth',2)

legend('MEC','LEC','ADJ')

X = [lags fliplr(lags)];
Y = [pyr_uc_dtrig fliplr(pyr_lc_dtrig)];

fill(X,Y,'b')

Y = [lec_uc_dtrig fliplr(lec_lc_dtrig)];
fill(X,Y,'g')

Y = [adj_uc_dtrig fliplr(adj_lc_dtrig)];
fill(X,Y,'k')

xlim([-2 2])
xlabel('Time Lag (s)')
ylabel('Amplitude (z)')
line([0 0],[-1 1.5],'Color','k')
line([-2 2],[0 0],'Color','k')