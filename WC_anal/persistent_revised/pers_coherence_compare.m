load C:\WC_Germany\persistent_revised\coherence\coherence_data

pyr_cells = 1:21;
lec_cells = 22:28;
lec_cells(4) = []; %get rid of 12-09-A
pyr_cells(18) = []; %get rid of 04-07

pyr_avg_coh = nanmean(Cmn(pyr_cells,:));
pyr_uc_coh = pyr_avg_coh+nanstd(Cmn(pyr_cells,:))/sqrt(length(pyr_cells));
pyr_lc_coh = pyr_avg_coh-nanstd(Cmn(pyr_cells,:))/sqrt(length(pyr_cells));

lec_avg_coh = nanmean(Cmn(lec_cells,:));
lec_uc_coh = lec_avg_coh+nanstd(Cmn(lec_cells,:))/sqrt(length(lec_cells));
lec_lc_coh = lec_avg_coh-nanstd(Cmn(lec_cells,:))/sqrt(length(lec_cells));

% adj_avg_coh = nanmean(Cmn(adj_cells,:));
% adj_uc_coh = adj_avg_coh+nanstd(Cmn(adj_cells,:))/sqrt(length(adj_cells));
% adj_lc_coh = adj_avg_coh-nanstd(Cmn(adj_cells,:))/sqrt(length(adj_cells));



figure
plot(f,pyr_avg_coh,'linewidth',2)
hold on
plot(f,lec_avg_coh,'g','linewidth',2)
% plot(f,adj_avg_coh,'k','linewidth',2)

legend('MEC','LEC')

X = [f fliplr(f)];
Y = [pyr_uc_coh fliplr(pyr_lc_coh)];

fill(X,Y,'b')

Y = [lec_uc_coh fliplr(lec_lc_coh)];
fill(X,Y,'g')

% Y = [adj_uc_coh fliplr(adj_lc_coh)];
% fill(X,Y,'k')

xlim([0 2])
xlabel('Frequency (Hz)')
ylabel('Coherence')