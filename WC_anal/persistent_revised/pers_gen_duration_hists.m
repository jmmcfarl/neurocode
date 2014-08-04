clear all
close all

load C:\WC_Germany\persistent_revised\pers_revised_dir
load C:\WC_Germany\persistent_revised\UDS_synch_state_dur\UDS_synch_state_dur_data
pyr_cells = 1:21;
lec_cells = 22:28;
lec_cells(4) = []; %get rid of 12-09-A
pyr_cells(18) = []; %get rid of 04-07

% adj_cells = (n_mec+1+n_lec):(n_mec+n_lec+n_adj);

%% for up state durations
pyr_avg_up_hist = nanmean(up_hist(pyr_cells,:));
pyr_uc_up_hist = pyr_avg_up_hist+nanstd(up_hist(pyr_cells,:))/sqrt(length(pyr_cells));
pyr_lc_up_hist = pyr_avg_up_hist-nanstd(up_hist(pyr_cells,:))/sqrt(length(pyr_cells));

pyr8_avg_up_hist = nanmean(up_hist8(pyr_cells,:));
pyr8_uc_up_hist = pyr8_avg_up_hist+nanstd(up_hist8(pyr_cells,:))/sqrt(length(pyr_cells));
pyr8_lc_up_hist = pyr8_avg_up_hist-nanstd(up_hist8(pyr_cells,:))/sqrt(length(pyr_cells));


lec_avg_up_hist = nanmean(up_hist(lec_cells,:));
lec_uc_up_hist = lec_avg_up_hist+nanstd(up_hist(lec_cells,:))/sqrt(length(lec_cells));
lec_lc_up_hist = lec_avg_up_hist-nanstd(up_hist(lec_cells,:))/sqrt(length(lec_cells));

% adj_avg_up_hist = nanmean(up_hist(adj_cells,:));
% adj_uc_up_hist = adj_avg_up_hist+nanstd(up_hist(adj_cells,:))/sqrt(length(adj_cells));
% adj_lc_up_hist = adj_avg_up_hist-nanstd(up_hist(adj_cells,:))/sqrt(length(adj_cells));

figure
stairs(up_grid,pyr_avg_up_hist,'linewidth',2)
hold on
stairs(up_grid,lec_avg_up_hist,'g','linewidth',2)
% stairs(up_grid,adj_avg_up_hist,'k','linewidth',2)
stairs(up_grid,pyr8_avg_up_hist,'r','linewidth',2)

legend('MEC','LEC','LF8')

for i = 1:length(up_grid)-1
    clear X Y
    bin_width = (up_grid(i+1)-up_grid(i));
    X = [up_grid(i) up_grid(i)+bin_width];
    Y = [pyr_uc_up_hist(i) pyr_uc_up_hist(i)];
    X = [X fliplr(X)];
    Y = [Y pyr_lc_up_hist(i) pyr_lc_up_hist(i)];
    fill(X,Y,'b','EdgeColor','none');
    hold on
end

for i = 1:length(up_grid)-1
    clear X Y
    bin_width = (up_grid(i+1)-up_grid(i));
    X = [up_grid(i) up_grid(i)+bin_width];
    Y = [pyr8_uc_up_hist(i) pyr8_uc_up_hist(i)];
    X = [X fliplr(X)];
    Y = [Y pyr8_lc_up_hist(i) pyr8_lc_up_hist(i)];
    fill(X,Y,'r','EdgeColor','none');
    hold on
end

for i = 1:length(up_grid)-1
    clear X Y
    bin_width = (up_grid(i+1)-up_grid(i));
    X = [up_grid(i) up_grid(i)+bin_width];
    Y = [lec_uc_up_hist(i) lec_uc_up_hist(i)];
    X = [X fliplr(X)];
    Y = [Y lec_lc_up_hist(i) lec_lc_up_hist(i)];
    fill(X,Y,'g','EdgeColor','none');
    hold on
end

% for i = 1:length(up_grid)-1
%     clear X Y
%     bin_width = (up_grid(i+1)-up_grid(i));
%     X = [up_grid(i) up_grid(i)+bin_width];
%     Y = [adj_uc_up_hist(i) adj_uc_up_hist(i)];
%     X = [X fliplr(X)];
%     Y = [Y adj_lc_up_hist(i) adj_lc_up_hist(i)];
%     fill(X,Y,'k','EdgeColor','none');
%     hold on
% end

set(gca,'yscale','log')
xlim([0 10])
ylim([1e-4 1e-1])
xlabel('Up State Duration (s)')
ylabel('Percent of Data')



%% for down state durations
pyr_avg_down_hist = nanmean(down_hist(pyr_cells,:));
pyr_uc_down_hist = pyr_avg_down_hist+nanstd(down_hist(pyr_cells,:))/sqrt(length(pyr_cells));
pyr_lc_down_hist = pyr_avg_down_hist-nanstd(down_hist(pyr_cells,:))/sqrt(length(pyr_cells));

pyr8_avg_down_hist = nanmean(down_hist8(pyr_cells,:));
pyr8_uc_down_hist = pyr8_avg_down_hist+nanstd(down_hist8(pyr_cells,:))/sqrt(length(pyr_cells));
pyr8_lc_down_hist = pyr8_avg_down_hist-nanstd(down_hist8(pyr_cells,:))/sqrt(length(pyr_cells));


lec_avg_down_hist = nanmean(down_hist(lec_cells,:));
lec_uc_down_hist = lec_avg_down_hist+nanstd(down_hist(lec_cells,:))/sqrt(length(lec_cells));
lec_lc_down_hist = lec_avg_down_hist-nanstd(down_hist(lec_cells,:))/sqrt(length(lec_cells));

adj_avg_down_hist = nanmean(down_hist(adj_cells,:));
adj_uc_down_hist = adj_avg_down_hist+nanstd(down_hist(adj_cells,:))/sqrt(length(adj_cells));
adj_lc_down_hist = adj_avg_down_hist-nanstd(down_hist(adj_cells,:))/sqrt(length(adj_cells));

figure
stairs(down_grid,pyr_avg_down_hist,'linewidth',2)
hold on
stairs(down_grid,lec_avg_down_hist,'g','linewidth',2)
stairs(down_grid,adj_avg_down_hist,'k','linewidth',2)
stairs(down_grid,pyr8_avg_down_hist,'r','linewidth',2)

legend('MEC','LEC','ADJ','LF8')

for i = 1:length(down_grid)-1
    clear X Y
    bin_width = (down_grid(i+1)-down_grid(i));
    X = [down_grid(i) down_grid(i)+bin_width];
    Y = [pyr_uc_down_hist(i) pyr_uc_down_hist(i)];
    X = [X fliplr(X)];
    Y = [Y pyr_lc_down_hist(i) pyr_lc_down_hist(i)];
    fill(X,Y,'b','EdgeColor','none');
    hold on
end

for i = 1:length(down_grid)-1
    clear X Y
    bin_width = (down_grid(i+1)-down_grid(i));
    X = [down_grid(i) down_grid(i)+bin_width];
    Y = [pyr8_uc_down_hist(i) pyr8_uc_down_hist(i)];
    X = [X fliplr(X)];
    Y = [Y pyr8_lc_down_hist(i) pyr8_lc_down_hist(i)];
    fill(X,Y,'r','EdgeColor','none');
    hold on
end

for i = 1:length(down_grid)-1
    clear X Y
    bin_width = (down_grid(i+1)-down_grid(i));
    X = [down_grid(i) down_grid(i)+bin_width];
    Y = [lec_uc_down_hist(i) lec_uc_down_hist(i)];
    X = [X fliplr(X)];
    Y = [Y lec_lc_down_hist(i) lec_lc_down_hist(i)];
    fill(X,Y,'g','EdgeColor','none');
    hold on
end

for i = 1:length(down_grid)-1
    clear X Y
    bin_width = (down_grid(i+1)-down_grid(i));
    X = [down_grid(i) down_grid(i)+bin_width];
    Y = [adj_uc_down_hist(i) adj_uc_down_hist(i)];
    X = [X fliplr(X)];
    Y = [Y adj_lc_down_hist(i) adj_lc_down_hist(i)];
    fill(X,Y,'k','EdgeColor','none');
    hold on
end

set(gca,'yscale','log')
xlim([0 10])
ylim([3e-4 1e-1])
xlabel('Down State Duration (s)')
ylabel('Percent of Data')
