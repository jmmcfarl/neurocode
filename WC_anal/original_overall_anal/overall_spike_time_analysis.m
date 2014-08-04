clear all
close all

load('C:\WC_Germany\overall_calcs\overall_dir.mat')
load C:\WC_Germany\overall_calcs\spk_locs\spk_id_data
load C:\WC_Germany\overall_calcs\UDS_dur_raw\UDS_raw_data
load C:\WC_Germany\overall_calcs\UDS_synch_state_dur\UDS_synch_state_dur_data

dsf = 16;
Fsd = 2016;
markerSize =4;
spk_win = 5;
sm_win = 10;

for d = 1:length(over_dir)

    cur_spike_ids = round(spk_ids{d}/dsf);

    cur_up_trans = up_trans{d}(synch_ups{d})*8;
    cur_down_trans = down_trans{d}(synch_ups{d})*8;
    cur_up_trans8 = up_trans8{d}(synch_ups8{d})*8;
    cur_down_trans8 = down_trans8{d}(synch_ups8{d})*8;



    up_spike_mat = cell(1,length(cur_up_trans8));
    burst_spike_mat = cell(1,length(cur_up_trans8));
    thi_spike_mat = cell(1,length(cur_up_trans8));
    down_spike_mat = cell(1,length(cur_up_trans8));
    dburst_spike_mat = cell(1,length(cur_up_trans8));
    dthi_spike_mat = cell(1,length(cur_up_trans8));

    num_good_spikes = 0;
    tot_good_up_dur = 0;
    num_good_bursts = 0;
    num_good_tris = 0;

    for i = 1:length(cur_up_trans8)

        spike_row = cur_spike_ids(abs(cur_spike_ids - cur_up_trans8(i)) < spk_win*Fsd);
        if ~isempty(spike_row)
            cur_isis = [0 diff(spike_row)]/Fsd;
            up_spike_mat{i} = spike_row - cur_up_trans8(i);
            burst_spike_mat{i} = spike_row(cur_isis < 0.05) - cur_up_trans8(i);

            num_good_spikes = num_good_spikes + length(up_spike_mat{i});
            tot_good_up_dur = tot_good_up_dur + (cur_down_trans8(i)-cur_up_trans8(i))/Fsd;
            num_good_bursts = num_good_bursts + length(burst_spike_mat{i});

            burst_spike_mat{i}(1) = [];
            if ~isempty(burst_spike_mat{i})
                cur_isis = [0 diff(burst_spike_mat{i})]/Fsd;
                thi_spike_mat{i} = burst_spike_mat{i}(cur_isis < 0.05);
                thi_spike_mat{i}(1) = [];
            else
                thi_spike_mat{i} = [];
            end

            num_good_tris = num_good_tris + length(thi_spike_mat{i});


            spike_row = cur_spike_ids(abs(cur_spike_ids - cur_down_trans8(i)) < spk_win*Fsd);
            if ~isempty(spike_row)
                cur_isis = [0 diff(spike_row)]/Fsd;

                down_spike_mat{i} = spike_row - cur_down_trans8(i);
                dburst_spike_mat{i} = spike_row(cur_isis < 0.05) - cur_down_trans8(i);
                dburst_spike_mat{i}(1) = [];
                if ~isempty(dburst_spike_mat{i})
                    cur_isis = [0 diff(dburst_spike_mat{i})]/Fsd;
                    dthi_spike_mat{i} = dburst_spike_mat{i}(cur_isis < 0.05);
                    dthi_spike_mat{i}(1) = [];
                else
                    dthi_spike_mat{i} = [];
                end
            end

        end

    end


    [up_dur,up_order] = sort((cur_down_trans8-cur_up_trans8)/Fsd);

    Fig = figure(1)
    clf
    set(Fig,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])

    for i = 1:length(cur_up_trans8)

        plot(up_spike_mat{up_order(i)}/Fsd,ones(size(up_spike_mat{up_order(i)}))*i,'.','MarkerSize',4)
        hold on
        plot(burst_spike_mat{up_order(i)}/Fsd,ones(size(burst_spike_mat{up_order(i)}))*i,'r.','MarkerSize',4)
        plot(thi_spike_mat{up_order(i)}/Fsd,ones(size(thi_spike_mat{up_order(i)}))*i,'g.','MarkerSize',4)


    end
    plot(up_dur,1:length(cur_up_trans8),'k')
    xlim([-4 4])
    line([0 0],[0 length(cur_up_trans8)],'Color','k')
    t_names = ['C:\WC_Germany\overall_calcs\spike_time\up_trig8_rast' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',t_names);
    close


    Fig = figure(1)
    clf
    set(Fig,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    for i = 1:length(cur_up_trans8)

        plot(down_spike_mat{up_order(i)}/Fsd,ones(size(down_spike_mat{up_order(i)}))*i,'.','MarkerSize',4)
        hold on
        plot(dburst_spike_mat{up_order(i)}/Fsd,ones(size(dburst_spike_mat{up_order(i)}))*i,'r.','MarkerSize',4)
        plot(dthi_spike_mat{up_order(i)}/Fsd,ones(size(dthi_spike_mat{up_order(i)}))*i,'g.','MarkerSize',4)

    end
    plot(-up_dur,1:length(cur_up_trans8),'k')
    xlim([-3 3])
    line([0 0],[0 length(cur_up_trans8)],'Color','k')

    t_names = ['C:\WC_Germany\overall_calcs\spike_time\down_trig8_rast' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',t_names);
    close


    %   mean_up_rate(d) = num_good_spikes/tot_good_up_dur;
    %   mean_burst_rate(d) = num_good_bursts/tot_good_up_dur;
    %   mean_tri_rate(d) = num_good_tris/tot_good_up_dur;


    rel_up_timing = -3:0.001:3;
    up_avg_mat = nan(length(cur_up_trans8),length(rel_up_timing));
    burst_avg_mat = nan(length(cur_up_trans8),length(rel_up_timing));
    thi_avg_mat = nan(length(cur_up_trans8),length(rel_up_timing));

    for i = 1:length(cur_up_trans8)
        temp_hist = hist(up_spike_mat{i}/Fsd,rel_up_timing);
        end_pt = find(rel_up_timing > (cur_down_trans8(i)-cur_up_trans8(i))/Fsd,1,'first');
        up_avg_mat(i,1:end_pt) = temp_hist(1:end_pt);
        temp_hist = hist(burst_spike_mat{i}/Fsd,rel_up_timing);
        burst_avg_mat(i,1:end_pt) = temp_hist(1:end_pt);
        temp_hist = hist(thi_spike_mat{i}/Fsd,rel_up_timing);
        thi_avg_mat(i,1:end_pt) = temp_hist(1:end_pt);
    end
    up_avg_rate = nanmean(up_avg_mat)*1000;
    up_avg_rate = jmm_smooth_1d(up_avg_rate,sm_win);
    burst_avg_rate = nanmean(burst_avg_mat)*1000;
    burst_avg_rate = jmm_smooth_1d(burst_avg_rate,2*sm_win);
    thi_avg_rate = nanmean(thi_avg_mat)*1000;
    thi_avg_rate = jmm_smooth_1d(thi_avg_rate,3*sm_win);



    rel_down_timing = -3:0.001:3;
    down_avg_mat = nan(length(cur_down_trans8),length(rel_down_timing));
    dburst_avg_mat = nan(length(cur_down_trans8),length(rel_down_timing));
    dthi_avg_mat = nan(length(cur_down_trans8),length(rel_down_timing));

    for i = 1:length(cur_down_trans8)
        temp_hist = hist(-down_spike_mat{i}/Fsd,rel_down_timing);
        end_pt = find(rel_down_timing > (cur_down_trans8(i)-cur_up_trans8(i))/Fsd,1,'first');
        down_avg_mat(i,1:end_pt) = temp_hist(1:end_pt);
        temp_hist = hist(-dburst_spike_mat{i}/Fsd,rel_down_timing);
        dburst_avg_mat(i,1:end_pt) = temp_hist(1:end_pt);
        temp_hist = hist(-dthi_spike_mat{i}/Fsd,rel_down_timing);
        dthi_avg_mat(i,1:end_pt) = temp_hist(1:end_pt);
    end
    down_avg_rate = nanmean(down_avg_mat)*1000;
    down_avg_rate = jmm_smooth_1d(down_avg_rate,sm_win);
    dburst_avg_rate = nanmean(dburst_avg_mat)*1000;
    dburst_avg_rate = jmm_smooth_1d(dburst_avg_rate,2*sm_win);
    dthi_avg_rate = nanmean(dthi_avg_mat)*1000;
    dthi_avg_rate = jmm_smooth_1d(dthi_avg_rate,3*sm_win);

    figure
    subplot(2,1,1)
    plot(rel_up_timing,up_avg_rate,'linewidth',2)
    xlim([-2.8 2.8])
    subplot(2,1,2)
    plot(-rel_down_timing,down_avg_rate,'linewidth',2)
    xlim([-2.8 2.8])
    t_names = ['C:\WC_Germany\overall_calcs\spike_time\avg_trig8_rate_' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',t_names);
    close

    %   figure
    %   subplot(2,1,1)
    %   plot(rel_up_timing,burst_avg_rate,'r','linewidth',2)
    %   hold on
    %   plot(rel_up_timing,thi_avg_rate/mean_tri_rate(d),'g','linewidth',2)
    %   xlim([0 3])
    %   subplot(2,1,2)
    %   plot(-rel_down_timing,dburst_avg_rate/mean_burst_rate(d),'r','linewidth',2)
    %   hold on
    %   plot(-rel_down_timing,dthi_avg_rate/mean_tri_rate(d),'g','linewidth',2)
    %   xlim([-3 0])
    %   t_names = ['C:\WC_Germany\overall_calcs\spike_time\avg_trig_burst_' num2str(cell_type(d)) '_' over_names{d}];
    %   print('-dpng',t_names);
    %   close


end