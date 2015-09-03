clear all

load F:\WC_Germany\overall_EC\overall_allcells_dir

addpath('F:\Code\WC_anal\general')
addpath('F:\WC_Germany\hsmm_state_detection\')
drive_letter = 'F';
Fsd = 2016/8;

% Calculate distributions of up and down states
up_range = [0.1 10];
down_range = [0.1 10];
numBins = 60;

lin_grid = linspace(up_range(1),up_range(2),numBins);
for d = 1:length(sess_data)
    
    fprintf('session %d\n',d)
    cdir = sess_data(d).directory;
    cdir(1) = drive_letter;
    cd(cdir)
    s_name = strcat(sess_data(d).region,'_l',sess_data(d).layer,'_',sess_data(d).name);
    
    load used_data lf8
    load ec_hmm_state_seq
    load ec_hmm_state_seq8
    
    mp_state_seq = hmm_bbstate_seq;
    lf8_state_seq = hmm_bbstate_seq8;
    
    [mp_state_durations{d}] = compute_state_durations_seg(mp_state_seq,Fsd);
    [lf8_state_durations{d}] = compute_state_durations_seg(lf8_state_seq,Fsd);
    
    if exist('ec_hmm_state_seq2r.mat')
        load ec_hmm_state_seq2r
        lf2r_state_seq = hmm_bbstate_seq2r;
        [lf2r_state_durations{d}] = compute_state_durations_seg(lf2r_state_seq,Fsd);
    end
    
%     [lup_hist(d,:),lup_grid] = log_histc_norm(mp_state_durations{d}{2},up_range,numBins);
%     [lup_hist8(d,:),lup_grid] = log_histc_norm(lfp_state_durations{d}{2},up_range,numBins);
%     
%     [ldown_hist(d,:),ldown_grid] = log_histc_norm(mp_state_durations{d}{1},down_range,numBins);
%     [ldown_hist8(d,:),ldown_grid] = log_histc_norm(lfp_state_durations{d}{1},down_range,numBins);
%  
%     up_hist(d,:) = histc(mp_state_durations{d}{2},lin_grid);
%     up_hist8(d,:) = histc(lfp_state_durations{d}{2},lin_grid);
%     up_hist(d,:) = up_hist(d,:)/sum(up_hist(d,:));
%     up_hist8(d,:) = up_hist8(d,:)/sum(up_hist8(d,:));
%     
%     down_hist(d,:) = histc(mp_state_durations{d}{1},lin_grid);
%     down_hist8(d,:) = histc(lfp_state_durations{d}{1},lin_grid);
%     down_hist(d,:) = down_hist(d,:)/sum(down_hist(d,:));
%     down_hist8(d,:) = down_hist8(d,:)/sum(down_hist8(d,:));
    
    %%
%     fig = figure('visible','off');
%     stairs(lup_grid,lup_hist(d,:),'linewidth',2)
%     set(gca,'yscale','log')
%     hold on
%     stairs(lup_grid,lup_hist8(d,:),'r','linewidth',2)
%     grid
%     title('Up State Duration','FontSize',14)
%     legend('MP','LFP')
%     ylim([1e-3 0.2])
%     xlim([0.2 10])
%     t_names = ['G:\WC_Germany\overall_EC\state_durs\uplog_' s_name];
%     print('-dpng',t_names);
%     close
%     
%      fig = figure('visible','off');
%         stairs(lin_grid,up_hist(d,:),'linewidth',2)
%     set(gca,'yscale','log')
%     hold on
%     stairs(lin_grid,up_hist8(d,:),'r','linewidth',2)
%     grid
%     title('Up State Duration','FontSize',14)
%     legend('MP','LFP')
%     ylim([1e-3 0.2])
%     xlim([0.2 10])
%     t_names = ['G:\WC_Germany\overall_EC\state_durs\uplin_' s_name];
%     print('-dpng',t_names);
%     close
% 
%      fig = figure('visible','off');
%     stairs(ldown_grid,ldown_hist(d,:),'linewidth',2)
%     set(gca,'yscale','log')
%     hold on
%     stairs(ldown_grid,ldown_hist8(d,:),'r','linewidth',2)
%     grid
%     title('Down State Duration','FontSize',14)
%     legend('MP','LFP')
%     ylim([1e-3 0.2])
%     xlim([0.2 10])
%     t_names = ['G:\WC_Germany\overall_EC\state_durs\downlog' s_name];
%     print('-dpng',t_names);
%     close
% 
%       fig = figure('visible','off');
%        stairs(lin_grid,down_hist(d,:),'linewidth',2)
%     set(gca,'yscale','log')
%     hold on
%     stairs(lin_grid,down_hist8(d,:),'r','linewidth',2)
%     grid
%     title('Down State Duration','FontSize',14)
%     legend('MP','LFP')
%     ylim([1e-3 0.2])
%     xlim([0.2 10])
%     t_names = ['G:\WC_Germany\overall_EC\state_durs\downlin' s_name];
%     print('-dpng',t_names);
%     close

    %%
%     net_up_time(d) = sum(mp_state_durations{d}{2});
%     net_down_time(d) = sum(mp_state_durations{d}{1});
%     net_up_time8(d) = sum(lfp_state_durations{d}{2});
%     net_down_time8(d) = sum(lfp_state_durations{d}{1});
%     
%     mean_up_dur(d) = mean(mp_state_durations{d}{2});
%     mean_down_dur(d) = mean(mp_state_durations{d}{1});
%     mean_up_dur8(d) = mean(lfp_state_durations{d}{2});
%     mean_down_dur8(d) = mean(lfp_state_durations{d}{1});
%  
%     max_up_dur(d) = max(mp_state_durations{d}{2});
%     max_down_dur(d) = max(mp_state_durations{d}{1});
%     max_up_dur8(d) = max(lfp_state_durations{d}{2});
%     max_down_dur8(d) = max(lfp_state_durations{d}{1});
% 
%     total_number_up_states(d) = length(mp_state_durations{d}{2});
%     total_number_up_states8(d) = length(lfp_state_durations{d}{2});
%     total_number_down_states(d) = length(mp_state_durations{d}{1});
%     total_number_down_states8(d) = length(lfp_state_durations{d}{1});
    
end


save F:\WC_Germany\overall_EC\state_dur_stats mp_* lf8_* lf2r_*