clear all
close all

thresholds = [-0.6 -0.58 -0.6 -0.55 -0.61 -0.54 -0.61 nan -0.59 nan nan nan nan -0.6 nan -0.58 -0.62 nan];


load ('C:\WC_Germany\Layer5\l5_heka_dir')
for cdir =16
    load(f_loc{cdir})
    dat_name = f_loc{cdir}(39:end);

    eval(['data = ' dat_name '_MP;'])
    eval(['time = ' dat_name '_sampletimes;'])
    Fs = 5000;
    dsf = 10;
    Fsd = Fs/dsf;
    niqf = Fs/2;
    [b,a] = butter(2,3/niqf,'low');
    [b2,a2] = butter(2,6/niqf,'low');

    sweep_t = [0;find(diff(time) < 0)];
    low = min(data);
    high = max(data)-0.3;
    dataseg = data(1:sweep_t(1));
    filtseg = filtfilt(b,a,dataseg);

    maxnoise = 0.04;
    mindowndur = 0.75; %min dur in sec
    minupdur = 1.0;
    minpeakheight = 0.2;

    down_cnt = 0;
    up_cnt = 0;

    for i = 1:length(sweep_t)-1

        dataseg = data(sweep_t(i)+1:sweep_t(i+1));
        dataseg(dataseg > -0.3) = -0.3;
        filtseg2 = filtfilt(b2,a2,dataseg);

        if isnan(thresholds(cdir))
            down_state_array{i} = dataseg;
        else

%             [cur_down_start, cur_down_stop] = down_state_loc_bim(dataseg,thresholds(cdir));
%                         [cur_up_start, cur_up_stop] = up_state_loc_bim(dataseg,thresholds(cdir));


            time_vec = time(sweep_t(i)+1:sweep_t(i+1));
            %
%             if ~isempty(cur_down_start) & ~isempty(cur_down_start)
%                 num_downs = length(cur_down_stop);
%                 if num_downs > 0
%                     for d = 1:num_downs
%                         down_cnt = down_cnt+1;
%                         down_state_array{down_cnt} = dataseg(cur_down_start(d):cur_down_stop(d));
%                     end
%                 end
%             end
            %             if ~isempty(cur_up_start) & ~isempty(cur_up_stop)
            %                 num_ups = length(cur_up_stop);
            %                 if num_ups > 0
            %                     for u = 1:num_ups
            %                         up_cnt = up_cnt + 1;
            %                         up_state_array{up_cnt} = dataseg(cur_up_start(u):cur_up_stop(u));
            %                     end
            %                 end
            %             end
            %
        end
        time_vec = time(sweep_t(i)+1:sweep_t(i+1));


                plot(time_vec,data(sweep_t(i)+1:sweep_t(i+1)))
                hold on
%                 plot(time_vec,filtseg2,'r','linewidth',2)
%                 hold on
%                 plot(time_vec(cur_down_start),dataseg(cur_down_start),'go','markersize',12,'linewidth',2)
%                 plot(time_vec(cur_down_stop),dataseg(cur_down_stop),'ro','markersize',12,'linewidth',2)
%                           plot(time_vec(cur_up_start),dataseg(cur_up_start),'g*','markersize',12,'linewidth',2)
%                           plot(time_vec(cur_up_stop),dataseg(cur_up_stop),'r*','markersize',12,'linewidth',2)
        
                ylim([low high])
        
                grid
                pause
                clf


    end


%     maxlag = 0.75*Fs;
%     lags = -maxlag:maxlag;
% 
%     %% now go through down states and create xcov matrix
%     niqf = Fs/2;
%     [bt,at] = butter(2,[2/niqf 40/niqf]);
%     xcov_mat = zeros(length(down_state_array),length(lags));
%     for i = 1:length(down_state_array)
%         mean_mp(i) = mean(down_state_array{i});
%         xcov_mat(i,:) = xcov(filtfilt(bt,at,down_state_array{i}),maxlag,'coeff');
%     end
%     %     [dummy,mp_order] = sort(mean_mp);
%     figure
%     pcolor(lags/Fs,1:length(down_state_array),xcov_mat);
%     shading flat
%     caxis([-0.3 0.4])
%     colorbar
%     xlim([0 maxlag/Fs])
%     t_names = ['C:\WC_Germany\down_state_cov_' f_names{cdir}];
%     print('-dpng',t_names);
%     close
%     %
%         clear down_state_array mean_mp
    % if ~isnan(thresholds(cdir))
    %     maxlag = 1.0*Fs;
    %     lags = -maxlag:maxlag;
    %
    %     %% now go through down states and create xcov matrix
    %     niqf = Fs/2;
    %     [bt,at] = butter(2,[1.5/niqf 40/niqf]);
    %     xcov_mat = zeros(length(up_state_array),length(lags));
    %     for i = 1:length(up_state_array)
    %         mean_mp(i) = mean(up_state_array{i});
    %         xcov_mat(i,:) = xcov(filtfilt(bt,at,up_state_array{i}),maxlag,'coeff');
    %     end
    %     [dummy,mp_order] = sort(mean_mp);
    %     figure
    %     pcolor(lags/Fs,1:length(up_state_array),xcov_mat);
    %     shading flat
    %     caxis([-0.3 0.4])
    %     colorbar
    %     xlim([0 maxlag/Fs])
    %     t_names = ['C:\WC_Germany\up_state_cov_' f_names{cdir}];
    %     print('-dpng',t_names);
    %     close
    % end
    clear up_state_array

end