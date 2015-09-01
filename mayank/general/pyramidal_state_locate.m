clear all
close all

% load('C:\WC_Germany\JMM_analysis_ste\stellate_heka_dir.mat')
load('C:\WC_Germany\JMM_analysis_pyr\pyr_heka_dir.mat')

thresholds = [-0.6 -0.55 -0.6 -0.61 -0.55 -0.6 -0.6 -0.55 -0.6 -0.65 -0.55 -0.61 -0.55 -0.5 -0.52 -0.55 -0.55];

hist_range = linspace(0,1000,100);

for cdir =1:17
    load(f_loc{cdir})
    dat_name = f_loc{cdir}(25:end);

    eval(['data = ' dat_name '_MP;'])
    eval(['time = ' dat_name '_sampletimes;'])
    Fs = 5000;
    dsf = 10;
    Fsd = Fs/dsf;
    niqf = Fs/2;
    [b,a] = butter(2,[3/niqf 20/niqf]);
    [b2,a2] = butter(2,[7/niqf 20/niqf]);

    sweep_t = [0;find(diff(time) < 0)];
    low = min(data);
    high = max(data)-0.3;
    dataseg = data(1:sweep_t(1));
    filtseg = filtfilt(b,a,dataseg);

    maxnoise = 0.1;
    mindowndur = 0.75; %min dur in sec
    minupdur = 1.0;
    minpeakheight = 0.2;

    threshold = 20;
    
    down_cnt = 0;
    up_cnt = 0;

    tot_ipspi{cdir} = [];
    tot_up_ipspi{cdir} = [];
    tot_down_ipspi{cdir} = [];
    
    
    for i = 1:length(sweep_t)-1

        dataseg = data(sweep_t(i)+1:sweep_t(i+1));
%         dataseg(dataseg > -0.3) = -0.3;
        filtseg = filter(b,a,dataseg);
        % filtseg2 = filtfilt(b2,a2,dataseg);
        [cur_down_start, cur_down_stop] = down_state_loc_bim(dataseg,thresholds(cdir));
        [cur_up_start, cur_up_stop] = up_state_loc_bim(dataseg,thresholds(cdir));
% 
        if ~isempty(cur_down_start) & ~isempty(cur_down_start)
            num_downs = length(cur_down_stop);
            if num_downs > 0
                for d = 1:num_downs
                    down_cnt = down_cnt+1;
                    down_state_array{down_cnt} = dataseg(cur_down_start(d):cur_down_stop(d));
                end
            end
        end
%                 if ~isempty(cur_up_start) & ~isempty(cur_up_stop)
%             num_ups = length(cur_up_stop);
%             if num_ups > 0
%                 for u = 1:num_ups
%                     up_cnt = up_cnt + 1;
%                     up_state_array{up_cnt} = dataseg(cur_up_start(u):cur_up_stop(u));
%                 end
%             end
%         end
        time_vec = time(sweep_t(i)+1:sweep_t(i+1));
        [psp_times,psp_sign,psp_det] = psp_detect(dataseg,threshold);
        tot_ipspi{cdir} = [tot_ipspi{cdir} diff(psp_times)];
        cur_up_psp_times = [];
        cur_down_psp_times = [];
        if ~isempty(cur_up_start) 
           for cu = 1:length(cur_up_start)
               [up_psp_times,temp_psp_sign,temp_psp_det] = ...
                   psp_detect(dataseg(cur_up_start(cu):cur_up_stop(cu)),threshold);
               cur_up_psp_times = [cur_up_psp_times (cur_up_start(cu)+up_psp_times)];
               tot_up_ipspi{cdir} = [tot_up_ipspi{cdir} diff(up_psp_times)];
           end
        end
        if ~isempty(cur_down_start) 
            for cd = 1:length(cur_down_start)
                [down_psp_times,temp_psp_sign,temp_psp_det] = ...
                    psp_detect(dataseg(cur_down_start(cd):cur_down_stop(cd)),threshold);
                cur_down_psp_times = [cur_down_psp_times (cur_down_start(cd)+down_psp_times)];
                tot_down_ipspi{cdir} = [tot_down_ipspi{cdir} diff(down_psp_times)];
            end
        end

%                plot(time_vec,dataseg)
%                hold on
% %                   plot(time_vec,filtseg+mean(dataseg),'r','linewidth',2)
%                   hold on
% %                   plot(time_vec(thresh_cross_up),filtseg(thresh_cross_up),'ko')
% %                   plot(time_vec(thresh_cross_down),filtseg(thresh_cross_down),'mo')
% %                   plot(time_vec(cur_down_start),dataseg(cur_down_start),'go')
% %                   plot(time_vec(cur_down_stop),dataseg(cur_down_stop),'ro')
%                   plot(time_vec(cur_up_psp_times),dataseg(cur_up_psp_times),'g.')
%                   plot(time_vec(cur_down_psp_times),dataseg(cur_down_psp_times),'r.')
% %                   plot(time_vec(cur_up_start),dataseg(cur_up_start),'g*')
% %                   plot(time_vec(cur_up_stop),dataseg(cur_up_stop),'r*')
%         
%                ylim([low high])
%             xlabel('Time (s)')
%             ylabel('MP (mV * 0.01)')
% 
%            grid
%            pause
%            clf

    end


    maxlag = 0.75*Fs;
    lags = -maxlag:maxlag;
    
    % now go through down states and create xcov matrix
%     niqf = Fs/2;
%     [bt,at] = butter(2,[2/niqf 40/niqf]);
%     xcov_mat = zeros(length(down_state_array),length(lags));
%     for i = 1:length(down_state_array)
%         mean_mp(i) = mean(down_state_array{i});
%         xcov_mat(i,:) = xcov(filtfilt(bt,at,down_state_array{i}),maxlag,'coeff');
%     end
%     [dummy,mp_order] = sort(mean_mp);
%     figure
%     pcolor(lags/Fs,1:length(down_state_array),xcov_mat);
%     shading flat
%     caxis([-0.3 0.4])
%     colorbar
%     xlim([0 maxlag/Fs])
%        t_names = ['C:\WC_Germany\down_state_cov_' f_names{cdir}];
%     print('-dpng',t_names);
%     close 
    
%      maxlag = 1.0*Fs;
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


    
clear up_state_array down_state_array mean_mp

%     %% now go through up states and create xcov matrix
%     niqf = Fs/2;
%     [b,a] = butter(2,2/niqf,'high');
%     xcov_mat_u = zeros(length(up_state_array),length(lags));
%     for i = 1:length(up_state_array)
%         mean_mp(i) = mean(up_state_array{i});
%         xcov_mat_u(i,:) = xcov(filtfilt(b,a,up_state_array{i}),maxlag,'coeff');
%     end
%     [dummy,mp_order] = sort(mean_mp);
%     figure
%     pcolor(lags/Fs,1:length(up_state_array),xcov_mat_u);
%     shading flat
%     colorbar
%     xlim([0 maxlag/Fs])
    
clear down_state_array cur_down_start cur_up_start mean_mp 

hist(tot_up_ipspi{cdir}/5,hist_range)
t_names = ['C:\WC_Germany\inter_psp\up_' f_names{cdir}];
xlim([0 1000])
print('-dpng',t_names);
close

hist(tot_down_ipspi{cdir}/5,hist_range)
t_names = ['C:\WC_Germany\inter_psp\down_' f_names{cdir}];
xlim([0 1000])
print('-dpng',t_names);
close

end


