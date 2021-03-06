clear all
close all

cd ~/Data/bruce/2_27_12/
load Blocks.mat


all_avg_lpos = [];
all_sac_times = [];
all_avg_rpos = [];
all_avg_disp = [];
all_std_disp = [];
all_std_lpos = [];
all_block_times = zeros(4,1);
cur_abs_offset = 0;
% block_id = 1;

thresh_eyespeed = 10;
for block_id = 1:4
    
    fprintf('Block %d of 4\n',block_id);
    block_times = Blocks{block_id}.blocktimes;
    stim_times = Blocks{block_id}.stimtime;
    stimids = Blocks{block_id}.stimids;
    is_stim_filtered = mod(stimids-1,4) >= 2;
    n_stims = length(stim_times);
    Pix2Deg = 1.1279 / 60;
    dsfrac = 2;
    
    cd ~/Data/bruce/2_27_12/saccades/
    load(sprintf('lemM232.5%d.em.sac.mat',block_id))
    sac_times = Expt.Trials.EyeMsac.sacT;
    
    
    %%
    n_units = length(Blocks{block_id}.spktimes);
    
    stimres = 0.03; %resolution of time bins for STA calc
    sac_win = [-0.1 0.2]; %time window around saccades to reject data
    sac_win_bins = round(sac_win/stimres);
    
    target_window = [-8 8; -8 8]; %acceptance window (in degrees) for gaze direction.
    target_window_pix = target_window/Pix2Deg/dsfrac;
    
    shift_time = -0.06;
    % delta_bin = round(-0.06/stimres); %number of time bins to shift stim relative to spike data
    
    %%
    reyepos = [Expt.Trials.Eyevals.rh Expt.Trials.Eyevals.rv];% right eye
    leyepos = [Expt.Trials.Eyevals.lh Expt.Trials.Eyevals.lv];%left eye
    
    clear disparity
    disparity(:,1) = leyepos(:,1)-reyepos(:,1);
    disparity(:,2) = leyepos(:,2)-reyepos(:,2);
    disparity_mag = sqrt(disparity(:,1).^2+disparity(:,2).^2);
    
    Eyedt = Expt.Header.CRrates(1); % resolution of eye signal (sec)
    eye_tx = Expt.Trials.Start/1e4:Eyedt:Expt.Trials.End/1e4;
    eye_tx(end) = []; %throw out last bin so that aligns to eye data vector
    
    clear sm_avg_eyepos eye_vel
    avg_eyepos = (reyepos + leyepos)/2;
    sm_avg_eyepos(:,1) = smooth(avg_eyepos(:,1),3);
    sm_avg_eyepos(:,2) = smooth(avg_eyepos(:,2),3);
    eye_vel(:,1) = [0; diff(sm_avg_eyepos(:,1))]/Eyedt;
    eye_vel(:,2) = [0; diff(sm_avg_eyepos(:,2))]/Eyedt;
    eye_speed = sqrt(eye_vel(:,1).^2+eye_vel(:,2).^2);
    
    time_axis = block_times(1,1):stimres:block_times(2,end);
    
    sac_start_times = eye_tx(eye_speed(1:end-1) < thresh_eyespeed & eye_speed(2:end) > thresh_eyespeed);
    sac_stop_times = eye_tx(eye_speed(1:end-1) > thresh_eyespeed & eye_speed(2:end) < thresh_eyespeed);
    
    sac_buffer_back = round(0.025/Eyedt);
    sac_buffer_forward = round(0.025/Eyedt);
    avg_disparity = nan(length(sac_times),2);
    std_disparity = nan(length(sac_times),2);
    avg_lpos = nan(length(sac_times),2);
    std_lpos = nan(length(sac_times),2);
    for i = 1:length(sac_times)-1
        cur_win = find(eye_tx >= sac_stop_times(i) & eye_tx < sac_start_times(i+1));
        if length(cur_win) > sac_buffer_back+sac_buffer_forward
            cur_win(1:sac_buffer_back) = [];
            cur_win(end-sac_buffer_forward+1:end) = [];
% cur_win(eye_speed(cur_win) > thresh_eyespeed) = [];
            avg_disparity(i,:) = mean(disparity(cur_win,:));
            std_disparity(i,:) = std(disparity(cur_win,:));
            avg_lpos(i,:) = mean(leyepos(cur_win,:));
            std_lpos(i,:) = std(leyepos(cur_win,:));
            avg_rpos(i,:) = mean(reyepos(cur_win,:));
            
            subplot(2,1,1)
            plot(eye_tx,leyepos)
            hold on
            plot(eye_tx(cur_win),leyepos(cur_win,:),'r.')
            xlim([eye_tx(cur_win(1)) - 1 eye_tx(cur_win(end)) + 1])
            subplot(2,1,2)
            plot(eye_tx,disparity)
            hold on
            plot(eye_tx(cur_win),disparity(cur_win,:),'r.')
            xlim([eye_tx(cur_win(1)) - 1 eye_tx(cur_win(end)) + 1])
            pause(1)
            clf
            
        end
    end
    
    used_sacs = find(avg_lpos(:,1) > target_window(1,1) & avg_lpos(:,1) < target_window(1,2) & ...
        avg_lpos(:,2) > target_window(2,1) & avg_lpos(:,2) < target_window(2,2));
    
    all_avg_lpos = [all_avg_lpos; avg_lpos(used_sacs,:)];
    all_sac_times = [all_sac_times sac_times(used_sacs)+cur_abs_offset];
    all_avg_rpos = [all_avg_rpos; avg_rpos(used_sacs,:)];
    all_avg_disp = [all_avg_disp; avg_disparity(used_sacs,:)];
    all_std_disp = [all_std_disp; std_disparity(used_sacs,:)];
    all_std_lpos = [all_std_lpos; std_lpos(used_sacs,:)];
    
    cur_abs_offset = cur_abs_offset + max(sac_times(used_sacs)) + 100;
    
end

%%
smooth_win = 5;
plot(all_sac_times,all_avg_disp(:,1),'.')
hold on
plot(all_sac_times,all_avg_disp(:,2),'r.')
ylim([-2 2]);
yl = ylim();
plot(all_sac_times,smooth(all_avg_disp(:,1),smooth_win,'rlowess'),'k','linewidth',2)
plot(all_sac_times,smooth(all_avg_disp(:,2),smooth_win,'rlowess'),'k','linewidth',2)

set(gca,'fontsize',14)
xlabel('Time (s)','fontsize',16)
ylabel('Disparity (degrees)','fontsize',16)

%%
smooth_win = 10;
plot(all_sac_times,all_std_disp(:,1),'.')
hold on
plot(all_sac_times,all_std_disp(:,2),'r.')
ylim([-2 2]);
yl = ylim();
plot(all_sac_times,smooth(all_std_disp(:,1),smooth_win,'rlowess'),'k','linewidth',2)
plot(all_sac_times,smooth(all_std_disp(:,2),smooth_win,'rlowess'),'k','linewidth',2)

set(gca,'fontsize',14)
xlabel('Time (s)','fontsize',16)
ylabel('Disparity (degrees)','fontsize',16)
%%
close all
smooth_win = 40;
[sorted_pos,pos_ord] = sort(all_avg_lpos(:,2));
plot(sorted_pos,all_avg_disp(pos_ord,2),'.')
hold on
plot(sorted_pos,smooth(all_avg_disp(pos_ord,2),smooth_win),'k','linewidth',2)

%%
abs_disp = sqrt(all_avg_disp(:,1).^2+all_avg_disp(:,2).^2);
% close all
smooth_win = 40;
[sorted_pos,pos_ord] = sort(all_avg_lpos(:,1));
plot(sorted_pos,abs_disp,'.')
hold on
plot(sorted_pos,smooth(abs_disp,smooth_win),'k','linewidth',2)

ylim([0 2])
% xlabel('Horizontal Eye position (degrees)','fontsize',16)
xlabel('Vertical Eye position (degrees)','fontsize',16)
ylabel('Total disparity (degrees)','fontsize',16)