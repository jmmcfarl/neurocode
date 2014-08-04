clear all

% load C:\WC_Germany\JMM_Analysis_pyr\UDS_dur_run_hist_v2\UDS_dur_data_over_smooth
load C:\WC_Germany\JMM_Analysis_pyr\UDS_dur_run_hist_v2\up_per_data_10_28
load C:\WC_Germany\JMM_Analysis_pyr\dir_tree_update
% load C:\WC_Germany\JMM_Analysis_pyr\time_domain_30\correl_data
% load C:\WC_Germany\JMM_Analysis_pyr\time_domain_30\lag_data

dsf = 8;
Fsd = 2016/dsf;
winsize = 30;
niqf = 2016/2;
hcf = 2/niqf;
lcf = 0.05/niqf;
[b,a] = butter(2,[lcf hcf]);
cnt = 0;
maxlag = round(10*Fsd);
tlags = [-maxlag:maxlag]/Fsd;

pers_thresh = 5;
back_time = 10;
forward_time = 30;
min_dur = 30;
no_pa = [];

for d = 1:length(dir_array)

    cd(dir_array{d})
    pwd

    load used_data lf8 wcv_minus_spike

    lf8 = filtfilt(b,a,lf8);
    wcv = filtfilt(b,a,wcv_minus_spike);

    lf8_d = downsample(lf8,dsf);
    wcv_d = downsample(wcv,dsf);

    lf8_d = zscore(lf8_d);
    wcv_d = zscore(wcv_d);

    total_dur = floor(length(lf8_d)/Fsd);

    %find all wcv up states lasting more than pers_dur seconds
    pers_ups = round(up_trans{d}(find(up_state_dur{d} > pers_thresh))/Fsd);
    pers_ups8 = round(up_trans8{d}(find(up_state_dur8{d} > pers_thresh))/Fsd);

    pers_begs = pers_ups - back_time;
    pers_ends = pers_ups + forward_time;
    pers_begs8= pers_ups8-back_time;
    pers_ends8 = pers_ups8+forward_time;

    pers_secs = [];
    pers_secs8 = [];
    for be = 1:length(pers_begs)
        pers_secs = [pers_secs pers_begs(be):pers_ends(be)];
    end
    for be = 1:length(pers_begs8)
        pers_secs8 = [pers_secs8 pers_begs8(be):pers_ends8(be)];
    end

    pers_secs = unique(pers_secs);
    pers_secs8 = unique(pers_secs8);
    
    pers_secs(pers_secs <= 0) = [];
    pers_secs(pers_secs > total_dur) = [];
    
    pers_secs8(pers_secs8 <= 0) = [];
    pers_secs8(pers_secs8 > total_dur) = [];

    pers_secs = setdiff(pers_secs,pers_secs8);
    npers_secs = setdiff([1:total_dur],pers_secs);
    npers_secs = setdiff(npers_secs,pers_secs8);
    
    pers_measure = zeros(1,total_dur);
    pers_measure(pers_secs) = 1;
    pers_diff = [0 diff(pers_measure)];
    pers_start = find(pers_diff==1);
    pers_stop = find(pers_diff == -1);
    
    npers_measure = zeros(1,total_dur);
    npers_measure(npers_secs) = 1;
    npers_diff = [0 diff(npers_measure)];
    npers_start = find(npers_diff ==1);
    npers_stop = find(npers_diff == -1);
    if npers_measure(1) == 1
        npers_start = [1 npers_start];
    elseif pers_measure(1) == 1
        pers_start = [1 pers_start];
    end
    if npers_measure(end) == 1
        npers_stop = [npers_stop total_dur];
    elseif pers_measure(end) == 1
        pers_stop = [pers_stop total_dur];
    end
    
    pers_dur = pers_stop - pers_start;
    npers_dur = npers_stop - npers_start;
    
    bad_pers = find(pers_dur < min_dur);
    bad_npers = find(npers_dur < min_dur);
    
    pers_start(bad_pers) = [];
    pers_stop(bad_pers) = [];
    npers_start(bad_npers) = [];
    npers_stop(bad_npers) = [];
    
    
% 
%     pers_measure8 = zeros(1,total_dur);
%     pers_measure8(pers_secs8) = 1;
%     pers_diff8 = [0 diff(pers_measure8)];
%     pers_start8 = find(pers_diff8==1);
%     pers_stop8 = find(pers_diff8 == -1);
%         
%     npers_start = pers_stop;
%     npers_stop = pers_start;
% 
%     if ~isempty(pers_start8)
%         npers_start(ismember(npers_start,pers_start8)) = [];
%         npers_stop(ismember(npers_stop,pers_stop8)) = [];
%         npers_start = [npers_start pers_stop8];
%         npers_stop = [npers_stop pers_start8];
%         npers_start = sort(npers_start);
%         npers_stop = sort(npers_stop);
%     end
    
    
    %if there is some persistent activity
    if length(pers_start)>0

%         if pers_measure(1) ==1
%             pers_start = [1 pers_start];
%         elseif pers_start(1) > 30
%             npers_start = [1 npers_start];
%         end
%         if pers_measure(end) == 1
%             pers_stop = [pers_stop total_dur];
%         elseif total_dur-pers_stop(end) > 30
%             npers_stop = [npers_stop total_dur];
%         else
%             npers_start(end) = [];
%         end
% 
%         pers_dur = pers_stop-pers_start;
%         npers_dur = npers_stop-npers_start;
% 
%         bad_npers = find(npers_dur < min_dur);
% %         bad_nstart = npers_start(bad_npers);
% %         bad_nstop = npers_stop(bad_npers);
%         npers_start(bad_npers) = [];
%         npers_stop(bad_npers) = [];
% %         pers_start(ismember(pers_start,bad_nstop)) = [];
% %         pers_stop(ismember(pers_stop,bad_nstart)) = [];
% 
%         bad_pers = find(pers_dur < min_dur);
% %         bad_start = pers_start(bad_pers);
% %         bad_stop = pers_stop(bad_pers);
%         pers_start(bad_pers) = [];
%         pers_stop(bad_pers) = [];
%         
%         %calculate percent persistent activity
        pers_dur = pers_stop - pers_start;
        npers_dur = npers_stop - npers_start;
        
        percent_pers(d) = sum(pers_dur)/total_dur;
        

        %convert to sample index
        pers_start = pers_start*Fsd;
        pers_stop = pers_stop*Fsd;
        npers_start = npers_start*Fsd;
        npers_stop = npers_stop*Fsd;


        %cycle through pers states
        %if there is some persistent activity


        num_pers_states = length(pers_start);
        num_npers_states = length(npers_start);

        %initialize
        pers_wcv = [];
        npers_wcv = [];
        pers_lf8 = [];
        npers_lf8 = [];
        
        %% Cycle through all persistent states and calculate the UTA and
        %% DTA

        for t = 1:num_pers_states

            state_beg = pers_start(t);
            state_end = pers_stop(t);
            pers_wcv = [pers_wcv;wcv_d(state_beg:state_end)];
            pers_lf8 = [pers_lf8;lf8_d(state_beg:state_end)];

        end



        %% Cycle through non-pers states

        %initialize

        for t = 1:num_npers_states

            state_beg = npers_start(t);
            state_end = npers_stop(t);

            npers_wcv = [npers_wcv;wcv_d(state_beg:state_end)];
            npers_lf8 = [npers_lf8;lf8_d(state_beg:state_end)];
            
        end


    else
        disp('no persistent activity')
        no_pa = [no_pa d]
        percent_pers(d) = 0;
    end

    [pers_wdist(d,:),gridv] = gpkde(pers_wcv,-3,[-4 4]);
    [pers_8dist(d,:),gridv] = gpkde(pers_lf8,-3,[-4 4]);
    [npers_wdist(d,:),gridv] = gpkde(npers_wcv,-3,[-4 4]);
    [npers_8dist(d,:),gridv] = gpkde(npers_lf8,-3,[-4 4]);

end

%get rid of cells without PApers_wdist(no_pa,:) = [];
pers_8dist(no_pa,:) = [];
npers_wdist(no_pa,:) = [];
npers_8dist(no_pa,:) = [];

%calculate and plot averages
m_pers_wdist = mean(pers_wdist);
u_pers_wdist = m_pers_wdist+2*std(pers_wdist)/sqrt(14);
l_pers_wdist = m_pers_wdist-2*std(pers_wdist)/sqrt(14);

m_pers_8dist = mean(pers_8dist);
u_pers_8dist = m_pers_8dist+2*std(pers_8dist)/sqrt(14);
l_pers_8dist = m_pers_8dist-2*std(pers_8dist)/sqrt(14);

m_npers_wdist = mean(npers_wdist);
u_npers_wdist = m_npers_wdist+2*std(npers_wdist)/sqrt(14);
l_npers_wdist = m_npers_wdist-2*std(npers_wdist)/sqrt(14);

m_npers_8dist = mean(npers_8dist);
u_npers_8dist = m_npers_8dist+2*std(npers_8dist)/sqrt(14);
l_npers_8dist = m_npers_8dist-2*std(npers_8dist)/sqrt(14);


plot(gridv,m_pers_wdist,'linewidth',2)
hold on
plot(gridv,m_npers_wdist,'r','linewidth',2)
legend('Persistent States','Non-Persistent States')
plot(gridv,u_pers_wdist,'--')
plot(gridv,l_pers_wdist,'--')
plot(gridv,m_npers_wdist,'r','linewidth',2)
plot(gridv,u_npers_wdist,'r--')
plot(gridv,l_npers_wdist,'r--')
xlim([-3.5 3.5])
xlabel('Amplitude (z-score)','FontSize',14)
ylabel('Probability','FontSize',14)


plot(gridv,m_pers_8dist,'linewidth',2)
hold on
plot(gridv,m_npers_8dist,'r','linewidth',2)
legend('Persistent States','Non-Persistent States')
plot(gridv,u_pers_8dist,'--')
plot(gridv,l_pers_8dist,'--')
plot(gridv,m_npers_8dist,'r','linewidth',2)
plot(gridv,u_npers_8dist,'r--')
plot(gridv,l_npers_8dist,'r--')
xlim([-3.5 3.5])
xlabel('Amplitude (z-score)','FontSize',14)
ylabel('Probability','FontSize',14)



%hist percent of time in persistent statee
percent_pers = percent_pers*100;
range = 0:5:100;
n = hist(percent_pers,range);
bar(range,n)
xlabel('Percent Persistent Activity','FontSize',14)
ylabel('Number of Cells','FontSize',14)
xlim([-5 100])