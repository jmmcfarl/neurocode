clear all

load C:\WC_Germany\JMM_Analysis_pyr\UDS_dur_run_hist_v2\UDS_dur_data_over_smooth
load C:\WC_Germany\JMM_Analysis_pyr\dir_tree_update
% load C:\WC_Germany\JMM_Analysis_pyr\time_domain_30\correl_data
% load C:\WC_Germany\JMM_Analysis_pyr\time_domain_30\lag_data

dsf = 8;
Fsd = 2016/dsf;
winsize = 30;
niqf = 2016/2;
hcf = 50/niqf;
[b,a] = butter(2,hcf,'low');

pers_thresh = 5;
back_time = 10;
forward_time = 60;
min_dur = 60;
window = [60 60];
params.Fs = Fsd;
params.fpass = [0 20];
params.tapers = [4 7];
params.err = [2 0.05];

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
    pers_measure = zeros(1,total_dur);
    pers_measure(pers_secs) = 1;
    pers_diff = [0 diff(pers_measure)];
    pers_start = find(pers_diff==1);
    pers_stop = find(pers_diff == -1);

    npers_start = pers_stop;
    npers_stop = pers_start;

    %if there is some persistent activity
    if length(pers_start)>0 & length(npers_start) > 0

        if pers_measure(1) ==1
            pers_start = [1 pers_start];
        elseif pers_start(1) > 60
            npers_start = [1 npers_start];
        else
            npers_stop(1) = [];
        end
        if pers_measure(end) == 1
            pers_stop = [pers_stop total_dur];
        elseif total_dur-pers_stop(end) > 60
            npers_stop = [npers_stop total_dur];
        else
            npers_start(end) = [];
        end

        
        pers_dur = pers_stop-pers_start;
        npers_dur = npers_stop-npers_start;

        bad_npers = find(npers_dur < min_dur);
        bad_nstart = npers_start(bad_npers);
        bad_nstop = npers_stop(bad_npers);
        npers_start(bad_npers) = [];
        npers_stop(bad_npers) = [];
        pers_start(ismember(pers_start,bad_nstop)) = [];
        pers_stop(ismember(pers_stop,bad_nstart)) = [];


        %convert to sample index
        pers_start = pers_start*Fsd;
        pers_stop = pers_stop*Fsd;
        npers_start = npers_start*Fsd;
        npers_stop = npers_stop*Fsd;

        if length(pers_start)>0 & length(npers_start) > 0
        
    [Sw_p(d,:), f, Spwerr]= mtspectrumc_unequal_length_trials(wcv_d, window, params, [pers_start' pers_stop']);        
    [S8_p(d,:), f, Sp8err]= mtspectrumc_unequal_length_trials(lf8_d, window, params, [pers_start' pers_stop']);        
    [Sw_n(d,:), f, Snwerr]= mtspectrumc_unequal_length_trials(wcv_d, window, params, [npers_start' npers_stop']);        
    [S8_n(d,:), f, Sn8err]= mtspectrumc_unequal_length_trials(lf8_d, window, params, [npers_start' npers_stop']);        
    [Cmn_p(d,:),Phimn_p,Smn,Smm,f,ConfC_p(d)] = coherencyc_unequal_length_trials( [wcv_d lf8_d], window, params, [pers_start' pers_stop']);
    [Cmn_n(d,:),Phimn_n,Smn,Smm,f,ConfC_n(d)] = coherencyc_unequal_length_trials( [wcv_d lf8_d], window, params, [npers_start' npers_stop']);
    

       Fig = figure(1)
       clf
       set(Fig,'PaperUnits','centimeters');
       set(gcf, 'PaperSize', [30 30]);% paper size is in [width height] format
        subplot(2,2,1)
        plot_vector(Sw_p(d,:),f,'l',Spwerr)
        hold on
        plot_vector(Sw_n(d,:),f,'l',Snwerr,'r')
        xlim([0 1])
        subplot(2,2,2)
        plot_vector(S8_p(d,:),f,'l',Sp8err)
        hold on
        plot_vector(S8_n(d,:),f,'l',Sn8err,'r')
        xlim([0 1])
        subplot(2,2,3)
        plot(f,Cmn_p(d,:),'linewidth',2)
        hold on
        plot(f,Cmn_n(d,:),'r','linewidth',2)
        xlim([0 2])
        line([0 2],[ConfC_p(d) ConfC_p(d)])
        line([0 2],[ConfC_n(d) ConfC_n(d)],'Color','r')
        subplot(2,2,4)
        plot(f,Phimn_p,'linewidth',2)
        hold on
        plot(f,Phimn_n,'r','linewidth',2)
        xlim([0 2])
        tname = ['C:\WC_Germany\JMM_Analysis_pyr\pers_state_anal\st_spec' f_names{d}];
        print('-dpng',tname);
        close all
        
    else
        disp('no persistent activity')
    end
    else
        disp('no persistent activity')
    end
end

%get rid of cells without PA
no_pa = [3 6 12];
Sw_p(no_pa,:) = [];
S8_p(no_pa,:) = [];
Sw_n(no_pa,:) = [];
S8_n(no_pa,:) = [];
Cmn_p(no_pa,:) = [];
Cmn_n(no_pa,:) = [];
Phimn_p(no_pa,:) = [];
Phimn_n(no_pa,:) = [];

%averages and plot
m_p_w = mean(10*log10(Sw_p));
u_p_w = m_p_w+2*std(10*log10(Sw_p))/sqrt(14);
l_p_w = m_p_w-2*std(10*log10(Sw_p))/sqrt(14);

m_n_w = mean(10*log10(Sw_n));
u_n_w = m_n_w+2*std(10*log10(Sw_n))/sqrt(14);
l_n_w = m_n_w-2*std(10*log10(Sw_n))/sqrt(14);

m_p_8 = mean(10*log10(S8_p));
u_p_8 = m_p_8+2*std(10*log10(S8_p))/sqrt(14);
l_p_8 = m_p_8-2*std(10*log10(S8_p))/sqrt(14);

m_n_8 = mean(10*log10(S8_n));
u_n_8 = m_n_8+2*std(10*log10(S8_n))/sqrt(14);
l_n_8 = m_n_8-2*std(10*log10(S8_n))/sqrt(14);

m_p_x = mean(Cmn_p);
u_p_x = m_p_x+2*std(Cmn_p)/sqrt(14);
l_p_x = m_p_x-2*std(Cmn_p)/sqrt(14);

m_n_x = mean(Cmn_n);
u_n_x = m_n_x+2*std(Cmn_n)/sqrt(14);
l_n_x = m_n_x-2*std(Cmn_n)/sqrt(14);

m_p_ph = mean(Phimn_p);
u_p_ph = m_p_x+2*std(Phimn_p)/sqrt(14);
l_p_ph = m_p_x-2*std(Phimn_p)/sqrt(14);

m_n_ph = mean(Phimn_n);
u_n_ph = m_n_x+2*std(Phimn_n)/sqrt(14);
l_n_ph = m_n_x-2*std(Phimn_n)/sqrt(14);



% plot(f,m_p_w,'linewidth',2)
% hold on
% plot(f,u_p_w,'--')
% plot(f,l_p_w,'--')
% plot(f,m_n_w,'r','linewidth',2)
% plot(f,u_n_w,'r--')
% plot(f,l_n_w,'r--')
% xlim([0 1])
% 
% plot(f,m_p_8,'linewidth',2)
% hold on
% plot(f,u_p_8,'--')
% plot(f,l_p_8,'--')
% plot(f,m_n_8,'r','linewidth',2)
% plot(f,u_n_8,'r--')
% plot(f,l_n_8,'r--')
% xlim([0 1])

plot(f,m_p_x,'linewidth',2)
hold on
plot(f,u_p_x,'--')
plot(f,l_p_x,'--')
plot(f,m_n_x,'r','linewidth',2)
plot(f,u_n_x,'r--')
plot(f,l_n_x,'r--')
xlim([0 10])

plot(f,m_p_ph,'linewidth',2)
hold on
plot(f,u_p_ph,'--')
plot(f,l_p_ph,'--')
plot(f,m_n_ph,'r','linewidth',2)
plot(f,u_n_ph,'r--')
plot(f,l_n_ph,'r--')
xlim([0 2])