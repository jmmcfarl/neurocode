clear all

load C:\WC_Germany\overall_EC\overall_EC_dir.mat
addpath('C:\WC_Germany\parietal_cortical_2010\')
addpath('C:\WC_Germany\hsmm_state_detection\')
addpath('C:\WC_Germany\persistent_2010\')
addpath('C:\Code\smoothing\software\')
addpath('C:\Code\general\')
addpath('C:\WT_KO_analysis_9_22_10\track_analysis\')

used_data = [l3mec_p l3lec_p];
sess_data = sess_data(used_data);

dsf = 8;
Fsd = 2016/dsf;
params.Fs = Fsd;
params.err = [2 .05];
params.fpass = [0 45];
movingwin = [50 50];
niqf = 2016/2;
hcf = 100/niqf;
% lcf = 0.05/niqf;
% [b1,a1] = butter(2,[lcf hcf]);
[b1,a1] = butter(2,hcf,'low');
params.tapers = [4 7];

hf1 = 15/niqf;
hf2 = 80/niqf;
hfsmooth = round(0.05*Fsd);
[b2,a2] = butter(2,[hf1 hf2]);

for d = 1:length(sess_data)
    
    cdir = sess_data(d).directory;
    cdir(1) = 'C';
    disp(sprintf('session %d',d))
    cd(cdir);
    s_name = strcat(sess_data(d).region,'_l',sess_data(d).layer,'_',sess_data(d).name);
    
    load ./used_data wcv_minus_spike lf8 lf3 lf5 lf2 lf6 lf7
    load ./desynch_times_lf8
    
    lf3_g = lf3 + lf5;
    
    %bandlimit signals
    down_w = filtfilt(b1,a1,wcv_minus_spike);
    down_8 = filtfilt(b1,a1,lf8);
    down_7 = filtfilt(b1,a1,lf7);
    down_6 = filtfilt(b1,a1,lf6);
    down_5 = filtfilt(b1,a1,lf5);
    down_3 = filtfilt(b1,a1,lf3);
    down_2 = filtfilt(b1,a1,lf2);
    down_3g = filtfilt(b1,a1,lf3_g);
   
    down_w = downsample(down_w,dsf);
    down_8 = downsample(down_8,dsf);
    down_7 = downsample(down_7,dsf);
    down_6 = downsample(down_6,dsf);
    down_5 = downsample(down_5,dsf);
    down_3 = downsample(down_3,dsf);
    down_2 = downsample(down_2,dsf);
    down_3g = downsample(down_3g,dsf);
        
    hf8 = downsample(filtfilt(b2,a2,lf8),dsf);
    hf8 = sqrt(jmm_smooth_1d_cor(hf8.^2,hfsmooth))';
    hf5 = downsample(filtfilt(b2,a2,lf5),dsf);
    hf5 = sqrt(jmm_smooth_1d_cor(hf5.^2,hfsmooth))';
    hf3 = downsample(filtfilt(b2,a2,lf3),dsf);
    hf3 = sqrt(jmm_smooth_1d_cor(hf3.^2,hfsmooth))';
    hf3g = downsample(filtfilt(b2,a2,lf3_g),dsf);
    hf3g = sqrt(jmm_smooth_1d_cor(hf3g.^2,hfsmooth))';
    hf2 = downsample(filtfilt(b2,a2,lf2),dsf);
    hf2 = sqrt(jmm_smooth_1d_cor(hf2.^2,hfsmooth))';
   
%     %zscore
%     zdown_w = zscore(down_w);
%     zdown_8 = zscore(down_8);
%     zdown_6 = zscore(down_6);
%     zdown_5 = zscore(down_5);
%     zdown_3 = zscore(down_3);
%     zdown_3g = zscore(down_3g);
%     zdown_2 = zscore(down_2);
    
%     [zdown_3hf] = get_hf_features(lf3,2016,252,[20 80],0.05);
    
    %compute markers indicating segments of data to be used
    t_axis = (1:length(down_w))/Fsd;
    if ~isempty(desynch_times_lf8)
        desynch_start = round(interp1(t_axis,1:length(t_axis),desynch_times_lf8(:,1)));
        desynch_stop = round(interp1(t_axis,1:length(t_axis),desynch_times_lf8(:,2)));
    else
        desynch_start = [];
        desynch_stop = [];
    end
    desynch_ind = zeros(size(down_w));
    for i = 1:length(desynch_start)
        desynch_ind(desynch_start(i):desynch_stop(i)) = 1;
    end
    synch_starts = find(desynch_ind(1:end-1)==1 & desynch_ind(2:end)==0)+1;
    if desynch_ind(1) == 0
        synch_starts = [1; synch_starts];
    end
    synch_stops = find(desynch_ind(1:end-1)==0 & desynch_ind(2:end)==1)+1;
    if desynch_ind(end) == 0
        synch_stops = [synch_stops; length(down_w)];
    end
    sMarkers = [synch_starts(:) synch_stops(:)];
    
    %
%     params.tapers = [5 9];
%     [Cmn(d,:),Phimn(d,:),Smn,Smm,f,ConfC(d),PhiStd,Cerr(d,:,:)] = ...
%         coherencyc_unequal_length_trials([down_w down_8],window, params, sMarkers);
    
%     [zPww(d,:), f, ~]= mtspectrumc_unequal_length_trials(zdown_w, movingwin, params, sMarkers );
%     [zP88(d,:), f, ~]= mtspectrumc_unequal_length_trials(zdown_8, movingwin, params, sMarkers );
%     [zP66(d,:), f, ~]= mtspectrumc_unequal_length_trials(zdown_6, movingwin, params, sMarkers );
%     [zP55(d,:), f, ~]= mtspectrumc_unequal_length_trials(zdown_5, movingwin, params, sMarkers );
%     [zP33(d,:), f, ~]= mtspectrumc_unequal_length_trials(zdown_3, movingwin, params, sMarkers );
%     [zP33g(d,:), f, ~]= mtspectrumc_unequal_length_trials(zdown_3g, movingwin, params, sMarkers );
%     [zP22(d,:), f, ~]= mtspectrumc_unequal_length_trials(zdown_2, movingwin, params, sMarkers );
%     zPww(d,:) = log10(zPww(d,:));
%     zP88(d,:) = log10(zP88(d,:));
%     zP66(d,:) = log10(zP66(d,:));
%     zP55(d,:) = log10(zP55(d,:));
%     zP22(d,:) = log10(zP22(d,:));
%     zP33(d,:) = log10(zP33(d,:));
%     zP33g(d,:) = log10(zP33g(d,:));
    
    [Pww(d,:), f, ~]= mtspectrumc_unequal_length_trials(down_w, movingwin, params, sMarkers );
    [P88(d,:), f, ~]= mtspectrumc_unequal_length_trials(down_8, movingwin, params, sMarkers );
    [P77(d,:), f, ~]= mtspectrumc_unequal_length_trials(down_7, movingwin, params, sMarkers );
    [P66(d,:), f, ~]= mtspectrumc_unequal_length_trials(down_6, movingwin, params, sMarkers );
    [P55(d,:), f, ~]= mtspectrumc_unequal_length_trials(down_5, movingwin, params, sMarkers );
    [P33(d,:), f, ~]= mtspectrumc_unequal_length_trials(down_3, movingwin, params, sMarkers );
    [P33g(d,:), f, ~]= mtspectrumc_unequal_length_trials(down_3g, movingwin, params, sMarkers );
    [P22(d,:), f, ~]= mtspectrumc_unequal_length_trials(down_2, movingwin, params, sMarkers );
    [hP88(d,:), f, ~]= mtspectrumc_unequal_length_trials(hf8, movingwin, params, sMarkers );
    [hP55(d,:), f, ~]= mtspectrumc_unequal_length_trials(hf5, movingwin, params, sMarkers );
    [hP33(d,:), f, ~]= mtspectrumc_unequal_length_trials(hf3, movingwin, params, sMarkers );
    [hP33g(d,:), f, ~]= mtspectrumc_unequal_length_trials(hf3g, movingwin, params, sMarkers );
    [hP22(d,:), f, ~]= mtspectrumc_unequal_length_trials(hf2, movingwin, params, sMarkers );
%     Pww(d,:) = log10(Pww(d,:));
%     P88(d,:) = log10(P88(d,:));
%     P66(d,:) = log10(P66(d,:));
%     P55(d,:) = log10(P55(d,:));
%     P22(d,:) = log10(P22(d,:));
%     P33(d,:) = log10(P33(d,:));
%     P33g(d,:) = log10(P33g(d,:));


%     zP33r(d,:) = log10(P33) - psi(0,dof/2) + log(dof/2);
%     zP33g(d,:) = log10(Pww_33g(:,2)) - psi(0,dof/2) + log(dof/2);
%     zP33h(d,:) = log10(Pww_33h(:,2)) - psi(0,dof/2) + log(dof/2);
    
%     correction(d) = (log(dof/2)-psi(0,dof/2));
    
    clear down_w down_8 wcv* lf8 lf3 lf3_g lf2 lf6
    
end

cd C:\WC_Germany\persistent_9_27_2010\
% save spectral_data_4 zPww zP88 zP33* f zP22* zP66 zP55
save spectral_data_12_10_2011 *P* f

%%
mec = 1:22;
lec = 23:36;
bad_lf3 = [7 9 11 12 13 16 17 25 31 32 35];
good_lf3 = setdiff(1:36,bad_lf3);
good_lf3_mec = good_lf3(ismember(good_lf3,mec));
good_lf3_lec = good_lf3(ismember(good_lf3,lec));
%%
uds_freqs = find(f > 0.2 & f < 1);
[peak_S8,peakf8] = max(P88(:,uds_freqs),[],2);
[peak_S7,peakf7] = max(P77(:,uds_freqs),[],2);
[peak_S6,peakf6] = max(P66(:,uds_freqs),[],2);
[peak_S5,peakf5] = max(P55(:,uds_freqs),[],2);
[peak_S3,peakf3] = max(P33(:,uds_freqs),[],2);
[peak_S3g,peakf3g] = max(P33g(:,uds_freqs),[],2);
[peak_S2,peakf2] = max(P22(:,uds_freqs),[],2);
 
[peak_S8h,peakf8h] = max(hP88(:,uds_freqs),[],2);
[peak_S5h,peakf5h] = max(hP55(:,uds_freqs),[],2);
[peak_S3gh,peakf3gh] = max(hP33g(:,uds_freqs),[],2);
[peak_S2h,peakf2h] = max(hP22(:,uds_freqs),[],2);

peak_f8 = f(uds_freqs(peakf8));

high_freqs = find(f >= 15);
hf_pow8 = trapz(P88(:,high_freqs),2);
hf_pow7 = trapz(P77(:,high_freqs),2);
hf_pow6 = trapz(P66(:,high_freqs),2);
hf_pow5 = trapz(P55(:,high_freqs),2);
hf_pow3 = trapz(P33(:,high_freqs),2);
hf_pow3g = trapz(P33g(:,high_freqs),2);
hf_pow2 = trapz(P22(:,high_freqs),2);

%%
cmap = colormap(jet(6));
for i = 1:36
    plot(f,P88(i,:),'color',cmap(1,:))
    hold on
    plot(f,P77(i,:),'color',cmap(2,:))
    plot(f,P66(i,:),'color',cmap(3,:))
    plot(f,P55(i,:),'color',cmap(4,:))
    plot(f,P33g(i,:),'color',cmap(5,:))
    plot(f,P22(i,:),'color',cmap(6,:))
    xlim([0 1])
   i
   pause
   clf
end

%%
cmap = colormap(jet(6));
for i = 1:36
    plot(f,hP88(i,:),'color',cmap(1,:))
    hold on
    plot(f,hP55(i,:),'color',cmap(4,:))
    plot(f,hP33g(i,:),'color',cmap(5,:))
    plot(f,hP22(i,:),'color',cmap(6,:))
    xlim([0 1])
   i
   pause
   clf
end


%%
figure
set(gca,'fontname','arial')
hold on
h = errorbar(f,nanmean(zP88),nanstd(zP88)/sqrt(36),'r');
errorbar_tick(h,.001,'units');
h = errorbar(f,nanmean(zP33(good_lf3,:)),nanstd(zP33(good_lf3,:))/sqrt(length(good_lf3)),'c');
errorbar_tick(h,.001,'units');
h = errorbar(f,nanmean(zPww(mec,:)),nanstd(zPww(mec,:))/sqrt(length(mec)),'b');
errorbar_tick(h,.001,'units');
h = errorbar(f,nanmean(zPww(lec,:)),nanstd(zPww(lec,:))/sqrt(length(lec)),'g');
errorbar_tick(h,.001,'units');
xlim([0 1])
legend('LF8','LF3','MEC','LEC')
xlabel('Frequency (Hz)','fontsize',14,'fontname','arial')
ylabel('Relative power','fontsize',14,'fontname','arial')


