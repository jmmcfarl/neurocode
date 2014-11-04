clear all
close all

cd C:\WC_Germany\sven_thomas_combined\
load ./combined_dir_nd.mat
addpath('C:\WC_Germany\parietal_cortical_2010\')
addpath('C:\WC_Germany\hsmm_state_detection\')
uset = sort([l3mec l3lec]);
all_cells = 1:61;
l3mec = find(ismember(all_cells(uset),l3mec));
l3lec = find(ismember(all_cells(uset),l3lec));
combined_dir = combined_dir(uset);
hpc_mua = hpc_mua(uset);
hpc_lfp = hpc_lfp(uset);
ctx_lfp = ctx_lfp(uset);

raw_Fs = 2016;
dsf = 8;
Fsd = raw_Fs/dsf;
minSegLength = 60;
maxLag = 10*Fsd;
lcf = .05;
hcf = 10;
rate_sm = round(Fsd*0.05);

%%
for d = 1:length(combined_dir)
    d
    cur_dir = combined_dir{d};
    cd(cur_dir)
    load ./used_data lf8 lf7 wcv_minus_spike     
    if ctx_lfp(d) == 7
        lf8 = lf7;
    end
    [desynch_times,desynch_inds,P_lf8,f,t] = locate_desynch_times_individual_v2(lf8);

    [lf8_lf,t_axis] = get_lf_features(lf8,raw_Fs,Fsd,[lcf hcf]);
    wcv_lf = get_lf_features(wcv_minus_spike,raw_Fs,Fsd,[lcf hcf]);          

    load ./sync_times.mat
    synct_d = downsample(synct,dsf);
    if ~isnan(hpc_mua(d))
        load ./mua_data2
        hpc_mua_times = mua_times{hpc_mua(d)};
        hpc_mua_times(hpc_mua_times > synct_d(end) | hpc_mua_times < synct_d(1)) = [];
        hpc_mua_rate =hist(hpc_mua_times,synct_d)*Fsd;
        hpc_mua_rate = jmm_smooth_1d_cor(hpc_mua_rate,rate_sm);
        if length(hpc_mua_rate) > length(synct_d)
            hpc_mua_rate = hpc_mua_rate(1:length(synct_d));
        end
    else
        hpc_mua_rate = nan(size(synct_d));
    end

    %compute markers indicating segments of data to be used
    if ~isempty(desynch_times)
        desynch_start = round(interp1(t_axis,1:length(t_axis),desynch_times(:,1)));
        desynch_stop = round(interp1(t_axis,1:length(t_axis),desynch_times(:,2)));
        synch_start = desynch_stop;
        synch_stop = desynch_start;
        if desynch_start(1) ~= 1
           synch_start = [1; synch_start]; 
        end
        if desynch_stop(end) ~= length(t_axis)
            synch_stop = [synch_stop; length(t_axis)];
        end
        seg_inds = [synch_start(:) synch_stop(:)];
    else
        desynch_start = [];
        desynch_stop = [];
        seg_inds = [1 length(t_axis)];
    end
    seg_durs = diff(seg_inds,[],2)/Fsd;
    too_short = find(seg_durs < minSegLength);
    seg_inds(too_short,:) = [];
    seg_durs(too_short) = [];
    
    cnt = 0;
    for i = 1:size(seg_inds,1)
        cnt = cnt+1;
        cur_seg_inds = seg_inds(i,1):seg_inds(i,2);
        [wcv_acorr(cnt,:),lags] = xcov(wcv_lf(cur_seg_inds),maxLag,'coeff');
        [lf8_acorr(cnt,:),lags] = xcov(lf8_lf(cur_seg_inds),maxLag,'coeff');
        
        [w8_x(cnt,:),lags] = xcov(wcv_lf(cur_seg_inds),lf8_lf(cur_seg_inds),maxLag,'coeff');
        
        if ~isnan(hpc_mua(d))
          [xm_w(cnt,:),lags] = xcov(wcv_lf(cur_seg_inds),hpc_mua_rate(cur_seg_inds),maxLag,'coeff');          
          [xm_8(cnt,:),lags] = xcov(lf8_lf(cur_seg_inds),hpc_mua_rate(cur_seg_inds),maxLag,'coeff');          
        end
    end

    %compute weighted averages (weighted by relative duration
    weighting = seg_durs/sum(seg_durs);
    tot_wcv_acorr(d,:) = sum(wcv_acorr.*repmat(weighting(:),1,length(lags)),1);
    tot_lf8_acorr(d,:) = sum(lf8_acorr.*repmat(weighting(:),1,length(lags)),1);
    tot_w8_x(d,:) = sum(w8_x.*repmat(weighting(:),1,length(lags)),1);   
    if ~isnan(hpc_mua(d))
       tot_wm_x(d,:) = sum(xm_w.*repmat(weighting(:),1,length(lags)),1);
       tot_8m_x(d,:) = sum(xm_8.*repmat(weighting(:),1,length(lags)),1);
    end
    
    clear wcv_acorr lf8_acorr w8_x xm_w xm_8
end

cd C:\WC_Germany\sven_thomas_combined\
save combined_time_domain_data_fin_nd tot* Fsd lags


%%
[peak_corr,peak_lag] = max(tot_w8_x,[],2);
peak_lag = lags(peak_lag)/Fsd;

%%
l3mec_m = l3mec(~isnan(hpc_mua(l3mec)));
[peak_corr8,peak_lag8] = max(tot_8m_x(l3mec_m,:),[],2);
peak_lag8 = lags(peak_lag8)/Fsd;
xl = [-1.5 1.5];
ul = find(lags/Fsd >= xl(1) & lags/Fsd <= xl(2));
l3mec_m = l3mec(~isnan(hpc_mua(l3mec)));
[peak_corrw,peak_lagw] = max(tot_wm_x(l3mec_m,:),[],2);
peak_lagw = lags(peak_lagw)/Fsd;
% tot_wm_x = fliplr(tot_wm_x);
figure
set(gca,'fontname','arial','fontsize',14)
plot(lags(ul)/Fsd,mean(tot_wm_x(l3mec_m,ul)),'r'); hold on
plot(lags(ul)/Fsd,mean(tot_8m_x(l3mec_m,ul)),'k'); hold on
legend('MECL3 - Hpc MUA','Ncx LFP - Hpc MUA')
shadedErrorBar(lags(ul)/Fsd,mean(tot_wm_x(l3mec_m,ul)),std(tot_wm_x(l3mec_m,ul))./sqrt(length(l3mec_m)),{'r'});
hold on
shadedErrorBar(lags(ul)/Fsd,mean(tot_8m_x(l3mec_m,ul)),std(tot_8m_x(l3mec_m,ul))./sqrt(length(l3mec_m)),{'k'});
xlim(xl)
ylim([-0.2 0.4])
yl = ylim();
line([0 0],yl,'color','k')
xlabel('Time lag (s)','fontsize',16,'fontname','arial')
ylabel('Correlation','fontsize',16,'fontname','arial')
% hold on
% plot(lags/Fsd,tot_wm_x(l3mec_m,:),'k','linewidth',0.5)

%%
figure; set(gca,'fontname','arial','fontsize',14)
plot(lags/Fsd,mean(tot_w8_x(l3mec,:)),'r','linewidth',2); hold on
plot(lags/Fsd,mean(tot_w8_x(l3lec,:)),'b','linewidth',2)
legend('MEC','LEC')
shadedErrorBar(lags/Fsd,mean(tot_w8_x(l3mec,:)),std(tot_w8_x(l3mec,:))./sqrt(length(l3mec)),{'r'});
shadedErrorBar(lags/Fsd,mean(tot_w8_x(l3lec,:)),std(tot_w8_x(l3lec,:))./sqrt(length(l3lec)),{'b'});
xlim([-2 2])
ylim([-0.4 0.8])
xlabel('Time lag (s)','fontsize',16,'fontname','arial')
ylabel('Correlation','fontsize',16,'fontname','arial')
yl = ylim();
line([0 0],yl,'color','k')

ul = find(lags/Fsd >= -1.5 & lags/Fsd <= 1.5);
[peak_corr,peak_lag] = max(tot_w8_x(:,ul),[],2);
peak_lag = lags(ul(peak_lag))/Fsd;

