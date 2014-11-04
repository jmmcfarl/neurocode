clear all
cd C:\WC_Germany\sven_thomas_combined\
load ./combined_core_analysis_fin_nd.mat
load ./combined_dir_nd.mat
uset = sort([l3mec l3lec]);
all_cells = 1:61;
l3mec = find(ismember(all_cells(uset),l3mec));
l3lec = find(ismember(all_cells(uset),l3lec));
combined_dir = combined_dir(uset);
hpc_mua = hpc_mua(uset);
hpc_lfp = hpc_lfp(uset);
ctx_lfp = ctx_lfp(uset);

dsf = 8;
raw_Fs = 2016;
Fsd = raw_Fs/dsf;
rate_sm = round(Fsd*0.05);
int_width_range = [2 8]; %range of spike widths (in 32kHz samples) for interneuron spikes
pyr_width_range = [10 18]; %range of spike widths for pyramidal spikes
min_corr = 0.4; %minimum correlation coefficient with the corresponding spike template


%%
% hpc_mua(~isnan(hpc_mua)) = 8;
for d = 1:length(combined_dir)
    if ~isnan(hpc_mua(d))
       d 
        cdir = combined_dir{d};
        disp(sprintf('session %d',d))
        cd(cdir);
        
        load ./spike_time_jmm
        load ./used_data wcv
        load ./sync_times.mat
%         s_fs = median(diff(synct));
%         s_fsd = s_fs*dsf;
        synct(length(wcv)+1:end) = [];
        synct_d = downsample(synct,dsf);
%         synct_d = synct(1):s_fsd:synct(end);
        load ./mua_data3
%         hpc_mua_times = mua_times{hpc_mua(d)};
        hpc_mua_times = mua_times{hpc_mua(d)+1}; %use MORE superficial channel for test

        %         good_mua = find(mua_corr_pyr{hpc_mua(d)} > min_corr | mua_corr_int{hpc_mua(d)} > min_corr);
%         hpc_mua_times_good = hpc_mua_times(good_mua);
%         pyr_mua = find(mua_corr_pyr{hpc_mua(d)} > min_corr & mua_widths{hpc_mua(d)}' > pyr_width_range(1) & mua_widths{hpc_mua(d)}' < pyr_width_range(2));
%         hpc_mua_times_pyr = hpc_mua_times(pyr_mua);
%         int_mua = find(mua_corr_int{hpc_mua(d)} > min_corr & mua_widths{hpc_mua(d)}' > int_width_range(1) & mua_widths{hpc_mua(d)}' < int_width_range(2));
%         hpc_mua_times_int = hpc_mua_times(int_mua);
        
        hpc_mua_times(hpc_mua_times > synct_d(end) | hpc_mua_times < synct_d(1)) = [];
%         hpc_mua_times_good(hpc_mua_times_good > synct_d(end) | hpc_mua_times_good < synct_d(1)) = [];
%         hpc_mua_times_pyr(hpc_mua_times_pyr > synct_d(end) | hpc_mua_times_pyr < synct_d(1)) = [];
%         hpc_mua_times_int(hpc_mua_times_int > synct_d(end) | hpc_mua_times_int < synct_d(1)) = [];
  
        hpc_mua_rate =hist(hpc_mua_times,synct_d)*Fsd;
        hpc_mua_rate = jmm_smooth_1d_cor(hpc_mua_rate,rate_sm);
        hpc_mua_rate = zscore(hpc_mua_rate');
 
%         hpc_mua_rate_good =hist(hpc_mua_times_good,synct_d)*Fsd;
%         hpc_mua_rate_good = jmm_smooth_1d_cor(hpc_mua_rate_good,rate_sm);
%         hpc_mua_rate_good = zscore(hpc_mua_rate_good');
% 
%         hpc_mua_rate_pyr =hist(hpc_mua_times_pyr,synct_d)*Fsd;
%         hpc_mua_rate_pyr = jmm_smooth_1d_cor(hpc_mua_rate_pyr,rate_sm);
%         hpc_mua_rate_pyr = zscore(hpc_mua_rate_pyr');

%         hpc_mua_rate_int =hist(hpc_mua_times_int,synct_d)*Fsd;
%         hpc_mua_rate_int = jmm_smooth_1d_cor(hpc_mua_rate_int,rate_sm);
%         hpc_mua_rate_int = zscore(hpc_mua_rate_int');
        
%                 [log_mua,offset(d)] = log_transform_sig(hpc_mua_rate);
%         hpc_mua_rate = zscore(log_mua(:));
        
        if ctx_lfp(d) == 8
            ctx_mua_times = mua_times{8};
            ctx_mua_times(ctx_mua_times > synct_d(end) | ctx_mua_times < synct_d(1)) = [];
        else
            ctx_mua_times = mua_times{7};
            ctx_mua_times(ctx_mua_times > synct_d(end) | ctx_mua_times < synct_d(1)) = [];
        end
        ctx_mua_rate =hist(ctx_mua_times,synct_d)*Fsd;
%         if length(hpc_mua_rate) > length(synct_d)
%             hpc_mua_rate = hpc_mua_rate(1:length(synct_d));
%             ctx_mua_rate = ctx_mua_rate(1:length(synct_d));
%         end
%         hpc_mua_rate = zscore(hpc_mua_rate(:));
       
        
        load ./pa_hsmm_state_seq_combined_fin_nd.mat
        load ./pa_hsmm_state_seq7_combined_fin_nd.mat
%         load ./pa_hsmm_state_seq_combined_fin
%         load ./pa_hsmm_state_seq7_combined_fin
        mp_state_seq = hsmm_bbstate_seq;
%         load ./pa_hsmm_state_seq7_combined
        lfp_state_seq = hsmm_bbstate_seq7;
        
        [new_seg_inds] = resample_uds_seg_inds_v2(hsmm.UDS_segs,hsmm.Fs,Fsd,mp_state_seq);
        seg_durs = new_seg_inds(:,2)-new_seg_inds(:,1)+1;
        
        [up_trans_inds,down_trans_inds] = compute_state_transitions_seg(new_seg_inds,mp_state_seq);
        up_state_durs = (down_trans_inds-up_trans_inds)/Fsd;
        
        wcv_up_log = nan(size(synct_d));
        lf8_up_log = nan(size(synct_d));
        
        for ns = 1:hmm.Nsegs
            cur_seg = new_seg_inds(ns,1):new_seg_inds(ns,2);
            wcv_up_log(cur_seg) = logical(mp_state_seq{ns}-1);
            lf8_up_log(cur_seg) = logical(lfp_state_seq{ns}-1);
        end
        
%         %only take persistent mp up states
%         for i = 1:length(up_trans_inds)
%             if ~ismember(i,rt2_ups{d})
%                 wcv_up_log(up_trans_inds(i):down_trans_inds(i)) = nan;
%             end
%         end
%         hpc_mua_rate(wcv_up_log==0) = nan;
 
%         pos_lags = round(Fsd*(-0.2:.05:0.2));
        pos_lags = 0;
        X = [ones(length(wcv_up_log),1) wcv_up_log(:) lf8_up_log(:)];
%         X = [ones(length(wcv_up_log),1) lf8_up_log(:)];
%          X = [ones(length(wcv_up_log),1) wcv_up_log(:)];
       uX = X;
        B = nan(length(pos_lags),length(pos_lags),3);
        R2 = nan(length(pos_lags),length(pos_lags));
        Stats = nan(length(pos_lags),length(pos_lags),3);
        for i = 1:length(pos_lags)
            fprintf('Lag %d of %d\n',i,length(pos_lags));
            for j = 1:length(pos_lags)
                if pos_lags(i) > 0
                    uX(:,2) = [X(pos_lags(i)+1:end,2); nan(pos_lags(i),1)];
                elseif pos_lags(i) < 0
                    uX(:,2) = [nan(-pos_lags(i),1); X(1:end+pos_lags(i),2)];
                end
                if pos_lags(j) > 0
                    uX(:,3) = [X(pos_lags(j)+1:end,3); nan(pos_lags(j),1)];
                elseif pos_lags(j) < 0
                    uX(:,3) = [nan(-pos_lags(j),1); X(1:end+pos_lags(j),3)];
                end
                temp_y = hpc_mua_rate;
                temp_y(any(isnan(uX),2)) = nan;
                [B(i,j,:),Bint,R,Rint,Stats] = regress(temp_y,uX);
%                 temp_y = hpc_mua_rate_good;
%                 temp_y(any(isnan(uX),2)) = nan;
%                 [B_good(i,j,:),Bint,R,Rint,Stats] = regress(temp_y,uX);
%                 temp_y = hpc_mua_rate_pyr;
%                 temp_y(any(isnan(uX),2)) = nan;
%                 [B_pyr(i,j,:),Bint,R,Rint,Stats] = regress(temp_y,uX);
%                 temp_y = hpc_mua_rate_int;
%                 temp_y(any(isnan(uX),2)) = nan;
%                 [B_int(i,j,:),Bint,R,Rint,Stats] = regress(temp_y,uX);

                
                % uX(:,1) = [];
% [B(i,j,:),Stats] = robustfit(uX,hpc_mua_rate,'fair');
% R2(i,j) = 1;
                R2(i,j) = Stats(1);
            end
        end
%         R2(i,j) = Stats(1);
        temp = find(R2==max(R2(:)),1,'first');
        [a,b] = ind2sub(size(R2),temp);
%         lags_w(d) = pos_lags(a)/Fsd;
%         lags_8(d) = pos_lags(b)/Fsd;
        alpha_w(d) = B(a,b,2);
        alpha_8(d) = B(a,b,3);
%          alpha_8(d) = B(a,b,2);
       const(d) = B(a,b,1);
 
%         alpha_w_good(d) = B_good(a,b,2);
%         alpha_8_good(d) = B_good(a,b,3);
%         const_good(d) = B_good(a,b,1);
% 
%         alpha_w_pyr(d) = B_pyr(a,b,2);
%         alpha_8_pyr(d) = B_pyr(a,b,3);
%         const_pyr(d) = B_pyr(a,b,1);
%         
%         alpha_w_int(d) = B_int(a,b,2);
%         alpha_8_int(d) = B_int(a,b,3);
%         const_int(d) = B_int(a,b,1);

%         if pos_lags(a) > 0
%             uX(:,2) = [X(pos_lags(a)+1:end,2); nan(pos_lags(a),1)];
%         elseif pos_lags(a) < 0
%             uX(:,2) = [nan(-pos_lags(a),1); X(1:end+pos_lags(a),2)];
%         end
%         if pos_lags(b) > 0
%             uX(:,3) = [X(pos_lags(b)+1:end,3); nan(pos_lags(b),1)];
%         elseif pos_lags(b) < 0
%             uX(:,3) = [nan(-pos_lags(b),1); X(1:end+pos_lags(b),3)];
%         end
%         uX(:,1) = [];
%         stats = regstats(hpc_mua_rate,uX,'linear',{'tstat'});
%         tstat_w(d) = stats.tstat.t(2);
%         tstat_8(d) = stats.tstat.t(3);
    end
end
l3mec_m = l3mec(~isnan(hpc_mua(l3mec)));

cd C:\WC_Germany\sven_thomas_combined\
save hpc_mua_modelsfin_nd_super *_w* *_8* const* pos_lags

%%

figure
set(gca,'fontname','arial','fontsize',14)
plot(alpha_8(l3mec_m),alpha_w(l3mec_m),'o','markersize',8,'linewidth',1)
hold on
% plot(alpha_8_good(l3mec_m),alpha_w_good(l3mec_m),'ko','markersize',8,'linewidth',2)
% plot(alpha_8_pyr(l3mec_m),alpha_w_pyr(l3mec_m),'ro','markersize',8,'linewidth',2)
% plot(alpha_8_int(l3mec_m),alpha_w_int(l3mec_m),'go','markersize',8,'linewidth',2)
% ylim([-1.4 1.4])
% xlim([-1.4 1.4])
ylim([0 1.5])
xlim([-0.75 0.75])
yl = ylim();xl = xlim();
line([-1.4 1.4],[-1.4 1.4],'color','k')
line([0 0],yl,'color','k','linestyle','--')
xlabel('Cortical modulation','fontsize',16,'fontname','arial')
ylabel('MECL3 modulation','fontsize',16,'fontname','arial')

%%
figure
cmap = colormap(jet(7));
load ./hpc_mua_models2
plot(alpha_8_good(l3mec_m),alpha_w_good(l3mec_m),'o','color',cmap(3,:),'markersize',8,'linewidth',2)
hold on
load ./hpc_mua_models2_el2
plot(alpha_8_good(l3mec_m),alpha_w_good(l3mec_m),'o','color',cmap(1,:),'markersize',8,'linewidth',2)
load ./hpc_mua_models2_el3
plot(alpha_8_good(l3mec_m),alpha_w_good(l3mec_m),'o','color',cmap(2,:),'markersize',8,'linewidth',2)
load ./hpc_mua_models2_el5
plot(alpha_8_good(l3mec_m),alpha_w_good(l3mec_m),'o','color',cmap(4,:),'markersize',8,'linewidth',2)
load ./hpc_mua_models2_el6
plot(alpha_8_good(l3mec_m),alpha_w_good(l3mec_m),'o','color',cmap(5,:),'markersize',8,'linewidth',2)
load ./hpc_mua_models2_el7
plot(alpha_8_good(l3mec_m),alpha_w_good(l3mec_m),'o','color',cmap(6,:),'markersize',8,'linewidth',2)
load ./hpc_mua_models2_el8
plot(alpha_8_good(l3mec_m),alpha_w_good(l3mec_m),'o','color',cmap(7,:),'markersize',8,'linewidth',2)
legend('MUA4','MUA2','MUA3','MUA5','MUA6','MUA7','MUA8')
xlabel('Cortical UDS','fontsize',14)
ylabel('MEC UDS','fontsize',14)

%%
load ./hpc_mua_modelsfin_nd.mat
alpha_8_all = alpha_8;
load ./hpc_mua_modelsfin_personly_nd.mat
alpha_8_pers = alpha_8;

figure
set(gca,'fontname','arial','fontsize',14)
plot(alpha_8_all(l3mec_m),alpha_8_pers(l3mec_m),'o','markersize',8,'linewidth',2)
hold on
ylim([-1 0.4])
xlim([-0.7 0.4])
yl = ylim();xl = xlim();
line([-1.4 1.4],[-1.4 1.4],'color','k')
line([0 0],yl,'color','k','linestyle','--')
line(xl,[0 0],'color','k','linestyle','--')
xlabel('Cortical modulation (all times)','fontsize',16,'fontname','arial')
ylabel('Cortical modulation (persistent MECL3 Up states)','fontsize',16,'fontname','arial')

%%
load ./spike_rate_data_fin_nd.mat
figure
set(gca,'fontname','arial','fontsize',14)
plot(hpc_npers_rate(l3mec_m),hpc_pers_rate(l3mec_m),'o','markersize',8,'linewidth',2)
hold on
ylim([-0.2 1])
xlim([-0.4 0.6])
yl = ylim();xl = xlim();
line([-1.4 1.4],[-1.4 1.4],'color','k')
line([0 0],yl,'color','k','linestyle','--')
line(xl,[0 0],'color','k','linestyle','--')
xlabel('Hpc MUA rate non-persistent MECL3 Up states (z)','fontsize',16,'fontname','arial')
ylabel('Hpc MUA rate persistent MECL3 Up states (z)','fontsize',16,'fontname','arial')
