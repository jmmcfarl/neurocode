clear all
cd C:\WC_Germany\sven_thomas_combined
load ./combined_dir_nd_dist.mat

raw_Fs = 2016;
dsf = 8;
Fsd = raw_Fs/dsf;
niqf = raw_Fs/2;
hcf = 2;
lcf = 0.05;
maxlag = round(Fsd*1);
mua_sm = round(Fsd*0.05);
for d = 1:length(combined_dir)
    cd(combined_dir{d})
    pwd
    
    if exist('./mua_data3.mat','file')
        load ./mua_data3
        load ./used_data lf7
        [lf8_lf,t_axis] = get_lf_features(lf7,raw_Fs,Fsd,[lcf hcf]);
        load ./sync_times
        synct_d = downsample(synct,dsf);
        
        mua_binned = hist(mua_times{8},synct_d);
        mua_binned([1 end]) = 0;
        lf8_lf = jmm_smooth_1d_cor(mua_binned,mua_sm);
        
        %%
        figure
        for ch = 2:8
            mua_binned = hist(mua_times{ch},synct_d);
            mua_binned([1 end]) = 0;
            rate(ch) = sum(mua_binned)/range(t_axis);
            mua_rate = jmm_smooth_1d_cor(mua_binned,mua_sm);
            [mua_lfp_xc(ch,:),lags] = xcov(lf8_lf,mua_rate,maxlag,'coeff');
            [xc_peakval(ch),xc_peakloc(ch)] = max(mua_lfp_xc(ch,:));
            subplot(7,2,(8-ch)*2+1)
            plot(lags/Fsd,mua_lfp_xc(ch,:),'linewidth',2)
            hold on
            plot(lags(xc_peakloc(ch))/Fsd,xc_peakval(ch),'ro','linewidth',2)
            xlim(lags([1 end])/Fsd)
            ylim([-0.4 1])
            grid on
         set(gca,'fontsize',20)
       end
        subplot(7,2,(1:7)*2)
        plot(rate(2:8),2:8,'o-')
        set(gca,'fontsize',20)
        fillPage(gcf,'Papersize',[15 20],'Margins',[0 0 0 0])

        temp = find(combined_dir{d} == '\',1,'last');
        cur_name = combined_dir{d}(temp+1:end);
        fname = ['C:\WC_Germany\sven_thomas_combined\MUA\depth2_prof_' cur_name];
        print('-djpeg',fname);close all

    end
end