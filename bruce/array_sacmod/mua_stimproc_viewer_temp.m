close all
n_mus = length(mua_data);
use_sua_hor = get_struct_data(mua_data,1:n_mus,'stc_hor_use');
use_sua_ver = get_struct_data(mua_data,1:n_mus,'stc_ver_use');
xl = [-0.2 0.4];
f1 = figure();
f2 = figure();
f3 = figure();
for ss = 1:n_mus
    if use_sua_hor(ss)
        figure(f1);
        subplot(6,3,1)
        cur_sta = mua_data(ss).sta_hor; ca = max(abs(cur_sta));
        imagesc(reshape(cur_sta,12,24)); caxis([-ca ca]);
        cur_stcs = mua_data(ss).stcs_hor; ca = max(abs(cur_stcs(:)));
        for ii = 1:3
        subplot(6,3,3+ii)
        imagesc(reshape(cur_stcs(:,ii),12,24)); caxis([-ca ca]);
        end
        for ii = 1:3
        subplot(6,3,6+ii)
        imagesc(reshape(cur_stcs(:,end-3+ii),12,24)); caxis([-ca ca]);
        end
        
        figure(f2);
        subplot(2,4,1);hold on
        plot(sac_bin_cents,mua_data(ss).gsac_ta_hor_sm);xlim(xl);
        subplot(2,4,2); hold on
        plot(sac_bin_cents,mua_data(ss).gsac_tkern_hor);xlim(xl);
        subplot(2,4,3); hold on
        plot(sac_bin_cents,mua_data(ss).gsac_skern_hor);xlim(xl);
        subplot(2,4,4); hold on
        plot(sac_bin_cents,mua_data(ss).gsacdep_info_hor);xlim(xl);
        line(xl,mua_data(ss).ov_info_hor([1 1]),'color','k','linestyle','--');
        
        figure(f3);
        subplot(2,4,1);hold on
        plot(sac_bin_cents,mua_data(ss).msac_ta_hor_sm,'r');xlim(xl);
        subplot(2,4,2); hold on
        plot(sac_bin_cents,mua_data(ss).msac_tkern_hor,'r');xlim(xl);
        subplot(2,4,3); hold on
        plot(sac_bin_cents,mua_data(ss).msac_skern_hor,'r');xlim(xl);
        subplot(2,4,4); hold on
        plot(sac_bin_cents,mua_data(ss).msacdep_info_hor,'r');xlim(xl);
        line(xl,mua_data(ss).ov_info_hor([1 1]),'color','k','linestyle','--');
    end
    if use_sua_ver(ss)
        figure(f1);
        subplot(6,3,9+1)
        cur_sta = mua_data(ss).sta_ver; ca = max(abs(cur_sta));
        imagesc(reshape(cur_sta,12,24)); caxis([-ca ca]);
        cur_stcs = mua_data(ss).stcs_ver; ca = max(abs(cur_stcs(:)));
        for ii = 1:3
        subplot(6,3,9+3+ii)
        imagesc(reshape(cur_stcs(:,ii),12,24)); caxis([-ca ca]);
        end
        for ii = 1:3
        subplot(6,3,9+6+ii)
        imagesc(reshape(cur_stcs(:,end-3+ii),12,24)); caxis([-ca ca]);
        end
        
        figure(f2);
        subplot(2,4,5);hold on
        plot(sac_bin_cents,mua_data(ss).gsac_ta_ver_sm);xlim(xl);
        subplot(2,4,6);hold on
        plot(sac_bin_cents,mua_data(ss).gsac_tkern_ver);xlim(xl);
        subplot(2,4,7);hold on
        plot(sac_bin_cents,mua_data(ss).gsac_skern_ver);xlim(xl);
        subplot(2,4,8);hold on
        plot(sac_bin_cents,mua_data(ss).gsacdep_info_ver);xlim(xl);
        line(xl,mua_data(ss).ov_info_ver([1 1]),'color','k','linestyle','--');
        
        figure(f3);
        subplot(2,4,5);hold on
        plot(sac_bin_cents,mua_data(ss).msac_ta_ver_sm,'r');xlim(xl);
        subplot(2,4,6);hold on
        plot(sac_bin_cents,mua_data(ss).msac_tkern_ver,'r');xlim(xl);
        subplot(2,4,7);hold on
        plot(sac_bin_cents,mua_data(ss).msac_skern_ver,'r');xlim(xl);
        subplot(2,4,8);hold on
        plot(sac_bin_cents,mua_data(ss).msacdep_info_ver,'r');xlim(xl);
        line(xl,mua_data(ss).ov_info_ver([1 1]),'color','k','linestyle','--');
    end
     
    pause
    figure(f1);clf;
    figure(f2);clf;
    figure(f3);clf;
end