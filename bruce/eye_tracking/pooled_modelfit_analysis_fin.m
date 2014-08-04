clear all
close all

expt_list = [85 86 87 88 89 93 95];
% expt_list = [86];
print_out = false;

load ~/Analysis/bruce/summary_analysis/su_data.mat

fig_dir = '~/Analysis/bruce/summary_analysis/eyetrack_figs/';

%%
tot_su_cnt = 1;
tot_mu_cnt = 1;
for ee = 1:length(expt_list)
    
    fprintf('Loading data for expt %d\n',expt_list(ee));
    anal_dir = ['~/Analysis/bruce/' sprintf('G0%d',expt_list(ee))];
    hor_data_name = [anal_dir '/monoc_eyecorr_hbar_finf.mat'];
    ver_data_name = [anal_dir '/monoc_eyecorr_vbar_finf.mat'];
    
    cur_expt_id = find(su_data.expt_nums == expt_list(ee));
    su_probes = find(su_data.is_su(cur_expt_id,:));
    
    %load in SU models for horizontal bars
    if exist(hor_data_name,'file')
        load(hor_data_name);
        hor_tr_set = et_tr_set;
        hor_su_inds = find(et_tr_set > 96);
        hor_su_probeinds = find(ismember(1:length(su_probes),hor_tr_set(hor_su_inds)-96));
        hor_mu_inds = find(et_tr_set <= 96);
        hor_mu_probeinds = find(ismember(1:96,hor_tr_set(hor_mu_inds)));
        
        hor_before_mods = it_mods{1}(hor_tr_set(hor_su_inds));
        hor_after_mods = dit_xv_mods{2}(hor_tr_set(hor_su_inds));
        hor_init_LL_imps = it_LLimp(1,hor_tr_set(hor_su_inds));
        hor_all_LL_imps = it_LLimp(:,hor_tr_set(hor_su_inds));
        hor_fin_LL_imps = dit_xv_LLimp(2,hor_tr_set(hor_su_inds));
        hor_init_true_xvLL = it_truexvLLimp(1,hor_tr_set(hor_su_inds));
        hor_fin_true_xvLL = dit_truexvLLimp(2,hor_tr_set(hor_su_inds));
                
        hor_MUbefore_mods = it_mods{1}(hor_tr_set(hor_mu_inds));
        hor_MUafter_mods = dit_xvmods{2}(hor_tr_set(hor_mu_inds));
        hor_MUinit_LL_imps = it_xvLLimp(1,hor_tr_set(hor_mu_inds));
        hor_MUall_LL_imps = it_xvLLimp(:,hor_tr_set(hor_mu_inds));
        hor_MUfin_LL_imps = dit_xvLLimp(2,hor_tr_set(hor_mu_inds));
        hor_MUinit_true_xvLL = it_truexvLLimp(1,hor_tr_set(hor_mu_inds));
        hor_MUfin_true_xvLL = dit_truexvLLimp(2,hor_tr_set(hor_mu_inds));
    else
        hor_tr_set = [];
        hor_su_inds = [];
        hor_su_probeinds = [];
        hor_mu_inds = [];
        hor_mu_probeinds = [];
    end
    
    %load in SU models for vertical bars
    if exist(ver_data_name,'file')
        load(ver_data_name);
        ver_tr_set = et_tr_set;
        ver_su_inds = find(et_tr_set > 96);
        ver_su_probeinds = find(ismember(1:length(su_probes),ver_tr_set(ver_su_inds)-96));
        ver_mu_inds = find(et_tr_set <= 96);
        ver_mu_probeinds = find(ismember(1:96,ver_tr_set(ver_mu_inds)));
        
        ver_before_mods = it_mods{1}(ver_tr_set(ver_su_inds));
        ver_after_mods = dit_xv_mods{2}(ver_tr_set(ver_su_inds));
        ver_init_LL_imps = it_LLimp(1,ver_tr_set(ver_su_inds));
        ver_all_LL_imps = it_LLimp(:,ver_tr_set(ver_su_inds));
        ver_fin_LL_imps = dit_xv_LLimp(2,ver_tr_set(ver_su_inds));
        ver_init_true_xvLL = it_truexvLLimp(1,ver_tr_set(ver_su_inds));
        ver_fin_true_xvLL = dit_truexvLLimp(2,ver_tr_set(ver_su_inds));
        
        ver_MUbefore_mods = it_mods{1}(ver_tr_set(ver_mu_inds));
        ver_MUafter_mods = dit_xvmods{2}(ver_tr_set(ver_mu_inds));
        ver_MUinit_LL_imps = it_xvLLimp(1,ver_tr_set(ver_mu_inds));
        ver_MUall_LL_imps = it_xvLLimp(:,ver_tr_set(ver_mu_inds));
        ver_MUfin_LL_imps = dit_xvLLimp(2,ver_tr_set(ver_mu_inds));
        ver_MUinit_true_xvLL = it_truexvLLimp(1,ver_tr_set(ver_mu_inds));
        ver_MUfin_true_xvLL = dit_truexvLLimp(2,ver_tr_set(ver_mu_inds));
    else
        ver_tr_set = [];
        ver_su_inds = [];
        ver_su_probeinds = [];
        ver_mu_inds = [];
        ver_mu_probeinds = [];
    end
    
    %% COMPUTE STATS ACROSS ALL FILTERS
    
    zpad_factor = 5;
    dx = 0.0565/et_params.spatial_usfac;
    dt = et_params.dt;
    flen = et_params.flen;
    use_nPix_us = et_params.use_nPix*et_params.spatial_usfac;
    Ft = linspace(-1/dt/2,1/dt/2,zpad_factor*flen);
    Fx = linspace(-1/dx/2,1/dx/2,zpad_factor*use_nPix_us);
    
    sig_Ft = 0.1; sig_Fx = 0.1;
    [FFx,FFt] = meshgrid(Fx,Ft);
    gauss_kern = FFt.^2/(2*sig_Ft^2) + FFx.^2/(2*sig_Fx^2);
    gauss_kern = exp(-gauss_kern);
    gauss_kern = gauss_kern/sum(gauss_kern(:));
    
    n_squared_filts = 2;
    
    %cycle through all SUs from this session
    for cc = 1:length(su_probes)
        fprintf('Processing stats for SU %d of %d\n',cc,length(su_probes));
        SU_DATA(tot_su_cnt).exptnum = expt_list(ee);
        SU_DATA(tot_su_cnt).probenum = su_probes(cc);
        
        cur = find(hor_su_probeinds == cc);
        if ~isempty(cur)
            SU_DATA(tot_su_cnt).hor_used = true;
            SU_DATA(tot_su_cnt).hor_init_LL_imp = hor_init_LL_imps(cur);
            SU_DATA(tot_su_cnt).hor_all_LL_imp = hor_all_LL_imps(:,cur);
            SU_DATA(tot_su_cnt).hor_fin_LL_imp = hor_fin_LL_imps(cur);
            SU_DATA(tot_su_cnt).hor_init_truexvLL = hor_init_true_xvLL(cur);
            SU_DATA(tot_su_cnt).hor_fin_truexvLL = hor_fin_true_xvLL(cur);
            
            Xtargs = [hor_before_mods(cur).mods(:).Xtarget];
            cor_filts = [hor_after_mods(cur).mods((Xtargs == 1)).filtK];
            uncor_filts = [hor_before_mods(cur).mods(Xtargs == 1).filtK];
            cor_zfilts = find(max(abs(cor_filts)) == 0);
            uncor_zfilts = find(max(abs(uncor_filts)) == 0);
            cor_filts = reshape(cor_filts,[flen use_nPix_us n_squared_filts+1]);
            uncor_filts = reshape(uncor_filts,[flen use_nPix_us n_squared_filts+1]);
            
            cor_spatial_profiles = squeeze(std(cor_filts));
            uncor_spatial_profiles = squeeze(std(uncor_filts));
            %gaussian fits to spatial profiles
            if max(cor_spatial_profiles(:,1)) > 0
                [fit_params,fit_z] = fitGaussianCurve((1:use_nPix_us)',cor_spatial_profiles(:,1));
                SU_DATA(tot_su_cnt).hor_cor_lin_mean = fit_params(1);
                SU_DATA(tot_su_cnt).hor_cor_lin_std = fit_params(2);
            else
                SU_DATA(tot_su_cnt).hor_cor_lin_mean = nan;
                SU_DATA(tot_su_cnt).hor_cor_lin_std = nan;
            end
            if max(uncor_spatial_profiles(:,1)) > 0
                [fit_params,fit_z] = fitGaussianCurve((1:use_nPix_us)',uncor_spatial_profiles(:,1));
                SU_DATA(tot_su_cnt).hor_uncor_lin_mean = fit_params(1);
                SU_DATA(tot_su_cnt).hor_uncor_lin_std = fit_params(2);
            else
                hor_uncor_lin_mean(cc) = nan;
                hor_uncor_lin_std(cc) = nan;
            end
            if max(max(cor_spatial_profiles(:,2:end))) > 0
                [fit_params,fit_z] = fitGaussianCurve((1:use_nPix_us)',mean(cor_spatial_profiles(:,2:end),2));
                SU_DATA(tot_su_cnt).hor_cor_sq_mean = fit_params(1);
                SU_DATA(tot_su_cnt).hor_cor_sq_std = fit_params(2);
            else
                SU_DATA(tot_su_cnt).hor_cor_sq_mean = nan;
                SU_DATA(tot_su_cnt).hor_cor_sq_std = nan;
            end
            if max(max(uncor_spatial_profiles(:,2:end))) > 0
                [fit_params,fit_z] = fitGaussianCurve((1:use_nPix_us)',mean(uncor_spatial_profiles(:,2:end),2));
                SU_DATA(tot_su_cnt).hor_uncor_sq_mean = fit_params(1);
                SU_DATA(tot_su_cnt).hor_uncor_sq_std = fit_params(2);
            else
                SU_DATA(tot_su_cnt).hor_uncor_sq_mean = nan;
                SU_DATA(tot_su_cnt).hor_uncor_sq_std = nan;
            end
            
            %FFT calculations
            cor_max_pow = nan(1,n_squared_filts+1);
            cor_max_ploc = nan(n_squared_filts+1,1);
            cor_tot_pow = nan(1,n_squared_filts+1);
            cor_filts_zpad = zeros(zpad_factor*flen,zpad_factor*use_nPix_us,n_squared_filts+1);
            cor_filts_zpad((flen+1):2*flen,(use_nPix_us+1):2*use_nPix_us,:) = cor_filts;
            cor_ffts = zeros(size(cor_filts_zpad));
            for ii = 1:3
                cur_ffts = abs(fftshift(fft2(cor_filts_zpad(:,:,ii))));
                cur_ffts = conv2(squeeze(cur_ffts),gauss_kern,'same');
                cor_max_pow(ii) = max(cur_ffts(:));
                cor_max_ploc(ii) = find(cur_ffts == cor_max_pow(ii),1);
                cor_tot_pow(ii) = sum(cur_ffts(:));
                cor_ffts(:,:,ii) = cur_ffts;
            end
            cor_FFx = abs(FFx(cor_max_ploc));
            cor_FFt = abs(FFt(cor_max_ploc));
            cor_max_pow(cor_zfilts) = nan;
            cor_tot_pow(cor_zfilts) = nan;
            cor_FFx(cor_zfilts) = nan; cor_FFt(cor_zfilts) = nan;
            
            SU_DATA(tot_su_cnt).hor_cor_lin_maxpow = cor_max_pow(1);
            SU_DATA(tot_su_cnt).hor_cor_sq_maxpow = nanmean(cor_max_pow(2:end));
            SU_DATA(tot_su_cnt).hor_cor_lin_totpow = cor_tot_pow(1);
            SU_DATA(tot_su_cnt).hor_cor_sq_totpow = nanmean(cor_tot_pow(2:end));
            SU_DATA(tot_su_cnt).hor_cor_lin_FFx = cor_FFx(1);
            SU_DATA(tot_su_cnt).hor_cor_sq_FFx = nanmean(cor_FFx(2:end));
            SU_DATA(tot_su_cnt).hor_cor_lin_FFt = cor_FFt(1);
            SU_DATA(tot_su_cnt).hor_cor_sq_FFt = nanmean(cor_FFt(2:end));
            
            uncor_max_pow = zeros(n_squared_filts+1,1);
            uncor_max_ploc = zeros(n_squared_filts+1,1);
            uncor_tot_pow = zeros(n_squared_filts+1,1);
            uncor_filts_zpad = zeros(zpad_factor*flen,zpad_factor*use_nPix_us,n_squared_filts+1);
            uncor_filts_zpad((flen+1):2*flen,(use_nPix_us+1):2*use_nPix_us,:) = uncor_filts;
            uncor_ffts = zeros(size(uncor_filts_zpad));
            for ii = 1:n_squared_filts+1
                cur_ffts = abs(fftshift(fft2(uncor_filts_zpad(:,:,ii))));
                cur_ffts = conv2(squeeze(cur_ffts),gauss_kern,'same');
                uncor_max_pow(ii) = max(cur_ffts(:));
                uncor_max_ploc(ii) = find(cur_ffts == uncor_max_pow(ii),1);
                uncor_tot_pow(ii) = sum(cur_ffts(:));
                uncor_ffts(:,:,ii) = cur_ffts;
            end
            uncor_FFx = abs(FFx(uncor_max_ploc));
            uncor_FFt = abs(FFt(uncor_max_ploc));
            uncor_max_pow(uncor_zfilts) = nan;
            uncor_tot_pow(uncor_zfilts) = nan;
            uncor_FFx(uncor_zfilts) = nan; uncor_FFt(uncor_zfilts) = nan;
            
            SU_DATA(tot_su_cnt).hor_uncor_lin_maxpow = uncor_max_pow(1);
            SU_DATA(tot_su_cnt).hor_uncor_sq_maxpow = nanmean(uncor_max_pow(2:end));
            SU_DATA(tot_su_cnt).hor_uncor_lin_totpow = uncor_tot_pow(1);
            SU_DATA(tot_su_cnt).hor_uncor_sq_totpow = nanmean(uncor_tot_pow(2:end));
            SU_DATA(tot_su_cnt).hor_uncor_lin_FFx = uncor_FFx(1);
            SU_DATA(tot_su_cnt).hor_uncor_sq_FFx = nanmean(uncor_FFx(2:end));
            SU_DATA(tot_su_cnt).hor_uncor_lin_FFt = uncor_FFt(1);
            SU_DATA(tot_su_cnt).hor_uncor_sq_FFt = nanmean(uncor_FFt(2:end));
        else
            SU_DATA(tot_su_cnt).hor_used = false;
            SU_DATA(tot_su_cnt).hor_init_LL_imp = nan;
            SU_DATA(tot_su_cnt).hor_init_truexvLL = nan;
            SU_DATA(tot_su_cnt).hor_fin_LL_imp = nan;
        end
        
        cur = find(ver_su_probeinds == cc);
        if ~isempty(cur)
            SU_DATA(tot_su_cnt).ver_used = true;
            SU_DATA(tot_su_cnt).ver_init_LL_imp = ver_init_LL_imps(cur);
            SU_DATA(tot_su_cnt).ver_all_LL_imp = ver_all_LL_imps(:,cur);
            SU_DATA(tot_su_cnt).ver_fin_LL_imp = ver_fin_LL_imps(cur);
            SU_DATA(tot_su_cnt).ver_init_truexvLL = ver_init_true_xvLL(cur);
            SU_DATA(tot_su_cnt).ver_fin_truexvLL = ver_fin_true_xvLL(cur);
            
            Xtargs = [ver_before_mods(cur).mods(:).Xtarget];
            cor_filts = [ver_after_mods(cur).mods((Xtargs == 1)).filtK];
            uncor_filts = [ver_before_mods(cur).mods(Xtargs == 1).filtK];
            cor_zfilts = find(max(abs(cor_filts)) == 0);
            uncor_zfilts = find(max(abs(uncor_filts)) == 0);
            cor_filts = reshape(cor_filts,[flen use_nPix_us n_squared_filts+1]);
            uncor_filts = reshape(uncor_filts,[flen use_nPix_us n_squared_filts+1]);
            
            cor_spatial_profiles = squeeze(std(cor_filts));
            uncor_spatial_profiles = squeeze(std(uncor_filts));
            %gaussian fits to spatial profiles
            if max(cor_spatial_profiles(:,1)) > 0
                [fit_params,fit_z] = fitGaussianCurve((1:use_nPix_us)',cor_spatial_profiles(:,1));
                SU_DATA(tot_su_cnt).ver_cor_lin_mean = fit_params(1);
                SU_DATA(tot_su_cnt).ver_cor_lin_std = fit_params(2);
            else
                SU_DATA(tot_su_cnt).ver_cor_lin_mean = nan;
                SU_DATA(tot_su_cnt).ver_cor_lin_std = nan;
            end
            if max(uncor_spatial_profiles(:,1)) > 0
                [fit_params,fit_z] = fitGaussianCurve((1:use_nPix_us)',uncor_spatial_profiles(:,1));
                SU_DATA(tot_su_cnt).ver_uncor_lin_mean = fit_params(1);
                SU_DATA(tot_su_cnt).ver_uncor_lin_std = fit_params(2);
            else
                ver_uncor_lin_mean(cc) = nan;
                ver_uncor_lin_std(cc) = nan;
            end
            if max(max(cor_spatial_profiles(:,2:end))) > 0
                [fit_params,fit_z] = fitGaussianCurve((1:use_nPix_us)',mean(cor_spatial_profiles(:,2:end),2));
                SU_DATA(tot_su_cnt).ver_cor_sq_mean = fit_params(1);
                SU_DATA(tot_su_cnt).ver_cor_sq_std = fit_params(2);
            else
                SU_DATA(tot_su_cnt).ver_cor_sq_mean = nan;
                SU_DATA(tot_su_cnt).ver_cor_sq_std = nan;
            end
            if max(max(uncor_spatial_profiles(:,2:end))) > 0
                [fit_params,fit_z] = fitGaussianCurve((1:use_nPix_us)',mean(uncor_spatial_profiles(:,2:end),2));
                SU_DATA(tot_su_cnt).ver_uncor_sq_mean = fit_params(1);
                SU_DATA(tot_su_cnt).ver_uncor_sq_std = fit_params(2);
            else
                SU_DATA(tot_su_cnt).ver_uncor_sq_mean = nan;
                SU_DATA(tot_su_cnt).ver_uncor_sq_std = nan;
            end
            
            %FFT calculations
            cor_max_pow = nan(1,n_squared_filts+1);
            cor_max_ploc = nan(n_squared_filts+1,1);
            cor_tot_pow = nan(1,n_squared_filts+1);
            cor_filts_zpad = zeros(zpad_factor*flen,zpad_factor*use_nPix_us,n_squared_filts+1);
            cor_filts_zpad((flen+1):2*flen,(use_nPix_us+1):2*use_nPix_us,:) = cor_filts;
            cor_ffts = zeros(size(cor_filts_zpad));
            for ii = 1:3
                cur_ffts = abs(fftshift(fft2(cor_filts_zpad(:,:,ii))));
                cur_ffts = conv2(squeeze(cur_ffts),gauss_kern,'same');
                cor_max_pow(ii) = max(cur_ffts(:));
                cor_max_ploc(ii) = find(cur_ffts == cor_max_pow(ii),1);
                cor_tot_pow(ii) = sum(cur_ffts(:));
                cor_ffts(:,:,ii) = cur_ffts;
            end
            cor_FFx = abs(FFx(cor_max_ploc));
            cor_FFt = abs(FFt(cor_max_ploc));
            cor_max_pow(cor_zfilts) = nan;
            cor_tot_pow(cor_zfilts) = nan;
            cor_FFx(cor_zfilts) = nan; cor_FFt(cor_zfilts) = nan;
            
            SU_DATA(tot_su_cnt).ver_cor_lin_maxpow = cor_max_pow(1);
            SU_DATA(tot_su_cnt).ver_cor_sq_maxpow = nanmean(cor_max_pow(2:end));
            SU_DATA(tot_su_cnt).ver_cor_lin_totpow = cor_tot_pow(1);
            SU_DATA(tot_su_cnt).ver_cor_sq_totpow = nanmean(cor_tot_pow(2:end));
            SU_DATA(tot_su_cnt).ver_cor_lin_FFx = cor_FFx(1);
            SU_DATA(tot_su_cnt).ver_cor_sq_FFx = nanmean(cor_FFx(2:end));
            SU_DATA(tot_su_cnt).ver_cor_lin_FFt = cor_FFt(1);
            SU_DATA(tot_su_cnt).ver_cor_sq_FFt = nanmean(cor_FFt(2:end));
            
            uncor_max_pow = zeros(n_squared_filts+1,1);
            uncor_max_ploc = zeros(n_squared_filts+1,1);
            uncor_tot_pow = zeros(n_squared_filts+1,1);
            uncor_filts_zpad = zeros(zpad_factor*flen,zpad_factor*use_nPix_us,n_squared_filts+1);
            uncor_filts_zpad((flen+1):2*flen,(use_nPix_us+1):2*use_nPix_us,:) = uncor_filts;
            uncor_ffts = zeros(size(uncor_filts_zpad));
            for ii = 1:n_squared_filts+1
                cur_ffts = abs(fftshift(fft2(uncor_filts_zpad(:,:,ii))));
                cur_ffts = conv2(squeeze(cur_ffts),gauss_kern,'same');
                uncor_max_pow(ii) = max(cur_ffts(:));
                uncor_max_ploc(ii) = find(cur_ffts == uncor_max_pow(ii),1);
                uncor_tot_pow(ii) = sum(cur_ffts(:));
                uncor_ffts(:,:,ii) = cur_ffts;
            end
            uncor_FFx = abs(FFx(uncor_max_ploc));
            uncor_FFt = abs(FFt(uncor_max_ploc));
            uncor_max_pow(uncor_zfilts) = nan;
            uncor_tot_pow(uncor_zfilts) = nan;
            uncor_FFx(uncor_zfilts) = nan; uncor_FFt(uncor_zfilts) = nan;
            
            SU_DATA(tot_su_cnt).ver_uncor_lin_maxpow = uncor_max_pow(1);
            SU_DATA(tot_su_cnt).ver_uncor_sq_maxpow = nanmean(uncor_max_pow(2:end));
            SU_DATA(tot_su_cnt).ver_uncor_lin_totpow = uncor_tot_pow(1);
            SU_DATA(tot_su_cnt).ver_uncor_sq_totpow = nanmean(uncor_tot_pow(2:end));
            SU_DATA(tot_su_cnt).ver_uncor_lin_FFx = uncor_FFx(1);
            SU_DATA(tot_su_cnt).ver_uncor_sq_FFx = nanmean(uncor_FFx(2:end));
            SU_DATA(tot_su_cnt).ver_uncor_lin_FFt = uncor_FFt(1);
            SU_DATA(tot_su_cnt).ver_uncor_sq_FFt = nanmean(uncor_FFt(2:end));
        else
            SU_DATA(tot_su_cnt).ver_used = false;
            SU_DATA(tot_su_cnt).ver_init_LL_imp = nan;
            SU_DATA(tot_su_cnt).ver_init_truexvLL = nan;
            SU_DATA(tot_su_cnt).ver_fin_LL_imp = nan;
        end
        
        %which orientation gives better LL for each SU
        hor_imp = SU_DATA(tot_su_cnt).hor_fin_LL_imp;
        ver_imp = SU_DATA(tot_su_cnt).ver_fin_LL_imp;
        if isnan(hor_imp) && isnan(ver_imp)
            SU_DATA(tot_su_cnt).better_ori = nan;
        else
            [a,b] = max([hor_imp ver_imp]);
            SU_DATA(tot_su_cnt).better_ori = b;
        end
        tot_su_cnt = tot_su_cnt + 1;
    end
    
    %% PROCESS MUs
    
    %cycle through all MUs from this session
    for cc = 1:96
        fprintf('Processing stats for MU %d of %d\n',cc,96);
        MU_DATA(tot_mu_cnt).exptnum = expt_list(ee);
        MU_DATA(tot_mu_cnt).probenum = cc;
        
        cur = find(hor_mu_probeinds == cc);
        if ~isempty(cur)
            MU_DATA(tot_mu_cnt).hor_used = true;
            MU_DATA(tot_mu_cnt).hor_init_LL_imp = hor_MUinit_LL_imps(cur);
            MU_DATA(tot_mu_cnt).hor_all_LL_imp = hor_MUall_LL_imps(:,cur);
            MU_DATA(tot_mu_cnt).hor_fin_LL_imp = hor_MUfin_LL_imps(cur);
            MU_DATA(tot_mu_cnt).hor_init_truexvLL = hor_MUinit_true_xvLL(cur);
            MU_DATA(tot_mu_cnt).hor_fin_truexvLL = hor_MUfin_true_xvLL(cur);
            
            Xtargs = [hor_MUbefore_mods(cur).mods(:).Xtarget];
            cor_filts = [hor_MUafter_mods(cur).mods(Xtargs == 1).filtK];
            uncor_filts = [hor_MUbefore_mods(cur).mods(Xtargs == 1).filtK];
            cor_zfilts = find(max(abs(cor_filts)) == 0);
            uncor_zfilts = find(max(abs(uncor_filts)) == 0);
            cor_filts = reshape(cor_filts,[flen use_nPix_us n_squared_filts+1]);
            uncor_filts = reshape(uncor_filts,[flen use_nPix_us n_squared_filts+1]);
            
            cor_spatial_profiles = squeeze(std(cor_filts));
            uncor_spatial_profiles = squeeze(std(uncor_filts));
            %gaussian fits to spatial profiles
            if max(cor_spatial_profiles(:,1)) > 0
                [fit_params,fit_z] = fitGaussianCurve((1:use_nPix_us)',cor_spatial_profiles(:,1));
                MU_DATA(tot_mu_cnt).hor_cor_lin_mean = fit_params(1);
                MU_DATA(tot_mu_cnt).hor_cor_lin_std = fit_params(2);
            else
                MU_DATA(tot_mu_cnt).hor_cor_lin_mean = nan;
                MU_DATA(tot_mu_cnt).hor_cor_lin_std = nan;
            end
            if max(uncor_spatial_profiles(:,1)) > 0
                [fit_params,fit_z] = fitGaussianCurve((1:use_nPix_us)',uncor_spatial_profiles(:,1));
                MU_DATA(tot_mu_cnt).hor_uncor_lin_mean = fit_params(1);
                MU_DATA(tot_mu_cnt).hor_uncor_lin_std = fit_params(2);
            else
                hor_uncor_lin_mean(cc) = nan;
                hor_uncor_lin_std(cc) = nan;
            end
            if max(max(cor_spatial_profiles(:,2:end))) > 0
                [fit_params,fit_z] = fitGaussianCurve((1:use_nPix_us)',mean(cor_spatial_profiles(:,2:end),2));
                MU_DATA(tot_mu_cnt).hor_cor_sq_mean = fit_params(1);
                MU_DATA(tot_mu_cnt).hor_cor_sq_std = fit_params(2);
            else
                MU_DATA(tot_mu_cnt).hor_cor_sq_mean = nan;
                MU_DATA(tot_mu_cnt).hor_cor_sq_std = nan;
            end
            if max(max(uncor_spatial_profiles(:,2:end))) > 0
                [fit_params,fit_z] = fitGaussianCurve((1:use_nPix_us)',mean(uncor_spatial_profiles(:,2:end),2));
                MU_DATA(tot_mu_cnt).hor_uncor_sq_mean = fit_params(1);
                MU_DATA(tot_mu_cnt).hor_uncor_sq_std = fit_params(2);
            else
                MU_DATA(tot_mu_cnt).hor_uncor_sq_mean = nan;
                MU_DATA(tot_mu_cnt).hor_uncor_sq_std = nan;
            end
            
            %FFT calculations
            cor_max_pow = nan(1,n_squared_filts+1);
            cor_max_ploc = nan(n_squared_filts+1,1);
            cor_tot_pow = nan(1,n_squared_filts+1);
            cor_filts_zpad = zeros(zpad_factor*flen,zpad_factor*use_nPix_us,n_squared_filts+1);
            cor_filts_zpad((flen+1):2*flen,(use_nPix_us+1):2*use_nPix_us,:) = cor_filts;
            cor_ffts = zeros(size(cor_filts_zpad));
            for ii = 1:3
                cur_ffts = abs(fftshift(fft2(cor_filts_zpad(:,:,ii))));
                cur_ffts = conv2(squeeze(cur_ffts),gauss_kern,'same');
                cor_max_pow(ii) = max(cur_ffts(:));
                cor_max_ploc(ii) = find(cur_ffts == cor_max_pow(ii),1);
                cor_tot_pow(ii) = sum(cur_ffts(:));
                cor_ffts(:,:,ii) = cur_ffts;
            end
            cor_FFx = abs(FFx(cor_max_ploc));
            cor_FFt = abs(FFt(cor_max_ploc));
            cor_max_pow(cor_zfilts) = nan;
            cor_tot_pow(cor_zfilts) = nan;
            cor_FFx(cor_zfilts) = nan; cor_FFt(cor_zfilts) = nan;
            
            MU_DATA(tot_mu_cnt).hor_cor_lin_maxpow = cor_max_pow(1);
            MU_DATA(tot_mu_cnt).hor_cor_sq_maxpow = nanmean(cor_max_pow(2:end));
            MU_DATA(tot_mu_cnt).hor_cor_lin_totpow = cor_tot_pow(1);
            MU_DATA(tot_mu_cnt).hor_cor_sq_totpow = nanmean(cor_tot_pow(2:end));
            MU_DATA(tot_mu_cnt).hor_cor_lin_FFx = cor_FFx(1);
            MU_DATA(tot_mu_cnt).hor_cor_sq_FFx = nanmean(cor_FFx(2:end));
            MU_DATA(tot_mu_cnt).hor_cor_lin_FFt = cor_FFt(1);
            MU_DATA(tot_mu_cnt).hor_cor_sq_FFt = nanmean(cor_FFt(2:end));
            
            uncor_max_pow = zeros(n_squared_filts+1,1);
            uncor_max_ploc = zeros(n_squared_filts+1,1);
            uncor_tot_pow = zeros(n_squared_filts+1,1);
            uncor_filts_zpad = zeros(zpad_factor*flen,zpad_factor*use_nPix_us,n_squared_filts+1);
            uncor_filts_zpad((flen+1):2*flen,(use_nPix_us+1):2*use_nPix_us,:) = uncor_filts;
            uncor_ffts = zeros(size(uncor_filts_zpad));
            for ii = 1:n_squared_filts+1
                cur_ffts = abs(fftshift(fft2(uncor_filts_zpad(:,:,ii))));
                cur_ffts = conv2(squeeze(cur_ffts),gauss_kern,'same');
                uncor_max_pow(ii) = max(cur_ffts(:));
                uncor_max_ploc(ii) = find(cur_ffts == uncor_max_pow(ii),1);
                uncor_tot_pow(ii) = sum(cur_ffts(:));
                uncor_ffts(:,:,ii) = cur_ffts;
            end
            uncor_FFx = abs(FFx(uncor_max_ploc));
            uncor_FFt = abs(FFt(uncor_max_ploc));
            uncor_max_pow(uncor_zfilts) = nan;
            uncor_tot_pow(uncor_zfilts) = nan;
            uncor_FFx(uncor_zfilts) = nan; uncor_FFt(uncor_zfilts) = nan;
            
            MU_DATA(tot_mu_cnt).hor_uncor_lin_maxpow = uncor_max_pow(1);
            MU_DATA(tot_mu_cnt).hor_uncor_sq_maxpow = nanmean(uncor_max_pow(2:end));
            MU_DATA(tot_mu_cnt).hor_uncor_lin_totpow = uncor_tot_pow(1);
            MU_DATA(tot_mu_cnt).hor_uncor_sq_totpow = nanmean(uncor_tot_pow(2:end));
            MU_DATA(tot_mu_cnt).hor_uncor_lin_FFx = uncor_FFx(1);
            MU_DATA(tot_mu_cnt).hor_uncor_sq_FFx = nanmean(uncor_FFx(2:end));
            MU_DATA(tot_mu_cnt).hor_uncor_lin_FFt = uncor_FFt(1);
            MU_DATA(tot_mu_cnt).hor_uncor_sq_FFt = nanmean(uncor_FFt(2:end));
        else
            MU_DATA(tot_mu_cnt).hor_used = false;
            MU_DATA(tot_mu_cnt).hor_init_truexvLL = nan;
            MU_DATA(tot_mu_cnt).hor_init_LL_imp = nan;
            MU_DATA(tot_mu_cnt).hor_fin_LL_imp = nan;
        end
        
        cur = find(ver_mu_probeinds == cc);
        if ~isempty(cur)
            MU_DATA(tot_mu_cnt).ver_used = true;
            MU_DATA(tot_mu_cnt).ver_init_LL_imp = ver_MUinit_LL_imps(cur);
            MU_DATA(tot_mu_cnt).ver_all_LL_imp = ver_MUall_LL_imps(:,cur);
            MU_DATA(tot_mu_cnt).ver_fin_LL_imp = ver_MUfin_LL_imps(cur);
            MU_DATA(tot_mu_cnt).ver_init_truexvLL = ver_MUinit_true_xvLL(cur);
            MU_DATA(tot_mu_cnt).ver_fin_truexvLL = ver_MUfin_true_xvLL(cur);
            
            Xtargs = [ver_MUbefore_mods(cur).mods(:).Xtarget];
            cor_filts = [ver_MUafter_mods(cur).mods((Xtargs == 1)).filtK];
            uncor_filts = [ver_MUbefore_mods(cur).mods(Xtargs == 1).filtK];
            cor_zfilts = find(max(abs(cor_filts)) == 0);
            uncor_zfilts = find(max(abs(uncor_filts)) == 0);
            cor_filts = reshape(cor_filts,[flen use_nPix_us n_squared_filts+1]);
            uncor_filts = reshape(uncor_filts,[flen use_nPix_us n_squared_filts+1]);
            
            cor_spatial_profiles = squeeze(std(cor_filts));
            uncor_spatial_profiles = squeeze(std(uncor_filts));
            %gaussian fits to spatial profiles
            if max(cor_spatial_profiles(:,1)) > 0
                [fit_params,fit_z] = fitGaussianCurve((1:use_nPix_us)',cor_spatial_profiles(:,1));
                MU_DATA(tot_mu_cnt).ver_cor_lin_mean = fit_params(1);
                MU_DATA(tot_mu_cnt).ver_cor_lin_std = fit_params(2);
            else
                MU_DATA(tot_mu_cnt).ver_cor_lin_mean = nan;
                MU_DATA(tot_mu_cnt).ver_cor_lin_std = nan;
            end
            if max(uncor_spatial_profiles(:,1)) > 0
                [fit_params,fit_z] = fitGaussianCurve((1:use_nPix_us)',uncor_spatial_profiles(:,1));
                MU_DATA(tot_mu_cnt).ver_uncor_lin_mean = fit_params(1);
                MU_DATA(tot_mu_cnt).ver_uncor_lin_std = fit_params(2);
            else
                ver_uncor_lin_mean(cc) = nan;
                ver_uncor_lin_std(cc) = nan;
            end
            if max(max(cor_spatial_profiles(:,2:end))) > 0
                [fit_params,fit_z] = fitGaussianCurve((1:use_nPix_us)',mean(cor_spatial_profiles(:,2:end),2));
                MU_DATA(tot_mu_cnt).ver_cor_sq_mean = fit_params(1);
                MU_DATA(tot_mu_cnt).ver_cor_sq_std = fit_params(2);
            else
                MU_DATA(tot_mu_cnt).ver_cor_sq_mean = nan;
                MU_DATA(tot_mu_cnt).ver_cor_sq_std = nan;
            end
            if max(max(uncor_spatial_profiles(:,2:end))) > 0
                [fit_params,fit_z] = fitGaussianCurve((1:use_nPix_us)',mean(uncor_spatial_profiles(:,2:end),2));
                MU_DATA(tot_mu_cnt).ver_uncor_sq_mean = fit_params(1);
                MU_DATA(tot_mu_cnt).ver_uncor_sq_std = fit_params(2);
            else
                MU_DATA(tot_mu_cnt).ver_uncor_sq_mean = nan;
                MU_DATA(tot_mu_cnt).ver_uncor_sq_std = nan;
            end
            
            %FFT calculations
            cor_max_pow = nan(1,n_squared_filts+1);
            cor_max_ploc = nan(n_squared_filts+1,1);
            cor_tot_pow = nan(1,n_squared_filts+1);
            cor_filts_zpad = zeros(zpad_factor*flen,zpad_factor*use_nPix_us,n_squared_filts+1);
            cor_filts_zpad((flen+1):2*flen,(use_nPix_us+1):2*use_nPix_us,:) = cor_filts;
            cor_ffts = zeros(size(cor_filts_zpad));
            for ii = 1:3
                cur_ffts = abs(fftshift(fft2(cor_filts_zpad(:,:,ii))));
                cur_ffts = conv2(squeeze(cur_ffts),gauss_kern,'same');
                cor_max_pow(ii) = max(cur_ffts(:));
                cor_max_ploc(ii) = find(cur_ffts == cor_max_pow(ii),1);
                cor_tot_pow(ii) = sum(cur_ffts(:));
                cor_ffts(:,:,ii) = cur_ffts;
            end
            cor_FFx = abs(FFx(cor_max_ploc));
            cor_FFt = abs(FFt(cor_max_ploc));
            cor_max_pow(cor_zfilts) = nan;
            cor_tot_pow(cor_zfilts) = nan;
            cor_FFx(cor_zfilts) = nan; cor_FFt(cor_zfilts) = nan;
            
            MU_DATA(tot_mu_cnt).ver_cor_lin_maxpow = cor_max_pow(1);
            MU_DATA(tot_mu_cnt).ver_cor_sq_maxpow = nanmean(cor_max_pow(2:end));
            MU_DATA(tot_mu_cnt).ver_cor_lin_totpow = cor_tot_pow(1);
            MU_DATA(tot_mu_cnt).ver_cor_sq_totpow = nanmean(cor_tot_pow(2:end));
            MU_DATA(tot_mu_cnt).ver_cor_lin_FFx = cor_FFx(1);
            MU_DATA(tot_mu_cnt).ver_cor_sq_FFx = nanmean(cor_FFx(2:end));
            MU_DATA(tot_mu_cnt).ver_cor_lin_FFt = cor_FFt(1);
            MU_DATA(tot_mu_cnt).ver_cor_sq_FFt = nanmean(cor_FFt(2:end));
            
            uncor_max_pow = zeros(n_squared_filts+1,1);
            uncor_max_ploc = zeros(n_squared_filts+1,1);
            uncor_tot_pow = zeros(n_squared_filts+1,1);
            uncor_filts_zpad = zeros(zpad_factor*flen,zpad_factor*use_nPix_us,n_squared_filts+1);
            uncor_filts_zpad((flen+1):2*flen,(use_nPix_us+1):2*use_nPix_us,:) = uncor_filts;
            uncor_ffts = zeros(size(uncor_filts_zpad));
            for ii = 1:n_squared_filts+1
                cur_ffts = abs(fftshift(fft2(uncor_filts_zpad(:,:,ii))));
                cur_ffts = conv2(squeeze(cur_ffts),gauss_kern,'same');
                uncor_max_pow(ii) = max(cur_ffts(:));
                uncor_max_ploc(ii) = find(cur_ffts == uncor_max_pow(ii),1);
                uncor_tot_pow(ii) = sum(cur_ffts(:));
                uncor_ffts(:,:,ii) = cur_ffts;
            end
            uncor_FFx = abs(FFx(uncor_max_ploc));
            uncor_FFt = abs(FFt(uncor_max_ploc));
            uncor_max_pow(uncor_zfilts) = nan;
            uncor_tot_pow(uncor_zfilts) = nan;
            uncor_FFx(uncor_zfilts) = nan; uncor_FFt(uncor_zfilts) = nan;
            
            MU_DATA(tot_mu_cnt).ver_uncor_lin_maxpow = uncor_max_pow(1);
            MU_DATA(tot_mu_cnt).ver_uncor_sq_maxpow = nanmean(uncor_max_pow(2:end));
            MU_DATA(tot_mu_cnt).ver_uncor_lin_totpow = uncor_tot_pow(1);
            MU_DATA(tot_mu_cnt).ver_uncor_sq_totpow = nanmean(uncor_tot_pow(2:end));
            MU_DATA(tot_mu_cnt).ver_uncor_lin_FFx = uncor_FFx(1);
            MU_DATA(tot_mu_cnt).ver_uncor_sq_FFx = nanmean(uncor_FFx(2:end));
            MU_DATA(tot_mu_cnt).ver_uncor_lin_FFt = uncor_FFt(1);
            MU_DATA(tot_mu_cnt).ver_uncor_sq_FFt = nanmean(uncor_FFt(2:end));
        else
            MU_DATA(tot_mu_cnt).ver_used = false;
            MU_DATA(tot_mu_cnt).ver_init_truexvLL = nan;
            MU_DATA(tot_mu_cnt).ver_init_LL_imp = nan;
            MU_DATA(tot_mu_cnt).ver_fin_LL_imp = nan;
        end
        
        %which orientation gives better LL for each SU
        hor_imp = MU_DATA(tot_mu_cnt).hor_fin_LL_imp;
        ver_imp = MU_DATA(tot_mu_cnt).ver_fin_LL_imp;
        if isnan(hor_imp) && isnan(ver_imp)
            MU_DATA(tot_mu_cnt).better_ori = nan;
        else
            [a,b] = max([hor_imp ver_imp]);
            MU_DATA(tot_mu_cnt).better_ori = b;
        end
        tot_mu_cnt = tot_mu_cnt + 1;
    end
    
end

%%
best_oris = [SU_DATA(:).better_ori];

su_LL_before = nan(size(best_oris));
su_LL_before(best_oris == 1) = [SU_DATA(best_oris == 1).hor_init_LL_imp];
su_LL_before(best_oris == 2) = [SU_DATA(best_oris == 2).ver_init_LL_imp];
su_LL_after = nan(size(best_oris));
su_LL_after(best_oris == 1) = [SU_DATA(best_oris == 1).hor_fin_LL_imp];
su_LL_after(best_oris == 2) = [SU_DATA(best_oris == 2).ver_fin_LL_imp];
su_trueLL_before(best_oris == 1) = [SU_DATA(best_oris == 1).hor_init_truexvLL];
su_trueLL_before(best_oris == 2) = [SU_DATA(best_oris == 2).ver_init_truexvLL];
su_trueLL_after(best_oris == 1) = [SU_DATA(best_oris == 1).hor_fin_truexvLL];
su_trueLL_after(best_oris == 2) = [SU_DATA(best_oris == 2).ver_fin_truexvLL];

su_used_hori = find([SU_DATA(:).hor_init_truexvLL] > 0);
su_used_ver = find([SU_DATA(:).ver_init_truexvLL] > 0);
su_all_LL_hor = nan(7,length(best_oris));
su_all_LL_hor(:,su_used_hori) = [SU_DATA(su_used_hori).hor_all_LL_imp; SU_DATA(su_used_hori).hor_fin_LL_imp];
su_all_LL_hor(:,su_used_hori) = bsxfun(@rdivide,su_all_LL_hor(:,su_used_hori),[SU_DATA(su_used_hori).hor_init_LL_imp]);

su_all_LL_ver = nan(7,length(best_oris));
su_all_LL_ver(:,su_used_ver) = [SU_DATA(su_used_ver).ver_all_LL_imp; SU_DATA(su_used_ver).ver_fin_LL_imp];
su_all_LL_ver(:,su_used_ver) = bsxfun(@rdivide,su_all_LL_ver(:,su_used_ver),[SU_DATA(su_used_ver).ver_init_LL_imp]);

su_RFmean_before = nan(length(best_oris),2);
su_RFmean_before(best_oris == 1,:) = [SU_DATA(best_oris == 1).hor_uncor_lin_mean; SU_DATA(best_oris == 1).hor_uncor_sq_mean]';
su_RFmean_before(best_oris == 2,:) = [SU_DATA(best_oris == 2).ver_uncor_lin_mean; SU_DATA(best_oris == 2).ver_uncor_sq_mean]';
su_RFmean_after = nan(length(best_oris),2);
su_RFmean_after(best_oris == 1,:) = [SU_DATA(best_oris == 1).hor_cor_lin_mean; SU_DATA(best_oris == 1).hor_cor_sq_mean]';
su_RFmean_after(best_oris == 2,:) = [SU_DATA(best_oris == 2).ver_cor_lin_mean; SU_DATA(best_oris == 2).ver_cor_sq_mean]';

su_RFstd_before = nan(length(best_oris),2);
su_RFstd_before(best_oris == 1,:) = [SU_DATA(best_oris == 1).hor_uncor_lin_std; SU_DATA(best_oris == 1).hor_uncor_sq_std]';
su_RFstd_before(best_oris == 2,:) = [SU_DATA(best_oris == 2).ver_uncor_lin_std; SU_DATA(best_oris == 2).ver_uncor_sq_std]';
su_RFstd_after = nan(length(best_oris),2);
su_RFstd_after(best_oris == 1,:) = [SU_DATA(best_oris == 1).hor_cor_lin_std; SU_DATA(best_oris == 1).hor_cor_sq_std]';
su_RFstd_after(best_oris == 2,:) = [SU_DATA(best_oris == 2).ver_cor_lin_std; SU_DATA(best_oris == 2).ver_cor_sq_std]';

su_tpow_before = nan(length(best_oris),2);
su_tpow_before(best_oris == 1,:) = [SU_DATA(best_oris == 1).hor_uncor_lin_totpow; SU_DATA(best_oris == 1).hor_uncor_sq_totpow]';
su_tpow_before(best_oris == 2,:) = [SU_DATA(best_oris == 2).ver_uncor_lin_totpow; SU_DATA(best_oris == 2).ver_uncor_sq_totpow]';
su_tpow_after = nan(length(best_oris),2);
su_tpow_after(best_oris == 1,:) = [SU_DATA(best_oris == 1).hor_cor_lin_totpow; SU_DATA(best_oris == 1).hor_cor_sq_totpow]';
su_tpow_after(best_oris == 2,:) = [SU_DATA(best_oris == 2).ver_cor_lin_totpow; SU_DATA(best_oris == 2).ver_cor_sq_totpow]';

su_FFx_before = nan(length(best_oris),2);
su_FFx_before(best_oris == 1,:) = [SU_DATA(best_oris == 1).hor_uncor_lin_FFx; SU_DATA(best_oris == 1).hor_uncor_sq_FFx]';
su_FFx_before(best_oris == 2,:) = [SU_DATA(best_oris == 2).ver_uncor_lin_FFx; SU_DATA(best_oris == 2).ver_uncor_sq_FFx]';
su_FFx_after = nan(length(best_oris),2);
su_FFx_after(best_oris == 1,:) = [SU_DATA(best_oris == 1).hor_cor_lin_FFx; SU_DATA(best_oris == 1).hor_cor_sq_FFx]';
su_FFx_after(best_oris == 2,:) = [SU_DATA(best_oris == 2).ver_cor_lin_FFx; SU_DATA(best_oris == 2).ver_cor_sq_FFx]';

su_FFt_before = nan(length(best_oris),2);
su_FFt_before(best_oris == 1,:) = [SU_DATA(best_oris == 1).hor_uncor_lin_FFt; SU_DATA(best_oris == 1).hor_uncor_sq_FFt]';
su_FFt_before(best_oris == 2,:) = [SU_DATA(best_oris == 2).ver_uncor_lin_FFt; SU_DATA(best_oris == 2).ver_uncor_sq_FFt]';
su_FFt_after = nan(length(best_oris),2);
su_FFt_after(best_oris == 1,:) = [SU_DATA(best_oris == 1).hor_cor_lin_FFt; SU_DATA(best_oris == 1).hor_cor_sq_FFt]';
su_FFt_after(best_oris == 2,:) = [SU_DATA(best_oris == 2).ver_cor_lin_FFt; SU_DATA(best_oris == 2).ver_cor_sq_FFt]';

su_used = find(su_trueLL_before > 0);
su_used_linfilts = su_used(su_tpow_after(su_used,1) > 1000);
su_used_sqfilts = su_used(su_tpow_after(su_used,2) > 1000);


%%
best_oris = [MU_DATA(:).better_ori];

mu_LL_before = nan(size(best_oris));
mu_LL_before(best_oris == 1) = [MU_DATA(best_oris == 1).hor_init_LL_imp];
mu_LL_before(best_oris == 2) = [MU_DATA(best_oris == 2).ver_init_LL_imp];
mu_LL_after = nan(size(best_oris));
mu_LL_after(best_oris == 1) = [MU_DATA(best_oris == 1).hor_fin_LL_imp];
mu_LL_after(best_oris == 2) = [MU_DATA(best_oris == 2).ver_fin_LL_imp];
mu_trueLL_before(best_oris == 1) = [MU_DATA(best_oris == 1).hor_init_truexvLL];
mu_trueLL_before(best_oris == 2) = [MU_DATA(best_oris == 2).ver_init_truexvLL];
mu_trueLL_after(best_oris == 1) = [MU_DATA(best_oris == 1).hor_fin_truexvLL];
mu_trueLL_after(best_oris == 2) = [MU_DATA(best_oris == 2).ver_fin_truexvLL];

used_hori = find([MU_DATA(:).hor_init_truexvLL] > 0);
used_ver = find([MU_DATA(:).ver_init_truexvLL] > 0);
mu_all_LL_hor = nan(7,length(best_oris));
mu_all_LL_hor(:,used_hori) = [MU_DATA(used_hori).hor_all_LL_imp; MU_DATA(used_hori).hor_fin_LL_imp];
mu_all_LL_hor(:,used_hori) = bsxfun(@rdivide,mu_all_LL_hor(:,used_hori),[MU_DATA(used_hori).hor_init_LL_imp]);

mu_all_LL_ver = nan(7,length(best_oris));
mu_all_LL_ver(:,used_ver) = [MU_DATA(used_ver).ver_all_LL_imp; MU_DATA(used_ver).ver_fin_LL_imp];
mu_all_LL_ver(:,used_ver) = bsxfun(@rdivide,mu_all_LL_ver(:,used_ver),[MU_DATA(used_ver).ver_init_LL_imp]);


mu_RFmean_before = nan(length(best_oris),2);
mu_RFmean_before(best_oris == 1,:) = [MU_DATA(best_oris == 1).hor_uncor_lin_mean; MU_DATA(best_oris == 1).hor_uncor_sq_mean]';
mu_RFmean_before(best_oris == 2,:) = [MU_DATA(best_oris == 2).ver_uncor_lin_mean; MU_DATA(best_oris == 2).ver_uncor_sq_mean]';
mu_RFmean_after = nan(length(best_oris),2);
mu_RFmean_after(best_oris == 1,:) = [MU_DATA(best_oris == 1).hor_cor_lin_mean; MU_DATA(best_oris == 1).hor_cor_sq_mean]';
mu_RFmean_after(best_oris == 2,:) = [MU_DATA(best_oris == 2).ver_cor_lin_mean; MU_DATA(best_oris == 2).ver_cor_sq_mean]';

mu_RFstd_before = nan(length(best_oris),2);
mu_RFstd_before(best_oris == 1,:) = [MU_DATA(best_oris == 1).hor_uncor_lin_std; MU_DATA(best_oris == 1).hor_uncor_sq_std]';
mu_RFstd_before(best_oris == 2,:) = [MU_DATA(best_oris == 2).ver_uncor_lin_std; MU_DATA(best_oris == 2).ver_uncor_sq_std]';
mu_RFstd_after = nan(length(best_oris),2);
mu_RFstd_after(best_oris == 1,:) = [MU_DATA(best_oris == 1).hor_cor_lin_std; MU_DATA(best_oris == 1).hor_cor_sq_std]';
mu_RFstd_after(best_oris == 2,:) = [MU_DATA(best_oris == 2).ver_cor_lin_std; MU_DATA(best_oris == 2).ver_cor_sq_std]';

mu_tpow_before = nan(length(best_oris),2);
mu_tpow_before(best_oris == 1,:) = [MU_DATA(best_oris == 1).hor_uncor_lin_totpow; MU_DATA(best_oris == 1).hor_uncor_sq_totpow]';
mu_tpow_before(best_oris == 2,:) = [MU_DATA(best_oris == 2).ver_uncor_lin_totpow; MU_DATA(best_oris == 2).ver_uncor_sq_totpow]';
mu_tpow_after = nan(length(best_oris),2);
mu_tpow_after(best_oris == 1,:) = [MU_DATA(best_oris == 1).hor_cor_lin_totpow; MU_DATA(best_oris == 1).hor_cor_sq_totpow]';
mu_tpow_after(best_oris == 2,:) = [MU_DATA(best_oris == 2).ver_cor_lin_totpow; MU_DATA(best_oris == 2).ver_cor_sq_totpow]';

mu_FFx_before = nan(length(best_oris),2);
mu_FFx_before(best_oris == 1,:) = [MU_DATA(best_oris == 1).hor_uncor_lin_FFx; MU_DATA(best_oris == 1).hor_uncor_sq_FFx]';
mu_FFx_before(best_oris == 2,:) = [MU_DATA(best_oris == 2).ver_uncor_lin_FFx; MU_DATA(best_oris == 2).ver_uncor_sq_FFx]';
mu_FFx_after = nan(length(best_oris),2);
mu_FFx_after(best_oris == 1,:) = [MU_DATA(best_oris == 1).hor_cor_lin_FFx; MU_DATA(best_oris == 1).hor_cor_sq_FFx]';
mu_FFx_after(best_oris == 2,:) = [MU_DATA(best_oris == 2).ver_cor_lin_FFx; MU_DATA(best_oris == 2).ver_cor_sq_FFx]';

mu_FFt_before = nan(length(best_oris),2);
mu_FFt_before(best_oris == 1,:) = [MU_DATA(best_oris == 1).hor_uncor_lin_FFt; MU_DATA(best_oris == 1).hor_uncor_sq_FFt]';
mu_FFt_before(best_oris == 2,:) = [MU_DATA(best_oris == 2).ver_uncor_lin_FFt; MU_DATA(best_oris == 2).ver_uncor_sq_FFt]';
mu_FFt_after = nan(length(best_oris),2);
mu_FFt_after(best_oris == 1,:) = [MU_DATA(best_oris == 1).hor_cor_lin_FFt; MU_DATA(best_oris == 1).hor_cor_sq_FFt]';
mu_FFt_after(best_oris == 2,:) = [MU_DATA(best_oris == 2).ver_cor_lin_FFt; MU_DATA(best_oris == 2).ver_cor_sq_FFt]';

mu_used = find(mu_trueLL_before > 0);
mu_used_linfilts = mu_used(mu_tpow_after(mu_used,1) > 1000);
mu_used_sqfilts = mu_used(mu_tpow_after(mu_used,2) > 1000);

%%
figure; hold on
errorbar(1:7,nanmean(su_all_LL_hor,2),nanstd(su_all_LL_hor,[],2)/sqrt(length(su_used_hori)),'o-','markersize',6,'linewidth',2);
errorbar(1:7,nanmean(su_all_LL_ver,2),nanstd(su_all_LL_ver,[],2)/sqrt(length(su_used_ver)),'ro-','markersize',6,'linewidth',2);

errorbar(1:7,nanmean(mu_all_LL_hor,2),nanstd(mu_all_LL_hor,[],2)/sqrt(length(used_hori)),'ko-','linewidth',2);
errorbar(1:7,nanmean(mu_all_LL_ver,2),nanstd(mu_all_LL_ver,[],2)/sqrt(length(used_ver)),'o-','color',[0.2 0.8 0.2],'linewidth',2);
legend('SU Vertical','SU Horizontal','MU Vertical','MU Horizontal');

 ylim([1 2.3])
xlim([0.5 7.5])
xlabel('Iteration number','fontsize',10);
ylabel('LL improvement','fontsize',10);
box off
set(gca,'fontsize',8,'fontname','arial');
fname = 'LL_imp_v_iter';
fillPage(gcf,'papersize',[4 4]);

%% LL SCATTER
mS = 5;
figure
% plot(mu_LL_before(mu_used),mu_LL_after(mu_used),'k.','markersize',mS);
hold on
plot(su_LL_before(su_used),su_LL_after(su_used),'k.','markersize',2*mS);
% b = robustfit(su_LL_before(su_used)',su_LL_after(su_used)',[],[],'off');
line([0 1.25],[0 1.25],'color','b');
% xx = linspace(0,1.25,100);
% plot(xx,b*xx,'k');
box off
xlabel('LL before (bits/spk)','fontsize',12);
ylabel('LL after (bits/spk)','fontsize',12);
set(gca,'fontsize',10,'fontname','arial');
xlim([0 1.2]); ylim([0 1.2]);
fillPage(gcf,'papersize',[4 4]);
fname = [fig_dir 'LL_scatter'];
if print_out
print(fname,'-dpdf');
close
end
%% RF WIDTH SCATTER
close all

mS = 4;
%lin filts
figure
% plot(mu_RFstd_before(mu_used_linfilts,1)*dx,mu_RFstd_after(mu_used_linfilts,1)*dx,'k.','markersize',mS);
hold on
plot(su_RFstd_before(su_used_linfilts,1)*dx,su_RFstd_after(su_used_linfilts,1)*dx,'r.','markersize',2*mS);
line([0 0.24],[0 0.24],'color','b');
% b = robustfit(su_RFstd_before(su_used_linfilts,1)*dx,su_RFstd_after(su_used_linfilts,1)*dx,[],[],'off');
% xx = linspace(0,0.24,100);
% plot(xx,b*xx,'k');
box off
xlabel('RF width before (deg)','fontsize',12);
ylabel('RF width after (deg)','fontsize',12);
set(gca,'fontsize',10,'fontname','arial');
xlim([0 0.24]); ylim([0 0.24]);
fillPage(gcf,'papersize',[4 4]);
% fname = [fig_dir 'width_scatter_lin'];
% if print_out
% print(fname,'-dpdf');
% close
% end

% figure
% plot(mu_RFstd_before(mu_used_sqfilts,2)*dx,mu_RFstd_after(mu_used_sqfilts,2)*dx,'k.','markersize',mS);
hold on
plot(su_RFstd_before(su_used_sqfilts,2)*dx,su_RFstd_after(su_used_sqfilts,2)*dx,'b.','markersize',2*mS);
% line([0 0.15],[0 0.15],'color','b');
% b = robustfit(su_RFstd_before(su_used_sqfilts,2)*dx,su_RFstd_after(su_used_sqfilts,2)*dx,[],[],'off');
% xx = linspace(0,0.15,100);
% plot(xx,b*xx,'k');
box off
xlabel('RF width before (deg)','fontsize',12);
ylabel('RF width after (deg)','fontsize',12);
set(gca,'fontsize',10,'fontname','arial');
% xlim([0 0.15]); ylim([0 0.15]);
fillPage(gcf,'papersize',[4 4]);
fname = [fig_dir 'width_scatter_both'];
if print_out
print(fname,'-dpdf');
close
end

%% RF FFx SCATTER
close all
FFx_jitter = 0.05;

mS = 4;
%lin filts
figure
% plot(mu_FFx_before(mu_used_linfilts,1) + randn(length(mu_used_linfilts),1)*FFx_jitter,...
%     mu_FFx_after(mu_used_linfilts,1) + randn(length(mu_used_linfilts),1)*FFx_jitter,'k.','markersize',mS);
hold on
plot(su_FFx_before(su_used_linfilts,1) + randn(length(su_used_linfilts),1)*FFx_jitter,...
    su_FFx_after(su_used_linfilts,1) + randn(length(su_used_linfilts),1)*FFx_jitter,'r.','markersize',2*mS);
% line([0 5],[0 5],'color','b');
% b = robustfit(su_FFx_before(su_used_linfilts,1),su_FFx_after(su_used_linfilts,1),[],[],'off');
% xx = linspace(0,5,100);
% plot(xx,b*xx,'k');
box off
xlabel('SF before (cyc/deg)','fontsize',12);
ylabel('SF after (cyc/deg)','fontsize',12);
set(gca,'fontsize',10,'fontname','arial');
% xlim([0 5]); ylim([0 5]);
fillPage(gcf,'papersize',[4 4]);
fname = [fig_dir 'FFx_scatter_lin'];
% if print_out
% print(fname,'-dpdf');
% close
% end

% figure
% plot(mu_FFx_before(mu_used_linfilts,2) + randn(length(mu_used_linfilts),1)*FFx_jitter,...
%     mu_FFx_after(mu_used_linfilts,2) + randn(length(mu_used_linfilts),1)*FFx_jitter,'k.','markersize',mS);
hold on
plot(su_FFx_before(su_used_linfilts,2) + randn(length(su_used_linfilts),1)*FFx_jitter,...
    su_FFx_after(su_used_linfilts,2) + randn(length(su_used_linfilts),1)*FFx_jitter,'b.','markersize',2*mS);
line([0 8],[0 8],'color','b');
box off
xlabel('SF before (cyc/deg)','fontsize',12);
ylabel('SF after (cyc/deg)','fontsize',12);
set(gca,'fontsize',10,'fontname','arial');
xlim([0 8]); ylim([0 8]);
fillPage(gcf,'papersize',[4 4]);
fname = [fig_dir 'FFx_scatter_both'];
if print_out
print(fname,'-dpdf');
close
end
%% TOT pow scatter
close all
mS = 4;
%lin filts
figure
% plot(mu_tpow_before(mu_used_linfilts,1)/1e4,mu_tpow_after(mu_used_linfilts,1)/1e4,'k.','markersize',mS);
hold on
plot(su_tpow_before(su_used_linfilts,1)/1e4,su_tpow_after(su_used_linfilts,1)/1e4,'r.','markersize',2*mS);
line([0 2],[0 2],'color','b');
box off
xlabel('RF Power before (AU)','fontsize',12);
ylabel('RF Power after (AU)','fontsize',12);
set(gca,'fontsize',10,'fontname','arial');
xlim([0 2]); ylim([0 2]);
fillPage(gcf,'papersize',[4 4]);
% fname = [fig_dir 'tpow_scatter_lin'];
% if print_out
% print(fname,'-dpdf');
% close
% end

% figure
% plot(mu_tpow_before(mu_used_linfilts,2)/1e4,mu_tpow_after(mu_used_linfilts,2)/1e4,'k.','markersize',mS);
hold on
plot(su_tpow_before(su_used_linfilts,2)/1e4,su_tpow_after(su_used_linfilts,2)/1e4,'b.','markersize',2*mS);
% line([0 2],[0 2],'color','b');
box off
xlabel('RF Power before (AU)','fontsize',12);
ylabel('RF Power after (AU)','fontsize',12);
set(gca,'fontsize',10,'fontname','arial');
% xlim([0 1.2]); ylim([0 1.2]);
fillPage(gcf,'papersize',[4 4]);
fname = [fig_dir 'tpow_scatter_both'];
if print_out
print(fname,'-dpdf');
close
end
%% RF WIDTH DISTRIBUTION 
xx = linspace(0,0.22,30);
mu_lin_nn = histc(mu_RFstd_before(mu_used_linfilts,1)*dx,xx)/length(mu_used_linfilts);
su_lin_nn = histc(su_RFstd_before(su_used_linfilts,1)*dx,xx)/length(su_used_linfilts);
mu_sq_nn = histc(mu_RFstd_before(mu_used_sqfilts,2)*dx,xx)/length(mu_used_sqfilts);
su_sq_nn = histc(su_RFstd_before(su_used_sqfilts,2)*dx,xx)/length(su_used_sqfilts);

figure; 
subplot(2,1,1);hold on
stairs(xx,mu_lin_nn,'k','linewidth',2);
stairs(xx,su_lin_nn,'r','linewidth',2);
box off
xlabel('RF width (deg)','fontsize',12);
ylabel('Relative frequency','fontsize',12);
set(gca,'fontsize',10,'fontname','arial');
xlim([0 0.22]);
title('Linear filters','fontsize',12);
legend('MUs','SUs');

subplot(2,1,2);hold on
stairs(xx,mu_sq_nn,'k','linewidth',2);
stairs(xx,su_sq_nn,'r','linewidth',2);
box off
xlabel('RF width (deg)','fontsize',12);
ylabel('Relative frequency','fontsize',12);
set(gca,'fontsize',10,'fontname','arial');
xlim([0 0.22]);
title('Squared filters','fontsize',12);
legend('MUs','SUs');

fillPage(gcf,'papersize',[4 6]);
fname = [fig_dir 'RFwidth_dists'];
if print_out
print(fname,'-dpdf');
close
end
%% SF DISTRIBUTION 
xx = linspace(0,8,30);
mu_lin_nn = histc(mu_FFx_before(mu_used_linfilts,1),xx)/length(mu_used_linfilts);
su_lin_nn = histc(su_FFx_before(su_used_linfilts,1),xx)/length(su_used_linfilts);
mu_sq_nn = histc(mu_FFx_before(mu_used_sqfilts,2),xx)/length(mu_used_sqfilts);
su_sq_nn = histc(su_FFx_before(su_used_sqfilts,2),xx)/length(su_used_sqfilts);

figure; 
subplot(2,1,1);hold on
stairs(xx,mu_lin_nn,'k','linewidth',2);
stairs(xx,su_lin_nn,'r','linewidth',2);
box off
xlabel('Spatial Freq (cyc/deg)','fontsize',12);
ylabel('Relative frequency','fontsize',12);
set(gca,'fontsize',10,'fontname','arial');
xlim([0 8]);
title('Linear filters','fontsize',12);
legend('MUs','SUs');

subplot(2,1,2);hold on
stairs(xx,mu_sq_nn,'k','linewidth',2);
stairs(xx,su_sq_nn,'r','linewidth',2);
box off
xlabel('Spatial Freq (cyc/deg)','fontsize',12);
ylabel('Relative frequency','fontsize',12);
set(gca,'fontsize',10,'fontname','arial');
xlim([0 8]);
title('Squared filters','fontsize',12);
legend('MUs','SUs');

fillPage(gcf,'papersize',[4 6]);
fname = [fig_dir 'FFx_dists'];
if print_out
print(fname,'-dpdf');
close
end
%% PLOT RF POSITIONS
load ~/Data/bruce/general_array_data/ArrayConfig.mat
X_pos = ArrayConfig.X;
Y_pos = ArrayConfig.Y;

MU_probenums = [MU_DATA(:).probenum];
SU_probenums = [SU_DATA(:).probenum];
all_probenums = [MU_probenums SU_probenums];
ALL_DATA = [MU_DATA SU_DATA];
for i = 1:96
   curset = find(all_probenums == i);
   cur_vpos = [ALL_DATA(curset).hor_uncor_sq_mean];
   cur_hpos = [ALL_DATA(curset).ver_uncor_sq_mean];
   all_sq_pos(i,:) = [nanmedian(cur_hpos) nanmedian(cur_vpos)];
   cur_vstd = [ALL_DATA(curset).hor_uncor_sq_std];
   cur_hstd = [ALL_DATA(curset).ver_uncor_sq_std];
   all_sq_std(i,:) = [nanmedian(cur_hstd) nanmedian(cur_vstd)];
end
all_sq_pos(:,1) = (all_sq_pos(:,1) - use_nPix_us/2)*dx + 0.34;
all_sq_pos(:,2) = (all_sq_pos(:,2) - use_nPix_us/2)*dx - 0.43;
all_sq_std = all_sq_std*dx;

for i = 1:96
    tempx(Y_pos(i),X_pos(i)) = all_sq_pos(i,1);
    tempy(Y_pos(i),X_pos(i)) = all_sq_pos(i,2);
    tempinds(Y_pos(i),X_pos(i)) = i;
end
weights = ones(10,10);
weights(tempx==0) = 0;
xpos_interp = smoothn(tempx,weights,'robust');
ypos_interp = smoothn(tempy,weights,'robust');

% figure;
% subplot(2,1,1); hold on
% for ii = 1:96
%     ellipse(2*all_sq_std(ii,1),2*all_sq_std(ii,2),0,all_sq_pos(ii,1),all_sq_pos(ii,2),'k','linewidth',0.5);
% end
% xlim([0 0.9]);ylim([-0.9 0]);
% xlabel('Horizontal position (deg)','fontsize',10);
% ylabel('Vertical position (deg)','fontsize',10);
% legend('MU','SU');
% set(gca,'fontsize',8,'fontname','arial');

% subplot(2,1,2); hold on
figure
hold on
plot(xpos_interp,ypos_interp,'ko','markersize',4,'linewidth',1.5)
% plot(all_sq_pos(:,1),all_sq_pos(:,2),'k.','markersize',10)
% plot(xpos_interp(su_inds),ypos_interp(su_inds),'ro','markersize',4);
xlim([0.1 0.75]);ylim([-0.75 -0.1]);
xlabel('Horizontal position (deg)','fontsize',10);
ylabel('Vertical position (deg)','fontsize',10);
set(gca,'fontsize',8,'fontname','arial');
fillPage(gcf,'papersize',[4 8]);
box off
fillPage(gcf,'papersize',[4 4]);
fname = [fig_dir 'RF_centers'];
if print_out
print(fname,'-dpdf');
close
end
