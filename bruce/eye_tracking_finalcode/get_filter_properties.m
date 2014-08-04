function filt_data = get_filter_properties(filtK,stim_params,dx)

zpad_factor = 5;
flen = stim_params.stim_dims(1);
SDIM = stim_params.stim_dims(2);

nfilts = size(filtK,2);
zero_filts = find(max(abs(filtK)) == 0);
filtK = reshape(filtK,[flen SDIM nfilts]);
spatial_profiles = squeeze(std(filtK));

% xdata = ((0.5:SDIM-0.5) - SDIM/2)*dx;

if nfilts > 1
    [a,best_sq_filter] = max(std(spatial_profiles(:,2:end)));
    best_sq_filter = 1 + best_sq_filter;
    worst_sq_filter = setdiff(2:3,best_sq_filter);
    
    filt_data.ford = [1 best_sq_filter worst_sq_filter];
    
    %gaussian fits to spatial profiles
    if max(spatial_profiles(:,1)) > 0
        [fit_params,lfit_z] = fitGaussianCurve((1:SDIM)',spatial_profiles(:,1));
        filt_data.lin_mean = fit_params(1)*dx;
        filt_data.lin_std = fit_params(2)*dx;
    else
        filt_data.lin_mean = nan;
        filt_data.lin_std = nan;
    end
    if max(spatial_profiles(:,best_sq_filter)) > 0
        %     [fit_params,fit_z] = fitGaussianCurve((1:SDIM)',mean(spatial_profiles(:,2:end),2));
        [fit_params,fit_z] = fitGaussianCurve((1:SDIM)',spatial_profiles(:,best_sq_filter));
        filt_data.sq_mean = fit_params(1)*dx;
        filt_data.sq_std = fit_params(2)*dx;
    else
        filt_data.sq_mean = nan;
        filt_data.sq_std = nan;
    end
    if max(spatial_profiles(:,worst_sq_filter)) > 0
        [fit_params,fit_z] = fitGaussianCurve((1:SDIM)',spatial_profiles(:,worst_sq_filter));
        filt_data.bsq_mean = fit_params(1)*dx;
        filt_data.bsq_std = fit_params(2)*dx;
    else
        filt_data.bsq_mean = nan;
        filt_data.bsq_std = nan;
    end
   
    lin_filt = squeeze(filtK(:,:,1));
    sq_filt = squeeze(filtK(:,:,best_sq_filter));
    bsq_filt = squeeze(filtK(:,:,worst_sq_filter));
    [~,best_tslice] = max(std(sq_filt,[],2));
    sq_best_slice = sq_filt(best_tslice,:);
    [~,best_tslice] = max(std(bsq_filt,[],2));
    bsq_best_slice = bsq_filt(best_tslice,:);
    [~,best_tslice] = max(std(lin_filt,[],2));
    lin_best_slice = lin_filt(best_tslice,:);
        
    %FFT calculations
    dt = stim_params.dt;
    Ft = linspace(-1/dt/2,1/dt/2,zpad_factor*flen);
    Fx = linspace(-1/dx/2,1/dx/2,zpad_factor*SDIM);
    sig_Ft = 0.1; sig_Fx = 0.1;
    [FFx,FFt] = meshgrid(Fx,Ft);
    gauss_kern = FFt.^2/(2*sig_Ft^2) + FFx.^2/(2*sig_Fx^2);
    gauss_kern = exp(-gauss_kern);
    gauss_kern = gauss_kern/sum(gauss_kern(:));
    
    max_pow = nan(nfilts,1);
    max_ploc = nan(nfilts,1);
    tot_pow = nan(nfilts,1);
    filts_zpad = zeros(zpad_factor*flen,zpad_factor*SDIM,nfilts);
    filts_zpad((flen+1):2*flen,(SDIM+1):2*SDIM,:) = filtK;
    filt_ffts = zeros(size(filts_zpad));
    for ii = 1:3
        cur_ffts = abs(fftshift(fft2(filts_zpad(:,:,ii))));
        cur_ffts = conv2(squeeze(cur_ffts),gauss_kern,'same');
        max_pow(ii) = max(cur_ffts(:));
        max_ploc(ii) = find(cur_ffts == max_pow(ii),1);
        tot_pow(ii) = sum(cur_ffts(:));
        filt_ffts(:,:,ii) = cur_ffts;
    end
    FFx = abs(FFx(max_ploc));
    FFt = abs(FFt(max_ploc));
    max_pow(zero_filts) = nan; tot_pow(zero_filts) = nan;
    FFx(zero_filts) = nan; FFt(zero_filts) = nan;
    
    filt_data.lin_FFx = FFx(1);
    filt_data.lin_FFt = FFt(1);
    % filt_data.sq_FFx = nanmean(FFx(2:end));
    % filt_data.sq_FFt = nanmean(FFt(2:end));
    filt_data.sq_FFx = FFx(best_sq_filter);
    filt_data.sq_FFt = FFt(best_sq_filter);
    filt_data.bsq_FFx = FFx(worst_sq_filter);
    filt_data.bsq_FFt = FFt(worst_sq_filter);
    
    filt_data.lin_max_pow = max_pow(1);
    filt_data.lin_tot_pow = tot_pow(1);
    % filt_data.sq_max_pow = nanmean(max_pow(2:end));
    % filt_data.sq_tot_pow = nanmean(tot_pow(2:end));
    filt_data.sq_max_pow = max_pow(best_sq_filter);
    filt_data.sq_tot_pow = tot_pow(best_sq_filter);
    filt_data.bsq_max_pow = max_pow(worst_sq_filter);
    filt_data.bsq_tot_pow = tot_pow(worst_sq_filter);
    
    filt_data.lin_fft = squeeze(filt_ffts(:,:,1));
    filt_data.sq_fft = squeeze(filt_ffts(:,:,best_sq_filter));
    filt_data.bsq_fft = squeeze(filt_ffts(:,:,worst_sq_filter));
    
    filt_data.simplicity = filt_data.lin_tot_pow/filt_data.sq_tot_pow;
    
    %% FIT GABOR FUNCTIONS
    % FOR BEST SQ SLICE
    opts = optimset('display','off');
    xdata = (1:SDIM)*dx;
    poss_phases = [0 pi/2 pi -pi/2];
    
    
%     LB = [dx 2*dx 1/2/(SDIM*dx) -2*pi 0 prctile(sq_best_slice,5)];
    LB = [dx dx 0 -2*pi 0 prctile(sq_best_slice,5)];
    UB = [SDIM*dx SDIM*dx/2 1/dx/2 2*pi range(sq_best_slice) prctile(sq_best_slice,95)];

    gparams(1) = filt_data.sq_mean;
    gparams(2) = filt_data.sq_std;
    gparams(3) = filt_data.sq_FFx;
    gparams(5) = prctile(sq_best_slice,95)-prctile(sq_best_slice,5);
    gparams(6) = prctile(sq_best_slice,50);
    
    if all(~isnan(gparams)) & max(abs(sq_best_slice)) > 0
    
    for ii = 1:length(poss_phases)
        gparams(4) = poss_phases(ii);
%         [gest(ii,:),resnorm(ii)] = lsqnonlin(@(x)gabor_1d_error(x,sq_best_slice,xx),gparams,[],[],opts);
        [gest(ii,:),resnorm(ii)] = lsqcurvefit(@gabor_1d,gparams,xdata,sq_best_slice,LB,UB,opts);
    end
    [best_sq_res,best_ind] = min(resnorm);
    sq_gabor_est = gest(best_ind,:);
    Gf = gabor_1d(sq_gabor_est,xdata);
    
    filt_data.sq_gest = sq_gabor_est;
    filt_data.sq_gabor_est = Gf;
    filt_data.sq_slice = sq_best_slice;
    filt_data.sq_est_R2 = 1-best_sq_res/sum(sq_best_slice.^2);
    
    else
    filt_data.sq_gest = nan(1,6);
    filt_data.sq_gabor_est = nan;
    filt_data.sq_slice = nan;
    filt_data.sq_est_R2 = nan;
    end
    
    gparams(1) = filt_data.sq_mean;
    gparams(2) = filt_data.sq_std;
    gparams(3) = filt_data.sq_FFx;
    gparams(5) = prctile(bsq_best_slice,95)-prctile(bsq_best_slice,5);
    gparams(6) = prctile(bsq_best_slice,50);
    
    if all(~isnan(gparams)) & max(abs(bsq_best_slice)) > 0
    
    for ii = 1:length(poss_phases)
        gparams(4) = poss_phases(ii);
        [gest(ii,:),resnorm(ii)] = lsqcurvefit(@gabor_1d,gparams,xdata,sq_best_slice,LB,UB,opts);
    end
    [best_bsq_res,best_ind] = min(resnorm);
    bsq_gabor_est = gest(best_ind,:);
    Gf = gabor_1d(bsq_gabor_est,xdata);
    
    filt_data.bsq_gest = bsq_gabor_est;
    filt_data.bsq_gabor_est = Gf;
    filt_data.bsq_slice = bsq_best_slice;
    filt_data.bsq_est_R2 = 1-best_bsq_res/sum(bsq_best_slice.^2);
    
    else
    filt_data.bsq_gest = nan(1,6);
    filt_data.bsq_gabor_est = nan;
    filt_data.bsq_slice = nan;
    filt_data.bsq_est_R2 = nan;
    end

    
%     LB = [dx 2*dx 1/2/(SDIM*dx) -2*pi 0 prctile(lin_best_slice,5)];
    LB = [dx dx 0 -2*pi 0 prctile(lin_best_slice,5)];
    UB = [SDIM*dx SDIM*dx/2 1/dx/2 2*pi range(lin_best_slice) prctile(lin_best_slice,95)];

    % FOR BEST LIN SLICE
    gparams(1) = filt_data.lin_mean;
    gparams(2) = filt_data.lin_std;
    gparams(3) = filt_data.lin_FFx;
    gparams(5) = prctile(lin_best_slice,95)-prctile(lin_best_slice,5);
    gparams(6) = prctile(lin_best_slice,50);
    
     if all(~isnan(gparams)) & max(abs(lin_best_slice)) > 0
   for ii = 1:length(poss_phases)
        gparams(4) = poss_phases(ii);
%         [gest(ii,:),resnorm(ii)] = lsqnonlin(@(x)gabor_1d_error(x,lin_best_slice,xx),gparams,[],[],opts);
        [gest(ii,:),resnorm(ii)] = lsqcurvefit(@gabor_1d,gparams,xdata,lin_best_slice,LB,UB,opts);
    end
    [best_lin_res,best_ind] = min(resnorm);
    lin_gabor_est = gest(best_ind,:);
    Gf = gabor_1d(lin_gabor_est,xdata);
    
    filt_data.lin_gest = lin_gabor_est;
    filt_data.lin_gabor_est = Gf;
    filt_data.lin_slice = lin_best_slice;
    filt_data.lin_est_R2 = 1-best_lin_res/sum(lin_best_slice.^2);

     else
    filt_data.lin_gest = nan(1,6);
    filt_data.lin_gabor_est = nan;
    filt_data.lin_slice = nan;
    filt_data.lin_est_R2 = nan;
         
     end
    %%
%         clf
%     subplot(2,2,1)
%     imagesc(lin_filt);
%     subplot(2,2,3)
%     hold on
%     plot((1:SDIM),lin_best_slice);
%     hold on
%     plot((1:SDIM),filt_data.lin_gabor_est,'r');
%     axis tight
%     subplot(2,2,2)
%     imagesc(sq_filt);
%     subplot(2,2,4)
%     plot((1:SDIM),sq_best_slice);
%     hold on
%     plot((1:SDIM),filt_data.sq_gabor_est,'r');
%     axis tight
%     fprintf('Lin R2: %.3f\n',filt_data.lin_est_R2);
%     fprintf('Sq R2: %.3f\n',filt_data.sq_est_R2);
%     pause

    
else
    filt_data.lin_mean = nan;
    filt_data.lin_std = nan;
    if max(spatial_profiles) > 0
        [fit_params,fit_z] = fitGaussianCurve((1:SDIM)',spatial_profiles(:));
        filt_data.sq_mean = fit_params(1)*dx;
        filt_data.sq_std = fit_params(2)*dx;
    else
        filt_data.sq_mean = nan;
        filt_data.sq_std = nan;
    end
    
    %FFT calculations
    dt = stim_params.dt;
    Ft = linspace(-1/dt/2,1/dt/2,zpad_factor*flen);
    Fx = linspace(-1/dx/2,1/dx/2,zpad_factor*SDIM);
    sig_Ft = 0.1; sig_Fx = 0.1;
    [FFx,FFt] = meshgrid(Fx,Ft);
    gauss_kern = FFt.^2/(2*sig_Ft^2) + FFx.^2/(2*sig_Fx^2);
    gauss_kern = exp(-gauss_kern);
    gauss_kern = gauss_kern/sum(gauss_kern(:));
    
    filts_zpad = zeros(zpad_factor*flen,zpad_factor*SDIM,nfilts);
    filts_zpad((flen+1):2*flen,(SDIM+1):2*SDIM,:) = filtK;
    filt_ffts = zeros(size(filts_zpad));
    cur_ffts = abs(fftshift(fft2(filts_zpad)));
    cur_ffts = conv2(squeeze(cur_ffts),gauss_kern,'same');
    max_pow = max(cur_ffts(:));
    max_ploc = find(cur_ffts == max_pow,1);
    tot_pow = sum(cur_ffts(:));
    filt_ffts = cur_ffts;
    
    FFx = abs(FFx(max_ploc));
    FFt = abs(FFt(max_ploc));
    max_pow(zero_filts) = nan; tot_pow(zero_filts) = nan;
    FFx(zero_filts) = nan; FFt(zero_filts) = nan;
    
    filt_data.lin_FFx = nan;
    filt_data.lin_FFt = nan;
    filt_data.sq_FFx = FFx;
    filt_data.sq_FFt = FFt;
    
    filt_data.lin_max_pow = nan;
    filt_data.lin_tot_pow = nan;
    filt_data.sq_max_pow = max_pow;
    filt_data.sq_tot_pow = tot_pow;
    
    filt_data.lin_fft = nan;
    filt_data.sq_fft = filt_ffts;
    
    filt_data.simplicity = filt_data.lin_tot_pow/filt_data.sq_tot_pow;
 
    %%

end
end


%%    
% 
% function efun = gabor_1d_error(X,y,k)
%     efun = (gabor_1d(X,k) - y);
% end

function G = gabor_1d(X,k)
    gauss = exp(-(k-X(1)).^2/(2*X(2)^2));
    G = gauss.*sin(2*pi*k*X(3) + X(4))*X(5) + X(6);
end

