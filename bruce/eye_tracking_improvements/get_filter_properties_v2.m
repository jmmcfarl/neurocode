function filt_data = get_filter_properties_v2(filtK,stim_params,dx)
% filt_data = get_filter_properties(filtK,stim_params,dx)
% compute stats on a set of spatiotemporal stimulus filters
% INPUTS:
% filtK, matrix of stim filters [KLEN x Nfilts]
% stim_params: struct of stimulus parameters
% dx: spatial resolution
% OUTPUTS:
% filt_data: struct with filter stats

zpad_factor = 5; %fold-expansion of filters with zero-padding to increase frequency resolution of FFTs
flen = stim_params.stim_dims(1); %n time lags
SDIM = stim_params.stim_dims(2); %n spatial coefs

nfilts = size(filtK,2); %number of filters
zero_filts = find(max(abs(filtK)) == 0); %filters that are all zeros
filtK = reshape(filtK,[flen SDIM nfilts]);
spatial_profiles = squeeze(std(filtK)); %SD across time lags gives the spatial profile

% %make column vecs
% if size(spatial_profiles,2) > size(spatial_profiles,1)
%     spatial_profiles = spatial_profiles';
% end

%FFT calculations
dt = stim_params.dt; %time res
Ft = linspace(-1/dt/2,1/dt/2,zpad_factor*flen); %temporal freq axis
Fx = linspace(-1/dx/2,1/dx/2,zpad_factor*SDIM); %spatial freq axis

%make gaussian smoothing kernel in 2d freq domain
sig_Ft = 0.1; sig_Fx = 0.1; %SD of gaussian smoothing kernel along each freq axis
[FFx,FFt] = meshgrid(Fx,Ft);
gauss_kern = FFt.^2/(2*sig_Ft^2) + FFx.^2/(2*sig_Fx^2);
gauss_kern = exp(-gauss_kern);
gauss_kern = gauss_kern/sum(gauss_kern(:));

uright_quad = find(FFx > 0 & FFt > 0); %indices of the upper-right quadrant
uleft_quad = find(FFx < 0 & FFt > 0); %indices of the upper left quadrant

%% FIT GAUSSIANS
for ii = 1:nfilts
    %gaussian fits to spatial profiles of each individual filter
    if max(spatial_profiles(:,ii)) > 0
        [fit_params,fit_z] = fitGaussianCurve((1:SDIM)',spatial_profiles(:,ii));
        filt_data.gauss_mean(ii) = fit_params(1)*dx;
        filt_data.gauss_std(ii) = fit_params(2)*dx;
    else
        filt_data.gauss_mean(ii) = nan;
        filt_data.gauss_std(ii) = nan;
    end
end

%now fit gaussian on the avg kernel spatial profile
avg_spatial_profile = squeeze(nanmean(spatial_profiles,2));
[fit_params,fit_z] = fitGaussianCurve((1:SDIM)',avg_spatial_profile);
filt_data.avg_gauss_mean = fit_params(1)*dx;
filt_data.avg_gauss_std = fit_params(2)*dx;

%% FOURIER ANALYSIS
[max_pow,max_ploc,tot_pow] = deal(nan(nfilts,1));

%make a zero-padded array with filters
filts_zpad = zeros(zpad_factor*flen,zpad_factor*SDIM,nfilts);
filts_zpad((flen+1):2*flen,(SDIM+1):2*SDIM,:) = filtK;
filt_ffts = zeros(size(filts_zpad));
for ii = 1:nfilts
    cur_ffts = abs(fftshift(fft2(filts_zpad(:,:,ii))));
    cur_ffts = conv2(squeeze(cur_ffts),gauss_kern,'same'); %smooth with gaussian kernel
    max_pow(ii) = max(cur_ffts(:)); %max of smoothed amp-spec
    max_ploc(ii) = find(cur_ffts == max_pow(ii),1); %location of max
    tot_pow(ii) = sum(cur_ffts(:)); %total power
    filt_ffts(:,:,ii) = cur_ffts;
end
max_pow(zero_filts) = nan; tot_pow(zero_filts) = nan; FFx(zero_filts) = nan; FFt(zero_filts) = nan; %for filters that were all zeros
filt_data.FFx = abs(FFx(max_ploc)); %spatial frequency at peak 
filt_data.FFt = abs(FFt(max_ploc)); %temporal frequency at peak
filt_data.max_pow = max_pow;
filt_data.tot_pow = tot_pow;

%estimate a measure of direction selectivity based on the assymetry of
%power in the upper-right and upper-left quadrants
cur_ffts = permute(filt_ffts,[3 1 2]);
tot_uright_pow = sum(cur_ffts(:,uright_quad),2);
tot_uleft_pow = sum(cur_ffts(:,uleft_quad),2);
filt_data.dir_selectivity = (tot_uright_pow - tot_uleft_pow)./(tot_uright_pow + tot_uleft_pow);
filt_data.filt_FFts = filt_ffts;
filt_data.ax_Fx = Fx; filt_data.ax_Ft = Ft;

%% FIT GABOR FUNCTIONS
opts = optimset('display','off');
xdata = (1:SDIM)*dx;
poss_phases = [0 pi/2 pi -pi/2];

for ii = 1:nfilts
    cur_filt = squeeze(filtK(:,:,ii));
    [~,best_tslice] = max(std(cur_filt,[],2)); %get index of best time slice (most power across space)
    best_slice = cur_filt(best_tslice,:); % %best time slice
    
    %[MU SDs SF spatial-phase amp offset]
    LB = [dx dx 0 -2*pi 0 prctile(best_slice,5)]; %lower bounds
    UB = [SDIM*dx SDIM*dx/2 1/dx/2 2*pi range(best_slice) prctile(best_slice,95)]; %upper bounds
    
    %initial estimates
    gparams(1) = filt_data.gauss_mean(ii);
    gparams(2) = filt_data.gauss_std(ii);
    gparams(3) = filt_data.FFx(ii);
    gparams(5) = prctile(best_slice,95)-prctile(best_slice,5);
    gparams(6) = prctile(best_slice,50);
    
    
    if all(~isnan(gparams)) && sum(best_slice ~= 0) > 3 %make sure there are at least three non-zero spatial coefs to work with
        gest = nan(length(poss_phases),6);resnorm = nan(length(poss_phases),1);
        for jj = 1:length(poss_phases) %cycle through a range of possible initial spatial phases
            gparams(4) = poss_phases(jj);
            [gest(jj,:),resnorm(jj)] = lsqcurvefit(@gabor_1d,gparams,xdata,best_slice,LB,UB,opts);
        end
        [best_res,best_ind] = min(resnorm); %find the best fit
        gabor_est = gest(best_ind,:); %best fit gabor params
        Gf = gabor_1d(gabor_est,xdata); %best-fit gabor function
        
        filt_data.gest(ii,:) = gabor_est;
        filt_data.gabor_est(ii,:) = Gf;
        filt_data.best_slice(ii,:) = best_slice;
        filt_data.gabor_R2(ii) = 1-best_res/sum(best_slice.^2);
        
    else
        filt_data.gest(ii,:) = nan(1,6);
        filt_data.gabor_est(ii,:) = nan;
        filt_data.best_slice(ii,:) = nan;
        filt_data.gabor_R2(ii) = nan;
    end
    
end

%%
% for ii = 1:nfilts
%     clf
%     subplot(2,1,1)
%     imagesc(squeeze(filtK(:,:,ii)));
%     subplot(2,1,2)
%     hold on
%     plot((1:SDIM),filt_data.best_slice(ii,:));
%     hold on
%     plot((1:SDIM),filt_data.gabor_est(ii,:),'r');
%     axis tight
%     fprintf('R2: %.3f\n',filt_data.gabor_R2(ii));
%     pause
% end

end


%%
function G = gabor_1d(X,k)
gauss = exp(-(k-X(1)).^2/(2*X(2)^2));
G = gauss.*sin(2*pi*k*X(3) + X(4))*X(5) + X(6);
end

