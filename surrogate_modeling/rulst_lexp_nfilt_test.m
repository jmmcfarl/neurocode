cd ~/James_scripts/surrogate_modeling/
addpath(genpath('~/James_scripts'));
clear all
% close all

dt = 0.01; %in s
max_rate = Inf; %in Hz
[X,Y] = meshgrid(-4:.02:4,-4:.02:4);

n_rep = 10;
poss_filts = 1:8;

%%
for rp = 1:n_rep
    
    NT = 200000; SDIM = 24; flen = 14;
    stim = round(2*rand(NT,SDIM))-1;
    xvNT = 100000;
    xvstim = round(2*rand(xvNT,SDIM))-1;
    
    % GENERATE GABOR FILTERSFILTER
    nfilts = 6;
    LAMBDA = 5; %5 4
    SIGMA = 1.5; %1.5
    x = repmat(1:SDIM,flen,1);
    sigma1 = repmat(SIGMA,flen,SDIM);
    lambda = repmat(LAMBDA,flen,SDIM);
    amp_vec = gampdf((1:flen-2)-1,3,1.3);
    amp_vec = [0 0 amp_vec];
    amp_vec = fliplr(amp_vec);
    b = repmat(amp_vec',1,SDIM);
    a = repmat(0,flen,SDIM);
    
    clear filt_mat
    desired_spacing = LAMBDA/5.75; %spacing between filter centers
    beg = SDIM/2-desired_spacing*floor(nfilts/2)+0.75;
    xovals = 0:desired_spacing:(desired_spacing*(nfilts-1));
    xovals = xovals - mean(xovals);
    xovals = xovals + SDIM/2+0.5;
    for i = 1:nfilts
        psi1 = repmat(((1:flen)'/6).^3+5,1,SDIM);
        
        x0 = repmat(xovals(i),flen,SDIM);
        temp = b.*exp(-((x-x0).^2./2./sigma1.^2)) .* (cos(2*pi.*(x-x0)./lambda+psi1))+a;
        filt(i,:,:) = temp/norm(temp(:));
        filt_mat(i,:) = temp(:)/norm(temp(:));
    end
    filt_mat = filt_mat';
    
    %define relative weight distribution
    c_cent = SDIM/2+0.5;
    c_std = 1.6; %1.6
    max_cval = 0.5; %0.5
    c_offset = 1;
    cvals = max_cval*exp(-(xovals-c_cent).^2/(2*c_std))+c_offset;
    
    % CREATE TIME-EMBEDDED STIMULUS
    stim_emb = makeStimRows(stim,flen);
    xvstim_emb = makeStimRows(xvstim,flen);
    
    % FILTER STIMULUS
    clear filt_stim xvfilt_stim
    for i = 1:nfilts
        temp = filt(i,:,:);
        filt_stim(i,:) = zscore(stim_emb*temp(:));
        xvfilt_stim(i,:) = zscore(xvstim_emb*temp(:));
    end
    
    
    %% create spike function
    %pass through internal NLs
    beta = 3; %3
    theta = 0;
    coefs = ones(1,nfilts).*cvals;
    g = zeros(1,NT);
    xvg = zeros(1,xvNT);
    clear lfilt_stim xvlfilt_stim
    for i = 1:nfilts
        lfilt_stim(i,:) = 1/beta*log(1+exp(beta*(filt_stim(i,:)-theta)));
        xvlfilt_stim(i,:) = 1/beta*log(1+exp(beta*(xvfilt_stim(i,:)-theta)));
        g = g + coefs(i)*lfilt_stim(i,:);
        xvg = xvg + coefs(i)*xvlfilt_stim(i,:);
    end
    sg = std(g);
    g = g/sg;
    xvg = xvg/sg;
    
    target_rate = 50; %in Hz (50)
    target_pspike = target_rate*dt;
    cur_theta = 1.75; %1.75  2.5
    cur_beta = 4; %4  4
    p_spike = 1/cur_beta*log(1+exp(cur_beta*(g-cur_theta)));
    xvp_spike = 1/cur_beta*log(1+exp(cur_beta*(xvg-cur_theta)));
    scale_f = target_rate/(mean(p_spike)/dt);
    p_spike = p_spike*scale_f;
    xvp_spike = xvp_spike*scale_f;
    
    %% GENERATE SPIKE DATA
    close all
    spikes = poissrnd(p_spike);
    spikebins = convert_to_spikebins(spikes);
    fprintf('Nspks: %d\n',length(spikebins));
    
    xvspikes = poissrnd(xvp_spike);
    xvspikebins = convert_to_spikebins(xvspikes);
    
    %% FIT QUAD MOD
    sdim = SDIM; flen = 14;
    stim_params.spatial_dims = 1;
    stim_params.sdim = sdim;
    stim_params.flen = flen;
    defmod.lambda_L1x = 50;
    defmod.lambda_d2XT = 0;
    
    for pp = 1:length(poss_filts)
        %     for pp = 8
        nmods = poss_filts(pp);
        
        init_kerns = randn(flen*sdim,nmods);
        init_kerns = bsxfun(@rdivide,init_kerns,sqrt(sum(init_kerns.^2)));
        init_signs = ones(1,nmods);
        clear kern_types
        kern_types{1} = 'lin';
        for j = 2:nmods
            kern_types{j} = 'quad';
        end
        quad_mod(rp,pp) = createGNM(init_kerns,init_signs,kern_types,defmod,stim_params);
        quad_mod(rp,pp) = fitGNM_filters(quad_mod(rp,pp),stim_emb,spikebins,'none',[],1e-4,1e-6);
        [~, ~, ~, ~, g,qint_g] = getLL_GNM(quad_mod(rp,pp),stim_emb,spikebins,'none');
        quad_mod(rp,pp) = fitGNM_spkNL(quad_mod(rp,pp),g,spikebins,0);
        quad_xvLL(rp,pp) = getLL_GNM(quad_mod(rp,pp),xvstim_emb,xvspikebins,'none');
        
        quad_modr(rp,pp,1) = quad_mod(rp,pp);
        quadr_xvLL(rp,pp,1) = quad_xvLL(rp,pp);
        n_iter = 3;
        for i = 2:n_iter+1
            quad_modr(rp,pp,i) = fitGNM_filters(quad_modr(rp,pp,i-1),stim_emb,spikebins,'none',[],1e-4,1e-6);
            [~, ~, ~, ~, g] = getLL_GNM(quad_modr(rp,pp,i),stim_emb,spikebins,'none');
            quad_modr(rp,pp,i) = fitGNM_spkNL(quad_modr(rp,pp,i),g,spikebins,0);
            quadr_xvLL(rp,pp,i) = getLL_GNM(quad_modr(rp,pp,i),xvstim_emb,xvspikebins,'none');
        end
    end
    
    %% FIT GNM MODEL IN QUAD SUBSPACE FOR INTITIALIZATION
    quad_basis = get_k_mat(quad_modr(rp,8,end));
    quad_out = stim_emb*quad_basis;
    quad_stim_params.spatial_dims = 1;
    quad_stim_params.sdim = size(quad_basis,2);
    quad_stim_params.flen = 1;
    n_quad_iter = 20;
    
    sdim = SDIM; flen = 14;
    stim_params.spatial_dims = 1;
    stim_params.sdim = sdim;
    stim_params.flen = flen;
    
    for pp = 1:length(poss_filts)
        nmods = poss_filts(pp);
        
        init_signs = ones(1,nmods);
        clear kern_types
        for i = 1:nmods
            kern_types{i} = 'threshlin';
        end
        
        clear temp_gnm temp_LL
        for jj = 1:n_quad_iter
            init_kerns = randn(size(quad_basis,2),nmods);
            temp_gnm(jj) = createGNM(init_kerns,init_signs,kern_types,[],quad_stim_params);
            temp_gnm(jj) = fitGNM_filters(temp_gnm(jj),quad_out,spikebins,'none',[],1e-4,1e-6);
            [~, ~, ~, ~, g] = getLL_GNM(temp_gnm(jj),quad_out,spikebins,'none');
            temp_gnm(jj) = fitGNM_spkNL(temp_gnm(jj),g,spikebins,0);
            temp_gnm(jj) = fitGNM_filters(temp_gnm(jj),quad_out,spikebins,'none',[],1e-4,1e-6);
            temp_LL(jj) = getLL_GNM(temp_gnm(jj),quad_out,spikebins,'none');
        end
        [~,best_loc] = min(temp_LL);
        init_mod = temp_gnm(best_loc);
        init_kerns = quad_basis*get_k_mat(init_mod);
        
        %%
        defmod.lambda_d2XT = 0;
        defmod.lambda_L1x = 50; %50
        gnm = createGNM(init_kerns,init_signs,kern_types,defmod,stim_params);
        gnm.spk_alpha = init_mod.spk_alpha;
        gnm.spk_beta = init_mod.spk_beta;
        gnm.spk_theta = init_mod.spk_theta;
        
        gnm = fitGNM_filters(gnm,stim_emb,spikebins,'none',[],1e-4,1e-6,0);
        [~, ~, ~, ~, g] = getLL_GNM(gnm,stim_emb,spikebins,'none');
        gnmu(rp,pp) = fitGNM_spkNL(gnm,g,spikebins,0);
        gnmu_xvLL(rp,pp) = getLL_GNM(gnmu(rp,pp),xvstim_emb,xvspikebins,'none')
        [~, ~, ~, prate] = getLL_GNM(gnmu(rp,pp),xvstim_emb,xvspikebins,'none');
        gnmu_dev(rp,pp) = compute_deviance(xvspikes,prate)
        
        %sort filters by COM
        [~,coms,peak_locs] = get_filter_coms_1d(gnmu(rp,pp));
        [~,ord] = sort(coms);
        gnmu(rp,pp).mods = gnmu(rp,pp).mods(ord);
        
        %%
        n_iter = 4;
        %Iteratively fit internal NLs
        gnmr(rp,pp,1) = adjust_all_reg(gnmu(rp,pp),'lnl2',500);
        gnmr(rp,pp,1) = setGNM_NLBFs(gnmr(rp,pp,1),stim_emb);
        gnmr(rp,pp,1) = adjust_all_reg(gnmr(rp,pp,1),'nltype','uncon');
        gnmr(rp,pp,1) = adjust_all_reg(gnmr(rp,pp,1),'nlmon',1);
        for j = 2:n_iter+1
            gnmr(rp,pp,j) = fitGNM_internal_NLs(gnmr(rp,pp,j-1),stim_emb,spikebins,1,2);
            gnmr(rp,pp,j) = fitGNM_filters(gnmr(rp,pp,j),stim_emb,spikebins,'none',[],1e-4,1e-6,0);
            [~, ~, ~, ~, g] = getLL_GNM(gnmr(rp,pp,j),stim_emb,spikebins,'none');
            gnmr(rp,pp,j) = fitGNM_spkNL(gnmr(rp,pp,j),g,spikebins,0);
            gnmr_xvLL(rp,pp,j) = getLL_GNM(gnmr(rp,pp,j),xvstim_emb,xvspikebins,'none');
        end
        
        %%
    end
    
    avg_rate = mean(spikes);
    null_prate = ones(size(xvspikes))*avg_rate;
    null_xvLL(rp) = -sum(xvspikes.*log(null_prate) - null_prate)/sum(xvspikes);
    
    %%
    cd ~/James_scripts/surrogate_modeling
    save lexp_sim_nfilt_fin_new3 *xvLL gnmr quad_mod poss_filts n_rep
end

%%
gnm_imp = -bsxfun(@minus,squeeze(gnmr_xvLL(:,:,end)),null_xvLL')/log(2);
quad_imp = -bsxfun(@minus,quadr_xvLL(:,:,end),null_xvLL')/log(2);

figure;
h = errorbar(poss_filts,mean(gnm_imp),std(gnm_imp),'b');
hold on
h = errorbar(poss_filts,mean(quad_imp),std(quad_imp),'r');
plot(poss_filts,mean(gnm_imp),'o')
plot(poss_filts,mean(quad_imp),'ro')
set(gca,'fontname','arial','fontsize',12)
xlabel('Number of subunits','fontsize',16);
ylabel('Log likelihood improvement (bits/spk)','fontsize',16)

%%
gnm_xvLLs = squeeze(gnmr_xvLL(:,:,end));
quad_xvLLs = squeeze(quadr_xvLL(:,:,end));
best_gnm_LL = min(gnm_xvLLs,[],2);
best_quad_LL = min(quad_xvLLs,[],2);
ov_best_LL = min([best_gnm_LL best_quad_LL],[],2);

gnm_perf = -bsxfun(@minus,gnm_xvLLs,ov_best_LL)/log(2);
quad_perf = -bsxfun(@minus,quad_xvLLs,ov_best_LL)/log(2);

figure;
h = errorbar(poss_filts,mean(gnm_perf),std(gnm_perf),'bo-','markersize',12);
hold on
h = errorbar(poss_filts,mean(quad_perf),std(quad_perf),'ro-','markersize',12);
set(gca,'fontname','arial','fontsize',12)
xlabel('Number of subunits','fontsize',16);
ylabel('Log likelihood improvement (bits/spk)','fontsize',16)
