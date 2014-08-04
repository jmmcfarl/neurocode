cd ~/James_scripts/surrogate_modeling/
addpath(genpath('~/James_scripts'));
clear all
close all

dt = 0.01; %in s
max_rate = Inf; %in Hz
[X,Y] = meshgrid(-4:.02:4,-4:.02:4);

n_ITER = 6;
for nn = 1:n_ITER
    NT = 200000; SDIM = 20; flen = 14;
    stim = round(2*rand(NT,SDIM))-1;
    xvNT = 30000;
    xvstim = round(2*rand(xvNT,SDIM))-1;
    
    % GENERATE A FILTER
    nfilts = 4;
    LAMBDA = 5; %4
    SIGMA = 1.5; %1.5
    x = repmat(1:SDIM,flen,1);
    sigma1 = repmat(SIGMA,flen,SDIM);
    lambda = repmat(LAMBDA,flen,SDIM);
    % amp_vec = ncx2pdf((1:flen-2)-1,4,1);
    % amp_vec = ncx2pdf((1:flen-2)-1,4,0);
    amp_vec = gampdf((1:flen-2)-1,3,1.3);
    amp_vec = [0 0 amp_vec];
    amp_vec = fliplr(amp_vec);
    b = repmat(amp_vec',1,SDIM);
    a = repmat(0,flen,SDIM);
    % psi1 = repmat(linspace(0,5*pi,flen)',1,SDIM);
    % psi1(1:5,:) = repmat(psi1(6,:),5,1);
    % psi1 = repmat(((1:flen)'/7).^3+4,1,SDIM);
    psi1 = repmat(((1:flen)'/6).^3+5,1,SDIM);
    
    desired_spacing = 0.69*LAMBDA;
    % desired_spacing = 3.25;
    xovals = 0:desired_spacing:(desired_spacing*(nfilts-1));
    xovals = xovals - mean(xovals);
    xovals = xovals + SDIM/2+0.5;
    clear filt filt_mat
    for i = 1:nfilts
        x0 = repmat(xovals(i),flen,SDIM);
        cur_psi = mod(psi1 + rand*2*pi,2*pi);
        temp = b.*exp(-((x-x0).^2./2./sigma1.^2)) .* (cos(2*pi.*(x-x0)./lambda+psi1))+a;
        filt(i,:,:) = temp/norm(temp(:));
        filt_mat(i,:) = temp(:)/norm(temp(:));
    end
    filt_mat = filt_mat';
    c_cent = 8.5;
    c_std = 5;
    max_cval = 1;
    c_offset = 1;
    cvals = max_cval*exp(-(xovals-c_cent).^2/(2*c_std))+c_offset;
    % cvals = [1 1.25 1.25 1];
    cvals = [1 1.25 1.25 1];
    
    % CREATE TIME-EMBEDDED STIMULUS
    stim_emb = makeStimRows(stim,flen);
    xvstim_emb = makeStimRows(xvstim,flen);
    
    % FILTER STIMULUS
    for i = 1:nfilts
        temp = filt(i,:,:);
        filt_stim(i,:) = zscore(stim_emb*temp(:));
        xvfilt_stim(i,:) = zscore(xvstim_emb*temp(:));
    end
    
    %% create spike function
    %pass through internal NLs
    beta = 4; theta = 0;
    coefs = ones(1,nfilts).*cvals;
    g = zeros(1,NT);
    xvg = zeros(1,xvNT);
    for i = 1:nfilts
        lfilt_stim(i,:) = filt_stim(i,:).^2;
        xvlfilt_stim(i,:) = xvfilt_stim(i,:).^2;
        g = g + coefs(i)*lfilt_stim(i,:);
        xvg = xvg + coefs(i)*xvlfilt_stim(i,:);
    end
    g = g/std(g);
    xvg = xvg/std(xvg);
    
    target_rate = 50; %in Hz
    target_pspike = target_rate*dt;
    cur_theta = 1.5; % 1.5 2.5
    cur_beta = 5; % 5 4
    p_spike = 1/cur_beta*log(1+exp(cur_beta*(g-cur_theta)));
    % p_spike = (g-cur_theta).^2; p_spike(p_spike<0) = 0;
    % p_spike = exp((g-1.5)/1);
    xvp_spike = 1/cur_beta*log(1+exp(cur_beta*(xvg-cur_theta)));
    scale_f = target_rate/(mean(p_spike)/dt);
    p_spike = p_spike*scale_f;
    % p_spike(p_spike/dt > 500) = 500*dt;
    xvp_spike = xvp_spike*scale_f;
    % xvp_spike(xvp_spike/dt > 500) = 500*dt;
    %%
    close all
    spikes = poissrnd(p_spike);
    rbins = (find(spikes>0.5));
    nsp = spikes(rbins);
    spk_vals = unique(spikes); spk_vals(spk_vals==0) = [];
    spikebins = [];
    for i = 1:length(spk_vals)
        cur_set = find(spikes == spk_vals(i));
        spikebins = [spikebins; repmat(cur_set(:),spk_vals(i),1)];
    end
    % spikebins = unique(spikebins); %max 1 spike per bin
    fprintf('Nspks: %d\n',length(spikebins));
    
    xvspikes = poissrnd(xvp_spike);
    xvrbins = (find(xvspikes>0.5));
    xvnsp = spikes(xvrbins);
    xvspk_vals = unique(xvspikes); xvspk_vals(xvspk_vals==0) = [];
    xvspikebins = [];
    for i = 1:length(xvspk_vals)
        cur_set = find(xvspikes == xvspk_vals(i));
        xvspikebins = [xvspikebins; repmat(cur_set(:),xvspk_vals(i),1)];
    end
    
    Robs_xv = convert_to_binned_spks(xvspikebins,xvNT);
    avg_rate = length(spikebins)/NT;
    prate = ones(xvNT,1)*avg_rate;
    null_xvLL(nn) = -sum(Robs_xv.*log(prate)-prate)/length(xvspikebins);
    
    %% QUAD MOD
    sdim = SDIM;
    stim_params.spatial_dims = 1;
    stim_params.sdim = sdim;
    stim_params.flen = flen;
    
    defmod.lambda_L1x = 0;
    defmod.lambda_d2XT = 0;
    
    poss_nmods = [1:6];
    for j = 1:length(poss_nmods)
        nmods = poss_nmods(j);
        init_kerns = randn(flen*sdim,nmods);
        init_kerns = bsxfun(@rdivide,init_kerns,sqrt(sum(init_kerns.^2)));
        init_signs = ones(1,nmods);
        clear kern_types
        for i = 1:nmods
            kern_types{i} = 'quad';
        end
        quad_mod(nn,j) = createGNM(init_kerns,init_signs,kern_types,defmod,stim_params);
        for k = 1:3
            quad_mod(nn,j) = fitGNM_filters(quad_mod(nn,j),stim_emb,spikebins,'none',[],1e-4,1e-6);
            [~, ~, ~, ~, g] = getLL_GNM(quad_mod(nn,j),stim_emb,spikebins,'none');
            quad_mod(nn,j) = fitGNM_spkNL(quad_mod(nn,j),g,spikebins,0);
        end
        [~,coms,peak_locs] = get_filter_coms_1d(quad_mod(nn,j));
        [~,ord] = sort(coms);
        quad_mod(nn,j).mods = quad_mod(nn,j).mods(ord);
        
        quad_xvLL(nn,j) = getLL_GNM(quad_mod(nn,j),xvstim_emb,xvspikebins,'none')
        [~, ~, ~, prate] = getLL_GNM(quad_mod(nn,j),xvstim_emb,xvspikebins,'none');
        quad_dev(nn,j) = compute_deviance(xvspikes,prate);
    end
    %%
    quad_basis = get_k_mat(quad_mod(nn,end));
    quad_out = stim_emb*quad_basis;
    quad_stim_params.spatial_dims = 1;
    quad_stim_params.sdim = size(quad_basis,2);
    quad_stim_params.flen = 1;
    n_quad_iter = 10; 
    
    sdim = SDIM; flen = 14;
    stim_params.spatial_dims = 1;
    stim_params.sdim = sdim;
    stim_params.flen = flen;
    
    defmod.lambda_L1x = 100;
    defmod.lambda_d2XT = 0;
    
    poss_nmods = [1:10];
    for j = 1:length(poss_nmods)
        nmods = poss_nmods(j);
%         init_kerns = randn(flen*sdim,nmods);
%         init_kerns = bsxfun(@rdivide,init_kerns,sqrt(sum(init_kerns.^2)));
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
        
        gnm = createGNM(init_kerns,init_signs,kern_types,defmod,stim_params);
        gnm.spk_alpha = init_mod.spk_alpha;
        gnm.spk_beta = init_mod.spk_beta;
        gnm.spk_theta = init_mod.spk_theta;
        gnm = fitGNM_filters(gnm,stim_emb,spikebins,'none',[],1e-4,1e-6);
        [~, ~, ~, ~, g] = getLL_GNM(gnm,stim_emb,spikebins,'none');
        gnm = fitGNM_spkNL(gnm,g,spikebins,0);
        %         gnm = fitGNM_filters(gnm,stim_emb,spikebins,'none',[],1e-4,1e-6);
        gnm_xvLL(nn,j) = getLL_GNM(gnm,xvstim_emb,xvspikebins,'none')
        [~, ~, ~, prate] = getLL_GNM(gnm,xvstim_emb,xvspikebins,'none');
        gnm_dev(nn,j) = compute_deviance(xvspikes,prate);
        [~,coms,peak_locs] = get_filter_coms_1d(gnm);
        [~,ord] = sort(coms);
        gnm.mods = gnm.mods(ord);
        gnmu(nn,j) = gnm;
        
        %for fitting internal NLs
        gnmr(nn,j,1) = adjust_all_reg(gnm,'lnl2',200);
        gnmr(nn,j,1) = setGNM_NLBFs(gnmr(nn,j,1),stim_emb);
        gnmr(nn,j,1) = adjust_all_reg(gnmr(nn,j,1),'nltype','uncon');
        gnmr(nn,j,1) = adjust_all_reg(gnmr(nn,j,1),'nlmon',1);
        gnmr_xvLL(nn,j,1) = getLL_GNM(gnmr(nn,j,1),xvstim_emb,xvspikebins,'none');
        %     [mod_seq,LL_seq,kscale_seq] = iterateGNM_filts_NLs(gnmr,stim_emb,spikebins);
        for i = 2:4
            gnmr(nn,j,i) = fitGNM_internal_NLs(gnmr(nn,j,i-1),stim_emb,spikebins,1,2);
            gnmr(nn,j,i) = fitGNM_filters(gnmr(nn,j,i),stim_emb,spikebins,'none',[],1e-4,1e-6);
            [~, ~, ~, ~, g] = getLL_GNM(gnmr(nn,j,i),stim_emb,spikebins,'none');
            gnmr(nn,j,i) = fitGNM_spkNL(gnmr(nn,j,i),g,spikebins,0);
            gnmr_xvLL(nn,j,i) = getLL_GNM(gnmr(nn,j,i),xvstim_emb,xvspikebins,'none');
        end
    end
end

%%


%%
cd ~/James_scripts/surrogate_modeling/sim_figs/
% save quad_sim_data2 gnm* quad* stc* filt* sdim flen
save quad_sim_data_nfilts_full gnm* quad* filt*

%%
gnm_use_xvLL = squeeze(gnmr_xvLL(:,:,4));
gnm_use_xvLL = -bsxfun(@minus,gnm_use_xvLL,null_xvLL');
quad_use_xvLL = -bsxfun(@minus,quad_xvLL,null_xvLL');

figure
errorbar(1:10,mean(gnm_use_xvLL),std(gnm_use_xvLL)/sqrt(n_ITER));
hold on
errorbar(1:6,mean(quad_use_xvLL),std(quad_use_xvLL)/sqrt(n_ITER),'r')
