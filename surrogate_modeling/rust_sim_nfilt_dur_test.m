cd ~/James_scripts/surrogate_modeling/
addpath(genpath('~/James_scripts'));
clear all
% close all

dt = 0.01; %in s
max_rate = Inf; %in Hz
[X,Y] = meshgrid(-4:.02:4,-4:.02:4);

NT = 50000; xvNT = 50000;
SDIM = 16; flen = 8;


n_iter = 10;
for n =1:n_iter
    stim = round(2*rand(NT,SDIM))-1;
    xvstim = round(2*rand(xvNT,SDIM))-1;
    % stim = trnd(10,NT,SDIM);
    % xvstim = trnd(10,xvNT,SDIM);
    
    % GENERATE A FILTER
    nfilts = 10;
    LAMBDA = 6; %5 4
    SIGMA = 1.6; %1.5
    x = repmat(1:SDIM,flen,1);
    sigma1 = repmat(SIGMA,flen,SDIM);
    lambda = repmat(LAMBDA,flen,SDIM);
    % amp_vec = ncx2pdf((1:flen-2)-1,4,1);
    % amp_vec = ncx2pdf((1:flen-2)-1,4,0);
    amp_vec = gampdf((1:flen)-1,3,1.3);
%     amp_vec = [0 0 amp_vec];
    amp_vec = fliplr(amp_vec);
    b = repmat(amp_vec',1,SDIM);
    a = repmat(0,flen,SDIM);
    
    desired_spacing = LAMBDA/5.75; %5.5/ 5 / 4.5n
    beg = SDIM/2-desired_spacing*floor(nfilts/2)+0.75;
    xovals = 0:desired_spacing:(desired_spacing*(nfilts-1));
    xovals = xovals - mean(xovals);
    xovals = xovals + SDIM/2+0.5;
    clear filt_mat
    for i = 1:nfilts
        psi1 = repmat(((1:flen)'/6).^3+5,1,SDIM);
        
        x0 = repmat(xovals(i),flen,SDIM);
            cur_psi = mod(psi1 + rand*pi,2*pi);
        temp = b.*exp(-((x-x0).^2./2./sigma1.^2)) .* (cos(2*pi.*(x-x0)./lambda+cur_psi))+a;
        filt(i,:,:) = temp/norm(temp(:));
        filt_mat(i,:) = temp(:)/norm(temp(:));
    end
    filt_mat = filt_mat';
    c_cent = SDIM/2+0.5;
    c_std = 1.6; %1.6
    max_cval = 0.5; %1
    c_offset = 1;
    cvals = max_cval*exp(-(xovals-c_cent).^2/(2*c_std))+c_offset;
    
    % CREATE TIME-EMBEDDED STIMULUS
    stim_emb = makeStimRows(stim,flen);
    xvstim_emb = makeStimRows(xvstim,flen);
    
    % FILTER STIMULUS
    for i = 1:nfilts
        temp = filt(i,:,:);
        %     filt_stim(i,:) = zscore(stim_emb*temp(:));
        %     xvfilt_stim(i,:) = zscore(xvstim_emb*temp(:));
        filt_stim(i,:) = stim_emb*temp(:);
        xvfilt_stim(i,:) = xvstim_emb*temp(:);
    end
    
    
    %% create spike function
    %pass through internal NLs
    beta = 3; %2.5
    theta = 0;
    coefs = ones(1,nfilts).*cvals;
    g = zeros(1,NT);
    xvg = zeros(1,xvNT);
    for i = 1:nfilts
        lfilt_stim(i,:) = 1/beta*log(1+exp(beta*(filt_stim(i,:)-theta)));
        xvlfilt_stim(i,:) = 1/beta*log(1+exp(beta*(xvfilt_stim(i,:)-theta)));
        g = g + coefs(i)*lfilt_stim(i,:);
        xvg = xvg + coefs(i)*xvlfilt_stim(i,:);
    end
    sg = std(g);
    g = g/sg;
    xvg = xvg/sg;
    
    target_rate = 25; %in Hz
    target_pspike = target_rate*dt;
    cur_theta = 1.75; %1.5  2.5
    cur_beta = 4; %4  4
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
    
    
    %%
    sdim = SDIM; %flen = 14;
    stim_params.spatial_dims = 1;
    stim_params.sdim = sdim;
    stim_params.flen = flen;
    
    defmod.lambda_d2XT = 20;
    
    clear gnmr*
    lam_vals = 5;
    defmod.lambda_L1x = lam_vals;
    poss_nmods = [1:10];
    for nm = 1:length(poss_nmods);
        nmods = poss_nmods(nm)
        init_kerns = 0.01*randn(flen*sdim,nmods);
        
        init_kerns = bsxfun(@rdivide,init_kerns,sqrt(sum(init_kerns.^2)));
        init_signs = ones(1,nmods);
        
        clear kern_types
        for j = 1:nmods
            kern_types{j} = 'threshlin';
        end
        gnm = createGNM(init_kerns,init_signs,kern_types,defmod,stim_params);
        tic
        gnm = fitGNM_filters(gnm,stim_emb,spikebins,'none',[],1e-4,1e-6);
        cur_dur_filt(nm,n) = toc;
        
        fo = gnm;
        fo = setGNM_NLBFs(fo,stim_emb);
        fo = adjust_all_reg(fo,'nlmon',1);
        fo = adjust_all_reg(fo,'lnl2',50);
        fo = adjust_all_reg(fo,'nltype','uncon');
        tic
        fo = fitGNM_internal_NLs(fo,stim_emb,spikebins,1,2);
        cur_dur_nl(nm,n) = toc;
        
    end
end

%%
cd ~/James_scripts/surrogate_modeling
save rust_nfilt_dur_test cur_dur* poss_nmods n_iter
%%
% figure
% shadedErrorBar(poss_nmods,mean(cur_dur_filt,2),std(cur_dur_filt,[],2)/sqrt(n_iter));
% hold on
% shadedErrorBar(poss_nmods,mean(cur_dur_nl,2),std(cur_dur_nl,[],2)/sqrt(n_iter),{'color','r'});
figure
shadedErrorBar(poss_nmods,mean(cur_dur_filt,2),std(cur_dur_filt,[],2));
figure
hold on
shadedErrorBar(poss_nmods,mean(cur_dur_nl,2),std(cur_dur_nl,[],2),{'color','r'});

