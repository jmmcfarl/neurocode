function [nll, pnll,lpen] = getLL_lexp_rot1d(kx,glmod,spkbs)
%% computes likelihood of generator potential
%% uses penalty terms for both nl and history parameters

%% pure ll based on model prediction (e.g. likelihood)
kx(kx > 50)    = 50; %saturate input to spiking NL
if strcmp(glmod.spk_nl,'logexp')
    rx = log(1+exp(kx)); %apply spiking NL tp get predicted rate
elseif strcmp(glmod.spk_nl,'exp')
    rx = exp(kx);
end
rs = rx(spkbs); %predicted rate at spike times
rs(rs < 1e-20) = 1e-20; %impose minimum on predicted rate

%% likelihood given by point process model
ll = sum(log(rs)) - sum(rx);

%% regularisaton penalty (e.g. enforcing smth like posterior)
lambdaW = glmod.lambdaW;
sdim = glmod.mods(1).SDIM;
nmods = length(glmod.mods);

nspks = length(spkbs);

%initialize penalty terms
lpen.loc = zeros(nmods,1);
lpen.w = zeros(nmods,1);
for i = 1:length(glmod.mods)
    mod     = glmod.mods(i); %for current module
    nlx = mod.nlx;
    nly = mod.nly;
    
    dy_dx2 = (nly(3:end)-nly(2:end-1))./(nlx(3:end)-nlx(2:end-1));
    dy_dx1 = (nly(2:end-1)-nly(1:end-2))./(nlx(2:end-1)-nlx(1:end-2));
    nlx_mids = (nlx(2:end)+nlx(1:end-1))/2;
    d2y_d2x = (dy_dx2-dy_dx1)./diff(nlx_mids);
    
    lpen.w(i) = lambdaW*abs(mod.w)/nspks;
    
    %kernel localization penalty
    %     lpen.loc(i) = mod.locLambda*kernel_std(glmod.STCbasis,mod.STCcf,sdim)/nspks;
    %     pll = pll - glmod.mods(i).locLambda*kernel_entropy(glmod.STCbasis,mod.STCcf,sdim);

    
end

loc_pens = arrayfun(@(x) x.locLambda,glmod.mods);
if max(loc_pens) > 0
    if strcmp(glmod.image_type,'2d')
        error('loc penalty not working for 2d yet')
    end
    Nmods = length(glmod.mods);
    kern_l = length(glmod.mods(1).k);
    sdim = glmod.mods(1).SDIM;
    kern_t = kern_l/sdim;
    loc_penalty = zeros(Nmods,1);
    [smooth_x,smooth_y] = meshgrid(-3:3,-3:3);
    smooth_sigma = 1;
    smooth_kern = exp(-(smooth_x.^2+smooth_y.^2)/(2*smooth_sigma^2));
    [X,T] = meshgrid(1:sdim,1:kern_t);
    for n = 1:Nmods
        cur_kmat = reshape(glmod.STCbasis*glmod.mods(n).STCcf,kern_t,sdim);
        cur_smooth_kmat = conv2(cur_kmat.^2,smooth_kern,'same');
        [~,peakloc] = max(cur_smooth_kmat(:));
        [peakt,peakx] = ind2sub(size(cur_smooth_kmat),peakloc);
        distmat = exp((X-peakx).^2/(2*glmod.mods(n).locSigmaX^2) + (T-peakt).^2/(2*glmod.mods(n).locSigmaT^2));
        l2_pen{n} = glmod.mods(n).locLambda*distmat(:);
        lpen.loc(n) = glmod.mods(n).locLambda*sum(cur_kmat(:).^2.*l2_pen{n})/nspks;
    end   
end

%% penalty terms for spiketrain moduls?
nll = -ll/nspks;
pnll = nll + sum(lpen.w) + sum(lpen.loc);

end