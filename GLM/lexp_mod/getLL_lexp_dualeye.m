function [nll, pnll,lpen] = getLL_lexp_dualeye(kx,glmod,spkbs)
%% computes likelihood of generator potential
%% uses penalty terms for both nl and history parameters

klen = length(glmod.mods(1).k);
flen = klen/glmod.mods(1).SDIM;

%% pure ll based on model prediction (e.g. likelihood)
kx(kx > 100)    = 100; %saturate input to spiking NL
if strcmp(glmod.spk_nl,'logexp')
    rx = log(1+exp(kx)); %apply spiking NL tp get predicted rate
elseif strcmp(glmod.spk_nl,'exp')
    rx = exp(kx);
end
rs = rx(spkbs); %predicted rate at spike times
rs(rs < 1e-10) = 1e-10; %impose minimum on predicted rate

%% likelihood given by point process model
ll = sum(log(rs(flen:end))) - sum(rx(flen:end));

%% regularisaton penalty (e.g. enforcing smth like posterior)
lambdaW = glmod.lambdaW;
fsdim = glmod.mods(1).fsdim;
nmods = length(glmod.mods);

nspks = length(spkbs);

%initialize penalty terms
lpen.l2x = zeros(nmods,1);
lpen.loc = zeros(nmods,1);
lpen.w = zeros(nmods,1);
for i = 1:length(glmod.mods)
    mod     = glmod.mods(i); %for current module
    lpen.w(i) = lambdaW*abs(mod.w)/nspks;
    
    if isfield(mod,'lambda_L2x')
        lpen.l2x(i) = 0.5*mod.lambda_L2x*sum(mod.pix.^2)/nspks;
    end
    
    %kernel localization penalty
    if mod.locLambda > 0 && isfield(mod,'filt_com')
        cur_pkern = mod.pix;
        kern_t = length(cur_pkern)/mod.fsdim;
        if strcmp(glmod.image_type,'2d')
            sdim = sqrt(fsdim);
            [xo,yo] = meshgrid(1:sdim,1:sdim);
            coords = [xo(:) yo(:)];
            dist_mat = repmat(coords,[1 1 kern_t]);
            com_dist_mat = squeeze(sum((dist_mat - repmat(mod.filt_com,[fsdim 1 kern_t])).^2,2))';
            pen_mat = exp(com_dist_mat/2/mod.locSigma^2);
            pen_mat(pen_mat > mod.maxLocPen) = mod.maxLocPen;
            lpen.loc(i) = cur_pkern.^2'*pen_mat(:)*mod.locLambda;
            lpen.loc(i) = lpen.loc(i)/nspks;
        elseif strcmp(glmod.image_type,'1d')
            sdim = fsdim;
            cur_dist_mat = repmat(1:sdim,kern_t,1);
            cur_pen_mat = exp((cur_dist_mat - mod.filt_com).^2/2/mod.locSigma^2);
            cur_pen_mat(cur_pen_mat > mod.maxLocPen) = mod.maxLocPen;
            kern_mat = reshape(cur_pkern,[kern_t sdim]);
            pen_mat = mod.locLambda*kern_mat.^2.*cur_pen_mat;
            lpen.loc(i) = sum(pen_mat(:));
            lpen.loc(i) = lpen.loc(i)/nspks;
        else
            error('Invalid image type')
        end
    else
        lpen.loc(i) = 0;
    end
    
    %kernel smoothness penalty
    if mod.lambda_dX > 0 | mod.lambda_dT > 0
        lpen.ksmooth(i) = 0;
        for eye = 1:2 %loop over eyes
            cur_pkern = mod.pix(glmod.eyeinds==eye);
            kern_t = length(cur_pkern)/mod.fsdim;
            if strcmp(glmod.image_type,'2d')
                sdim = sqrt(mod.fsdim);
                if sdim ~= round(sdim)
                    error('non-square pixel array')
                end
                kern_mat = reshape(cur_pkern,[kern_t,sdim,sdim]);
                x_derivs = (kern_mat(:,2:end,:) - kern_mat(:,1:end-1,:));
                y_derivs = (kern_mat(:,:,2:end) - kern_mat(:,:,1:end-1));
                t_derivs = (kern_mat(2:end,:,:) - kern_mat(1:end-1,:,:));
                lpen.ksmooth(i) = lpen.ksmooth(i) + (mod.lambda_dX*sum(x_derivs(:).^2) + mod.lambda_dX*sum(y_derivs(:).^2) + mod.lambda_dT*sum(t_derivs(:).^2));
%                 lpen.ksmooth(i) = lpen.ksmooth(i)/nspks;
            elseif strcmp(glmod.image_type,'1d')
                sdim = fsdim;
                kern_mat = reshape(cur_pkern,[kern_t sdim]);
                x_derivs = (kern_mat(:,2:end) - kern_mat(:,1:end-1));
                t_derivs = (kern_mat(2:end,:) - kern_mat(1:end-1,:));
                lpen.ksmooth(i) = lpen.ksmooth(i) + (mod.lambda_dX*sum(x_derivs(:).^2) + mod.lambda_dT*sum(t_derivs(:).^2));
%                 lpen.ksmooth(i) = lpen.ksmooth(i)/nspks;
            end
        end
        lpen.ksmooth(i) = lpen.ksmooth(i)/nspks;
    else
        lpen.ksmooth(i) = 0;
    end
    
    %L1 kernel penalty
    if mod.lambda_L1x > 0
        cur_pkern = mod.pix;
        lpen.l1x(i) = sum(abs(cur_pkern).*mod.lambda_L1x/nspks);
    else
        lpen.l1x(i) = 0;
    end
    
end

%% penalty terms for spiketrain moduls?
nll = -ll/nspks;
pnll = nll + sum(lpen.w) + sum(lpen.loc) + sum(lpen.ksmooth) + sum(lpen.l1x) + sum(lpen.l2x);

end