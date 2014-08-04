function [nll, pnll,lpen,rx] = getLL_lexp(kx,glmod,spkbs)
%% computes likelihood of generator potential
%% uses penalty terms for both nl and history parameters

if strcmp(glmod.basis,'pix')
    klen = length(glmod.mods(1).k);
    flen = klen/glmod.mods(1).SDIM;
elseif strcmp(glmod.basis,'white')
    klen = length(glmod.mods(1).k);
    flen = 1;
end
%% pure ll based on model prediction (e.g. likelihood)
kx(kx > 100)    = 100; %saturate input to spiking NL
if strcmp(glmod.spk_nl,'logexp')
    rx = log(1+exp(kx)); %apply spiking NL tp get predicted rate
elseif strcmp(glmod.spk_nl,'exp')
    rx = exp(kx);
else
    error('Not accepted spiking NL')
end
rs = rx(spkbs); %predicted rate at spike times
rs(rs < 1e-10) = 1e-10; %impose minimum on predicted rate

%% likelihood given by point process model
% flen = 1; %if you want to include the initial bins
ll = sum(log(rs(flen:end))) - sum(rx(flen:end));

%% regularisaton penalty (e.g. enforcing smth like posterior)
lambdaW = glmod.lambdaW;
fsdim = glmod.mods(1).fsdim;
nmods = length(glmod.mods);

nspks = length(spkbs);

%initialize penalty terms
lpen.l2x = zeros(nmods,1);
lpen.ksmooth = zeros(nmods,1);
lpen.lapl = zeros(nmods,1);
lpen.laplXT = zeros(nmods,1);
lpen.loc = zeros(nmods,1);
lpen.w = zeros(nmods,1);
lpen.laplT = zeros(nmods,1);
for i = 1:length(glmod.mods)
    mod = glmod.mods(i); %for current module
    lpen.w(i) = lambdaW*abs(mod.w)/nspks;
    
    if isfield(mod,'lambda_L2x')
        lpen.l2x(i) = 0.5*mod.lambda_L2x*sum((mod.w*mod.pix).^2)/nspks; %scale reintroduced to k
    end
    
    %kernel localization penalty
    if mod.locLambda > 0 && isfield(mod,'filt_com')
        %         cur_pkern = mod.pix;
        cur_pkern = mod.pix*mod.kscale; %reintroduce scale into k
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
        %         cur_pkern = mod.pix;
        cur_pkern = mod.pix*mod.kscale; %reintroduce scale into k
        kern_t = length(cur_pkern)/mod.fsdim;
        if strcmp(glmod.image_type,'2d')
            sdim = sqrt(mod.fsdim);
            if sdim ~= round(sdim)
                error('non-square pixel array')
            end
            kern_mat = reshape(cur_pkern,[kern_t,sdim,sdim]);
            
            %zero pad kernel matrix
            kern_mat = cat(1,zeros(1,size(kern_mat,2),size(kern_mat,3)),kern_mat,zeros(1,size(kern_mat,2),size(kern_mat,3)));
            kern_mat = cat(2,zeros(size(kern_mat,1),1,size(kern_mat,3)),kern_mat,zeros(size(kern_mat,1),1,size(kern_mat,3)));
            kern_mat = cat(3,zeros(size(kern_mat,1),size(kern_mat,2),1),kern_mat,zeros(size(kern_mat,1),size(kern_mat,2),1));
            
            x_derivs = (kern_mat(:,2:end,:) - kern_mat(:,1:end-1,:));
            y_derivs = (kern_mat(:,:,2:end) - kern_mat(:,:,1:end-1));
            t_derivs = (kern_mat(2:end,:,:) - kern_mat(1:end-1,:,:));
            lpen.ksmooth(i) = (mod.lambda_dX*sum(x_derivs(:).^2) + mod.lambda_dX*sum(y_derivs(:).^2) + mod.lambda_dT*sum(t_derivs(:).^2));
            lpen.ksmooth(i) = lpen.ksmooth(i)/nspks;
        elseif strcmp(glmod.image_type,'1d') || strcmp(glmod.image_type,'0d')
            sdim = fsdim;
            kern_mat = reshape(cur_pkern,[kern_t sdim]);
            
            % buffer with zeros
            kern_mat = cat(1,zeros(1,size(kern_mat,2)),kern_mat,zeros(1,size(kern_mat,2)));
            kern_mat = cat(2,zeros(size(kern_mat,1)),kern_mat,zeros(size(kern_mat,1),1));
            
            x_derivs = (kern_mat(:,2:end) - kern_mat(:,1:end-1));
            t_derivs = (kern_mat(2:end,:) - kern_mat(1:end-1,:));
            lpen.ksmooth(i) = (mod.lambda_dX*sum(x_derivs(:).^2) + mod.lambda_dT*sum(t_derivs(:).^2));
            lpen.ksmooth(i) = lpen.ksmooth(i)/nspks;
        end
    else
        lpen.ksmooth(i) = 0;
    end
    if mod.lambda_d2X > 0
        %         cur_pkern = mod.pix;
        cur_pkern = mod.pix*mod.kscale; %reintroduce scale into k
        kern_t = length(cur_pkern)/mod.fsdim;
        if strcmp(glmod.image_type,'2d')
            sdim = sqrt(mod.fsdim);
            if sdim ~= round(sdim)
                error('non-square pixel array')
            end
            kern_mat = reshape(cur_pkern,[kern_t,sdim,sdim]);
            
            %zero pad kernel matrix
            kern_mat = cat(1,zeros(1,size(kern_mat,2),size(kern_mat,3)),kern_mat,zeros(1,size(kern_mat,2),size(kern_mat,3)));
            kern_mat = cat(2,zeros(size(kern_mat,1),1,size(kern_mat,3)),kern_mat,zeros(size(kern_mat,1),1,size(kern_mat,3)));
            kern_mat = cat(3,zeros(size(kern_mat,1),size(kern_mat,2),1),kern_mat,zeros(size(kern_mat,1),size(kern_mat,2),1));
            
            lapl_oper = [0 1 0;1 -4 1;0 1 0];
            lapl = zeros(size(kern_mat,1)-2,size(kern_mat,2),size(kern_mat,3));
            for ii = 2:(size(kern_mat,1)-1)
                lapl(ii-1,:,:) = conv2(squeeze(kern_mat(ii,:,:)),lapl_oper,'same');
            end
            lapl = lapl(:,2:end-1,2:end-1);
            lpen.lapl(i) = mod.lambda_d2X*sum(lapl(:).^2);
            lpen.lapl(i) = lpen.lapl(i)/nspks;
        elseif strcmp(glmod.image_type,'1d')
            sdim = fsdim;
            kern_mat = reshape(cur_pkern,[kern_t sdim]);
            
            % buffer with zeros
            kern_mat = cat(1,zeros(1,size(kern_mat,2)),kern_mat,zeros(1,size(kern_mat,2)));
            kern_mat = cat(2,zeros(size(kern_mat,1)),kern_mat,zeros(size(kern_mat,1),1));
            
            lapl_oper = [1 -2 1];
            lapl = zeros(size(kern_mat,1)-2,size(kern_mat,2));
            for ii = 2:(size(kern_mat,1)-1)
                lapl(ii-1,:) = conv(squeeze(kern_mat(ii,:)),lapl_oper,'same');
            end
            lapl = lapl(:,2:end-1);
            lpen.lapl(i) = mod.lambda_d2X*sum(lapl(:).^2);
            lpen.lapl(i) = lpen.lapl(i)/nspks;
        end
    else
        lpen.lapl(i) = 0;
    end
    if mod.lambda_d2T > 0
        if strcmp(glmod.image_type,'0d')
            cur_pkern = mod.pix*mod.kscale; %reintroduce scale into k
            kern_t = length(cur_pkern)/mod.fsdim;
            
            kern_mat = reshape(cur_pkern,[kern_t 1]);
            % buffer with zeros
            kern_mat = cat(1,zeros(1,size(kern_mat,2)),kern_mat,zeros(1,size(kern_mat,2)));
            
            lapl_oper = [1 -2 1];
            lapl = zeros(size(kern_mat,1),1);
            lapl = conv(kern_mat,lapl_oper,'same');
            lapl = lapl(2:end-1);
            lpen.laplT(i) = mod.lambda_d2T*sum(lapl(:).^2);
            lpen.laplT(i) = lpen.lapl(i)/nspks;
        elseif strcmp(glmod.image_type,'2d')
            cur_pkern = mod.pix*mod.kscale;
            kern_t = length(cur_pkern)/mod.fsdim;
            sdim = sqrt(mod.fsdim);
            kern_mat = reshape(cur_pkern,[kern_t sdim sdim]);
            lapl = zeros(size(kern_mat));
            lapl(2:end-1,:,:) = kern_mat(1:end-2,:,:) + kern_mat(3:end,:,:) - 2*kern_mat(2:end-1,:,:);
            temp = lapl(2:end-1,2:end-1,2:end-1);
            lpen.laplT(i) = mod.lambda_d2T*sum(temp(:).^2);
            lpen.laplT(i) = lpen.laplT(i)/nspks;
        else
            error('Unsupported');
        end
    end
    %for space-time laplacian
    if mod.lambda_d2XT > 0
        %         cur_pkern = mod.pix;
        cur_pkern = mod.pix*mod.kscale; %reintroduce scale into k
        kern_t = length(cur_pkern)/mod.fsdim;
        if strcmp(glmod.image_type,'2d')
            error('Cant use space-time Laplacian with 2 spacial dims yet')
        end
        sdim = fsdim;
        kern_mat = reshape(cur_pkern,[kern_t sdim]);
        
        % buffer with zeros
        kern_mat = cat(1,zeros(1,size(kern_mat,2)),kern_mat,zeros(1,size(kern_mat,2)));
        kern_mat = cat(2,zeros(size(kern_mat,1)),kern_mat,zeros(size(kern_mat,1),1));
        
        lapl_oper = [0 1 0;1 -4 1;0 1 0];
        lapl = conv2(kern_mat,lapl_oper,'same');
        lapl = lapl(2:end-1,2:end-1);
        lpen.laplXT(i) = mod.lambda_d2XT*sum(lapl(:).^2);
        lpen.laplXT(i) = lpen.laplXT(i)/nspks;
    else
        lpen.laplXT(i) = 0;
    end
    %L1 kernel penalty
    if mod.lambda_L1x > 0
        %         cur_pkern = mod.pix;
        cur_pkern = mod.pix*mod.kscale; %reintroduce scale into k
        lpen.l1x(i) = sum(abs(cur_pkern).*mod.lambda_L1x/nspks);
    else
        lpen.l1x(i) = 0;
    end
    
end


%% penalty terms for spiketrain moduls?
nll = -ll/nspks;
pnll = nll + sum(lpen.w) + sum(lpen.loc) + sum(lpen.ksmooth) + sum(lpen.l1x) + sum(lpen.l2x) + sum(lpen.lapl) + sum(lpen.laplXT) + sum(lpen.laplT);

end