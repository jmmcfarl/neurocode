function [LL, LLgrad] = FULLBF_lexp_LLinternal(params, Robs, X, model, fprimes)
%
% Usage: [LL, LLgrad] = RG_LLinternal_elog( params, Robs, X, Skb, model, targets, lamrange, fprimes )
%

Nmods = length(model.mods);
[NT,wklen] = size(X);
if strcmp(model.basis,'white')
    klen = size(model.pix_conv_mat,2);
elseif strcmp(model.basis,'pix')
    klen = size(X,2);
end
stimlen = size(X,1);

hlen = 1; %no psc term
flen = klen/model.mods(1).SDIM;

b = params(end);
g = zeros(stimlen,1);
min_x = min(model.mods(1).nlx);
for n = 1:Nmods
    % Calculate convolution with internal receptive field
    if strcmp(model.basis,'white')
        gint{n} = X * (params((n-1)*klen + (1:klen))' * model.kern_conv_mat)';
    elseif strcmp(model.basis,'pix')
        gint{n} = X * params((n-1)*klen + (1:klen));
    end
    
    if strcmp(model.mods(n).nltype,'lexp')
        bgint = (gint{n}-model.mods(n).theta)*model.mods(n).beta;
        too_large = bgint > 50;
        fgint = 1/model.mods(n).beta*log(1+exp(bgint)) - 1/model.mods(n).beta*log(1+exp(model.mods(n).beta*min_x));
        fgint(too_large) = 1/model.mods(n).beta*bgint(too_large) - 1/model.mods(n).beta*min_x;
    elseif strcmp(model.mods(n).nltype,'uncon')
        fgint = nlin_proc_stim(gint{n},model.mods(n).nly,model.mods(n).nlx);
    elseif strcmp(model.mods(n).nltype,'quad')
        fgint = gint{n}.^2;
    elseif strcmp(model.mods(n).nltype,'lin')
        fgint = gint{n};
    elseif strcmp(model.mods(n).nltype,'threshlin')
        fgint = gint{n};
        fgint(fgint < 0) = 0;
    elseif strcmp(model.mods(n).nltype,'rquad')
        fgint = gint{n}.^2;
        fgint(gint{n} < 0) = 0;
    else
        error('Invalid internal NL');
    end
    
    %multiply by weight and add to generating function
    g = g + fgint*model.mods(n).w;
end
g(g > 100) = 100;
expg = exp(g + b)';
if strcmp(model.spk_nl,'logexp')
    r = log(1+expg);
elseif strcmp(model.spk_nl,'exp')
    r = expg;
else
    error('invalid spk nl');
end
r(r < 1e-20) = 1e-20; %minimum predicted rate

LL = sum(Robs(flen:end) .* log(r(flen:end)) - r(flen:end));

if strcmp(model.spk_nl,'logexp')
    residual = (Robs./r - 1) .* expg ./ (1+expg);
else
    residual = Robs - r;
end

%initialize LL gradient
LLgrad = zeros(length(params),1);

% Calculate derivatives with respect to constant
LLgrad(end) = sum(residual);

if model.mods(1).locLambda > 0
    if strcmp(model.image_type,'1d')
        [model,COMs] = get_filter_coms_1d(model);
    end
end

% Calculate output of derivative module
chunk_size = 1000; %maximum chunk size for processing
for n = 1:Nmods
    
    if strcmp(model.mods(n).nltype,'lexp')
        %internal gen fun processed by f'
        %     g = exp(model.mods(n).beta*gint{n})./(1+exp(model.mods(n).beta*gint{n}));
        bgint = (gint{n}-model.mods(n).theta)*model.mods(n).beta;
        g = exp(bgint)./(1+exp(bgint));
        too_large = bgint > 50;
        g(too_large) = 1;
    elseif strcmp(model.mods(n).nltype,'uncon')
        g = fprime(gint{n},fprimes{n}, model.mods(n).nlx);
    elseif strcmp(model.mods(n).nltype,'quad')
        g = 2*gint{n};
    elseif strcmp(model.mods(n).nltype,'lin')
        g = ones(size(gint{n}));
    elseif strcmp(model.mods(n).nltype,'threshlin')
        g = ones(size(gint{n}));
        g(gint{n} < 0) = 0;
    elseif strcmp(model.mods(n).nltype,'rquad')
        g = 2*gint{n};
        g(gint{n} < 0) = 0;
    end
    
    if strcmp(model.basis,'white')
        temp = X.*repmat(g,1,wklen)*model.mods(n).w;
        LLgrad(((n-1)*klen+1):n*klen) = (residual * temp)*model.pix_conv_mat;
    elseif strcmp(model.basis,'pix')
        
        if klen <= chunk_size
            %         temp = bsxfun(@times,X,g)*model.mods(n).w;
            temp = X.*repmat(g,1,wklen)*model.mods(n).w;
            %         LLgrad(((n-1)*klen+1):(n*klen)) = (residual(flen:end) * temp(flen:end,:));
            LLgrad(((n-1)*klen+1):(n*klen)) = residual*temp;
            
        else
            NChunks = ceil(klen/chunk_size);
            for cc = 1:NChunks
                cs = (cc-1)*chunk_size + 1;
                ce = min(klen,cs+chunk_size);
                cur_chunklen = ce-cs+1;
                temp = X(:,cs:ce).*repmat(g,1,cur_chunklen)*model.mods(n).w;
                %                 LLgrad(((n-1)*klen+cs):((n-1)*klen+ce)) =
                %                 residual(flen:end) * temp(flen:end,:); %this is
                %                 technically correct, but much slower
                LLgrad(((n-1)*klen+cs):((n-1)*klen+ce)) = residual*temp;
            end
        end
    end
end


%%********Add penalty terms for smoothness*******
fsdim = model.mods(1).fsdim;
kern_t = klen/fsdim;
if strcmp(model.image_type,'2d')
    sdim = sqrt(fsdim);
elseif strcmp(model.image_type,'1d')
    sdim = fsdim;
end

smooth_penalty = zeros(Nmods,1);
lapl_penalty = zeros(Nmods,1);
laplT_penalty = zeros(Nmods,1);
laplXT_penalty = zeros(Nmods,1);
loc_penalty = zeros(Nmods,1);
l2_penalty = zeros(Nmods,1);
if strcmp(model.image_type,'1d')
    dist_mat = repmat(1:sdim,kern_t,1);
elseif strcmp(model.image_type,'2d')
    [xo,yo] = meshgrid(1:sdim,1:sdim);
    coords = [xo(:) yo(:)];
    dist_mat = repmat(coords,[1 1 kern_t]);
elseif strcmp(model.image_type,'0d')
    dist_mat = nan;
else
    error('Invalid image type')
end
for n = 1:Nmods
    cur_kern = params((n-1)*klen + (1:klen));
    if model.mods(n).lambda_dX > 0 || model.mods(n).locLambda > 0 || model.mods(n).lambda_dT > 0 || model.mods(n).lambda_d2X > 0 || model.mods(n).lambda_d2XT > 0 || model.mods(n).lambda_d2T > 0
        if strcmp(model.image_type,'2d')
            %             kern_mat{n} = reshape(cur_kern,[kern_t,sdim,sdim]);
            kern_mat{n} = model.mods(n).kscale*reshape(cur_kern,[kern_t,sdim,sdim]); %reintroduce w into scale of k
            
            %buffer with zeros
            kern_mat{n} = cat(1,zeros(1,size(kern_mat{n},2),size(kern_mat{n},3)),kern_mat{n},zeros(1,size(kern_mat{n},2),size(kern_mat{n},3)));
            kern_mat{n} = cat(2,zeros(size(kern_mat{n},1),1,size(kern_mat{n},3)),kern_mat{n},zeros(size(kern_mat{n},1),1,size(kern_mat{n},3)));
            kern_mat{n} = cat(3,zeros(size(kern_mat{n},1),size(kern_mat{n},2),1),kern_mat{n},zeros(size(kern_mat{n},1),size(kern_mat{n},2),1));
            
            %first spatial derivative (L2 norm of gradient)
            if model.mods(n).lambda_dX > 0
                x_derivs = (kern_mat{n}(:,2:end,:) - kern_mat{n}(:,1:end-1,:));
                y_derivs = (kern_mat{n}(:,:,2:end) - kern_mat{n}(:,:,1:end-1));
                t_derivs = (kern_mat{n}(2:end,:,:) - kern_mat{n}(1:end-1,:,:));
                smooth_penalty(n) = (model.mods(n).lambda_dX*sum(x_derivs(:).^2) + model.mods(n).lambda_dX*sum(y_derivs(:).^2) + model.mods(n).lambda_dT*sum(t_derivs(:).^2));
            end
            % Discrete Laplacian
            if model.mods(n).lambda_d2X > 0
                lapl_oper = [0 1 0;1 -4 1;0 1 0];
                lapl{n} = zeros(size(kern_mat{n},1)-2,size(kern_mat{n},2),size(kern_mat{n},3));
                for i = 2:(size(kern_mat{n},1)-1)
                    lapl{n}(i-1,:,:) = conv2(squeeze(kern_mat{n}(i,:,:)),lapl_oper,'same');
                end
                temp = lapl{n}(:,2:end-1,2:end-1);
                lapl_penalty(n) = model.mods(n).lambda_d2X*sum(temp(:).^2);
            end
            if model.mods(n).lambda_d2T > 0
                lapl_oper = [1 -2 1];
                lapl_T{n} = zeros(size(kern_mat{n}));
                lapl_T{n}(2:end-1,:,:) = kern_mat{n}(1:end-2,:,:) + kern_mat{n}(3:end,:,:) - 2*kern_mat{n}(2:end-1,:,:);
                temp = lapl_T{n}(2:end-1,2:end-1,2:end-1);
                laplT_penalty(n) = model.mods(n).lambda_d2T*sum(temp(:).^2);
            end
            if model.mods(n).locLambda > 0
                com_dist_mat = squeeze(sum((dist_mat - repmat(model.mods(n).filt_com,[fsdim 1 kern_t])).^2,2))';
                pen_mat = exp(com_dist_mat/2/model.mods(n).locSigma^2);
                pen_mat(pen_mat > model.mods(n).maxLocPen) = model.mods(n).maxLocPen;
                cur_kern = kern_mat{n}(2:end-1,2:end-1,2:end-1);
                pen_mat = reshape(cur_kern.^2,[kern_t fsdim]).*pen_mat*model.mods(n).locLambda;
                loc_penalty(n) = sum(pen_mat(:));
            end
        elseif strcmp(model.image_type,'1d')
            %             kern_mat{n} = reshape(cur_kern,[kern_t,sdim]);
            kern_mat{n} = model.mods(n).kscale*reshape(cur_kern,[kern_t,sdim]); %reintroduce w into the scale of k
            
            %buffer with zeros
            kern_mat{n} = cat(1,zeros(1,size(kern_mat{n},2)),kern_mat{n},zeros(1,size(kern_mat{n},2)));
            kern_mat{n} = cat(2,zeros(size(kern_mat{n},1),1),kern_mat{n},zeros(size(kern_mat{n},1),1));
            if model.mods(n).lambda_dX > 0
                x_derivs = (kern_mat{n}(:,2:end) - kern_mat{n}(:,1:end-1));
                t_derivs = (kern_mat{n}(2:end,:) - kern_mat{n}(1:end-1,:));
                
                smooth_penalty(n) = (model.mods(n).lambda_dX*sum(x_derivs(:).^2) + model.mods(n).lambda_dT*sum(t_derivs(:).^2));
            end
            % Discrete 1d Laplacian
            if model.mods(n).lambda_d2X > 0
                lapl_oper = [1 -2 1];
                lapl{n} = zeros(size(kern_mat{n},1)-2,size(kern_mat{n},2));
                for i = 2:(size(kern_mat{n},1)-1)
                    lapl{n}(i-1,:) = conv(squeeze(kern_mat{n}(i,:)),lapl_oper,'same');
                end
                temp = lapl{n}(:,2:end-1);
                lapl_penalty(n) = model.mods(n).lambda_d2X*sum(temp(:).^2);
            end
            % Space-time Laplacian
            if model.mods(n).lambda_d2XT > 0
                lapl_oper_XT = [0 1 0;1 -4 1;0 1 0];
                lapl_XT{n} = conv2(kern_mat{n},lapl_oper_XT,'same');
                temp = lapl_XT{n}(2:end-1,2:end-1);
                laplXT_penalty(n) = model.mods(n).lambda_d2XT*sum(temp(:).^2);
            end
            %localization pen
            if model.mods(n).locLambda > 0
                pen_mat = exp((dist_mat - model.mods(n).filt_com).^2/2/model.mods(n).locSigma^2);
                pen_mat(pen_mat > model.mods(n).maxLocPen) = model.mods(n).maxLocPen;
                cur_kern = kern_mat{n}(2:end-1,2:end-1);
                pen_mat = cur_kern.^2.*pen_mat*model.mods(n).locLambda;
                loc_penalty(n) = sum(pen_mat(:));
            end
        end
    end
    if isfield(model.mods(n),'lambda_L2x')
        l2_penalty(n) = 0.5*model.mods(n).lambda_L2x*sum(cur_kern.^2);
    end
    if strcmp(model.image_type,'0d')
        kern_mat{n} = model.mods(n).kscale*reshape(cur_kern,[kern_t,1]);
        %buffer with zeros
        kern_mat{n} = cat(1,zeros(1,size(kern_mat{n},2)),kern_mat{n},zeros(1,size(kern_mat{n},2)));
        if model.mods(n).lambda_dT > 0 || model.mods(n).lambda_d2T > 0
            %             kern_mat{n} = model.mods(n).kscale*reshape(cur_kern,[kern_t,1]);
            
            t_derivs = kern_mat{n}(2:end,:) - kern_mat{n}(1:end-1,:);
            smooth_penalty(n) = model.mods(n).lambda_dT*sum(t_derivs(:).^2);
            
            if model.mods(n).lambda_d2T > 0
                lapl_oper = [1 -2 1];
                lapl{n} = conv(kern_mat{n},lapl_oper,'same');
                temp = lapl{n}(:,2:end-1);
                laplT_penalty(n) = model.mods(n).lambda_d2T*sum(temp(:).^2);
            end
        end
    end
end

LL = LL - sum(smooth_penalty) - sum(loc_penalty) - sum(l2_penalty) - sum(lapl_penalty) - ...
    sum(laplXT_penalty) - sum(laplT_penalty);

% zerovec = zeros(1,klen)';
LLgrad_pen = zeros(size(LLgrad));
for n = 1:Nmods
    if strcmp(model.image_type,'2d') || strcmp(model.image_type,'1d')
        if model.mods(n).lambda_dX > 0 || model.mods(n).lambda_dT > 0
            if strcmp(model.image_type,'2d')
                
                [t_gradmat,x_gradmat,y_gradmat] = deal(zeros(size(kern_mat{n})));
                x_gradmat(:,2:end,:) = kern_mat{n}(:,2:end,:) - kern_mat{n}(:,1:end-1,:);
                y_gradmat(:,:,2:end) = kern_mat{n}(:,:,2:end) - kern_mat{n}(:,:,1:end-1);
                
                t_gradmat(2:end,:,:) = kern_mat{n}(2:end,:,:) - kern_mat{n}(1:end-1,:,:);
                gradmat = model.mods(n).lambda_dX*x_gradmat(2:end-1,2:end-1,2:end-1) + model.mods(n).lambda_dX*y_gradmat(2:end-1,2:end-1,2:end-1)...
                    + model.mods(n).lambda_dT*t_gradmat(2:end-1,2:end-1,2:end-1);
            elseif strcmp(model.image_type,'1d')
                
                [t_gradmat,x_gradmat] = deal(zeros(size(kern_mat{n})));
                x_gradmat(:,2:end) = kern_mat{n}(:,2:end) - kern_mat{n}(:,1:end-1);
                
                t_gradmat(2:end,:) = kern_mat{n}(2:end,:) - kern_mat{n}(1:end-1,:);
                gradmat = model.mods(n).lambda_dX*x_gradmat(2:end-1,2:end-1) + model.mods(n).lambda_dT*t_gradmat(2:end-1,2:end-1);
            end
            LLgrad_pen((n-1)*klen + (1:klen)) = 2*gradmat(:);
        end
    end
    if model.mods(n).lambda_d2X > 0
        if strcmp(model.image_type,'2d')
            temp = zeros(size(lapl{n}));
            for i = 1:size(lapl{n},1)
                temp(i,:,:) = conv2(squeeze(lapl{n}(i,:,:)),lapl_oper,'same');
            end
            temp = temp(:,2:end-1,2:end-1);
            LLgrad_pen((n-1)*klen + (1:klen)) = LLgrad_pen((n-1)*klen + (1:klen)) + 2*model.mods(n).lambda_d2X*temp(:);
        elseif strcmp(model.image_type,'1d')
            temp = zeros(size(lapl{n}));
            for i = 1:size(lapl{n},1)
                temp(i,:) = conv(squeeze(lapl{n}(i,:)),lapl_oper,'same');
            end
            temp = temp(:,2:end-1);
            LLgrad_pen((n-1)*klen + (1:klen)) = LLgrad_pen((n-1)*klen + (1:klen)) + 2*model.mods(n).lambda_d2X*temp(:);
        end
    end
    if model.mods(n).lambda_d2T > 0
        if strcmp(model.image_type,'0d')
            t_gradmat = zeros(size(kern_mat{n}));
            t_gradmat(2:end,:) = kern_mat{n}(2:end) - kern_mat{n}(1:end-1);
            gradmat = model.mods(n).lambda_dT*t_gradmat(2:end-1);
            LLgrad_pen((n-1)*klen + (1:klen)) = 2*gradmat(:);
            
            if model.mods(n).lambda_d2T > 0
                temp = conv(squeeze(lapl{n}),lapl_oper,'same');
                temp = temp(2:end-1);
                LLgrad_pen((n-1)*klen + (1:klen)) = LLgrad_pen((n-1)*klen + (1:klen)) + 2*model.mods(n).lambda_d2T*temp(:);
            end
        elseif strcmp(model.image_type,'2d')
            temp = zeros(size(lapl_T{n}));
            temp(2:end-1,:,:) = lapl_T{n}(1:end-2,:,:) + lapl_T{n}(3:end,:,:) - 2*lapl_T{n}(2:end-1,:,:);
            temp = temp(2:end-1,2:end-1,2:end-1);
            LLgrad_pen((n-1)*klen + (1:klen)) = LLgrad_pen((n-1)*klen + (1:klen)) + 2*model.mods(n).lambda_d2T*temp(:);
        else
            error('Not supported');
        end
    end
    if model.mods(n).lambda_d2XT > 0
        if strcmp(model.image_type,'2d')
            error('Cant use space-time laplacian in 2 spatial dims yet');
        elseif strcmp(model.image_type,'1d')
            temp = conv2(squeeze(lapl_XT{n}),lapl_oper_XT,'same');
            temp = temp(2:end-1,2:end-1);
            LLgrad_pen((n-1)*klen + (1:klen)) = LLgrad_pen((n-1)*klen + (1:klen)) + 2*model.mods(n).lambda_d2XT*temp(:);
        end
    end
    if model.mods(n).locLambda > 0
        if strcmp(model.image_type,'2d')
            com_dist_mat = squeeze(sum((dist_mat - repmat(model.mods(n).filt_com,[fsdim 1 kern_t])).^2,2))';
            pen_mat = exp(com_dist_mat/2/model.mods(n).locSigma^2);
            pen_mat(pen_mat > model.mods(n).maxLocPen) = model.mods(n).maxLocPen;
            temp_kern = kern_mat{n}(2:end-1,2:end-1,2:end-1);
            gradmat = 2*model.mods(n).locLambda*pen_mat(:).*temp_kern(:);
        elseif strcmp(model.image_type,'1d')
            pen_mat = exp((dist_mat - model.mods(n).filt_com).^2/2/model.mods(n).locSigma^2);
            pen_mat(pen_mat > model.mods(n).maxLocPen) = model.mods(n).maxLocPen;
            temp_kern = kern_mat{n}(2:end-1,2:end-1);
            gradmat = 2*model.mods(n).locLambda*pen_mat(:).*temp_kern(:);
        else
            error('Invalid image type')
        end
        LLgrad_pen((n-1)*klen+(1:klen)) = LLgrad_pen((n-1)*klen+(1:klen)) + gradmat;
    end
    if isfield(model.mods(n),'lambda_L2x')
        cur_kern = params((n-1)*klen + (1:klen));
        LLgrad_pen((n-1)*klen+(1:klen)) = LLgrad_pen((n-1)*klen+(1:klen)) + model.mods(n).lambda_L2x*cur_kern;
    end
end
% eps = 1e-8;
% LLgrad_locpen = zeros(size(LLgrad));
% for n = 1:Nmods
%     cur_cfs = params((n-1)*klen + (1:klen));
%     %     for numerical gradient calculation
%     for j = 1:model.klen
%         cur_cfs_eps = cur_cfs;
%         cur_cfs_eps(j) = cur_cfs_eps(j) + eps;
%         ent_dx = kernel_std(model.STCbasis,cur_cfs_eps,sdim);
%         ent = kernel_std(model.STCbasis,cur_cfs,sdim);
%         LLgrad_locpen((n-1)*klen+j) = model.mods(n).locLambda*(ent_dx-ent)/eps;
%     end
%
% end

LLgrad = LLgrad - LLgrad_pen;
%%****************************

Nspks = sum(Robs);
LL = -LL/Nspks;
LLgrad = -LLgrad/Nspks;

