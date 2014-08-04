function [LL, LLgrad] = FULLBF_lexp_LLinternal_const(params, Robs, X, const, model)
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

b = const;
g = zeros(stimlen,1);
for n = 1:Nmods
    % Calculate convolution with internal receptive field
    if strcmp(model.basis,'white')
        gint{n} = X * (params((n-1)*klen + (1:klen))' * model.kern_conv_mat)';
    elseif strcmp(model.basis,'pix')
        gint{n} = X * params((n-1)*klen + (1:klen));
    end
    
    %process by lexp internal NL
    fgint = log(1+exp(model.mods(n).beta*gint{n}));
      
    %multiply by weight and add to generating function
    g = g + fgint*model.mods(n).w;
end
g(g > 100) = 100;
expg = exp(g + b)';
if strcmp(model.spk_nl,'logexp')
    r = log(1+expg);
elseif ~strcmp(model.spk_nl,'exp')
    r = expg;
    error('invalid spk nl');
end
r(r < 1e-20) = 1e-20; %minimum predicted rate

LL = sum(Robs .* log(r) - r);

if strcmp(model.spk_nl,'logexp')
    residual = (Robs./r - 1) .* expg ./ (1+expg);
else
    residual = Robs - r;
end

%initialize LL gradient
LLgrad = zeros(length(params),1);

% Calculate output of derivative module
chunk_size = 1000; %maximum chunk size for processing
for n = 1:Nmods
    
    %internal gen fun processed by f'
    g = model.mods(n).beta*exp(model.mods(n).beta*gint{n})./(1+exp(model.mods(n).beta*gint{n}));
    
    if strcmp(model.basis,'white')
        temp = X.*repmat(g,1,wklen)*model.mods(n).w;
        LLgrad(((n-1)*klen+1):n*klen) = (residual * temp)*model.pix_conv_mat;
    elseif strcmp(model.basis,'pix')
        NChunks = ceil(klen/chunk_size);
        for cc = 1:NChunks
            cs = (cc-1)*chunk_size + 1;
            ce = min(klen,cs+chunk_size);
            cur_chunklen = ce-cs+1;
            temp = X(:,cs:ce).*repmat(g,1,cur_chunklen)*model.mods(n).w;
            LLgrad(((n-1)*klen+cs):((n-1)*klen+ce)) = residual * temp;
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
loc_penalty = zeros(Nmods,1);
if strcmp(model.image_type,'1d')
    dist_mat = repmat(1:sdim,kern_t,1);
elseif strcmp(model.image_type,'2d')
    [xo,yo] = meshgrid(1:sdim,1:sdim);
    coords = [xo(:) yo(:)];
    dist_mat = repmat(coords,[1 1 kern_t]);
else
    error('Invalid image type')
end
for n = 1:Nmods
    cur_kern = params((n-1)*klen + (1:klen));
    if model.mods(n).lambda_dX > 0 || model.mods(n).locLambda > 0
        if strcmp(model.image_type,'2d')
            kern_mat{n} = reshape(cur_kern,[kern_t,sdim,sdim]);
            if model.mods(n).lambda_dX > 0
                x_derivs = (kern_mat{n}(:,2:end,:) - kern_mat{n}(:,1:end-1,:));
                y_derivs = (kern_mat{n}(:,:,2:end) - kern_mat{n}(:,:,1:end-1));
                t_derivs = (kern_mat{n}(2:end,:,:) - kern_mat{n}(1:end-1,:,:));
                smooth_penalty(n) = (model.mods(n).lambda_dX*sum(x_derivs(:).^2) + model.mods(n).lambda_dX*sum(y_derivs(:).^2) + model.mods(n).lambda_dT*sum(t_derivs(:).^2));
            end
            if model.mods(n).locLambda > 0
                com_dist_mat = squeeze(sum((dist_mat - repmat(model.mods(n).filt_com,[fsdim 1 kern_t])).^2,2))';
                pen_mat = exp(com_dist_mat/2/model.mods(n).locSigma^2);
                pen_mat(pen_mat > model.mods(n).maxLocPen) = model.mods(n).maxLocPen;
                pen_mat = reshape(kern_mat{n}.^2,[kern_t fsdim]).*pen_mat*model.mods(n).locLambda;
                loc_penalty(n) = sum(pen_mat(:));
            end
        elseif strcmp(model.image_type,'1d')
            kern_mat{n} = reshape(cur_kern,[kern_t,sdim]);
            if model.mods(n).lambda_dX > 0
                x_derivs = (kern_mat{n}(:,2:end) - kern_mat{n}(:,1:end-1));
                t_derivs = (kern_mat{n}(2:end,:) - kern_mat{n}(1:end-1,:));
                smooth_penalty(n) = (model.mods(n).lambda_dX*sum(x_derivs(:).^2) + model.mods(n).lambda_dT*sum(t_derivs(:).^2));
            end
            if model.mods(n).locLambda > 0
                pen_mat = exp((dist_mat - model.mods(n).filt_com).^2/2/model.mods(n).locSigma^2);
                pen_mat(pen_mat > model.mods(n).maxLocPen) = model.mods(n).maxLocPen;
                pen_mat = kern_mat{n}.^2.*pen_mat*model.mods(n).locLambda;
                loc_penalty(n) = sum(pen_mat(:));
            end
        end
    end
end

LL = LL - sum(smooth_penalty) - sum(loc_penalty);

% zerovec = zeros(1,klen)';
LLgrad_pen = zeros(size(LLgrad));
for n = 1:Nmods
    if model.mods(n).lambda_dX > 0
        if strcmp(model.image_type,'2d')
            
            [t_gradmat,x_gradmat,y_gradmat] = deal(zeros(size(kern_mat{n})));
            x_gradmat(:,2:end,:) = kern_mat{n}(:,2:end,:) - kern_mat{n}(:,1:end-1,:);
            x_gradmat(:,1:end-1,:) = x_gradmat(:,1:end-1,:) + kern_mat{n}(:,1:end-1,:) - kern_mat{n}(:,2:end,:);
            
            y_gradmat(:,:,2:end) = kern_mat{n}(:,:,2:end) - kern_mat{n}(:,:,1:end-1);
            y_gradmat(:,:,1:end-1) = y_gradmat(:,:,1:end-1) + kern_mat{n}(:,:,1:end-1) - kern_mat{n}(:,:,2:end);
            
            t_gradmat(2:end,:,:) = kern_mat{n}(2:end,:,:) - kern_mat{n}(1:end-1,:,:);
            t_gradmat(1:end-1,:,:) = t_gradmat(1:end-1,:,:) + kern_mat{n}(1:end-1,:,:) - kern_mat{n}(2:end,:,:);
            
            gradmat = model.mods(n).lambda_dX*x_gradmat + model.mods(n).lambda_dX*y_gradmat + model.mods(n).lambda_dT*t_gradmat;
        elseif strcmp(model.image_type,'1d')
            
            [t_gradmat,x_gradmat] = deal(zeros(size(kern_mat{n})));
            x_gradmat(:,2:end) = kern_mat{n}(:,2:end) - kern_mat{n}(:,1:end-1);
            x_gradmat(:,1:end-1) = x_gradmat(:,1:end-1) + kern_mat{n}(:,1:end-1) - kern_mat{n}(:,2:end);
                        
            t_gradmat(2:end,:) = kern_mat{n}(2:end,:) - kern_mat{n}(1:end-1,:);
            t_gradmat(1:end-1,:) = t_gradmat(1:end-1,:) + kern_mat{n}(1:end-1,:) - kern_mat{n}(2:end,:);
            
            gradmat = model.mods(n).lambda_dX*x_gradmat + model.mods(n).lambda_dT*t_gradmat;
        end
        LLgrad_pen((n-1)*klen + (1:klen)) = 2*gradmat(:);
    end
    if model.mods(n).locLambda > 0
        if strcmp(model.image_type,'2d')
            com_dist_mat = squeeze(sum((dist_mat - repmat(model.mods(n).filt_com,[fsdim 1 kern_t])).^2,2))';
            pen_mat = exp(com_dist_mat/2/model.mods(n).locSigma^2);
            pen_mat(pen_mat > model.mods(n).maxLocPen) = model.mods(n).maxLocPen;
            gradmat = 2*model.mods(n).locLambda*pen_mat(:).*kern_mat{n}(:);
        elseif strcmp(model.image_type,'1d')
           pen_mat = exp((dist_mat - model.mods(n).filt_com).^2/2/model.mods(n).locSigma^2);
           pen_mat(pen_mat > model.mods(n).maxLocPen) = model.mods(n).maxLocPen;
           gradmat = 2*model.mods(n).locLambda*pen_mat(:).*kern_mat{n}(:);
       else
           error('Invalid image type')
       end
       LLgrad_pen((n-1)*klen+(1:klen)) = LLgrad_pen((n-1)*klen+(1:klen)) + gradmat;
    end
end

LLgrad = LLgrad - LLgrad_pen;
%%****************************

Nspks = sum(Robs);
LL = -LL/Nspks;
LLgrad = -LLgrad/Nspks;

