function [LL, LLgrad] = fitGNM_filters_internal(params, Robs, X, model, fprimes)

Nmods = length(model.mods);
[stimlen,klen] = size(X);
flen = model.stim_params.flen;
fsdim = model.stim_params.fsdim;
sdim = model.stim_params.sdim;

b = params(end);
g = zeros(stimlen,1);
for n = 1:Nmods
    % Calculate convolution with internal receptive field
    gint{n} = X * params((n-1)*klen + (1:klen));
    
    if strcmp(model.mods(n).nltype,'lexp')
        min_x = min(model.mods(1).nlx);
        bgint = (gint{n}-model.mods(n).lexp_theta)*model.mods(n).lexp_beta;
        too_large = bgint > 50;
        fgint = 1/model.mods(n).lexp_beta*log(1+exp(bgint)) - 1/model.mods(n).lexp_beta*log(1+exp(model.mods(n).lexp_beta*min_x));
        fgint(too_large) = 1/model.mods(n).lexp_beta*bgint(too_large) - 1/model.mods(n).lexp_beta*min_x;
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
g = g + b;

if strcmp(model.spk_nl,'logexp')
    internal = model.spk_beta*g;
    too_big = find(internal > 60);
    expg = exp(internal)';
    r = model.spk_alpha*log(1+expg);
    r(too_big) = model.spk_alpha*internal(too_big);
elseif strcmp(model.spk_nl,'exp')
    if max(g) > 100
        g(g > 100) = 100;%max g value to prevent overflow
        disp('Warning, overflow in g')
    end
    
    expg = exp(g)';
    r = expg;
elseif strcmp(model.spk_nl,'gauss')
    r = g;
else
    error('invalid spk nl');
end

if ~strcmp(model.spk_nl,'gauss')
    if min(r) < 1e-20
        r(r < 1e-20) = 1e-20; %minimum predicted rate
        disp('Warning, underflow in predicted r')
    end
    LL = sum(Robs.* log(r) - r); %point process likelihood
else
    LL = -sum((Robs - r).^2); %SSE for Gaussian likelihood
end


if strcmp(model.spk_nl,'logexp')
    residual = model.spk_alpha*model.spk_beta*(Robs./r - 1) .* expg ./ (1+expg);
    residual(too_big) = model.spk_alpha*model.spk_beta*(Robs(too_big)./r(too_big)-1);
elseif strcmp(model.spk_nl,'exp')
    residual = Robs - r;
elseif strcmp(model.spk_nl,'gauss')
    residual = 2*(Robs-r)';
else
    error('Unsupported spiking NL')
end

%initialize LL gradient
LLgrad = zeros(length(params),1);

% Calculate derivatives with respect to constant
LLgrad(end) = sum(residual);

% Calculate output of derivative module
chunk_size = 1000; %maximum chunk size for processing
for n = 1:Nmods
    
    if strcmp(model.mods(n).nltype,'lexp')
        %internal gen fun processed by f'
        bgint = (gint{n}-model.mods(n).lexp_theta)*model.mods(n).lexp_beta;
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
    
    if klen <= chunk_size
        temp = X.*repmat(g,1,klen)*model.mods(n).w;
        LLgrad(((n-1)*klen+1):(n*klen)) = residual*temp;
    else
        NChunks = ceil(klen/chunk_size);
        for cc = 1:NChunks
            cs = (cc-1)*chunk_size + 1;
            ce = min(klen,cs+chunk_size);
            cur_chunklen = ce-cs+1;
            temp = X(:,cs:ce).*repmat(g,1,cur_chunklen)*model.mods(n).w;
            LLgrad(((n-1)*klen+cs):((n-1)*klen+ce)) = residual*temp;
        end
    end
end


%%********Add penalty terms*******
smooth_penalty = zeros(Nmods,1);
lapl_penalty = zeros(Nmods,1);
laplT_penalty = zeros(Nmods,1);
laplXT_penalty = zeros(Nmods,1);
l2_penalty = zeros(Nmods,1);

for n = 1:Nmods
    cur_kern = params((n-1)*klen + (1:klen));
    if model.stim_params.spatial_dims == 2
        kern_mat{n} = model.mods(n).kscale*reshape(cur_kern,[flen,sdim,sdim]);
        %buffer with zeros along each dimension
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
    elseif model.stim_params.spatial_dims == 1
        kern_mat{n} = model.mods(n).kscale*reshape(cur_kern,[flen,sdim]); %reintroduce w into the scale of k
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
        if model.mods(n).lambda_d2T > 0
            lapl_oper = [1 -2 1];
            lapl_T{n} = zeros(size(kern_mat{n}));
            lapl_T{n}(2:end-1,:) = kern_mat{n}(1:end-2,:) + kern_mat{n}(3:end,:) - 2*kern_mat{n}(2:end-1,:);
            temp = lapl_T{n}(2:end-1,2:end-1);
            laplT_penalty(n) = model.mods(n).lambda_d2T*sum(temp(:).^2);
        end
        
    elseif model.stim_params.spatial_dims == 0
        kern_mat{n} = model.mods(n).kscale*reshape(cur_kern,[flen,1]);
        %buffer with zeros
        kern_mat{n} = cat(1,zeros(1,size(kern_mat{n},2)),kern_mat{n},zeros(1,size(kern_mat{n},2)));
        
        if model.mods(n).lambda_dT > 0
            t_derivs = kern_mat{n}(2:end,:) - kern_mat{n}(1:end-1,:);
            smooth_penalty(n) = model.mods(n).lambda_dT*sum(t_derivs(:).^2);
        end
        if model.mods(n).lambda_d2T > 0
            lapl_oper = [1 -2 1];
            lapl{n} = conv(kern_mat{n},lapl_oper,'same');
            temp = lapl{n}(:,2:end-1);
            laplT_penalty(n) = model.mods(n).lambda_d2T*sum(temp(:).^2);
        end
    end
end
if isfield(model.mods(n),'lambda_L2x')
    l2_penalty(n) = 0.5*model.mods(n).lambda_L2x*sum(cur_kern.^2);
end


LL = LL - sum(smooth_penalty) - sum(l2_penalty) - sum(lapl_penalty) - ...
    sum(laplXT_penalty) - sum(laplT_penalty);

LLgrad_pen = zeros(size(LLgrad));
for n = 1:Nmods
    if model.stim_params.spatial_dims == 2
        if model.mods(n).lambda_dX > 0 || model.mods(n).lambda_dT > 0
            
            [t_gradmat,x_gradmat,y_gradmat] = deal(zeros(size(kern_mat{n})));
            x_gradmat(:,2:end,:) = kern_mat{n}(:,2:end,:) - kern_mat{n}(:,1:end-1,:);
            y_gradmat(:,:,2:end) = kern_mat{n}(:,:,2:end) - kern_mat{n}(:,:,1:end-1);
            
            t_gradmat(2:end,:,:) = kern_mat{n}(2:end,:,:) - kern_mat{n}(1:end-1,:,:);
            gradmat = model.mods(n).lambda_dX*x_gradmat(2:end-1,2:end-1,2:end-1) + model.mods(n).lambda_dX*y_gradmat(2:end-1,2:end-1,2:end-1)...
                + model.mods(n).lambda_dT*t_gradmat(2:end-1,2:end-1,2:end-1);
            
            LLgrad_pen((n-1)*klen + (1:klen)) = 2*gradmat(:);
        end
        if model.mods(n).lambda_d2X > 0
            temp = zeros(size(lapl{n}));
            for i = 1:size(lapl{n},1)
                temp(i,:,:) = conv2(squeeze(lapl{n}(i,:,:)),lapl_oper,'same');
            end
            temp = temp(:,2:end-1,2:end-1);
            LLgrad_pen((n-1)*klen + (1:klen)) = LLgrad_pen((n-1)*klen + (1:klen)) + 2*model.mods(n).lambda_d2X*temp(:);
        end
        if model.mods(n).lambda_d2T > 0
            temp = zeros(size(lapl_T{n}));
            temp(2:end-1,:,:) = lapl_T{n}(1:end-2,:,:) + lapl_T{n}(3:end,:,:) - 2*lapl_T{n}(2:end-1,:,:);
            temp = temp(2:end-1,2:end-1,2:end-1);
            LLgrad_pen((n-1)*klen + (1:klen)) = LLgrad_pen((n-1)*klen + (1:klen)) + 2*model.mods(n).lambda_d2T*temp(:);
        end
    elseif model.stim_params.spatial_dims == 1
        if model.mods(n).lambda_dX > 0 || model.mods(n).lambda_dT > 0
            
            [t_gradmat,x_gradmat] = deal(zeros(size(kern_mat{n})));
            x_gradmat(:,2:end) = kern_mat{n}(:,2:end) - kern_mat{n}(:,1:end-1);
            
            t_gradmat(2:end,:) = kern_mat{n}(2:end,:) - kern_mat{n}(1:end-1,:);
            gradmat = model.mods(n).lambda_dX*x_gradmat(2:end-1,2:end-1) + model.mods(n).lambda_dT*t_gradmat(2:end-1,2:end-1);
            
            LLgrad_pen((n-1)*klen + (1:klen)) = 2*gradmat(:);
        end
        if model.mods(n).lambda_d2X > 0
            temp = zeros(size(lapl{n}));
            for i = 1:size(lapl{n},1)
                temp(i,:) = conv(squeeze(lapl{n}(i,:)),lapl_oper,'same');
            end
            temp = temp(:,2:end-1);
            LLgrad_pen((n-1)*klen + (1:klen)) = LLgrad_pen((n-1)*klen + (1:klen)) + 2*model.mods(n).lambda_d2X*temp(:);
        end
        if model.mods(n).lambda_d2XT > 0
            temp = conv2(squeeze(lapl_XT{n}),lapl_oper_XT,'same');
            temp = temp(2:end-1,2:end-1);
            LLgrad_pen((n-1)*klen + (1:klen)) = LLgrad_pen((n-1)*klen + (1:klen)) + 2*model.mods(n).lambda_d2XT*temp(:);
        end
        if model.mods(n).lambda_d2T > 0
            temp = conv2(squeeze(lapl_T{n}),lapl_oper,'same');
            temp = temp(2:end-1,2:end-1);
            LLgrad_pen((n-1)*klen + (1:klen)) = LLgrad_pen((n-1)*klen + (1:klen)) + 2*model.mods(n).lambda_d2T*temp(:);
        end
    elseif model.stim_params.spatial_dims == 0
        if model.mods(n).lambda_dT > 0
            t_gradmat = zeros(size(kern_mat{n}));
            t_gradmat(2:end) = kern_mat{n}(2:end) - kern_mat{n}(1:end-1);
            gradmat = model.mods(n).lambda_dT*t_gradmat(2:end-1);
            LLgrad_pen((n-1)*klen + (1:klen)) = 2*gradmat(:);
        end
        if model.mods(n).lambda_d2T > 0
            %             temp = conv(squeeze(lapl{n}),lapl_oper,'same');
            %     temp = temp(2:end-1);
            temp = zeros(size(lapl{n}));
            temp(2:end-1) = -2*lapl{n}(2:end-1) + lapl{n}(1:end-2) + lapl{n}(3:end);
            temp(1) = -lapl{n}(2);
            temp(end) = -lapl{n}(end-1);
            temp = temp(2:end-1);
            LLgrad_pen((n-1)*klen + (1:klen)) = LLgrad_pen((n-1)*klen + (1:klen)) + 2*model.mods(n).lambda_d2T*temp(:);        end
    end
end

if isfield(model.mods(n),'lambda_L2x')
    cur_kern = params((n-1)*klen + (1:klen));
    LLgrad_pen((n-1)*klen+(1:klen)) = LLgrad_pen((n-1)*klen+(1:klen)) + model.mods(n).lambda_L2x*cur_kern;
end

LLgrad = LLgrad - LLgrad_pen;
%%****************************

if ~strcmp(model.spk_nl,'gauss')
    Nspks = sum(Robs);
else
    Nspks = length(Robs);
end
LL = -LL/Nspks;
LLgrad = -LLgrad/Nspks;

