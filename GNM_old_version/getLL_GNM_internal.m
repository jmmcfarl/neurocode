function [nll, pnll,lpen,rx] = getLL_GNM_internal(g,glmod,spkbs)
%% computes likelihood of generator potential
%% uses penalty terms for both nl and history parameters

nmods = length(glmod.mods);
if nmods > 0
    klen = length(glmod.mods(1).k);
end
flen = glmod.stim_params.flen;
fsdim = glmod.stim_params.fsdim;
sdim = glmod.stim_params.sdim;

nspks = length(spkbs);

%% pure ll based on model prediction (e.g. likelihood)
if strcmp(glmod.spk_nl,'logexp')
    %     rx = glmod.spk_alpha*log(1+exp(glmod.spk_beta*(g - glmod.spk_theta))); %apply spiking NL tp get predicted rate
    internal = glmod.spk_beta*g;
    too_big = find(internal > 60);
    rx = glmod.spk_alpha*log(1+exp(internal)); %apply spiking NL tp get predicted rate
    rx(too_big) = glmod.spk_alpha*internal(too_big);
elseif strcmp(glmod.spk_nl,'exp')
    if max(g) > 100
        g(g > 100) = 100; %saturate input to spiking NL
        disp('Overflow detected')
    end
    rx = exp(g);
elseif strcmp(glmod.spk_nl,'gauss')
    rx = g;
else
    error('Not accepted spiking NL')
end
if ~strcmp(glmod.spk_nl,'gauss')
    rs = rx(spkbs); %predicted rate at spike times
    if min(rx) < 1e-20
        disp('Underflow detected')
        rs(rs < 1e-20) = 1e-20; %impose minimum on predicted rate
        rx(rx < 1e-20) = 1e-20;
    end
    %% likelihood given by point process model
    ll = sum(log(rs)) - sum(rx);
else
   ll = -sum((spkbs - rx).^2);
end

%% regularisaton penalty (e.g. enforcing smth like posterior)
%initialize penalty terms
lpen.l2x = zeros(nmods,1);
lpen.grad = zeros(nmods,1);
lpen.lapl = zeros(nmods,1);
lpen.laplXT = zeros(nmods,1);
lpen.l1x = zeros(nmods,1);
for i = 1:length(glmod.mods)
    mod = glmod.mods(i); %for current module
    if isfield(mod,'lambda_L2x')
        lpen.l2x(i) = 0.5*mod.lambda_L2x*sum((mod.w*mod.k).^2)/nspks; %scale reintroduced to k
    end
    %kernel smoothness penalty
    if mod.lambda_dX > 0 || mod.lambda_dT > 0
        cur_kern = mod.k*mod.kscale; %reintroduce scale into k
        if glmod.stim_params.spatial_dims == 2
            kern_mat = reshape(cur_kern,[flen,sdim,sdim]);
            %zero pad kernel matrix
            kern_mat = cat(1,zeros(1,size(kern_mat,2),size(kern_mat,3)),kern_mat,zeros(1,size(kern_mat,2),size(kern_mat,3)));
            kern_mat = cat(2,zeros(size(kern_mat,1),1,size(kern_mat,3)),kern_mat,zeros(size(kern_mat,1),1,size(kern_mat,3)));
            kern_mat = cat(3,zeros(size(kern_mat,1),size(kern_mat,2),1),kern_mat,zeros(size(kern_mat,1),size(kern_mat,2),1));
            x_derivs = (kern_mat(:,2:end,:) - kern_mat(:,1:end-1,:));
            y_derivs = (kern_mat(:,:,2:end) - kern_mat(:,:,1:end-1));
            t_derivs = (kern_mat(2:end,:,:) - kern_mat(1:end-1,:,:));
            lpen.grad(i) = (mod.lambda_dX*sum(x_derivs(:).^2) + mod.lambda_dX*sum(y_derivs(:).^2) + mod.lambda_dT*sum(t_derivs(:).^2));
            lpen.grad(i) = lpen.grad(i)/nspks;
        elseif glmod.stim_params.spatial_dims == 1
            kern_mat = reshape(cur_pkern,[flen sdim]);
            % buffer with zeros
            kern_mat = cat(1,zeros(1,size(kern_mat,2)),kern_mat,zeros(1,size(kern_mat,2)));
            kern_mat = cat(2,zeros(size(kern_mat,1)),kern_mat,zeros(size(kern_mat,1),1));
            x_derivs = (kern_mat(:,2:end) - kern_mat(:,1:end-1));
            t_derivs = (kern_mat(2:end,:) - kern_mat(1:end-1,:));
            lpen.grad(i) = (mod.lambda_dX*sum(x_derivs(:).^2) + mod.lambda_dT*sum(t_derivs(:).^2));
            lpen.grad(i) = lpen.grad(i)/nspks;
        elseif glmod.stim_params.spatial_dims == 0
            kern_mat = reshape(cur_pkern,[flen sdim]);
            % buffer with zeros
            kern_mat = [0; kern_mat; 0];
            t_derivs = (kern_mat(2:end) - kern_mat(1:end-1));
            lpen.grad(i) = mod.lambda_dT*sum(t_derivs(:).^2);
            lpen.grad(i) = lpen.grad(i)/nspks;
        end
    end
    if mod.lambda_d2X > 0
        cur_kern = mod.k*mod.kscale; %reintroduce scale into k
        if glmod.stim_params.spatial_dims == 2
            kern_mat = reshape(cur_kern,[flen,sdim,sdim]);
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
        elseif glmod.stim_params.spatial_dims == 1
            kern_mat = reshape(cur_kern,[flen sdim]);
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
    end
    if mod.lambda_d2T > 0
        if glmod.stim_params.spatial_dims == 0
            cur_kern = mod.k*mod.kscale; %reintroduce scale into k
            kern_mat = reshape(cur_kern,[flen 1]);
            % buffer with zeros
            kern_mat = cat(1,zeros(1,size(kern_mat,2)),kern_mat,zeros(1,size(kern_mat,2)));
            
            lapl_oper = [1 -2 1];
            lapl = zeros(size(kern_mat,1),1);
            lapl = conv(kern_mat,lapl_oper,'same');
            lapl = lapl(2:end-1);
            lpen.laplT(i) = mod.lambda_d2T*sum(lapl(:).^2);
            lpen.laplT(i) = lpen.lapl(i)/nspks;
        elseif glmod.stim_params.spatial_dims == 1
            cur_kern = mod.k*mod.kscale;
            kern_mat = reshape(cur_kern,[flen sdim]);
            lapl = zeros(size(kern_mat));
            lapl(2:end-1,:) = kern_mat(1:end-2,:) + kern_mat(3:end,:) - 2*kern_mat(2:end-1,:);
            temp = lapl(2:end-1,2:end-1);
            lpen.laplT(i) = mod.lambda_d2T*sum(temp(:).^2);
            lpen.laplT(i) = lpen.laplT(i)/nspks;
        elseif glmod.stim_params.spatial_dims == 2
            cur_kern = mod.k*mod.kscale;
            kern_mat = reshape(cur_kern,[flen sdim sdim]);
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
        if glmod.stim_params.spatial_dims == 1
            cur_kern = mod.k*mod.kscale; %reintroduce scale into k
            kern_mat = reshape(cur_kern,[flen sdim]);
            % buffer with zeros
            kern_mat = cat(1,zeros(1,size(kern_mat,2)),kern_mat,zeros(1,size(kern_mat,2)));
            kern_mat = cat(2,zeros(size(kern_mat,1)),kern_mat,zeros(size(kern_mat,1),1));
            lapl_oper = [0 1 0;1 -4 1;0 1 0];
            lapl = conv2(kern_mat,lapl_oper,'same');
            lapl = lapl(2:end-1,2:end-1);
            lpen.laplXT(i) = mod.lambda_d2XT*sum(lapl(:).^2);
            lpen.laplXT(i) = lpen.laplXT(i)/nspks;
        else
            error('Unsupported')
        end
    end
    %L1 kernel penalty
    if mod.lambda_L1x > 0
        cur_kern = mod.k*mod.kscale; %reintroduce scale into k
        lpen.l1x(i) = sum(abs(cur_kern).*mod.lambda_L1x/nspks);
    else
        lpen.l1x(i) = 0;
    end
    
end


%% penalty terms for spiketrain moduls?
nll = -ll/nspks;
pnll = nll + sum(lpen.grad) + sum(lpen.l1x) + sum(lpen.l2x) + sum(lpen.lapl);

end