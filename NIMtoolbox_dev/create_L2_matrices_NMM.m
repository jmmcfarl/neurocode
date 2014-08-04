function L2_mats = create_L2_matrices_NMM(nim)
%
% L2_matrices = create_L2_matrices_NMM(nim)
%
% Creates a set of matrices A specifying different L2 penalties as ||Ak||^2
% INPUTS:
%     nim: Input model to which we're adding regularization
% OUTPUTS:
%     L2_mats: struct containing all the needed regularization matrices
%
% The method of computing sparse differencing matrices used here is adapted from
% Bryan C. Smith's and Andrew V. Knyazev's function "laplacian", available
% here: http://www.mathworks.com/matlabcentral/fileexchange/27279-laplacian-in-1d-2d-or-3d


%% EXTRACT ALL REG PARAMS

stim_params = nim.stim_params;
nXmats = length(stim_params);
Xtargs = [nim.mods(:).Xtarget];
nmods = length(nim.mods);
all_d2T = [];
all_dT = [];
all_dX = [];
all_d2X = [];
all_d2XT = [];
for ii = 1:nmods
    if nim.mods(ii).reg_params.lambda_d2T > 0
        all_d2T = cat(1,all_d2T,Xtargs(ii));
    end
    if nim.mods(ii).reg_params.lambda_dT > 0
        all_dT = cat(1,all_dT,Xtargs(ii));
    end
    if nim.mods(ii).reg_params.lambda_d2X > 0
        all_d2X = cat(1,all_d2X,Xtargs(ii));
    end
    if nim.mods(ii).reg_params.lambda_dX > 0
        all_dX = cat(1,all_dX,Xtargs(ii));
    end
    if nim.mods(ii).reg_params.lambda_d2XT > 0
        all_d2XT = cat(1,all_d2XT,Xtargs(ii));
    end
end
all_d2T = unique(all_d2T);
all_dT = unique(all_dT);
all_d2X = unique(all_d2X);
all_dX = unique(all_dX);
all_d2XT = unique(all_d2XT);

bound_conds = zeros(nXmats,3);
mixing_props = zeros(nXmats,3);
for ii = 1:nXmats
    cur_mod = find(Xtargs == ii,1);
    if ~isempty(cur_mod)
       bound_conds(ii,:) = nim.mods(cur_mod).reg_params.boundary_conds; 
       mixing_props(ii,:) = nim.mods(cur_mod).reg_params.mixing_prop;
    end
end

for ii = 1:nXmats
    
    nLags = stim_params(ii).stim_dims(1);
    nPix = squeeze(stim_params(ii).stim_dims(2:3));
    
    %for 0-dimensional stimuli
    if nPix == 1
        
        %check for temporal laplacian
        if ismember(ii,all_d2T)
            et = ones(nLags,1);
            if isinf(bound_conds(ii,1)) %if 'free boundary'
                et([1 end]) = 0;
            end
            L2_mats.L2_d2T{ii} = spdiags([et -2*et et], [-1 0 1], nLags, nLags)';
        else
           L2_mats.L2_d2T{ii} = []; 
        end
        
        %check for temporal laplacian
        if ismember(ii,all_dT)
            et = ones(nLags,1);
            if isinf(bound_conds(ii,1)) %if free boundary
                et([end]) = 0;
            end
            L2_mats.L2_dT{ii} = spdiags([et -et], [0 1], nLags, nLags)';
        else
           L2_mats.L2_dT{ii} = []; 
        end
        
        %for 1-dimensional stimuli
    elseif nPix(2) == 1
        
        %check for temporal derivative
        if ismember(ii,all_dT)
            et = ones(nLags,1);
            if isinf(bound_conds(ii,1)) %if free boundary
                et([end]) = 0;
            end
            D1t = spdiags([et -et], [0 1], nLags, nLags)';
            
            Ix = speye(nPix(1));
            L2_mats.L2_dT{ii} = kron(Ix,D1t);
        else
           L2_mats.L2_dT{ii} = []; 
        end
        
        %check for spatial derivative
        if ismember(ii,all_dX)
            ex = ones(nPix(1),1);
            if isinf(bound_conds(ii,2)) %if free boundary
                ex([end]) = 0;
            end
            D1x = spdiags([ex -ex], [0 1], nPix(1), nPix(1))';
            
            It = speye(nLags);
            L2_mats.L2_dX{ii} = kron(D1x,It);
        else
           L2_mats.L2_dX{ii} = []; 
        end
        
        %check for temporal laplacian
        if ismember(ii,all_d2T)
            et = ones(nLags,1);
            if isinf(bound_conds(ii,1)) %if free boundary
                et([1 end]) = 0;
            end
            D1t = spdiags([et -2*et et], [-1 0 1], nLags, nLags)';
            
            Ix = speye(nPix(1));
            L2_mats.L2_d2T{ii} = kron(Ix,D1t);
        else
           L2_mats.L2_d2T{ii} = []; 
        end
        
        %check for temporal laplacian
        if ismember(ii,all_d2X)
            ex = ones(nPix(1),1);
            if isinf(bound_conds(ii,2)) %if free boundary
                ex([1 end]) = 0;
            end
            D1x = spdiags([ex -2*ex ex], [-1 0 1], nPix(1), nPix(1))';
            
            It = speye(nLags);
            L2_mats.L2_d2X{ii} = kron(D1x,It);
        else
           L2_mats.L2_d2X{ii} = []; 
        end
        
        %check for spatio-temporal laplacian
        if ismember(ii,all_d2XT)
            et = ones(nLags,1)*sqrt(mixing_props(ii,1));
            if isinf(bound_conds(ii,1)) %if free boundary
                et([1 end]) = 0;
            end
            ex = ones(nPix(1),1)*sqrt(mixing_props(ii,2));
            if isinf(bound_conds(ii,2)) %if free boundary
                ex([1 end]) = 0;
            end
            D1t = spdiags([et -2*et et], [-1 0 1], nLags, nLags)';
            D1x = spdiags([ex -2*ex ex], [-1 0 1], nPix(1), nPix(1))';
            
            It = speye(nLags);
            Ix = speye(nPix(1));
            L2_mats.L2_d2XT{ii} = kron(Ix,D1t) + kron(D1x,It);
        else
           L2_mats.L2_d2XT{ii} = []; 
        end
    else %for 2d stim
        %check for temporal derivative
        if ismember(ii,all_dT)
            et = ones(nLags,1);
            if isinf(bound_conds(ii,1)) %if free boundary
                et([end]) = 0;
            end
            D1t = spdiags([et -et], [0 1], nLags, nLags)';
            
            Ix = speye(nPix(1));
            Iy = speye(nPix(2));
            L2_mats.L2_dT{ii} = kron(Iy, kron(Ix, D1t));
        else
            L2_mats.L2_dT{ii} = [];
        end
        
        %check for temporal laplacian
        if ismember(ii,all_d2T)
            et = ones(nLags,1);
            if isinf(bound_conds(ii,1)) %if free boundary
                et([1 end]) = 0;
            end
            D1t = spdiags([et -2*et et], [-1 0 1], nLags, nLags)';
            
            Ix = speye(nPix(1));
            Iy = speye(nPix(2));
            L2_mats.L2_d2T{ii} = kron(Iy, kron(Ix, D1t));
        else
            L2_mats.L2_d2T{ii} = [];
        end
        
        %check for spatial derivative
            if ismember(ii,all_dX)
            ex = ones(nPix(1),1);
            ey = ones(nPix(2),1);
            if isinf(bound_conds(ii,2)) %if free boundary
                ex([end]) = 0;
                ey([end]) = 0;
            end
            D1x = spdiags([ex -ex], [0 1], nPix(1), nPix(1))';
            D1y = spdiags([ey -ey], [0 1], nPix(2), nPix(2))';
            
            It = speye(nLags); Ix = speye(nPix(1)); Iy = speye(nPix(2));
            L2_mats.L2_dX{ii} = kron(Iy, kron(D1x, It)) + kron(kron(D1y,Ix),It);
            else
               L2_mats.L2_dX{ii} = []; 
            end
        
        %check for temporal laplacian
        if ismember(ii,all_d2X)
            ex = ones(nPix(1),1)*sqrt(mixing_props(ii,2));
            ey = ones(nPix(2),1)*sqrt(mixing_props(ii,3));
            if isinf(bound_conds(ii,2)) %if free boundary
                ex([1 end]) = 0;
            end
            if isinf(bound_conds(ii,3))
                ey([1 end]) = 0;
            end
            D1x = spdiags([ex -2*ex ex], [-1 0 1], nPix(1), nPix(1))';
            D1y = spdiags([ey -2*ey ey], [-1 0 1], nPix(2), nPix(2))';
            
            It = speye(nLags); Ix = speye(nPix(1)); Iy = speye(nPix(2));
            L2_mats.L2_d2X{ii} = kron(Iy, kron(D1x, It)) + kron(kron(D1y,Ix),It);
        else
           L2_mats.L2_d2X{ii} = []; 
        end
        
        %check for spatio-temporal laplacian
        if ismember(ii,all_d2XT)
            et = ones(nLags,1)*sqrt(mixing_props(ii,1));
            ex = ones(nPix(1),1)*sqrt(mixing_props(ii,2));
            ey = ones(nPix(2),1)*sqrt(mixing_props(ii,3));
            if isinf(bound_conds(ii,1)) %if free boundary
                et([1 end]) = 0;
            end
            if isinf(bound_conds(ii,2)) %if free boundary
                ex([1 end]) = 0;
            end
            if isinf(bound_conds(ii,3))
                ey([1 end]) = 0;
            end
            D1t = spdiags([et -2*et et], [-1 0 1], nLags, nLags)';
            D1x = spdiags([ex -2*ex ex], [-1 0 1], nPix(1), nPix(1))';
            D1y = spdiags([ey -2*ey ey], [-1 0 1], nPix(2), nPix(2))';
            
            It = speye(nLags); Ix = speye(nPix(1)); Iy = speye(nPix(2));
            L2_mats.L2_d2XT{ii} = kron(Iy, kron(Ix, D1t)) + kron(Iy, kron(D1x, It))...
                + kron(kron(D1y,Ix),It);
        else
            L2_mats.L2_d2XT{ii} = [];
        end
    end
end