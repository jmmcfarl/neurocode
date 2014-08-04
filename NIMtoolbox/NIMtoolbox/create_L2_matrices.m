function L2_mats = create_L2_matrices(nim)
%
% L2_matrices = create_L2_matrices(nim)
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
nmods = length(nim.mods);
all_d2T = zeros(nmods,1);
all_dT = zeros(nmods,1);
all_dX = zeros(nmods,1);
all_d2X = zeros(nmods,1);
all_d2XT = zeros(nmods,1);
t_bound = zeros(nmods,1); %[1 for free-boundaries, 0 for 0-boundaries]
x_bound = zeros(nmods,1);
for ii = 1:nmods
    all_d2T(ii) = nim.mods(ii).reg_params.lambda_d2T;
    all_dT(ii) = nim.mods(ii).reg_params.lambda_dT;
    all_d2X(ii) = nim.mods(ii).reg_params.lambda_d2X;
    all_dX(ii) = nim.mods(ii).reg_params.lambda_dX;
    all_d2XT(ii) = nim.mods(ii).reg_params.lambda_d2XT;
    if strcmp(nim.mods(ii).reg_params.spatial_boundaries,'free')
        x_bound(ii) = 1;
    end
    if strcmp(nim.mods(ii).reg_params.temporal_boundaries,'free')
        t_bound(ii) = 1;
    end
end
x_bound = unique(x_bound); t_bound = unique(t_bound);
if length(x_bound) > 1 || length(t_bound) > 1
    error('Regularization boundary conditions must be the same for all subunits');
end

%% CREATE L2 REG MATRICES
L2_mats = [];
nLags = stim_params.stim_dims(1);
nPix = squeeze(stim_params.stim_dims(2:3));

%for 0-dimensional stimuli
if nPix == 1
    
    %check for temporal laplacian
    if max(all_d2T) > 0
        et = ones(nLags,1);
        if t_bound == 1 %if 'free boundary'
            et([1 end]) = 0;
        end
        L2_mats.L2_d2T = spdiags([et -2*et et], [-1 0 1], nLags, nLags)';
    end
    
    %check for temporal laplacian
    if max(all_dT) > 0
        et = ones(nLags,1);
        if t_bound == 1 %if free boundary
            et([end]) = 0;
        end
        L2_mats.L2_dT = spdiags([et -et], [0 1], nLags, nLags)';
    end
    
    %for 1-dimensional stimuli
elseif nPix(2) == 1
    
    %check for temporal derivative
    if max(all_dT) > 0
        et = ones(nLags,1);
        if t_bound == 1 %if free boundary
            et([end]) = 0;
        end
        D1t = spdiags([et -et], [0 1], nLags, nLags)';
        
        Ix = speye(nPix(1));
        L2_mats.L2_dX = kron(Ix,D1t);
    end
    
    %check for spatial derivative
    if max(all_dX) > 0
        ex = ones(nPix(1),1);
        if x_bound == 1 %if free boundary
            ex([end]) = 0;
        end
        D1x = spdiags([ex -ex], [0 1], nPix(1), nPix(1))';
        
        It = speye(nLags);
        L2_mats.L2_dX = kron(D1x,It);
    end
    
    %check for temporal laplacian
    if max(all_d2T) > 0
        et = ones(nLags,1);
        if t_bound == 1 %if free boundary
            et([1 end]) = 0;
        end
        D1t = spdiags([et -2*et et], [-1 0 1], nLags, nLags)';
        
        Ix = speye(nPix(1));
        L2_mats.L2_d2T = kron(Ix,D1t);
    end
    
    %check for temporal laplacian
    if max(all_d2X) > 0
        ex = ones(nPix(1),1);
        if x_bound == 1 %if free boundary
            ex([1 end]) = 0;
        end
        D1x = spdiags([ex -2*ex ex], [-1 0 1], nPix(1), nPix(1))';
        
        It = speye(nLags);
        L2_mats.L2_d2X = kron(D1x,It);
    end
    
    %check for spatio-temporal laplacian
    if max(all_d2XT) > 0
        et = ones(nLags,1);
        if t_bound == 1 %if free boundary
            et([1 end]) = 0;
        end
        ex = ones(nPix(1),1);
        if x_bound == 1 %if free boundary
            ex([1 end]) = 0;
        end
        D1t = spdiags([et -2*et et], [-1 0 1], nLags, nLags)';
        D1x = spdiags([ex -2*ex ex], [-1 0 1], nPix(1), nPix(1))';
        
        It = speye(nLags);
        Ix = speye(nPix(1));
        L2_mats.L2_d2XT = kron(Ix,D1t) + kron(D1x,It);
    end
else %for 2d stim
    %check for temporal derivative
    if max(all_dT) > 0
        et = ones(nLags,1);
        if t_bound == 1 %if free boundary
            et([end]) = 0;
        end
        D1t = spdiags([et -et], [0 1], nLags, nLags)';
        
        Ix = speye(nPix(1));
        Iy = speye(nPix(2));
        L2_mats.L2_dX = kron(Iy, kron(Ix, D1t));
    end
    
    %check for temporal laplacian
    if max(all_d2T) > 0
        et = ones(nLags,1);
        if t_bound == 1 %if free boundary
            et([1 end]) = 0;
        end
        D1t = spdiags([et -2*et et], [-1 0 1], nLags, nLags)';
        
        Ix = speye(nPix(1));
        Iy = speye(nPix(2));
        L2_mats.L2_d2T = kron(Iy, kron(Ix, D1t));
    end
    
    %check for spatial derivative
    if max(all_dX) > 0
        ex = ones(nPix(1),1);
        ey = ones(nPix(2),1);
        if x_bound == 1 %if free boundary
            ex([end]) = 0;
            ey([end]) = 0;
        end
        D1x = spdiags([ex -ex], [0 1], nPix(1), nPix(1))';
        D1y = spdiags([ey -ey], [0 1], nPix(2), nPix(2))';
        
        It = speye(nLags); Ix = speye(nPix(1)); Iy = speye(nPix(2));
        L2_mats.L2_dX = kron(Iy, kron(D1x, It)) + kron(kron(D1y,Ix),It);
    end

    %check for temporal laplacian
    if max(all_d2X) > 0
        ex = ones(nPix(1),1);
        ey = ones(nPix(2),1);
        if x_bound == 1 %if free boundary
            ex([1 end]) = 0;
            ey([1 end]) = 0;
        end
        D1x = spdiags([ex -2*ex ex], [-1 0 1], nPix(1), nPix(1))';
        D1y = spdiags([ey -2*ey ey], [-1 0 1], nPix(2), nPix(2))';
        
        It = speye(nLags); Ix = speye(nPix(1)); Iy = speye(nPix(2));
        L2_mats.L2_d2X = kron(Iy, kron(D1x, It)) + kron(kron(D1y,Ix),It);
    end
    
    %check for spatio-temporal laplacian
    if max(all_d2XT) > 0
        et = ones(nLags,1);
        ex = ones(nPix(1),1);
        ey = ones(nPix(2),1);
        if t_bound == 1 %if free boundary
            et([1 end]) = 0;
        end
        if x_bound == 1 %if free boundary
            ex([1 end]) = 0;
            ey([1 end]) = 0;
        end
        D1t = spdiags([et -2*et et], [-1 0 1], nLags, nLags)';
        D1x = spdiags([ex -2*ex ex], [-1 0 1], nPix(1), nPix(1))';
        D1y = spdiags([ey -2*ey ey], [-1 0 1], nPix(2), nPix(2))';
        
        It = speye(nLags); Ix = speye(nPix(1)); Iy = speye(nPix(2));
        L2_mats.L2_d2XT = kron(Iy, kron(Ix, D1t)) + kron(Iy, kron(D1x, It))...
            + kron(kron(D1y,Ix),It);
    end
end