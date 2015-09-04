function Tmat = create_Tikhonov_matrix(stim_params, reg_type)
%
% Tmat = create_Tikhonov_matrix(stim_params, direction, order)
%
% Creates a matrix specifying a form of L2-regularization of the form
% ||T*k||^2. Currently only supports second derivative/Laplacian operations
% 
% INPUTS:
%     stim_params: parameter struct associated with the target stimulus.
%             must contain .dims field specifying the number of stimulus elements along each dimension
%             <.boundary_conds> specifies boundary conditions: Inf is free boundary, 0 is tied to 0, and -1 is periodic
%       direction: direction of the derivative relative to the stimulus dimensions. e.g. 1 is along the first dim, 2 is along the second, [1 2] is a laplacian
%
% OUTPUTS:
%     Tmat: sparse matrix specifying the desired Tikhonov operation
%
% The method of computing sparse differencing matrices used here is adapted from
% Bryan C. Smith's and Andrew V. Knyazev's function "laplacian", available
% here: http://www.mathworks.com/matlabcentral/fileexchange/27279-laplacian-in-1d-2d-or-3d


nLags = stim_params.dims(1); %first dimension is assumed to represent time
nPix = squeeze(stim_params.dims(2:3)); %additional dimensions are treated as spatial
allowed_reg_types = {'d2xt','d2x','d2t'};
assert(isstr(reg_type) && ismember(reg_type,allowed_reg_types),'not an allowed regularization type');

if ismember(reg_type,{'d2xt','d2t'}) %if there is a temporal component
    et = ones(nLags,1);
    if isinf(stim_params.boundary_conds(1)) %if temporal dim has free boundary
        et([1 end]) = 0;
    end
end
if ismember(reg_type,{'d2xt','d2x'}) %if there is a spatial component
    ex = ones(nPix(1),1); 
    if isinf(stim_params.boundary_conds(2)) %if first spatial dim has free boundary
        ex([1 end]) = 0;
    end
    ey = ones(nPix(2),1);
    if isinf(stim_params.boundary_conds(3)); %if second spatial dim has free boundary
        ey([1 end]) = 0;
    end
end

if nPix == 1 %for 0-spatial-dimensional stimuli can only do temporal
    assert(ismember(reg_type,{'d2t'}),'can only do temporal reg for stimuli without spatial dims');
    Tmat = spdiags([et -2*et et],[-1 0 1], nLags, nLags);
    if stim_params.boundary_conds(1) == -1 %if periodic boundary cond
        Tmat(end,1) = 1; Tmat(1,end) = 1;
    end
    
elseif nPix(2) == 1 %for 1-spatial dimensional stimuli
    if strcmp(reg_type,'d2t') %if temporal deriv
        D1t = spdiags([et -2*et et],[-1 0 1],nLags,nLags)';
        if stim_params.boundary_conds(1) == -1 %if periodic boundary cond
            D1t(end,1) = 1; D1t(1,end) = 1;
        end
        Ix = speye(nPix(1));
        Tmat = kron(Ix,D1t);
    elseif strcmp(reg_type,'d2x') %if spatial deriv
        It = speye(nLags);
        D1x = spdiags([ex -2*ex ex], [-1 0 1], nPix(1), nPix(1))';
        if stim_params.boundary_conds(2) == -1 %if periodic boundary cond
           D1x(end,1) = 1; D1x(1,end) = 1; 
        end
        Tmat = kron(D1x,It);
    elseif strcmp(reg_type,'d2xt') %if spatiotemporal laplacian
        D1t = spdiags([et -2*et et], [-1 0 1], nLags, nLags)';
        if stim_params.boundary_conds(1) == -1 %if periodic boundary cond
           D1t(end,1) = 1; D1t(1,end) = 1; 
        end
        D1x = spdiags([ex -2*ex ex], [-1 0 1], nPix(1), nPix(1))';
        if stim_params.boundary_conds(2) == -1 %if periodic boundary cond
           D1x(end,1) = 1; D1x(1,end) = 1; 
        end
        It = speye(nLags);
        Ix = speye(nPix(1));
        Tmat = kron(Ix,D1t) + kron(D1x,It);
    end
else %for stimuli with 2-spatial dimensions
    if strcmp(reg_type,'d2t') %temporal deriv
        D1t = spdiags([et -2*et et], [-1 0 1], nLags, nLags)';
        if stim_params.boundary_conds(1) == -1 %if periodic boundary cond
           D1t(end,1) = 1; D1t(1,end) = 1; 
        end
        Ix = speye(nPix(1));
        Iy = speye(nPix(2));
        Tmat = kron(Iy, kron(Ix, D1t));
    elseif strcmp(reg_type,'d2x') %spatial laplacian
        It = speye(nLags);
        Ix = speye(nPix(1));
        Iy = speye(nPix(2));
        D1x = spdiags([ex -2*ex ex], [-1 0 1], nPix(1), nPix(1))';
        if stim_params.boundary_conds(2) == -1 %if periodic boundary cond
           D1x(end,1) = 1; D1x(1,end) = 1; 
        end
        D1y = spdiags([ey -2*ey ey], [-1 0 1], nPix(2), nPix(2))';
        if stim_params.boundary_conds(3) == -1 %if periodic boundary cond
           D1y(end,1) = 1; D1y(1,end) = 1; 
        end
        Tmat = kron(Iy,kron(D1x,It)) + kron(D1y, kron(Ix,It));
    else
        error('unsupported regularization type');
    end
end

Tmat = Tmat';
