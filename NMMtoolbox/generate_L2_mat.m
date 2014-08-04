function full_mat = generate_L2_mat(params,klen)
%
% full_mat = generate_L2_mat(params,klen)
%
% INPUTS: 
% params: [rx1] struct array specifying r different regularization matrices to apply
% klen: total number of predictor dimensions.
% OUTPUTS:
% full_mat: cell array of L2 matrices
%%

krange = params.krange;
kdims = params.kdims;
deriv_order = params.deriv_order;
direction = params.direction;
boundaries = params.boundaries;
mix_prop = params.mix_prop;

%%

if deriv_order > 0
%for derivatives along first dimension
if direction == 1
    et = ones(kdims(1),1);
    if deriv_order == 1
        if boundaries == Inf
            et(end) = 0;
        end
        D1t = spdiags([et -et], [0 1], kdims(1), kdims(1))';
    elseif deriv_order == 2
        if boundaries == Inf
            et([1 end]) = 0;
        end
        D1t = spdiags([et -2*et et], [-1 0 1], kdims(1), kdims(1))';
    else
        error('');
    end
    
    if length(kdims) > 1
        Ix = speye(kdims(2));
        L2_mat = kron(Ix,D1t);
    else
        L2_mat = D1t;
    end
    
elseif direction == 2   %for derivative along second dimension
    
    ex = ones(kdims(2),1);
    if deriv_order == 1
        if boundaries == Inf
            ex(end) = 0;
        end
        D1x = spdiags([ex -ex], [0 1], kdims(2), kdims(2))';
    elseif deriv_order == 2
        if boundaries == Inf
            ex([1 end]) = 0;
        end
        D1x = spdiags([ex -2*ex ex], [-1 0 1], kdims(2), kdims(2))';
    else
        error('');
    end
    It = speye(kdims(1));
    L2_mat = kron(D1x,It);
    
elseif direction == 3   %for mixed derivative
    et = ones(kdims(1),1)*sqrt(mix_prop(1));
    if deriv_order == 1
        if boundaries(1) == Inf
            et(end) = 0;
        end
        D1t = spdiags([et -et], [0 1], kdims(1), kdims(1))';
    elseif deriv_order == 2
        if boundaries(1) == Inf
            et([1 end]) = 0;
        end
        D1t = spdiags([et -2*et et], [-1 0 1], kdims(1), kdims(1))';
    else
        error('');
    end
    
    ex = ones(kdims(2),1)*sqrt(mix_prop(2));
    if deriv_order == 1
        if boundaries(2) == Inf
            ex(end) = 0;
        end
        D1x = spdiags([ex -ex], [0 1], kdims(2), kdims(2))';
    elseif deriv_order == 2
        if boundaries(2) == Inf
            ex([1 end]) = 0;
        end
        D1x = spdiags([ex -2*ex ex], [-1 0 1], kdims(2), kdims(2))';
    else
        error('');
    end
    
    It = speye(kdims(1));
    Ix = speye(kdims(2));
    L2_mat = kron(Ix,D1t) + kron(D1x,It);
    
else
    error('');
end
else
    L2_mat = speye(diff(krange)+1);
end
%%
full_mat = sparse(klen,klen,0); %initialize sparse zero mat
full_mat(krange(1):krange(2),krange(1):krange(2)) = L2_mat;
