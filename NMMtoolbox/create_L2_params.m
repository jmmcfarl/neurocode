function L2_params = create_L2_params(prev_struct,krange,kdims,deriv_order,direction,boundaries,mix_prop)
%
% L2_params = create_L2_params(prev_struct,krange,<kdims>,<deriv_order>,<direction>,<boundaries>,<mix_prop>)
%
%INPUTS: 
% prev_struct: struct array of regularization parameters to add to (null struct [] if creating a new param struct).
% krange: 2x1 vector of parameter indices specifying the first and last parameters to apply the current regularization to
% kdims: specifies size of the predictor matrix
% deriv_order: order of derivative to regularize
% direction: dimension along which to regularize [1 = first dim, 2 = second dim, 3 = mixed]
% boundaries: boundary conditions on derivative reg. (0 for smoothing towards 0 and Inf for free boundaries)
% mix_prop: for mixed derivatives, this determines the relative weight of regularization along each dimension
%OUTPUTS:
% L2_params: structure storing specified regularization params.
% 
% EXAMPLE:
% % create param struct for 2nd derivative regularization. Spatiotemporal smoothing with flen time lags and nPIx spatial dims. 
% % Use mixing proportion [1 0.5] so twice as much smoothness along the temporal dim
% L2_params = create_L2_params([],1:flen*nPix,[flen nPix],2,3,[1 0.5]);
% % Add another element to the param struct array specifying regularization
% % of a 1-dimensional filter of length tlen. Use first derivative.
% L2_params = create_L2_params(L2_params,(flen*nPix+1):(flen*nPix+tlen),[tlen 1],1);

%%
if length(krange) > 2
    error('Input Krange as 2-element vector [kbeg kend]');
end
if nargin < 3 || isempty(kdims)
    kdims = diff(krange)+1;
end
if nargin < 4 || isempty(deriv_order)
    deriv_order = 2; 
end
if nargin < 5 || isempty(direction)
    direction = 1; %[1= first dim, 2 = second dim, 3 = mixed]
end
if nargin < 6 || isempty(boundaries)
    if length(kdims) == 1
        boundaries = [0];
    else
        boundaries = [0 0];
    end
end
if nargin < 7 || isempty(mix_prop)
    mix_prop = ones(length(kdims));
end

cur_L2_params = struct('krange',krange,'kdims',kdims,'deriv_order',deriv_order,'direction',direction,'boundaries',boundaries,'mix_prop',mix_prop);
L2_params = cat(1,prev_struct,cur_L2_params);
