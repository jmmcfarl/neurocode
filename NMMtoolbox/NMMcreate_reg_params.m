function reg_params = NMMcreate_reg_params(varargin)
%
% reg_params = NMMcreate_reg_params('lambda_type1',lambda1_val,'lambda_type2',lambda2_val',...)
%
% INPUTS: 
%     Lambda_type: String that specifies type of regularization
%     Lambda_val: Value of regularization hyperparameters
%     CAN SET ANY NUMBER OF PARAMETERS BY INPUTTING PAIRS OF THESE INPUTS
% OUTPUTS: 
%     reg_params: struct of regularization parameters
%
% Different regularization types:
%    lambda_NLd2: second derivative of tent basis coefs
%    lambda_dX: first spatial deriv
%    lambda_dT: first temporal deriv
%    lambda_d2XT: spatiotemporal laplacian
%    lambda_d2X: 2nd spatial deriv
%    lambda_d2T: 2nd temporal deriv
%    lambda_L2: L2 on filter coefs
%    lambda_L1: L1 on filter coefs
%    lambda_custom: for custom L2 reg matrices
%
%NOTE: lambda_vals can either be input in column-vector format, or as
%scalars (in which case it is assumed that all subunits have the same
%regularization strength). If using column-vector format, you need to do
%this for all input values.

%% INITIALIZE REG_PARAMS WITH DEFAULT VALUES
reg_params.lambda_NLd2 = 0; %second derivative of tent basis coefs
reg_params.lambda_dX = 0; %first spatial deriv
reg_params.lambda_dT = 0; %first temporal deriv
reg_params.lambda_d2XT = 0; %spatiotemporal laplacian
reg_params.lambda_d2X = 0; %2nd spatial deriv
reg_params.lambda_d2T = 0; %2nd temporal deriv
reg_params.lambda_L2 = 0; %L2 on filter coefs
reg_params.lambda_L1 = 0; %L1 on filter coefs
reg_params.lambda_custom = 0; % whatever

%boundary conditions (Inf is free boundary and 0 for smoothing towards 0 at
%boundaries
reg_params.boundary_conds = [Inf 0 0]; %default to free boundary on first dim (usually time), and zero-boundaries in other dims (usually space)
reg_params.mixing_prop = [1 1 1]; %default mixing ratio


%% PARSE INPUTS AND ADD TO REG_PARAMS
if mod(length(varargin),2) ~= 0
    error('Inputs must be in the form: "lambda_name",lambda_val"');
end
n_inputs = length(varargin)/2;

if n_inputs > 0
    nmods = length(varargin{2});
    reg_params = repmat(reg_params,nmods,1);
end

for ii = 1:n_inputs
    input_name = varargin{(ii-1)*2+1};
    input_val = varargin{(ii-1)*2 + 2};
    if size(input_val,1) ~= nmods
        error('All lambda vectors must have same length');
    end
    if ~isfield(reg_params,input_name)
        error('Invalid regularization type');
    else
        cur_d = size(input_val,2);
        for ii = 1:nmods
            reg_params(ii).(input_name)(1:cur_d) = input_val(ii,:);
        end
    end
end