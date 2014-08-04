function reg_params = NIMcreate_reg_params(varargin)
%
% reg_params = NIMcreate_reg_params('lambda_type1',lambda1_val,'lambda_type2',lambda2_val',...)
%
% INPUTS: 
%     Lambda_type: String that specifies type of regularization
%     Lambda_val: Value of regularization hyperparameters
%     CAN SET ANY NUMBER OF PARAMETERS BY INPUTTING PAIRS OF THESE INPUTS
% OUTPUTS: 
%     reg_params: struct of regularization parameters

%% INITIALIZE REG_PARAMS WITH DEFAULT VALUES
reg_params.lambda_NLd2 = 0; %second derivative of tent basis coefs
reg_params.lambda_dX = 0; %first spatial deriv
reg_params.lambda_dT = 0; %first temporal deriv
reg_params.lambda_d2XT = 0; %spatiotemporal laplacian
reg_params.lambda_d2X = 0; %2nd spatial deriv
reg_params.lambda_d2T = 0; %2nd temporal deriv
reg_params.lambda_L2 = 0; %L2 on filter coefs
reg_params.lambda_L1 = 0; %L1 on filter coefs

%boundary conditions
reg_params.spatial_boundaries = 'zero'; %default is assume zeros-boundaries in spatial dims
reg_params.temporal_boundaries = 'free'; %assume 'free' boundaries for temporal dims

%% PARSE INPUTS AND ADD TO REG_PARAMS
if mod(length(varargin),2) ~= 0
    error('Inputs must be in the form: "lambda_name",lambda_val"');
end
n_inputs = length(varargin)/2;
for ii = 1:n_inputs
   input_name = varargin{(ii-1)*2+1};
   input_val = varargin{(ii-1)*2 + 2};
   if ~isfield(reg_params,input_name)
       error('Invalid regularization type');
   else
      reg_params = setfield(reg_params,input_name,input_val); 
   end
end