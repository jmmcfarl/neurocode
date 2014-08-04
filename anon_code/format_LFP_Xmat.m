function [X_LFP,L2_mat,stim_params] = format_LFP_Xmat(ampgram,phasegram,mix_prop,dt)
%
% [X_LFP,L2_mat,stim_params] = format_LFP_Xmat(ampgram,phasegram,<mix_prop>,<dt>)
% 
% INPUTS: 
%   ampgram: amplitude at each scale and channel must be in form: [NT
% n_freqs n_chs] 
%   phasegram: phase at each scale and channel, must be same
% dims as ampgram 
%   <mix_prop>: if using mixed-derivative smoothing (assumed
% if there are more than 1 freqs and chs used), this gives the relative
% strength of smoothing along the freq and ch dims [p_freq p_ch]
%   <dt>: dt for use in stim_params struct
%
% OUTPUTS:
%     X_LFP: [NT 2*n_freqs*n_chs] predictor mat [a*cos(phi) a*sin(phi)]
%     L2_mat: [2*n_freqs*n_chs 2*n_freqs*n_chs] L2 matrix
%     stim_params: param struct for use with NMMfit functions
%%

if nargin < 3 || isempty(mix_prop)
    mix_prop = [1 1];
end
if nargin < 4 || isempty(dt)
    dt = 1;
end

[NT,n_freqs,n_chs] = size(ampgram);
tot_kLen = n_freqs*n_chs;

%z-score normalize amplitudes at each freq/ch
ampgram = bsxfun(@minus,ampgram,mean(ampgram));
ampgram = bsxfun(@rdivide,ampgram,std(ampgram));

X_LFP = [reshape(ampgram,NT,tot_kLen).*reshape(cos(phasegram),NT,tot_kLen) ...
    reshape(ampgram,NT,tot_kLen).*reshape(sin(phasegram),NT,tot_kLen)];

%create stim params struct
up_samp_fac = 1; tent_spacing = 1;
stim_params = struct('stim_dims',[n_freqs n_chs],'dt',dt,'up_fac',up_samp_fac,...
    'tent_spacing',tent_spacing);

%% CREATE L2 MAT
if n_freqs > 1 && n_chs > 1 % IF USING MULTIPLE FREQS AND CHS
    bound_types = [Inf Inf]; %free boundaries in both freq and ch dimensions
    
    %create L2 mat for cosine term
    L2_params = create_L2_params([],[1 tot_kLen],[n_freqs n_chs],2,3,bound_types,mix_prop);
    L2_mat = generate_L2_mat(L2_params(1),2*tot_kLen);
    
    %add L2 mat for sine term
    L2_params = create_L2_params(L2_params,[tot_kLen+1 2*tot_kLen],[n_freqs n_chs],2,3,bound_types,mix_prop);
    L2_mat = L2_mat + generate_L2_mat(L2_params(2),2*tot_kLen);
    
elseif n_freqs > 1 && n_chs == 1 % IF USING MULTIPLE FREQS AND ONE CH
    bound_types = [Inf];
    
    %create L2 mat for cosine term
    L2_params = create_L2_params([],[1 tot_kLen],[n_freqs],2,1,bound_types);
    L2_mat = generate_L2_mat(L2_params(1),2*tot_kLen);
    
    %add L2 mat for sine term
    L2_params = create_L2_params(L2_params,[tot_kLen+1 2*tot_kLen],[n_freqs],2,1,bound_types);
    L2_mat = L2_mat + generate_L2_mat(L2_params(2),2*tot_kLen);
  
elseif n_freqs == 1 && n_chs > 1 % IF USING MULTIPLE CHS AND ONE FREQ
    bound_types = [Inf];
    
    %create L2 mat for cosine term
    L2_params = create_L2_params([],[1 tot_kLen],[n_chs],2,1,bound_types);
    L2_mat = generate_L2_mat(L2_params(1),2*tot_kLen);
    
    %add L2 mat for sine term
    L2_params = create_L2_params(L2_params,[tot_kLen+1 2*tot_kLen],[n_chs],2,1,bound_types);
    L2_mat = L2_mat + generate_L2_mat(L2_params(2),2*tot_kLen);

elseif n_freqs == 1 && n_chs == 1 % IF USING ONE CH AND ONE FREQ DONT NEED SMOOTHING
    L2_mat = [];
else
    error('Check format of inputs');
end