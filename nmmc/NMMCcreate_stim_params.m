function stim_params = NMMCcreate_stim_params( stim_dims, stim_dt, up_samp_fac, tent_spacing )
%
% stim_params = NMMCcreate_stim_params(stim_dims,<stim_dt>,<up_samp_fac>,<tent_spacing> )
%
% Creates a struct containing stimulus parameters 
% INPUTS:
%     stim_dims: dimensionality of the (time-embedded) stimulus, in the
%         form: [nLags nXPix nYPix]. For 1 spatial dimension use only nXPix
%     <stim_dt>: time resolution (in ms) of Xmatrix (used only for plotting)
%     <up_samp_fac>: optional up-sampling of the stimulus from its raw form
%     <tent_spacing>: optional spacing of tent-basis functions when using a tent-basis 
%         representaiton of the stimulus (allows for the stimulus filters to be 
%         represented at a lower time resolution than other model
%         components). 
% OUTPUTS:
%     stim_params: struct of stimulus parameters
    
%%
%defaults
if nargin < 2 || isempty(stim_dt)
	stim_dt = 1;
end
if nargin < 3 || isempty(up_samp_fac)
	up_samp_fac = 1;
end
if (nargin < 4) || isempty(tent_spacing)
	tent_spacing = [];
end

%make sure stim_dims input has form [nLags nXPix nYPix] and concatenate 1's
%if necessary
if length(stim_dims) == 1
    stim_dims = [stim_dims 1 1];
elseif length(stim_dims) == 2
    stim_dims = [stim_dims 1];
elseif length(stim_dims) > 3
    error('Stimulus dimensions must be in format: [nLags XPix YPix]');
end

dt = stim_dt/up_samp_fac; %model fitting dt

% Set PSC parameters
%PSC_dt = 1;
%PSC_NT = 0;
%if (nargin > 4)
%	PSC_NT = PSC_params(1);
%	if length(PSC_params) > 1
%		PSC_dt = PSC_params(2);
%	end
%end

stim_params = struct('stim_dims',stim_dims,'dt',dt,'up_fac',up_samp_fac,...
    'tent_spacing',tent_spacing );
%    'tent_spacing',tent_spacing, 'PSC_NT',0, 'PSC_dt', 1, 'SPconv',[] );


