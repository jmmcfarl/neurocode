function nim = NMMinitialize_upstreamNLs( nim, Xstims, target_mods, lambda_smooth, NLmon, edge_p, n_bfs, space_type )
%
% Usage: nim = NMMinitialize_upstreamNLs( nim, Xstims, <target_mods>, <lambda_smooth>, <NLmon>, <edge_p>, <n_bfs>, <space_type>)
%
% Initializes the specified model subunits to have nonparametric
% (tent-basis) upstream NLs
%
% INPUTS:
%     nim: input model structure
%     Xstim: time-embedded stim matrix
%     <target_mods>: vector specifying which subunits will be made to have nonparametric upstream NLs
%     <lambda_smooth>: specifies strength of smoothness regularization for the tent-basis coefs
%     <NLmon>: Set to +1 to constrain NL coefs to be monotonic increasing and -1 to make mono. 
%        decreasing. 0 means no constraint. Default here is +1 (monotonic increasing)
%     <edge_p>: Determines the locations of the outermost tent-bases relative to the underlying generating distribution
%     <n_bfs>: Number of tent-basis functions to use
%     <space_type>: Use either 'equispace' for uniform bin spacing, or 'equipop' for 'equipopulated bins'
% 
% OUTPUTS:
%     nim: output model

%%
Nmods = length(nim.mods);

% Process Xstims (in case multiple Xstims)
if ~iscell(Xstims)
	tmp = Xstims;
	clear Xstims
	Xstims{1} = tmp;
end

if (nargin < 3) || isempty(target_mods)
	target_mods = 1:Nmods;
end
if (nargin < 4) || isempty(lambda_smooth)
	lambda_smooth = 0;
end
if (nargin < 5) || isempty(NLmon)
	NLmon = 1; % monotonic increasing
end
if (nargin < 6) || isempty(edge_p)
	edge_p=0.05; %relative to the generating distribution (pth percentile) where to put the outermost tent bases
end
if (nargin < 7) || isempty(edge_p)
	n_bfs = 25; %default number of tent bases
end
if nargin < 8
	space_type = 'equispace'; %default uninimrm tent basis spacing
end

target_mods(target_mods <= 0) = [];

%store NL tent-basis parameters
NL_tb_params.edge_p = edge_p;
NL_tb_params.n_bfs = n_bfs;
NL_tb_params.space_type = space_type;
nim.NL_tb_params = NL_tb_params;

% Compute internal generating functions
%Kmat = [nim.mods(:).filtK];
%gint = Xstim*Kmat;
gint = nan(size(Xstims{1},1),Nmods);
for imod = 1:Nmods		
	gint(:,imod) = Xstims{nim.mods(imod).Xtarget} * nim.mods(imod).filtK;
end

for imod = target_mods 	%nimr each module
    
    prev_NL_type = nim.mods(imod).NLtype;
    nim.mods(imod).NLtype = 'nonpar'; %set the subunit NL type to nonpar
    
    if ~isempty(n_bfs)
        n_dxs = n_bfs; %if specifying the number of basis functions
    else
        %update existing tent basis X-axis
        if ~isfield(nim.mods(imod).NLx)
            error('No existing tent basis, need to specify params to create one');
        else
            n_dxs = length(nim.mods(imod).NLx);
        end
    end
    if strcmp(space_type,'equispace')
        %adjust edges so that one of the bins lands on 0
        left_edge = my_prctile(gint(:,imod),edge_p);
        right_edge = my_prctile(gint(:,imod),100-edge_p);
        if left_edge == right_edge
          % make unit length
          left_edge = right_edge - 0.5;
          right_edge = right_edge + 0.5;
        end
        spacing = (right_edge - left_edge)/n_dxs;
        left_edge = ceil(left_edge/spacing)*spacing;
        right_edge = floor(right_edge/spacing)*spacing;
        
        nim.mods(imod).NLx = linspace(left_edge,right_edge,n_dxs); %equispacing
    elseif strcmp(space_type,'equipop')
        if std(gint(:,imod)) == 0  % constant term
          nim.mods(imod).NLx = mean(gint(:,imod)) + linspace(-0.5,0.5,n_dxs);
        else  
          nim.mods(imod).NLx = my_prctile(gint(:,imod),linspace(edge_p,100-edge_p,n_dxs)); %equipopulated
        end
    else
        error('Not supported')
    end
    %set nearest tent basis to 0
    [~,nearest] = min(abs(nim.mods(imod).NLx));
    nim.mods(imod).NLx(nearest) = 0;
    
    %initalize tent basis coefs
    if strcmp(prev_NL_type,'lin');
        nim.mods(imod).NLy = nim.mods(imod).NLx;
    elseif strcmp(prev_NL_type,'threshlin')
        nim.mods(imod).NLy = nim.mods(imod).NLx;
        nim.mods(imod).NLy(nim.mods(imod).NLy < 0) = 0;
    elseif strcmp(prev_NL_type,'quad')
        nim.mods(imod).NLy = nim.mods(imod).NLx.^2;
    elseif strcmp(prev_NL_type,'nonpar')
        %disp('TB Upstream NL already initialized. Adjusting regularization...');
    else
        error('Unsupported NL type');
    end
    if ~isnan(lambda_smooth) %set lambda_smooth to nan to avoid resetting lambda_NLd2
    nim.mods(imod).reg_params.lambda_NLd2 = lambda_smooth;
    end
    if ~isnan(NLmon) %set NLmon to nan to avoid resetting NLmon
    nim.mods(imod).NLmon = NLmon;
    end
end
