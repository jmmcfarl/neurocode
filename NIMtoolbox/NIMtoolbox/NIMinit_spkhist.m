function nim = NIMinit_spkhist(nim,n_bins,init_spacing,doubling_time,negCon)
%
% Usage: nim = NIMinit_spkhist(nim,n_bins,init_spacing,doubling_time,negCon)
%
% Adds a spike history term with specified parameters to an existing NIM.
%
% INPUTS:
%     nim: input model structure
%     n_bins: number of coefficients in spike history filter
%     [init_spacing]: Initial spacing (in time bins) of piecewise constant 'basis functions'
%     [doubling_time]: Make bin spacing logarithmically increasing, with given doubling time
%     [negCon]: If == 1, constrain spike history filter coefs to be non-positive
% OUTPUTS:
%     nim: output model


%%
% default inputs
if nargin < 3
    init_spacing = 1;
end
if nargin < 4
    doubling_time = n_bins;
end
if nargin < 5
    negCon = 0;
end

%% COMPUTE RECTANGULAR BASIS FUNCTIONS
bin_edges = zeros(n_bins+1,1);
inc = init_spacing;
pos = 1; count = 0;
for n = 1:n_bins+1
    bin_edges(n) = pos;
    pos = pos + inc;  count = count + 1;
    if count >= doubling_time
        count = 0; inc = inc * 2;
    end
end

%% LOAD INTO A SPK_HIST SUBSTRUCT IN NIM
coefs = zeros(n_bins,1);
nim.spk_hist.bin_edges = bin_edges;
nim.spk_hist.coefs = coefs;
nim.spk_hist.negCon = negCon;
nim.spk_hist.spkhstlen = n_bins;
