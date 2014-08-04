vars.BrainBound = 3;
vars.topCond = []; %100;
vars.diam = 0.45; % in mm
vars.BadChannels = [];%M.BadChannels;
vars.Fs = 1000; %M.DataG.(GroupName).DOA{Cond}.Fs;
vars.useVaknin = true;
vars.useHamming = true;

type = {...
%         'standard'...
%         'delta'...
        'step'...
%         'spline'...
    };
type = type{1};

Data = data;
% Data = pot1;

var = [1 2 3 4 5 6]; % BrainBound
% var = [0.3 0.45 0.6 1 5 10]; % diam3 
% var = [0 0.1 0.276 0.5 1 2 5 10 1000]; % topCond
if ~exist('var', 'var'); var = 1; end
for i = 1:length(var);
    vars.BrainBound = var(i);
%     vars.diam = var(i);
%     vars.topCond = var(i);
    
    CSD = PettersenCSD(Data, type, vars);
    
    % Plot
    TimeVec = 0:1/vars.Fs:(size(CSD,2)-1)/vars.Fs;
    Amp = max(nanmax(abs(CSD)));
    figure; imagesc(TimeVec, 1:size(CSD,1), CSD, [-Amp Amp]);colorbar;
    xlabel('Time (s)'); ylabel('Channels');
    if ~isfield(vars,'topCond') || isempty(vars.topCond)
        topCond = 1.654;
    else
        topCond = vars.topCond;
    end
    title(['CSD "' type '" with diam: ' num2str(vars.diam) ' mm  and topCond: ' num2str(topCond) ' S/m']);
end
clear('var', 'topCond', 'diams', 'TimeVec', 'Amp', 'i', 'vars');
% - BadChannels (optional, can also be empty): default is none
% - Fs (required): sampling frequency
% - BrainBound (required for all except standard): The number of the first
% contact in the brain. All the new CSD methods can only use Data from
% channels within the brain.
% - ChanSep (optional, can also be empty): distance in mm between adjacent
% channels [default: 0.1 mm]
% - electrodePos (optional, can also be empty): the position in mm of the
% electrodes on the probe [default: 0:ChanSep:ChanSep*(NumChan -1)]
% - exCond (optional, can also be empty): conductivity in S/m of extracellular
% medium [default: 0.276 S/m]
% - topCond (optional, can also be empty): conductivity in S/m on top
% of the brain [default: 1.654 (CSF-Ajay)]
% - diam (required for all except standard): diameter of the source
% cylinders in mm
% - useVaknin (required only standard): either true or false
% - useHamming (required only standard and delta): use a Hamming filter, either true or false. This
% will reduce the number of estimated contacts by 2
% - b0 (optional, can also be empty): center weight of hamming filter
% [default: 0.56];
% - b1 (optional, can also be empty): neighbor weight of hamming filter
% [default: 0.23];
% - dev (optional only step and spline - can also be empty): standard deviation of the gaussian
% filter used to filter the CSD. [default: 0.1]