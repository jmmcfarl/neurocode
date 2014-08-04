function CSD = PettersenCSD(Data, type, vars)
% Questions:
% -is it possible to put negative channel positions? (I think no?
% - how about the change in conductivity that occurs at the white matter?
% There is much more fat there, so conductivity should be lower...

%% General
% => This function computes the current-source density matrix by using the
% functions provided in Petterson's iCSD toolbox.
% => If the data has more than 1 trial, the CSD is computed separately for
% each trial
% => For all methods except the standard CSD method, it is required that all
% the channels are inside cortex!!!
% => The step and spline methods predict a continuous function and the
% standard and delta a discrete CSD
% => The  Bad Channels are handled differently depending on the method:
% - In the standard method, the LFP at the bad channel is interpolated and
% then the CSD is formed
% - In the other methods, the corresponding row of the data matrix is
% simply removed
% - In the delta method, the resulting CSD is then interpolated to retrieve
% the same number of channels than originally.
% - In the spline and step method there is no need for subsequent
% interpolation because the retrieved CSD is continuous.

%% Inputs
% Data: in the form of channels x time x trials
% type: 'standard', 'delta', 'step', 'spline'
% vars: structure of variables necessary for each method

% Variables:
% - BadChannels (optional, can also be empty): default is none
% - Fs (required): sampling frequency
% - BrainBound (required for all except standard): The number of the first
% contact in the brain. All the new CSD methods can only use Data from
% channels within the brain.
% - offset (optional, used for all except standard, can also be empty): The distance between the
% first contact and the top of the brain. This needs to be smaller than the
% Channel Separation. [default: ChanSep/2]
% - ChanSep (optional, can also be empty): distance in mm between adjacent
% channels [default: 0.1 mm]
% - electrodePos (optional, can also be empty): the position in mm of the
% electrodes on the probe [default: 0:ChanSep:ChanSep*(NumChan -1)]
% - exCond (optional, can also be empty): conductivity of extracellular
% medium in S/m [default: 0.276 S/m]
% - topCond (optional, used for all except standard, can also be empty): conductivity in S/m on top
% of the brain [default: 1.654]
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


%% Outputs
% CSD: current source density matrix of the same size as the Data matrix.

%% Parse the different variables
dt = 1/vars.Fs;
NumChan = size(Data,1);

if ~isfield(vars,'ChanSep') || isempty(vars.ChanSep)
    ChanSep = 0.1;
else
    ChanSep = vars.ChanSep;
end

if ~isfield(vars,'offset') || isempty(vars.offset)
    offset = ChanSep/2;
else
    offset = vars.offset;
    if offset > ChanSep
        warning('The offset between the first channel in the brain and the surface is larger than the Channel separation!!! This is certainly wrong! Change the offset or adapt the BrainBound.');
    end
end

if ~isfield(vars, 'electrodePos') || isempty(vars.electrodePos)
    electrodePos = 0:ChanSep:(ChanSep*(NumChan-1));
else
    electrodePos = vars.electrodePos;
end

if ~isfield(vars,'exCond') || isempty(vars.exCond)
    exCond = 0.276;
else
    exCond = vars.exCond;
end

if ~isfield(vars,'topCond') || isempty(vars.topCond)
    topCond = 1.654;
else
    topCond = vars.topCond;
end

if ismember(lower(type), {'standard', 'delta'})
    if vars.useHamming
        if ~isfield(vars,'b0') || isempty(vars.b0)
            b0 = 0.56;
        else
            b0 = vars.b0;
        end
        if ~isfield(vars,'b1') || isempty(vars.b1)
            b1 = 0.23;
        else
            b1 = vars.b1;
        end
    else
        b1 = 0;
        b0 = 0;
    end
end

if ismember(lower(type), {'step', 'spline'})
    if ~isfield(vars,'dev') || isempty(vars.dev)
        dev = 0.1;
    else
        dev = vars.dev;
    end
end


%% Run the computation
CSD = nan(size(Data));

if ismember(lower(type), {'delta', 'step', 'spline'})
    electrodePos_noBad = electrodePos; % electrode positions without the bad channels removed
    chans = 1:NumChan;
    if isfield(vars,'BadChannels') && ~isempty(vars.BadChannels)
        BadChan = vars.BadChannels;
        BadChan(BadChan < vars.BrainBound) = [];
        if ~isempty(BadChan)
            
            chans(BadChan) = []; % This is for the interpolation of the CSD
            electrodePos(BadChan) = [];
            Data(BadChan,:,:) = []; % remove the data that is not inside the brain
        end
    end
    chans = chans(vars.BrainBound:end);
    Data(1:vars.BrainBound-1,:,:) = []; % remove the data that is not inside the brain
    electrodePos_noBad = electrodePos_noBad(vars.BrainBound:end) - (vars.BrainBound - 1)*ChanSep + offset;
    electrodePos = electrodePos(vars.BrainBound:end) - (vars.BrainBound - 1)*ChanSep + offset; % only use the electrodes in the brain
end

switch lower(type)
    case 'standard'
        if vars.useVaknin
            VakIdx = 0;
        else
            VakIdx = 1;
        end
        if isfield(vars,'BadChannels') && ~isempty(vars.BadChannels)
            % Interpolate Data across BadChannels
            DataTemp = Data;
            DataTemp(vars.BadChannels,:,:) = [];
            chans = 1:NumChan;
            chans(vars.BadChannels) = [];
            Data = interp1(chans, DataTemp, 1:NumChan, 'spline');
        end
        
        for i = 1:size(Data,3)
            CSDTemp = newStandardCSD( squeeze(Data(:,:,i)), electrodePos, dt, exCond, vars.useVaknin, b0, b1 );
            CSD(1+VakIdx:end-VakIdx,:,i) = CSDTemp;
        end
        
    case 'delta'
        
        for i = 1:size(Data,3)
            CSDTemp = newDeltaCSD( squeeze(Data(:,:,i)), electrodePos, dt, b0, b1, exCond, topCond, vars.diam );
            if ~isempty(setdiff(chans, vars.BrainBound:NumChan));
                CSDTemp = interp1(chans, CSDTemp, vars.BrainBound:NumChan, 'spline');
            end
            CSD(vars.BrainBound:end,:,i) = CSDTemp;
        end
        
    case 'step'
        
        for i = 1:size(Data,3)
            [CSDTemp, zs] = newStepCSD(squeeze(Data(:,:,i)), electrodePos, dt, dev, exCond, topCond, vars.diam);
            if i == 1
                chan_idx = nan(1, length(electrodePos_noBad));
                for j = 1:length(electrodePos_noBad)
                    chan_idx(j) = find(electrodePos_noBad(j)*1e-3 <= zs, 1, 'first');
                end
            end
            CSD(vars.BrainBound:end,:,i) = CSDTemp(chan_idx,:);
        end
        
    case 'spline'
        
        for i = 1:size(Data,3)
            [CSDTemp, zs] = newSplineCSD(squeeze(Data(:,:,i)), electrodePos, dt, dev, exCond, topCond, vars.diam);
            if i == 1
                chan_idx = nan(1, length(electrodePos_noBad));
                for j = 1:length(electrodePos_noBad)
                    chan_idx(j) = find(electrodePos_noBad(j)*1e-3 <= zs, 1, 'first');
                end
            end
            CSD(vars.BrainBound:end,:,i) = CSDTemp(chan_idx,:);
        end
        
end

CSD = squeeze(CSD);
end