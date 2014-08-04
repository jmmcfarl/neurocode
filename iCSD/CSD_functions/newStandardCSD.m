function [ CSD ] = newStandardCSD( data, electrodePos, dt, conductivity, useVaknin, b0, b1 )
%NEWSTANDARDCSD Summary of this function goes here
%   Detailed explanation goes here
electrodePos = electrodePos * 1E-3;
switch(nargin)
    case 0
        error('Too few input parameters')
    case 1
        % TODO
        % set default values
end

if( useVaknin == 0)
    useElectrodes = electrodePos(2:end-1);
    N = length(electrodePos);
    h = mean(diff(electrodePos));
    [numElectrodes, numData] = size(data);
else
    useElectrodes = electrodePos;
    buf = data;
    data(1,:) = buf(1,:);
    data(2:end+1,:) = buf;
    data(end+1,:) = buf(end,:);
    clear buf;
    N = length(electrodePos);
    h = mean(diff(electrodePos));
    [numElectrodes, numData] = size(data);
end

CSD = -conductivity*D1(numElectrodes,h)*data;
if(b1 ~= 0)
    [n1, n2] = size(CSD);
    CSD_add(1,:) = zeros(1,n2);
    CSD_add(n1+2,:) = zeros(1,n2);
    CSD_add(2:n1+1,:) = CSD;
    CSD = S_general(n1+2,b0,b1)*CSD_add;
end

% plot_CSD(CSD,useElectrodes,dt,1,0);

end