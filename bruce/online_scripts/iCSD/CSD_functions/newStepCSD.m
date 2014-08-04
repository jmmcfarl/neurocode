function [CSD, zs] = newStepCSD(data, electrodePos, dt, dev, exCond, topCond, diam)

electrodePos = electrodePos * 1E-3;
dev = dev * 1E-3;
diam = diam * 1E-3;
filter_range = 5*dev;

CSD = F_const(electrodePos,diam,exCond, topCond)^(-1) * data;
numElectrodes = length(electrodePos);
h = mean( diff(electrodePos));

first_z = electrodePos(1) - h/2;
mfactor = ceil(200/numElectrodes);
minizs = 0:h/mfactor:(mfactor-1)*h/mfactor;

for i=1:size(CSD,1) % all rows
    zs((1:mfactor)+mfactor*(i-1)) = first_z+(i-1)*h+minizs;
    new_CSD_matrix((1:mfactor)+mfactor*(i-1),:) = repmat(CSD(i,:),mfactor,1);
end;

if(dev ~= 0)
    [zs,new_CSD_matrix]=gaussian_filtering(zs,new_CSD_matrix,dev,filter_range);
end

% plot_CSD2(new_CSD_matrix,zs,dt,1,0);
CSD = new_CSD_matrix;

end