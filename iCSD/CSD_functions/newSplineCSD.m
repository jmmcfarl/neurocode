function [ CSD, zs ] = newSplineCSD(data, electrodePos, dt, dev, exCond, topCond, diam)
electrodePos = electrodePos * 1E-3;
dev = dev * 1E-3;
diam = diam * 1E-3;
filter_range = 5*dev;

Fcs = F_cubic_spline(electrodePos,diam,exCond,topCond);
[zs, CSD_cs] = make_cubic_splines(electrodePos,data,Fcs);

if(dev ~= 0)
    [zs,CSD_cs] = gaussian_filtering(zs,CSD_cs,dev,filter_range);
end

% plot_CSD(CSD_cs,zs,dt,1,0);
CSD = CSD_cs;

end