function [] = plotfilterbank(allks,sdim,spixids)
% USAGE: [] = plotfilterbank(allks,sdim,spixids
%   plot filters in longk-format

%figure;
nslices = size(allks,1)/length(spixids);
nfilts  = size(allks,2);
zmin    = min(allks(:)); zmax=max(allks(:));
% disp('*** plotting normalized over all filters & times')
for ifilt = 1:nfilts;
    zmin    = min(allks(:,ifilt)); zmax=max(allks(:,ifilt));
%     largest = max(abs([zmin zmax]));
%     zmin = -largest; zmax = largest;
    if zmin == zmax
        zmin = -1; zmax = 1;
    end
    tcell = pad2dslices(allks(:,ifilt),spixids,sdim);
    % zmin    = min(allks(:,ifilt)); zmax=max(allks(:,ifilt));
    for islice = 1:nslices;
        subplot(nfilts,nslices,(ifilt-1)*nslices+islice);
        plot2drfmat(flipud(tcell{islice}),[zmin,zmax]);
    end
end

