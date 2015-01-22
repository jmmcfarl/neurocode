function h=PlotLFPphase(weight, phase, usechannel, bandnames)

Xmin=bandnames(1); Xmax=bandnames(end); % Min and max of Y-axis
logx=linspace(log(Xmin),log(Xmax),length(bandnames));
%%
weight = weight';
% weight = quadraticscaling(weight);
phase = phase';
%%

[Nch,Nfreq] = size(weight);
phaplot = ones(Nch,Nfreq,3);

cmap = colormap(hsv(128));
% phac = linspace(-pi,pi,129);

% phase value is in the range of [0, 2 * pi]
phac = linspace(0,2*pi,129);
phase(phase<0) = phase(phase<0)+2*pi;
weight = weight/max(weight(:));

for y = 1:Nfreq
    for x = 1:Nch
        cxy = cmap(find(phase(x,y)>phac,1,'last'),:);
        cxy = cxy*weight(x,y);
        phaplot(x,y,1:3) = 1-cxy; %
    end
end
phaplot(phaplot<0) = 0;
phaplot(phaplot>1) = 1;

image(logx,usechannel,phaplot);
set(gca, 'YDir', 'normal');

axis tight;axis square;
Xtick=[0.1 0.5 1 2 5 10 20  40 70];

set(gca,'XTick', log(Xtick), ...
    'XTickLabel', arrayfun(@num2str, Xtick(:), 'UniformOutput', false));
xlabel(' freq (Hz) '); ylabel('chan #');xlim([logx(1) logx(end)])