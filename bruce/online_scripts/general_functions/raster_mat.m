function raster_mat( spmat, range, rep_bound_times, mk_size, colors)
%
% Usage: raster( splist, <range>, <color> )
%
% appended sequential spike lists - if next spike time
% is smaller than previous, new raster
%
% ns returns the number of spikes in each column
%
% if 'color' is used, assume plot raster on same figure
%
% Created 17 Apr 2003 DAB


y=1;
ncells = length(spmat);
mint = range(1);
maxt = range(2);

spacing = 1;
y = 0;

figure;set(gcf,'Position',[500 1000 800 1000])
hold on
if nargin < 5
    colors = zeros(ncells,3);
    colors([2:2:ncells],1) = 1;
%     colors([2:2:ncells],1) = 1;
%     colors([2:2:ncells],2) = 0.5;
    colors([1:2:ncells],3) = 1;
end

nreps = size(rep_bound_times,1);

for cc = 1:ncells
    for n = 1:nreps
        cur_spks = spmat{cc}(spmat{cc} > rep_bound_times(n,1) & spmat{cc} < rep_bound_times(n,2));
        cur_spks = cur_spks - rep_bound_times(n,1);
        cur_spks = cur_spks(cur_spks > range(1) & cur_spks < range(2));
        plot(cur_spks,ones(size(cur_spks))*y,'*','color',colors(cc,:),'markersize',mk_size);
        y = y + spacing;
    end
    y = y + spacing;
    line([mint maxt],[y y],'color','k')
    y = y + spacing;
end

axis([mint maxt 0 (y+1)])
set(gca,'LineWidth',0.6,'FontSize',8,'YTickLabel',[])
