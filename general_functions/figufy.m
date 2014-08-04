function figufy(fig_handle)

allAxes = findall(fig_handle,'type','axes');
allText = findall(fig_handle,'type','text');

set(allAxes,'fontsize',12,'fontname','arial');
set(allText,'fontsize',16,'fontname','times');

box off
ch = get(fig_handle,'children');
for ii = 1:length(ch)
    if strcmp(get(ch(ii),'type'),'axes')
    set(ch(ii),'box','off');
    end
end
