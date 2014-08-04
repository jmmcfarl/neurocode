function sigout = piececonst_process(sigin, ypts, xedges)
%
% sigout = piececonst_process(sigin, ypts, xedges)
%
% Process input sigin by the piecewise constant function defined by edge
% points xedges and values in ypts

%%
sigout = zeros(length(sigin),1);
for n = 1:length(xedges)-1
  sigout((sigin >= xedges(n)) & (sigin < xedges(n+1))) = ypts(n);
end

