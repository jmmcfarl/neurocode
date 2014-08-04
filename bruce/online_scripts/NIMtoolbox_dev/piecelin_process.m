function gout = piecelin_process(gin, NLy, NLx)
%
% gout = piecelin_process(gin, NLy, NLx)
%
% Process input gin by the piecewise linear function specified by
% bound-points in vector NLx and center-heights NLy. Equivalent to the
% summed outputs of the tent basis functions specified by NLx and NLy

%%
[NT NX] = size(gin);
if NX > 1
    error('Must input a column vector')
end
gout = zeros(size(gin));

%Data where X < NLx(1) are determined entirely by the first tent basis
left_edge = find(gin < NLx(1));
gout(left_edge) = gout(left_edge) + NLy(1);

%similarly for the right edge
right_edge = find(gin > NLx(end));
gout(right_edge) = gout(right_edge) + NLy(end);

slopes = diff(NLy)./diff(NLx);
for j = 1:length(NLy)-1
   cur_set = find(gin >= NLx(j) & gin < NLx(j+1));
   gout(cur_set) = gout(cur_set) + NLy(j) + slopes(j)*(gin(cur_set) - NLx(j));
end

