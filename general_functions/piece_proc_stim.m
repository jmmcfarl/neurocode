function proc_s = piece_proc_stim( s, center, di )
%
% Usage: proc_s = piece_proc_stim( s, center, [left_pt right_pt] )

%computes the contribution of a given NL tent basis function to the overall
%processed stimulus.

if length(di) == 2
  rl = center-di(1);
  rr = di(2)-center;
else %if only given one di input
  if di < center %if di is left of the center
    rl = center-di;
    rr = [];
  else %if di is right of center
    rl = []; %no left step
    rr = di-center; %rightwards x step size
  end
end

proc_s = zeros(length(s),1); %initialize NL processed stimulus

if ~isempty(rl) %if there is a left boundary
  a = find((s > center-rl) & (s <= center)); %find all points left of center and right of boundary
  proc_s(a) = 1-(center-s(a))/rl;%contribution of this basis function to the processed stimulus
else %if no left boundary
  %a = find(s <= center);
  %proc_s(a) = 1;
  proc_s(s <= center) = 1; %set everything left of the center to 1 (outside range of NL function)
end
if ~isempty(rr) %if there is a right boundary
  a = find((s >= center) & (s < center+rr)); %find all points right of center and left of boundary
  proc_s(a) = 1-(s(a)-center)/rr; %contribution of this basis function to the processed stimulus
else
  %a = find(s >= center);
  %proc_s(a) = 1;
  proc_s(s >= center) = 1;
end
