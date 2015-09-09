function tent_out = get_tentbasis_output( gin, tent_cent, tent_edges )
% 
% tent_out = get_tentbasis_output( gin, tent_cent, tent_edges )
% 
% Takes an input vector and passes it through the tent basis function
% specified by center location tent_cent and the 2-element vector tent_edges = [left_edge right_edge]
% specifying the tent bases 'edges'

%%
tent_out = zeros(size(gin)); %initialize NL processed stimulus

%for left side
if ~isinf(tent_edges(1)) %if there is a left boundary
  cur_set = (gin > tent_edges(1)) & (gin <= tent_cent); %find all points left of center and right of boundary
  tent_out(cur_set) = 1-(tent_cent-gin(cur_set))/(tent_cent - tent_edges(1));%contribution of this basis function to the processed stimulus    
else
  cur_set = gin <= tent_cent;
  tent_out(cur_set) = 1;
end

%for right side
if ~isinf(tent_edges(2)) %if there is a left boundary
  cur_set = (gin >= tent_cent) & (gin < tent_edges(2)); %find all points left of center and right of boundary
  tent_out(cur_set) = 1-(gin(cur_set)-tent_cent)/(tent_edges(2)-tent_cent);%contribution of this basis function to the processed stimulus    
else
  cur_set = gin >= tent_cent;
  tent_out(cur_set) = 1;
end

