function out = circ_pert_params(cur_params)

n_filts = (length(cur_params)-1)/2;
thetas = rand(n_filts,1)*2*pi;
out = zeros(length(cur_params),1);
for i = 1:n_filts
    cur_pt = [cur_params(2*(i-1)+1:2*i)];
    cur_rmat = [cos(thetas(i)) -sin(thetas(i)); sin(thetas(i)) cos(thetas(i))];
   out(2*(i-1)+1:2*i) = cur_rmat*cur_pt;
end
