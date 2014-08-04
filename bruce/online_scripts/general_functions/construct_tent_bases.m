function tbmat = construct_tent_bases(tent_centers,dx)

n_tents = length(tent_centers);
NX = ceil(range(tent_centers)/dx);

tbmat = zeros(n_tents,NX);
for i = 1:n_tents
    if i > 1
        cur_lset = tent_centers(i-1):dx:tent_centers(i);
        cur_linds = 1+round(cur_lset/dx);
        cur_x = 0:dx:(tent_centers(i)-tent_centers(i-1));
        cur_a = 1/(tent_centers(i)-tent_centers(i-1));
        tbmat(i,cur_linds) = cur_x*cur_a;
    end
    if i < n_tents
        cur_rset = tent_centers(i):dx:tent_centers(i+1);
        cur_rinds = 1+round(cur_rset/dx);
        cur_x = 0:dx:(tent_centers(i+1)-tent_centers(i));
        cur_a = 1/(tent_centers(i+1)-tent_centers(i));
        tbmat(i,cur_rinds) = fliplr(cur_x*cur_a);
    end
end
