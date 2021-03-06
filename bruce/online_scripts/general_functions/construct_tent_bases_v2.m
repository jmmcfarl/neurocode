function tbmat = construct_tent_bases_v2(tent_centers,x_ax)

n_tents = length(tent_centers);
NX = length(x_ax);
dx = median(diff(x_ax));
tbmat = zeros(n_tents,NX);
for i = 1:n_tents
    if i == 1
%         cur_xloc = round(tent_centers(i)/dx);
%         cur_inds = 0:2*cur_xloc;
%         cur_x = cur_inds;
%         cur_x((cur_xloc+2):end) = (cur_xloc-dx):-dx:0;
%         bad = find(cur_inds < 1);
%         cur_inds(bad) = [];
%         cur_x(bad) = [];
%         tbmat(i,cur_inds) = cur_x/max(cur_x);

            cur_rset = tent_centers(i):dx:tent_centers(i+1);
        cur_rinds = round(cur_rset/dx);
        cur_x = 0:dx:(tent_centers(i+1)-tent_centers(i));
        cur_a = 1/(tent_centers(i+1)-tent_centers(i));
        tbmat(i,cur_rinds) = fliplr(cur_x*cur_a);
end
    if i == n_tents
        cur_lset = tent_centers(i-1):dx:tent_centers(i);
        cur_linds = round(cur_lset/dx);
        cur_x = 0:dx:(tent_centers(i)-tent_centers(i-1));
        cur_a = 1/(tent_centers(i)-tent_centers(i-1));
        tbmat(i,cur_linds) = cur_x*cur_a;
    end
    if i > 1 && i < n_tents
        cur_lset = tent_centers(i-1):dx:tent_centers(i);
        cur_linds = round(cur_lset/dx);
        cur_x = 0:dx:(tent_centers(i)-tent_centers(i-1));
        cur_a = 1/(tent_centers(i)-tent_centers(i-1));
        tbmat(i,cur_linds) = cur_x*cur_a;
        cur_rset = tent_centers(i):dx:tent_centers(i+1);
        cur_rinds = round(cur_rset/dx);
        cur_x = 0:dx:(tent_centers(i+1)-tent_centers(i));
        cur_a = 1/(tent_centers(i+1)-tent_centers(i));
        tbmat(i,cur_rinds) = fliplr(cur_x*cur_a);
    end
end
