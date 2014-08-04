function d2 = dist_shift3d( d1, sh, dim)

sh = floor(sh+0.5);
L = size(d1);

if sh > L(dim)
    sh = L(dim);
end
if sh < -L(dim)
    sh = -L(dim);
end

mapper = mod((0:(L(dim)-1))-sh,L(dim))+1;
if dim == 1
  d2 = d1(mapper,:,:);
elseif dim == 2
  d2 = d1(:,mapper,:);
else
  d2 = d1(:,:,mapper);
end

if dim == 1
    if sh > 0
        d2(1:sh,:,:) = 0;
    else
        d2(L(1)-(0:(-sh-1)),:,:) = 0;
    end
elseif dim == 2
    if sh > 0
        d2(:,1:sh,:) = 0;
    else
        d2(:,L(2)-(0:(-sh-1)),:) = 0;
    end
else
   if sh > 0
      d2(:,:,1:sh) = 0; 
   else
       d2(:,:,L(3)-(0:(-sh-1))) = 0;
   end    
end

