function new = shift_mat2(old,x,y,wrap,bval)


if wrap ~= 0 & wrap ~= 1
    error('Invalid wrap value')
end

new = old;
[sy,sx] = size(old);
if x > 0
    if wrap == 0
        rep_buf = bval*ones(sy,x);
    else
        rep_buf = new(:,end-x+1:end);
    end
    new = [rep_buf new(:,1:end-x)]; 
elseif x < 0
    x = abs(x);
    if wrap == 0
        rep_buf = bval*ones(sy,x);
    else
        rep_buf = new(:,1:x);
    end
    new = [new(:,x+1:end) rep_buf];
end

if y > 0
    if wrap == 0
       rep_buf = bval*ones(y,sx); 
    else
       rep_buf = new(end-y+1:end,:);
    end
    new = [rep_buf; new(1:end-y,:)]; 
elseif y < 0
    y = abs(y);
    if wrap == 0
       rep_buf = bval*ones(y,sx); 
    else
       rep_buf = new(1:y,:);
    end
    new = [new(y+1:end,:); rep_buf];
end