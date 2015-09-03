function [data,time] = convert_ascii_heka_sweeps(fname)

fid = fopen(fname);

data = [];
% block_num = 1;
while ~feof(fid)
    tline = fgetl(fid);
    next_vals = sscanf(tline,'%g %g %g');
    while isempty(next_vals)
        tline = fgetl(fid);
        next_vals = sscanf(tline,'%g %g %g');
    end
    n_cols = length(next_vals);
    next_data = fscanf(fid,'%g',[n_cols Inf]);
    data = [data next_vals next_data];
%     block_num = block_num + 1
end

time = data(2,:);
data = data(3,:);