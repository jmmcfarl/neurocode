function [data,Fs] = convert_ascii_heka_sweeps(fname)

fid = fopen(fname);

data = [];
while ~feof(fid)
    tline = fgetl(fid);
    next_vals = sscanf(tline,'%g %g %g');
    while isempty(next_vals)
        tline = fgetl(fid);
        next_vals = sscanf(tline,'%g %g %g');
    end
    next_data = fscanf(fid,'%g %g %g',[3 Inf]);
    data = [data next_vals next_data];
end

time = data(2,:);
Fs = 1/median(diff(time));