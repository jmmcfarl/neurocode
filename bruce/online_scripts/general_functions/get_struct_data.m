function data = get_struct_data(dstruct,sids,fname)

data = [];
for ss = 1:length(sids)
    cur_data = getfield(dstruct(sids(ss)),fname);
    data = cat(2,data,cur_data(:));
end
data = data';
