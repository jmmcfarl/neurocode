function locations = find_in_struct_vec(cur_struct,field,target)

struct_length = length(cur_struct);

locations = [];
for i = 1:struct_length
    if isstr(target)
        if strcmp(target, getfield(cur_struct,{i},field))
            locations = [locations i];
        end
    else
        if target ==  getfield(cur_struct,{i},field);
              locations = [locations i];
        end          
    end
end