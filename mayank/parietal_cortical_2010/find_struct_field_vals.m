function [indices] = find_struct_field_vals(struct_in,field_in,match)

indices = [];
for i = 1:length(struct_in)
    cur_val = getfield(struct_in(i),field_in);
    if isa(match,'char')
        if strcmp(cur_val,match)
            indices = [indices i];
        end
    else
        if match==cur_val
            indices = [indices i];
        end
    end
end

if isempty(indices) && strcmp(field_in,'name')
    if match(5) == '_'
        match(5) = '-';
    else
        match(5) = '_';
    end
    for i = 1:length(struct_in)
        cur_val = getfield(struct_in(i),field_in);
        if isa(match,'char')
            if strcmp(cur_val,match)
                indices = [indices i];
            end
        else
            if match==cur_val
                indices = [indices i];
            end
        end
    end
end