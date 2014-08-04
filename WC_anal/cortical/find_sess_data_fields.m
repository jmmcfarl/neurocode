function [indices, issame] = find_sess_data_fields(sess_data, field, target_val)

issame = zeros(1,length(sess_data));

for i = 1:length(sess_data)
    
   eval(['issame(i) = strcmp(sess_data(i).' field ',target_val);'])
    
end

indices = find(issame);