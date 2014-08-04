function [new_Expts] = add_field_from_online(online_fname,Expts,expt_nums,fieldname)

fid = fopen(online_fname,'r');
if fid < 0
    mycprintf('errors','Cant Read %s\n',online_fname);
    return;
end
a = textscan(fid,'%s','delimiter','\n','bufsize',2^14-1);
fclose(fid);
txt = a{1};

%%
idid = find(strncmp('id',txt,2));
trialid = zeros(length(idid),1);
for ii = 1:length(idid)
    [trialid(ii)] = sscanf(txt{idid(ii)},'id%d');
end
[trialid,ia] = unique(trialid);
idid = idid(ia);

%%
all_field_lines = find(strncmp([fieldname '='],txt,length(fieldname)+1));

%%
new_Expts = Expts;
%%
for ee = 1:length(expt_nums)
    cur_expt_num = expt_nums(ee);
    expt_tids = [Expts{cur_expt_num}.Trials(:).id];
    
    n_trials = length(Expts{cur_expt_num}.Trials);
    field_values = cell(n_trials,1);
    for nn = 1:n_trials
       cur_line = idid(trialid == expt_tids(nn)); 
       if expt_tids(nn) == max(trialid)
          search_to = length(txt); 
       else
           search_to = idid(find(trialid > expt_tids(nn),1));
       end
       cur_range = cur_line:search_to;
       target = cur_range(ismember(cur_range,all_field_lines));
       if length(target) > 1
           target = target(1);
       elseif isempty(target)
           target = nan; 
       end
       if ~isnan(target)
       field_values{nn} = txt{target}((length(fieldname)+2):end);
       else
          field_values{nn} = nan; 
          fprintf('Warning, missing values in online file\n');
       end
       new_Expts{cur_expt_num}.Trials(nn).(fieldname) = field_values{nn};
    end
end
%%
