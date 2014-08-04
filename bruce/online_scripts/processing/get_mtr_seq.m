function [mtr_seq,mtr_is_match] = get_mtr_seq(online_fname,Expts,expt_num)

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
expt_tids = [Expts{expt_num}.Trials(:).id];
[use_ids,corresp_ids] = ismember(expt_tids,trialid);
%%
mtrPid = find(strncmp('mtrP',txt,2));
mtr_is_match = zeros(length(expt_tids),1);
mtr_seq = cell(length(expt_tids),1);
for ii = 1:length(expt_tids)
    next_mtr = find(mtrPid > idid(corresp_ids(ii)),1,'first');
    if corresp_ids(ii) < length(idid)
        if mtrPid(next_mtr) < idid(corresp_ids(ii)+1)
            mtr_is_match(ii) = 1;
        end
    elseif ~isempty(next_mtr)
        mtr_is_match(ii) = 1;
    end
    mtr_seq{ii} = sscanf(txt{mtrPid(next_mtr)},'mtrP=%s');
end