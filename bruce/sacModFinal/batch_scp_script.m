clear all

fit_unCor = true;

base_fname = 'sacStimProcFin_noXV';
% base_fname = 'sac_info_timing';
% base_fname = 'sac_info_timing';
% base_fname = 'corrected_models2';

Expt_list = {'G085','G086','G087','G088','G089','G091','G093','G095'};
ori_list = [0 90; 0 90; 0 90; 0 90; 0 90; 0 90; 0 90; 0 nan];
% ori_list = [0 nan; 0 nan; 0 nan; 0 nan; 0 nan; 0 nan; 0 nan; 0 nan];

copy_to = 'local';

results = zeros(length(Expt_list),2);
for ee = 1:length(Expt_list)
    Expt_name = Expt_list{ee};
    Expt_num = str2num(Expt_name(2:end));
%     source_dir = ['/home/james/Analysis/bruce/' Expt_name '/models/'];
%     dest_dir = ['/home/james/Analysis/bruce/' Expt_name '/models/'];
    source_dir = ['/home/james/Analysis/bruce/' Expt_name '/FINsac_mod/'];
    dest_dir = ['/home/james/Analysis/bruce/' Expt_name '/FINsac_mod/'];
    
    if strcmp(copy_to,'local')
        if ~exist(dest_dir,'dir')
            r = system(sprintf('mkdir -p %s',dest_dir));
            if r
                error('couldnt make directory');
            end
            fprintf('made directory %s\n',dest_dir);
        end
    end
    
    for oo = 1:2
        cur_ori = ori_list(ee,oo);
        if ~isnan(cur_ori)
            cur_fname = strcat(source_dir,base_fname,sprintf('_ori%d',cur_ori));
            if fit_unCor
                cur_fname = strcat(cur_fname,'_unCor');
            end
            cur_fname = strcat(cur_fname,'.mat');
            if strcmp(copy_to,'local')
                fprintf('Copying %s from remote to %s\n',cur_fname,dest_dir);
                results(ee,oo) = system(sprintf('scp james@Retina:%s %s',cur_fname,dest_dir));
            elseif strcmp(copy_to,'remote')
                fprintf('Copying %s to remote %s\n',cur_fname,dest_dir);
                results(ee,oo) = system(sprintf('scp %s james@Retina:%s',cur_fname,dest_dir));
            else
                error('invalid copy to');
            end
            %             results(ee,oo) = system(sprintf('scp james@Retina:%s %s',cur_fname,sac_dir));
            %             results(ee,oo) = system(sprintf('scp james@CA1:%s %s',cur_fname,targ_dir));
        end
    end
end

if any(results==1)
    fprintf('ERRORS DETECTED IN SOME TRANSFERS!\n');
end


%%
Expt_list = {'M266','M270','M275','M277','M281','M287','M289','M294','M296'};
ori_list = [80 nan; 60 nan; 135 nan; 70 nan; 140 nan; 90 nan; 160 nan; 40 nan; 45 nan; 0 90];

results = zeros(length(Expt_list),2);
for ee = 1:length(Expt_list)
    Expt_name = Expt_list{ee};
    Expt_num = str2num(Expt_name(2:end));
   
    source_dir = ['/home/james/Analysis/bruce/' Expt_name '/models/'];
    dest_dir = ['/home/james/Analysis/bruce/' Expt_name '/models/'];
   
    if strcmp(copy_to,'local')
        if ~exist(dest_dir,'dir')
            r = system(sprintf('mkdir -p %s',dest_dir));
            if r
                error('couldnt make directory');
            end
            fprintf('made directory %s\n',dest_dir);
        end
    end
    
    for oo = 1:2
        cur_ori = ori_list(ee,oo);
        if ~isnan(cur_ori)
            cur_fname = strcat(source_dir,base_fname,sprintf('_ori%d',cur_ori));
            if fit_unCor
                cur_fname = strcat(cur_fname,'_unCor');
            end
            cur_fname = strcat(cur_fname,'.mat');
            if strcmp(copy_to,'local')
                fprintf('Copying %s to remote %s\n',cur_fname,dest_dir);
                results(ee,oo) = system(sprintf('scp james@Retina:%s %s',cur_fname,dest_dir));
            elseif strcmp(copy_to,'remote')
                fprintf('Copying %s from remote to %s\n',cur_fname,dest_dir);
                results(ee,oo) = system(sprintf('scp %s james@Retina:%s',cur_fname,dest_dir));
            else
                error('invalid copy to');
            end
            %             results(ee,oo) = system(sprintf('scp james@Retina:%s %s',cur_fname,sac_dir));
            %             results(ee,oo) = system(sprintf('scp james@CA1:%s %s',cur_fname,targ_dir));
        end
    end
end

if any(results==1)
    fprintf('ERRORS DETECTED IN SOME TRANSFERS!');
end

