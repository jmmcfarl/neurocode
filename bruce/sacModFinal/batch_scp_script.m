base_fname = 'sacStimProc3';

Expt_list = {'G085','G086','G087','G088','G089','G091','G093','G095'};
ori_list = [0 90; 0 90; 0 90; 0 90; 0 90; 0 90; 0 90; 0 nan];

results = zeros(length(Expt_list),2);
for ee = 1:length(Expt_list)
    Expt_name = Expt_list{ee};
    Expt_num = str2num(Expt_name(2:end));
    sac_dir = ['/home/james/Analysis/bruce/' Expt_name '/FINsac_mod/'];
    
    for oo = 1:2
        cur_ori = ori_list(ee,oo);
        if ~isnan(cur_ori)
            cur_fname = strcat(sac_dir,base_fname,sprintf('_ori%d.mat',cur_ori));
            fprintf('Copying from %s to %s\n',cur_fname,sac_dir);
            results(ee,oo) = system(sprintf('scp james@Retina:%s %s',cur_fname,sac_dir));
        end
    end
end

if any(results==1)
    fprintf('ERRORS DETECTED IN SOME TRANSFERS!');
end


%%
% Expt_list = {'M266','M270','M275','M277','M281','M287','M289','M294','M296'};
% ori_list = [80 nan; 60 nan; 135 nan; 70 nan; 140 nan; 90 nan; 160 nan; 40 nan; 45 nan; 0 90];
% 
% results = zeros(length(Expt_list),2);
% for ee = 1:length(Expt_list)
%     Expt_name = Expt_list{ee};
%     Expt_num = str2num(Expt_name(2:end));
%     sac_dir = ['/home/james/Analysis/bruce/' Expt_name '/FINsac_mod/'];
%     
%     for oo = 1:2
%         cur_ori = ori_list(ee,oo);
%         if ~isnan(cur_ori)
%             cur_fname = strcat(sac_dir,base_fname,sprintf('_ori%d.mat',cur_ori));
%             fprintf('Copying from %s to %s\n',cur_fname,sac_dir);
%             results(ee,oo) = system(sprintf('scp james@Retina:%s %s',cur_fname,sac_dir));
%         end
%     end
% end
% 
% if any(results==1)
%     fprintf('ERRORS DETECTED IN SOME TRANSFERS!');
% end
