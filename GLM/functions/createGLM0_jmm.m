function [glm0] = createGLM0_jmm(ks,defmod,mname)
%%USAGE: [glm0] = createGLM0(ks,defmod,mname)
%   takes the filters and default module structure to create a glm
mods = [];
for imod = 1:size(ks,2);
	tmod  = defmod; %initialize module param structure
	tmod.k = ks(:,imod); %store current filter coefs
	[matchprof,sliceid] = getMatchProf(ks(:,imod),defmod.SDIM); 
	tmod.sliceid        = sliceid;
	tmod.matchprof      = matchprof;
	mods                = [mods,tmod]; 
end;
glm0   = struct('mods',mods,'const',0,'LL',0,'LLnadj',0,'mname',mname,'lambdaW',0);


