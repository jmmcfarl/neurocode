function [r,gen_functions] = get_predicted_rate_stcbf(glmod,kern_output)

%extract matrix of STCcfs across modules
STCcf_mat = get_STCcf_mat(glmod);
%internal filter outputs of each module
gen_functions = kern_output*STCcf_mat;

fsdim = glmod.mods(1).fsdim;
hlen = 1;

kx    = glmod.const; %initialize kx with model constant term
nmods = length(glmod.mods);

NT = size(kern_output,1);

for imod = 1:nmods %loop across NL modules
    
    mod = glmod.mods(imod);
    fg = gen_functions(:,imod);
    fg(fg < 0) = 0;
    
    kx = kx + fg*mod.w;
    
end

kx = kx(hlen:end);

kx(kx > 20)    = 20; %saturate input to spiking NL
r             = log(1+exp(kx)); %apply spiking NL tp get predicted rate

