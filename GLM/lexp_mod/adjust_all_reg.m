function new_mod = adjust_all_reg(fo,regname,newval)

nmods = length(fo.mods);

new_mod = fo;
for i = 1:nmods    
    new_mod.mods(i) = setfield(fo.mods(i),regname,newval);
end
