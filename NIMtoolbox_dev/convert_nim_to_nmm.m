function nmm = convert_nim_to_nmm(nim)

nmm = nim;
nmods = length(nim.mods);
for ii = 1:nmods
  nmm.mods(ii).Xtarget = 1;
end