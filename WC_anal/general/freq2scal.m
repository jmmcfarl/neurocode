function scales = freq2scal(freqs,wname,delta)

Fc = centfrq(wname);
scales = Fc./(freqs*delta);