function znngamma=znn_gamma(nerrZnn,nerrEC)
global gamma0
fis_a=fis2();
fis_out=evalfis([nerrZnn,nerrEC],fis_a);
v=norm(fis_out,2);
znngamma = gamma0 +v;