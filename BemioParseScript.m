PlateBEM.FexcRe = h5read('D:\Mitacs Globalink 2022\Sys-MoDEL\Modelica Codes\rm3.h5','/body1/hydro_coeffs/excitation/components/re/3_1');
PlateBEM.FexcRe = PlateBEM.FexcRe(1,:);
PlateBEM.FexcIm = h5read('D:\Mitacs Globalink 2022\Sys-MoDEL\Modelica Codes\rm3.h5','/body1/hydro_coeffs/excitation/components/im/3_1');
PlateBEM.FexcIm = PlateBEM.FexcIm(1,:);
PlateBEM.w = h5read('D:\Mitacs Globalink 2022\Sys-MoDEL\Modelica Codes\rm3.h5','/simulation_parameters/w');
save('PlateBEM.mat')