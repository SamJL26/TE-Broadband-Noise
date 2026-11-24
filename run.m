clear

addpath("Flat Plat Models\");
addpath("Function Classes\")

w     =  CommonConst.kX * CommonConst.Mc * CommonConst.c0;
theta = pi/4;
alpha = pi/4;


Sqq = Source_Term_Sqq();
Dl  = Polar_Directivity();
Lte = Turbulent_Lift();

sqq_result = Sqq.get(w);
dl_result  = Dl.get(theta, alpha, w, CommonConst.azi_mode(11));
lte_result = Lte.get(theta, alpha, w);



disp("Sqq Result: " + sqq_result);
disp("Dl Result: " + dl_result);
disp("Lte Result: " + double(lte_result));
