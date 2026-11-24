clear

addpath("Flat Plat Models\");
addpath("Function Classes\")

w     =  linspace(100,10000,100);
theta = [pi/6,4*pi/6,5*pi/6];
alpha = 4*pi/9;

Sqq = Source_Term_Sqq();
Dl  = Polar_Directivity();
Lte = Turbulent_Lift();

spp_result = zeros(1,length(w));

for a = 1:length(theta)
    
    for i = 1:length(w)
    
        spp_summation = 0;
    
        for L = 1: length(CommonConst.azi_mode)
           
           wl  = w(i) + (CommonConst.azi_mode(L) * CommonConst.Omega);

           KXl = wl/CommonConst.Uc;
         
           spp_summation = spp_summation + ( (Dl.get(theta(a), alpha, w(i), CommonConst.azi_mode(L))) * (Sqq.get(KXl)) * (abs(double(Lte.get(theta(a), alpha, wl,  CommonConst.azi_mode(L), KXl)))^2));
        end
    
        spp_result(i) = (CommonConst.B / (2*pi)) * ((CommonConst.get_k0(w(i)) * CommonConst.b) / (CommonConst.r0))^2 * CommonConst.deltar * spp_summation;
    
    end
    
    
    semilogx(w,10*log10(spp_result));
    hold on
end
