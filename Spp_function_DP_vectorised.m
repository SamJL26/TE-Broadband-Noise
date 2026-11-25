clear clc

addpath("Flat Plat Models\"); 
addpath("Function Classes\")

w     = 10*(2*pi);
theta = linspace(0, pi, 100);

Sqq = Source_Term_Sqq();
Dl  = Polar_Directivity();
Lte = Turbulent_Lift();

l   = CommonConst.azi_mode(:);
wl  = w + l*CommonConst.Omega;
KXl = wl/CommonConst.Uc;
k0  = CommonConst.get_k0(w);


S = Sqq.get(KXl);                         
D = Dl.get_all(theta, w, l);  

spp_result = zeros(length(theta),1);

for i = 1:length(theta)
    L = abs(Lte.get(theta(i), w, KXl)).^2;
    spp_result(i) = (CommonConst.B/(2*pi)) * ((k0*CommonConst.b)/CommonConst.r0)^2 ...
                    * CommonConst.deltar * sum(D(:,i) .* S(:) .* L(:));
end

polarplot(theta, 10*log10(2*(4*pi)*spp_result) + 94);