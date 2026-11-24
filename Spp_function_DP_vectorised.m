clear clc

addpath("Flat Plat Models\"); 
addpath("Function Classes\")

w     = 100*(2*pi);
theta = linspace(0, pi, 50);

Sqq = Source_Term_Sqq();
Dl  = Polar_Directivity();
Lte = Turbulent_Lift();

l   = CommonConst.azi_mode(:);
wl  = w + l*CommonConst.Omega;
KXl = wl/CommonConst.Uc;
k0  = CommonConst.get_k0(w);

tic
B = Sqq.get(KXl);                         % once
D = Dl.get_all(theta, w, l);              % once, all angles at once
fprintf('D and B: %.1f s\n', toc);

spp_result = zeros(length(theta),1);
for i = 1:length(theta)
    C = abs(Lte.get(theta(i), w, KXl)).^2;
    spp_result(i) = (CommonConst.B/(2*pi)) * ((k0*CommonConst.b)/CommonConst.r0)^2 ...
                    * CommonConst.deltar * sum(D(:,i) .* B(:) .* C(:));
end

polarplot(theta, 10*log10(2*(4*pi)*spp_result) + 94);
thetalim([0 180]);