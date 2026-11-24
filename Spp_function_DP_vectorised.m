clear

addpath("Flat Plat Models\");
addpath("Function Classes\")


w     = [10*(2*pi)];%,100*(2*pi),10000*(2*pi)];
theta = linspace(0, pi, 50);


Sqq = Source_Term_Sqq();
Dl  = Polar_Directivity();
Lte = Turbulent_Lift();

spp_result = zeros(length(theta),length(w));

tiledlayout(1,length(w));

p = polaraxes;

for a = 1:length(w)

    for i = 1:length(theta)

        spp_summation = 0;

        

            wl  = w(a) + (CommonConst.azi_mode * CommonConst.Omega);


            KXl = wl/CommonConst.Uc;

            spp_summation = spp_summation + ( (Dl.get(theta(i), alpha, w(a), CommonConst.azi_mode)) * (Sqq.get(KXl)) * (abs(double(Lte.get(theta(i), alpha, w(a),  CommonConst.azi_mode, KXl)))^2));

       

        spp_result(i,a) = (CommonConst.B / (2*pi)) * ((CommonConst.get_k0(w(a)) * CommonConst.b) / (CommonConst.r0))^2 * CommonConst.deltar * spp_summation;

    end

    nexttile
    polarplot(theta,10*log10(2*(4*pi)*spp_result(:,a))+94);
end
