clear

addpath("Flat Plat Models\");
addpath("Function Classes\")


w     = [10000*(2*pi)];%,100*(2*pi),10000*(2*pi)];
theta = linspace(0, pi, 30);


Sqq = Source_Term_Sqq();
Dl  = Polar_Directivity();
Lte = Turbulent_Lift();

spp_result = zeros(length(theta),length(w));

for i = 1:length(theta)
        tic
        wl  = w + (CommonConst.azi_mode .* CommonConst.Omega);


        KXl = wl/CommonConst.Uc;


        A = Dl.get(theta(i), w);
        
        B = Sqq.get(KXl);

        C = (abs(Lte.get(theta(i), w, KXl))).^2;


        product = A .* B(:) .* C(:);
        spp_summation = sum(product);
        

        spp_result(i) = (CommonConst.B / (2*pi)) * ((CommonConst.get_k0(w) * CommonConst.b) / (CommonConst.r0))^2 * CommonConst.deltar * spp_summation;
        toc
end


polarplot(theta,10*log10(2*(4*pi)*spp_result(:))+94);


