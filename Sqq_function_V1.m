clear clc
function Sqq = Sqq_function_V1(w)

    function lr = lr(w,kr)
    MX = 0.6;
    Mc = 0.8*MX;
    Uc = 343*Mc;
    xi2 = 1.6;
    l2 = xi2*Uc/w;
    lr = (l2/pi)*(1/(1+l2^2*kr^2));
    
    end
    
    MX = 0.6;
    Mc = 0.8*MX;
    UX = 343*MX;
    Uc = 343*Mc;
    kX = 10;
    kr = 0;
    
    
    function phi_pp = phi_pp(w)
    MX = 0.6;
    UX = 343*MX;
    rho0 = 1.225;
    delta_star = 1;
    w_bar = w .* delta_star ./ UX;
    phi_bar = (2e-5) ./(1 + w_bar + 0.217.*w_bar.^2 + 0.00562.*w_bar.^4);
    phi_pp = phi_bar .* ((0.5*rho0*UX^2)^2 .* delta_star ./ UX);
    
    end
    
    
    Sqq = (1/pi)*lr(kX*Uc,kr)*phi_pp(w);

end