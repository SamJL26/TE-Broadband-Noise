function phi_pp = Goody(w)
    tau_w = CommonConst.Cf * 0.5 * CommonConst.rho0 * CommonConst.UX^2;          % wall shear stress
    u_tau = sqrt(tau_w / CommonConst.rho0);              % friction velocity
    RT    = (CommonConst.delta/CommonConst.UX) / (CommonConst.nu/u_tau^2);       % = (u_tau*delta/nu)*sqrt(Cf/2)
    
    wb = abs(w) * CommonConst.delta / CommonConst.UX;                % dimensionless frequency
    
    num = 3.0 * wb.^2;
    den = (wb.^0.75 + 0.5).^3.7 + (1.1 * RT^(-0.57) * wb).^7;
    
    phi_pp = 0.5 * (tau_w^2 * CommonConst.delta / CommonConst.UX) .* num ./ den;
end