function phi_pp = Willmarth_Roos_Amiet(w)
    w_bar   = w * CommonConst.delta_star / CommonConst.UX;
    phi_bar = (2e-5) / (1 + w_bar + 0.217 * w_bar.^2 + 0.00562 * w_bar.^4);
    phi_pp  = phi_bar * ((0.5*CommonConst.rho0 * CommonConst.UX^2)^2 * CommonConst.delta_star / CommonConst.UX);
end