classdef Turbulent_Lift
    methods
        function obj = Turbulent_Lift()
           
        end
        function result = get(obj, theta, w, SWWN)
            SWWN    = SWWN(:);
            l       = CommonConst.azi_mode(:);
            kappa_l = (l ./ CommonConst.rbar) * sin(CommonConst.alpha) ...
                      - CommonConst.get_k0(w) * cos(CommonConst.alpha) * cos(theta);
            tA     = obj.thetaA(SWWN);
            tB     = obj.thetaB(kappa_l, SWWN);
            tB_abs = obj.thetaB_abs(kappa_l, SWWN);
        
            result = (exp(2i.*tB) ./ (1i.*tB_abs)) .* ( ...
                        exp(-2i.*tB) .* sqrt(tA ./ (tA - tB)) .* erfz(sqrt(2i.*(tA - tB))) ...
                        - erfz(sqrt(2i.*tA)) + 1 );
        end

        function thetaA = thetaA(~, SWWN)
            thetaA = CommonConst.b * (SWWN + CommonConst.get_muainf(SWWN) + CommonConst.get_mua(SWWN) * CommonConst.MX);
        end

        function thetaB = thetaB(~, kappa, SWWN)
            thetaB = CommonConst.b * (SWWN + kappa);
        end

        function thetaB_abs = thetaB_abs(~, kappa, SWWN)
            thetaB_abs = CommonConst.b * (abs(SWWN) + abs(kappa));
        end
    end

    
end