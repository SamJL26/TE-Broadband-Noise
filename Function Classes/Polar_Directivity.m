classdef Polar_Directivity

    properties
       
    end

    methods
        function obj = Polar_Directivity()
            
        end

        function result = get(obj, theta, w)
            result = (1 / CommonConst.deltar) * integral(obj.integral_function(theta, w), CommonConst.rbar - CommonConst.deltar/2, CommonConst.rbar + CommonConst.deltar/2);
        end

        function intfunc = integral_function(obj, theta, w)
            intfunc    = @(r) ( (CommonConst.azi_mode ./(CommonConst.get_k0(w) .* r)) .* cos(CommonConst.alpha) + cos(theta) .* sin(CommonConst.alpha) ).^2 .* (obj.Jn(CommonConst.get_k0(w) * r * sin(theta))) .^ 2;
        end

        function Jn = Jn(~, bessel_arg)
            Jn = besselj(CommonConst.azi_mode, bessel_arg);
        end

    end
end