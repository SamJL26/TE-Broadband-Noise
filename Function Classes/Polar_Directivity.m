classdef Polar_Directivity

    properties
       
    end

    methods
        function obj = Polar_Directivity()
            
        end

        function result = get(obj, theta, alpha, w, azi_mode)
            result = (1 / CommonConst.deltar) * integral(obj.integral_function(theta, alpha, w, azi_mode), CommonConst.rbar - CommonConst.deltar/2, CommonConst.rbar + CommonConst.deltar/2);
        end

        function intfunc = integral_function(obj, theta, alpha, w, azi_mode)
            intfunc    = @(r) ( (azi_mode ./(CommonConst.get_k0(w) .* r)) .* cos(alpha) + cos(theta) .* sin(alpha) ).^2 .* (obj.Jn(CommonConst.get_k0(w) * r * sin(theta),azi_mode)) .^ 2;
        end

        function Jn = Jn(~, bessel_arg,azi_mode)
            Jn = besselj(azi_mode, bessel_arg);
        end

    end
end