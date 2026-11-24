classdef Polar_Directivity

    properties
       
    end

    methods
        function obj = Polar_Directivity()
            
        end

        function result = get(obj, theta, w)
            result = (1 / CommonConst.deltar) * integral(obj.integral_function(theta, w), ...
                      CommonConst.rbar - CommonConst.deltar/2, ...
                      CommonConst.rbar + CommonConst.deltar/2, ...
                      'ArrayValued', true);
        end

        function intfunc = integral_function(obj, theta, w)
            azi_cof = @(r) cos(CommonConst.alpha) ./ (CommonConst.get_k0(w) .* r);
            azi_dep = @(r) (azi_cof(r) .* CommonConst.azi_mode(:) + cos(theta)*sin(CommonConst.alpha)).^2;
            bessel  = @(r) (obj.Jn(CommonConst.get_k0(w) * r * sin(theta))).^2;

            intfunc = @(r) (azi_dep(r) .* bessel(r));
        end

        function Jn = Jn(~, bessel_arg)
            Jn = besselj(CommonConst.azi_mode(:), bessel_arg);
        end

    end
end