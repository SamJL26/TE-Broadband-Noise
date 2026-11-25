classdef Polar_Directivity

    methods
        function obj = Polar_Directivity()
        end

        function D = get_all(obj, theta, w, l)

            alpha = CommonConst.alpha;
            k0    = CommonConst.get_k0(w);
            dr    = CommonConst.deltar;
            rmin  = CommonConst.rbar - dr/2;
            rmax  = CommonConst.rbar + dr/2;

            theta = theta(:).';
            l     = l(:);
            l_abs = abs(l);

            thm         = min(theta, pi - theta);
            [ut, ~, ic] = uniquetol(thm, 1e-12);
            ic          = ic(:).';

            D = zeros(numel(l), numel(theta));

            for j = 1:numel(ut)
                s    = sin(ut(j));                          
                xm   = k0*rmax*s;
                lcut = ceil(xm + 10*xm^(1/3) + 10);         % modes beyond this are negligible
               
                if lcut > max(l_abs)
                    warning('azi_mode too small at %.1f deg: need |l| up to %d', ...
                            rad2deg(ut(j)), lcut);
                end

                idx = find(l_abs <= lcut);
                lj  = unique(l_abs(idx));                      

                % Gauss-Legendre nodes over the strip
                N       = ceil(2.2*k0*dr*s/pi) + 20;        
                [x, wq] = Polar_Directivity.gauss_legendre(N);
                
                r  = 0.5*(rmax - rmin)*x + 0.5*(rmax + rmin);   
                wq = 0.5*(rmax - rmin)*wq;                     

                
                [LL, XX] = ndgrid(lj, k0*s*r(:).');
                J2 = besselj(LL, XX).^2;

                I0 = J2*wq;                 
                I1 = J2*(wq ./ r);          
                I2 = J2*(wq ./ r.^2);       

                [~, loc] = ismember(l_abs(idx), lj);
                lv = l(idx);


                for c = find(ic == j)
                    cth = cos(theta(c));
                    D(idx, c) = ( cos(alpha)^2 * lv.^2 .* I2(loc) / k0^2 ...
                                + 2*cos(alpha)*sin(alpha)*cth * lv .* I1(loc) / k0 ...
                                + sin(alpha)^2 * cth^2 * I0(loc) ) / dr;
                end
            end
        end
    end

    methods (Static)
        function [x, wq] = gauss_legendre(N)
            a    = 1:N-1;
            beta = a ./ sqrt(4*a.^2 - 1);
            [V, E] = eig(diag(beta,1) + diag(beta,-1));
            [x, k] = sort(diag(E));
            wq = 2*(V(1,k).').^2;
        end
    end
end