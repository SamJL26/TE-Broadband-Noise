classdef Polar_Directivity

    methods
        function obj = Polar_Directivity()
        end

        function D = get_all(obj, theta, w, l)
            % theta : vector of observer angles [rad]
            % w     : observer frequency [rad/s] (scalar)
            % l     : azimuthal modes (vector, e.g. CommonConst.azi_mode)
            % D     : Nl x Nt matrix, D(m,c) = D_l(m)(theta(c))

            alpha = CommonConst.alpha;
            k0    = CommonConst.get_k0(w);
            dr    = CommonConst.deltar;
            rmin  = CommonConst.rbar - dr/2;
            rmax  = CommonConst.rbar + dr/2;

            theta = theta(:).';
            l     = l(:);
            la    = abs(l);

            % theta and pi-theta share sin(theta), so they share the Bessel work
            thm = min(theta, pi - theta);
            [ut, ~, ic] = uniquetol(thm, 1e-12);
            ic = ic(:).';

            D = zeros(numel(l), numel(theta));

            for j = 1:numel(ut)
                s    = sin(ut(j));                          % scalar
                xm   = k0*rmax*s;
                lcut = ceil(xm + 10*xm^(1/3) + 10);         % modes beyond this are negligible
                if lcut > max(la)
                    warning('azi_mode too small at %.1f deg: need |l| up to %d', ...
                            rad2deg(ut(j)), lcut);
                end
                idx = find(la <= lcut);
                lj  = unique(la(idx));                      % distinct |l| >= 0

                % Gauss-Legendre nodes over the strip
                N       = ceil(2.2*k0*dr*s/pi) + 20;        % scalar
                [x, wq] = Polar_Directivity.gauss_legendre(N);
                r  = 0.5*(rmax - rmin)*x + 0.5*(rmax + rmin);   % N x 1
                wq = 0.5*(rmax - rmin)*wq;                      % N x 1

                % One Bessel matrix per unique angle: numel(lj) x N
                [LL, XX] = ndgrid(lj, k0*s*r(:).');
                J2 = besselj(LL, XX).^2;

                % Three moments of J_l^2 over the strip
                I0 = J2*wq;                 % integral of J^2
                I1 = J2*(wq ./ r);          % integral of J^2 / r
                I2 = J2*(wq ./ r.^2);       % integral of J^2 / r^2

                [~, loc] = ismember(la(idx), lj);
                lv = l(idx);

                % Assemble D_l for every angle that shares this sin(theta)
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