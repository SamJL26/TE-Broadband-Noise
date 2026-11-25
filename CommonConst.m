classdef CommonConst

    properties (Constant)
        
        %{
          The following are set up parameters for each type of propeller
          used in the validation from VB's thesis. In the thesis, he uses
          values of delta and Cf from XFOIL calculations, the following are 
          chosen to best match his results across the frequency range.

            Wind Turbine: 
                        - delta   = 0.045
                        - Cf      = 0.0025
                        - r_blade = 29
                        - alpha   = 4pi/9
                        - omega   = 2.618
                        - b       = 1
                        - B       = 3
            
            Cooling Fan: 
                        - Boundary Layer = 0.0105
                        - Cf = 0.0016 
                        - r_blade = 0.4
                        - alpha   = 14pi/45
                        - omega   = 62.83
                        - b       = 0.065
                        - B       = 2
         
        %}

        


        rho0 = 1.225;   
        c0   = 343;
        nu   = 1.5*10^-5;


        Omega   = 62.83;
        B       = 2;
        alpha   = 14*pi/45;

        r0      = 1.0;
        r_blade = 0.4;
        rbar    = 0.75*CommonConst.r_blade;
        deltar  = CommonConst.rbar/3;

        b       = 0.065;


        UX = 0.75*CommonConst.Omega*CommonConst.r_blade;
        Uc = 0.8*CommonConst.UX; 
        
        MX = CommonConst.UX/CommonConst.c0;
        Mc = 0.8*CommonConst.MX;
        
        Cf = 0.0016;

        
        kr = 0.0;

        delta      = 0.0105;
        delta_star = CommonConst.delta/8;


        xi2 = 1.6;
        

        betaX  = sqrt(1 - CommonConst.MX^2);
        
        

       
        azi_mode = [-5000:5000];
        
    end

    methods(Static)
        function k0 = get_k0(w)
            k0 = w / CommonConst.c0;
        end
        
        function mua = get_mua(SWWN)
            mua = (SWWN * CommonConst.Mc) / (CommonConst.betaX)^2;
        end

        function muainf = get_muainf(SWWN)
            muainf = sqrt((CommonConst.get_mua(SWWN)).^2 - ...
                ((CommonConst.kr) / CommonConst.betaX)^2);
        end
    end
end