% From McCoy1974

imp = 'Fortran';

projectile = 'BRL-1';
switch projectile
    case 'BRL-1'
        DREF    =   5.7;	% Bullet diameter (mm)
        LT      =   5.480;	% Bullet length (calibers)
        LN      =   3.000;  % Nose length (calibers)
        RTR     =   0.50;	% Headshape parameter. 0 for cone, 1 for tangent ogive nose.
        LBT     =   1.000;	% Boat-tail length (calibers)
        DB      =   0.754; 	% Base diameter (calibers)
        DM      =   0.000; 	% Meplat diameter (calibers)
        DBND   	=   1.000; 	% Rotating band diameter (calibers)
        XCG   	=   3.34;   % Center of gravity (calibers from nose)
        BLC     =   "L/T"; 	% "L/L" or "L/T" or "T/T"
        IDENT   =   "BRL-1"; % Projectile Name
    case 'MIN/MAN'
        DREF    =   55.6;	% Bullet diameter (mm)
        LT      =   3.250;	% Bullet length (calibers)
        LN      =   0.967;  % Nose length (calibers)
        RTR     =   0.00;	% Headshape parameter. 0 for cone, 1 for tangent ogive nose.
        LBT     =   1.180;	% Boat-tail length (calibers)
        DB      =   1.630; 	% Base diameter (calibers)
        DM      =   0.200; 	% Meplat diameter (calibers)
        DBND   	=   1.000; 	% Rotating band diameter (calibers)
        XCG   	=   1.76;   % Center of gravity (calibers from nose)
        BLC     =   "T/T"; 	% "L/L" or "L/T" or "T/T"
        IDENT   =   "MIN/MAN"; % Projectile Name
    case 'M549'
        DREF    =   155;	% Bullet diameter (mm)
        LT      =   5.65;	% Bullet length (calibers)
        LN      =   3.010;  % Nose length (calibers)
        RTR     =   0.5;	% Headshape parameter. 0 for cone, 1 for tangent ogive nose.
        LBT     =   0.58;	% Boat-tail length (calibers)
        DB      =   0.848; 	% Base diameter (calibers)
        DM      =   0.09; 	% Meplat diameter (calibers)
        DBND   	=   1.020; 	% Rotating band diameter (calibers)
        XCG   	=   3.53;   % Center of gravity (calibers from nose)
        BLC     =   "T/T"; 	% "L/L" or "L/T" or "T/T"
        IDENT   =   "M549"; % Projectile Name
end

M = [0.5:0.1:0.8 0.85 0.9 0.925 0.95 0.975 1:0.1:1.8 2 2.2 2.5:0.5:4].';   % Free stream Mach number
CDH_Fortran = zeros(size(M));
CDSF_Fortran = zeros(size(M));
CDBND_Fortran = zeros(size(M));
CDBT_Fortran = zeros(size(M));
CDB_Fortran = zeros(size(M));
CDO_Fortran = zeros(size(M));
PBP1_Fortran = zeros(size(M));

fprintf('%s \t %s \t %s \t %s \t %s \t %s \t %s \t %s\n', 'M', 'CDO', 'CDH', 'CDSF', 'CDBND', 'CDBT', 'CDB', 'PB/P1');

if(strcmp(imp,'Fortran'))
    % From Fortran:
    for ii=1:length(M)
        TA = (1-DM)/LN;
        M2 = M(ii)^2;
        
        RE = 23296.3*M(ii)*LT*DREF;
        RET = 0.4343 * log(RE);
        CFT = (0.455/(RET^2.58))*((1 + 0.21*M2)^(-0.32));
        DUM = 1 + ((0.333 + (0.02/(LN^2)))*RTR);
        SWM = 1.5708*LN*DUM*(1 + 1/(8*LN^2));
        SWCYL = 3.1416 * (LT-LN);
        SW = SWM+SWCYL;
        if(strcmp(BLC,'L/T'))
            CFL=(1.328/(sqrt(RE))) * ((1+0.12*M2)^(-0.12));
        end
        if(strcmp(BLC, 'T/T'))
            CFL = CFT;
        end
        
        CDSFL = 1.2732*SW*CFL;
        CDSFT = 1.2732*SW*CFT;
        CDSF_Fortran(ii) = (CDSFL * SWM + CDSFT * SWCYL)/SW;
        
        CHI = (M2-1)/(2.4*M2);
        
        if(M(ii)<=1)
            PTP = (1+0.2*M2)^3.5;
        end
        if(M(ii) > 1)
            PTP = ((1.2*M2)^3.5)*((6/(7*M2-1))^2.5);
        end
        CMEP = (1.122*(PTP-1)*(DM^2))/M2;
        
        if(M(ii) <= 0.91)
            CDHM = 0;
        end
        if(M(ii) >= 1.41)
            CDHM = 0.85 * CMEP;
        end
        if(M(ii) > 0.91 && M(ii) < 1.41)
            CDHM = (0.254+2.88*CHI)*CMEP;
        end
        
        if(M(ii) < 1)
            PB2 = 1/(1 + 0.1875*M2 + 0.0531*M2^2);
        end
        if(M(ii) >= 1)
            PB2 = 1/(1 + 0.2477*M2 + 0.0345*M2^2);
        end
        
        PB4 = (1 + 0.09*M2*(1-exp(LN-LT)))*(1+ 0.25*M2*(1-DB));
        PBP1_Fortran(ii) = PB2*PB4;
        CDB_Fortran(ii) = (1.4286*(1-PBP1_Fortran(ii))*DB^2)/M2;
        
        if(M(ii)<0.95)
            CDBND_Fortran(ii) = (M(ii)^12.5)*(DBND-1);
        end
        if(M(ii) >= 0.95)
            CDBND_Fortran(ii) = (0.21+0.28/M2)*(DBND-1);
        end
        
        if(M(ii)<=1)
            % Subsonic and transonic speeds
            if(LBT <=0)
                CDBT_Fortran(ii) = 0;
            else
                if(M(ii) <= 0.85)
                    CDBT_Fortran(ii) = 0;
                else
                    TB = (1-DB)/(2*LBT);
                    TB23 = 2*TB^2 + TB^3;
                    EBT = exp(-2*LBT);
                    BBT = 1- EBT +2*TB*((EBT*(LBT+0.5)) - 0.5);
                    CDBT_Fortran(ii) = 2*TB23 *BBT*(1/(0.564 + 1250*CHI^2));
                end
            end
            XMC =  1/sqrt(1+ 0.552*TA^0.8);
            if(M(ii) <= XMC)
                CDHT = 0;
            else
                CDHT = 0.368 * (TA^1.8) + 1.6 * TA * CHI;
            end
        else
            % Supersonic speeds
            BE2 = M2 - 1;
            BE = sqrt(BE2);
            ZE = BE;
            SSMC = 1+ 0.368 * (TA^1.85);
            if(M(ii) < SSMC)
                ZE = sqrt(SSMC^2-1);
            end
            
            C1 = 0.7156 - 0.5313*RTR + 0.595*RTR^2;
            C2 = 0.0796+0.0779*RTR;
            C3 = 1.587 + 0.049*RTR;
            C4 = 0.1122 + 0.1658*RTR;
            RZ2 = 1./ZE^2;
            
            CDHT = (C1 - C2*TA^2)*RZ2*((TA*ZE)^(C3+C4*TA));
            
            if(LBT > 0)
                TB = (1-DB)/(2*LBT);
                if(M(ii)<=1.1)
                    TB23=2*TB^2+TB^3;
                    EBT=exp(-2*LBT);
                    BBT=1 - EBT + 2*TB*((EBT*(LBT+0.5))-0.5);
                    CDBT_Fortran(ii) = 2*TB23*BBT*(1.774-9.3*CHI);
                end
                if(M(ii) > 1.1)
                    BB = 0.85 / BE;
                    AA2 = (5 * TA)/(6*BE) + (0.5*TA)^2 - (0.7435 / M2) * ((TA*M(ii))^1.6);
                    AA1 = (1 - ((0.6*RTR)/ M(ii))) * AA2;
                    EXL = exp(((-1.1952)/ M(ii)) * (LT - LN -LBT));
                    XXM = ((2.4 * M2^2 - 4*BE2) * (TB^2)) / (2*BE2^2);
                    AA = AA1 * EXL - XXM + ((2*TB)/ BE);
                    RB = 1/BB;
                    EXBT = exp((-BB)*LBT);
                    AAB = 1 - EXBT + (2*TB * (EXBT * (LBT + RB) - RB));
                    CDBT_Fortran(ii) = 4 * AA * TB * AAB * RB;
                end
            else
                CDBT_Fortran(ii) = 0;
            end
        end
        
        CDH_Fortran(ii) = CDHT + CDHM;
        
        CDO_Fortran(ii) = CDH_Fortran(ii) + CDSF_Fortran(ii) + CDBND_Fortran(ii) + CDBT_Fortran(ii) + CDB_Fortran(ii);
        
        fprintf('%3.3f \t %3.3f \t %3.3f \t %3.3f \t %3.3f \t %3.3f \t %3.3f \t %3.3f\n', M(ii), CDO_Fortran(ii), CDH_Fortran(ii), CDSF_Fortran(ii), CDBND_Fortran(ii), CDBT_Fortran(ii), CDB_Fortran(ii), PBP1_Fortran(ii));
    end
end

plot(M, [CDO_Fortran CDH_Fortran CDSF_Fortran CDB_Fortran CDBND_Fortran CDBT_Fortran CDB_Fortran]);
legend('CDO', 'CDH', 'CDSF', 'CDB', 'CDBND', 'CDBT', 'CDB');
grid on;
xlabel('Mach Number');
ylabel('Drag Coefficients');
title(IDENT);