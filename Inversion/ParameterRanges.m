function PmRg = ParameterRanges%(Pm)

%parameters and ranges varied during MCMC.
%note that this assumes uniform distributions (only specifying end points
%for prior)


PmRg(1,1)=25; PmRg(1,2)=125; %Tbg, background T range

PmRg(2,1)=.03; PmRg(2,2)=10; %Tw, dike shut off ramp duration scale (yrs, although not directly)

PmRg(3,1)=0.1; PmRg(3,2)=10; %Tc, dike longevity (yrs)

PmRg(4,1)=2; PmRg(4,2)=10; %wall rock thermal conductivity

%if Pm.DataFlag == 1|| Pm.DataFlag == 3|| Pm.DataFlag == 4 || Pm.DataFlag == 6

PmRg(5,1)=120; PmRg(5,2)=145; %EaA, activation energy Ap, kJ/mol
%end
%if Pm.DataFlag == 1|| Pm.DataFlag == 3|| Pm.DataFlag == 4 || Pm.DataFlag == 5

PmRg(6,1)=160; PmRg(6,2)=175; %EaZ, activation energy Zr, kJ/mol
%end

%if Pm.DataFlag == 3|| Pm.DataFlag == 4 || Pm.DataFlag == 5
PmRg(7,1)=190; PmRg(7,2)=250; %activation energy Bt, kJ/mol
%end
%if Pm.DataFlag == 3|| Pm.DataFlag == 4 || Pm.DataFlag == 6
PmRg(8,1)=.8; PmRg(8,2)=.9; %rmr0 Ap fission track

    %PmRg(9,1)=0.1; PmRg(9,2)=1; %rmr0 Zr fission track
%end



