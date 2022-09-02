function [fbasalt,ftonalite,dfbdT,dftdT,T]=meltfraction(Pm)

T=700:1300;

% Pm.Tliq=1100; %tonalite liquidus (C)
% Pm.Tsol=725;
% Pm.Tliqd=1165; %basalt liquidus (C)
% Pm.Tsold=1015;

ftonalite=zeros(size(T));
fbasalt=zeros(size(T));

for i=1:length(T)
    if T(i)>=Pm.Tsol && T(i)<=Pm.Tliq
        %model tonalite as piecewise linear, based on fit to P&D2005 
        if T(i)<=800
            ftonalite(i) = (0.2/(800-Pm.Tsol))*(T(i)-Pm.Tsol);
        elseif T(i)> 800 && T(i) <= 895
            ftonalite(i) = 0.2;
        elseif T(i) > 895 && T(i) <= 950
            ftonalite(i) = 0.2+(0.9-0.2)/(950-895)*(T(i)-895);
        elseif T(i)>950
            ftonalite(i) = 0.9 + (1.0-0.9)/(Pm.Tliq-950) *(T(i)-950);
        end
    elseif T(i)>Pm.Tliq
        ftonalite(i)=1;
    end
    
    if T(i)>=Pm.Tsold && T(i)<=Pm.Tliqd
        %model basalt as powerlaw, based on fit to P&D2005 
        fbasalt(i) = ((T(i)-Pm.Tsold)/(Pm.Tliqd-Pm.Tsold))^Pm.bd;
    elseif T(i)>Pm.Tliqd
        fbasalt(i) = 1;
    end    
end

dfbdT=gradient(fbasalt,T(2)-T(1));
dftdT=gradient(ftonalite,T(2)-T(1));

dftdT(dftdT<0)=0;
dfbdT(dfbdT<0)=0;

% plot(T(ftonalite>0),ftonalite(ftonalite>0),T(ftonalite>0),dftdT(ftonalite>0))
% xlim([700 1200])


