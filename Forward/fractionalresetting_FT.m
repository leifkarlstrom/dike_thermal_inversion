function f = fractionalresetting_FT(T,t,Tm,rmr0,year)
%implements parameterization of the Ketcham (1999) fanning curvilinear model, 
%based on FTIndex c code by Ketcham

%Approach: given a t-T path, split up into isothermal intervals and compute
%an equivalent time in the final empirical equation for fractional track annealing
finegrid=true;
%first subdivide into equal time chunks
% maybe play with this to see whether we can use coarser time steps
% (smaller numbers numt)
if finegrid
    numt=5e4;

    teq = linspace(0,max(t),numt); %time in seconds
    dt = teq(2)-teq(1);
%disp(dt/year * 365)
%interpolate temperature onto this
    Teq = interp1(t,T,teq);

%%%%to compare with HeFTy
%AgeSt = 0.5e6*year; %startting model age (0.5 Ma) 
%note that "resetting" is only really set by the center time, which is AgeSt-350 yr 
%teqMa = AgeSt - teq;
%save('FT_tTime_finegrid.mat','teq','Teq','teqMa')
else
    teq = t; %time in seconds
    dt = teq(2)-teq(1);
%disp(dt/year * 365)
%interpolate temperature onto this
    Teq = T;
end

%%%check to see whether isothermal hold ("fading") is same as equivalent time
fading=false; %or true

if fading
    Tiso = 200; %in C
    Teq = Tiso*ones(size(teq));
end
%%%%%%%%%%%%%%

%now compute resetting
c0 = Tm.c0;
c1 = Tm.c1;
c2 = Tm.c2;
c3 = Tm.c3;

c4 = Tm.alpha; %alpha
c5 = Tm.beta; %beta

%rmr0 =  0.8204;%0.7815; %0.7507;

%NOTE: Ketcham 1999 also suggests that one can use %kappa = 1 - rmr0;
k = 1 - rmr0;%0.2753; %0.2803;

%initialize
Rhos = 0;
currLength =0.99999;

currSlope = (Teq(2)-Teq(1))/dt;
tempStep = currSlope * dt;
currTemp = Teq(1) + 273.15;

currPt = 1;
currTime = teq(1);

while currTime < teq(end) - dt
    
    %update avg temperature via local slope
    Tavg = currTemp + tempStep/2.0;
    currSlope = (Teq(currPt+1)-Teq(currPt))/dt;
    tempStep = currSlope * dt;
    currTemp = currTemp + tempStep;% Teq(currPt) + 273.15;   
    %above approach is equivalent to 
    %Tavg = (Teq(currPt+1)+Teq(currPt))/2 + 273.15;
    
    %iterative method (Ketcham approach, fanning curvilinear)
    equivTime = (currLength)^(1/k) *(1-rmr0) + rmr0;    
    equivTime = (((1.0-equivTime^c5)/c5)^c4 - 1.0)/c4 - c0  ;
    
    equivTime = equivTime/c1;    
    equivTime = exp(equivTime*(log(1.0/Tavg)-c3)+c2);
        
    timeInt = equivTime+dt; %increment time forward
    currTime = currTime + dt;    
    Tpt(currPt) = Tavg;
    
    %iterative method (Ketcham approach)
    currLength = 1.0+c4*(c0+c1*(log(timeInt)-c2)/(log(1.0/Tavg)-c3));
    currLength = 1.0-c5*(currLength)^(1/c4);
    
    
    if currLength < 0
        currLength = 0.0 ;
    else
        currLength = (currLength)^(1/c5);
    end
    
    if imag(((currLength - rmr0)/(1 - rmr0))^k)>0
        currLength = 0;
        break
    else
        currLength = ((currLength - rmr0)/(1 - rmr0))^k;
    end    
    %currTemp = (Teq(currPt+1)+Teq(currPt))/2; 
    
    currPt = currPt + 1;    
end

r = currLength;

%%% check to see whether isothermal hold is same as equivalent time
if fading 
        %purely isothermal 
    for j = 1:length(teq)    
        t = teq(j);
        T=Teq(j)+273.15;

        lambda = 1+c4*(c0+c1*(log(t)-c2)/(log(1/T)-c3));
        rstar = (1 - c5*lambda^(1/c4))^(1/c5);
        if imag(((rstar - rmr0)/(1-rmr0))^k)>0
            riso = 0;
            break
        else
            riso = ((rstar - rmr0)/(1-rmr0))^k;
        end
    end
        disp(['Residual: (Isothermal eqn - equiv time) = ' num2str(currLength-riso)])
        if riso == 0
            disp(['Fading time for T=' num2str(Tiso) ' is ' num2str(t/year) ' years'])
        else
            disp(['Fading time is greater than timevec of ' num2str(t/year) ' years'])
        end
disp('Returning without computing resetting')
%keyboard
return
end

%compute density from reduced length
if r>0.757
        Rhos = r*1.6-0.6;
elseif r>=0.53
        Rhos = 9.205*r^2 - 9.157*r + 2.269;
else
        Rhos = 0;
end

%compute reset fraction    
f = 1 - Rhos;

%disp('reduced length, density, reset fraction =')
%disp([r,Rhos,f])


