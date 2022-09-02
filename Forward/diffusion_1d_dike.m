function [x,time,theta] = diffusion_1d_dike(Pm)
%1D diffusion model of host rock heating by dike
%used centered in space finite differences for spatial
%derivatives, ode45 matlab solver for time evolution

%written by Leif Karlstrom, October 2017
%TotT,DikeT,BackgndT,L,b,RightBC
%inputs: structure Pm with fields
%TotT - total time of simulation in years
%DikeT - total time the dike is active (total time to keep boundary
%BackgndT - farfield temperature in C (would want to calculate this based on depth in geotherm!!)
%condition at left side equal to dike temperature
%L - total length of domain in meters

%outputs
%x - distance vector in meters
%time - time vector in years
%theta - array of temperature in C at each spatial and time grid point
%f - array of host rock melt fraction at each node

%update, November 2017
%time varying dike temperature and seperate melt-fraction temp relations
%for dike, specifying half-width
%flag for non-uniform grid, to improve efficiency. Arctan
%coordinate transform from Erickson et al 2017 concentrates points near
%dike

%update, December 2017
%simplified for efficiency - removed variable coefficients
%added melt fraction temperature relations from Petcoviv and Dufek 2005
%no longer outputting melt fraction explicitly

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define parameters
year=Pm.year;%year = 3600*24*365;  %one year in seconds
T=Pm.TotT*year; %total time of simulation

Pm.N = Pm.L/Pm.h ; %number of nodes = domain length in km/spacing

if Pm.dononuniform
    dZ=1/(Pm.N);
    dikeThickZ=atan(Pm.DikeThick/Pm.ell)/atan(Pm.L/Pm.ell);
    dikepts = 1:round(dikeThickZ/dZ);
    Z = 0:dZ:1; %distance vector in scaled coordinates
    x = Pm.ell*tan(atan(Pm.L/Pm.ell)*Z);
else
    dikepts = 1:round(Pm.DikeThick/Pm.h);
    x = 0:Pm.h:Pm.L; %distance vector
end
    
%calculate meltfraction, only using tonalite parameterization now
[~,~,~,Pm.dftdT,Pm.Teval]=meltfraction(Pm);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initial condition
T0=ones(length(x),1)*Pm.Tbackground;

T0(dikepts)=Pm.Tliqd;

%solve the problem
TSPAN=[0 T];
%change BOTH RelTol and AbsTol if you want to modify runtime of ODE solver,
%e.g., down to 1e-3

options = odeset('RelTol',5e-6,'AbsTol',5e-6);
sol=ode45(@(t,y) RHS(t,y,Pm),TSPAN,T0,options); 

%process the output
theta=sol.y;
time=sol.x;

MaxT=max(theta,[],2);
if floor(max(MaxT))>Pm.Tliqd+1
    %disp(t/Pm.year)
    keyboard
end

%keyboard
%f=zeros(size(theta)); 
% for i=1:length(time)
% meltcr=find(theta(dikepts(end)+1:end,i)>Pm.Tsol&theta(dikepts(end)+1:end,i)<=Pm.Tliq);
% meltdk=find(theta(1:dikepts(end),i)>Pm.Tsold&theta(1:dikepts(end),i)<=Pm.Tliqd);
%  if sum(isfinite(meltcr))>0
%  f(meltcr,i) = ((theta(meltcr,i)-Pm.Tsol)./(Pm.Tliq-Pm.Tsol)).^Pm.b;
%  end
%  if sum(isfinite(meltdk))>0
%  f(meltdk,i) = ((theta(meltdk,i)-Pm.Tsold)./(Pm.Tliqd-Pm.Tsold)).^Pm.bd;      
%  end
%  if sum(find(f(:,i)>1))>0
%      f(f(:,i)>1,i)=1;
%  end
% end
%keyboard
end

function theta1=RHS(t,theta0,Pm)
%RHS for ode solver
    if Pm.dononuniform
        dZ=1/(Pm.N);
        dikeThickZ=atan(Pm.DikeThick/Pm.ell)/atan(Pm.L/Pm.ell);
        dikepts = 1:round(dikeThickZ/dZ);

        Z =0:dZ:1;Z=Z';%atan(x/Pm.ell)/atan(Pm.L/Pm.ell); %transformed coordinates

        %dZdx=Pm.ell./((Pm.ell^2+x.^2).*atan(Pm.L/Pm.ell));
        dxdZinv = 1./(Pm.ell*atan(Pm.L/Pm.ell).*sec(Z.*atan(Pm.L/Pm.ell)).^2);    
    else
        dikepts = 1:round(Pm.DikeThick/Pm.h);
    end


    N=Pm.N;
    T=theta0+273; %convert to Kelvin
    %temp-dependent material parameters from Whittington 2009
    %1200*ones(size(T));%

    cp = Pm.cp*ones(size(T));%(199.50+0.0857.*T-5*1e-6 .*T.^(-2)) / 0.22178 ; %specific heat capacity in J/(kg*K)
    %1e-6 *ones(size(T));%
    kappa = Pm.kappa *ones(size(T));%(567.3./T - 0.062)*1e-6; %thermal diffusivity in m^2/s
        
    dfdT=zeros(size(theta0));
    
    %define the grid points not associated with dike
    nondikepts = length(dikepts)+1:length(theta0); 
    
     meltcr=find(theta0(nondikepts)>Pm.Tsol);%&theta0(dikepts(end)+1:end)<=Pm.Tliq);
     meltdk=find(theta0(dikepts)>Pm.Tsold&theta0(dikepts)<=Pm.Tliqd);

 if sum(isfinite(meltcr))>0
     %keyboard
     %f(melt) = ((theta0(melt)-Pm.Tsol)./(Pm.Tliq-Pm.Tsol)).^Pm.b;
     dfdT(nondikepts(meltcr)) = interp1(Pm.Teval',Pm.dftdT',theta0(nondikepts(meltcr)));
    %Pm.b*((theta0(meltcr)-Pm.Tsol)./(Pm.Tliq-Pm.Tsol)).^(Pm.b-1);

     %meltcrS=find(theta0(dikepts(end)+1:end)>Pm.Tliq);
     dfdT(theta0(dikepts(end)+1:end)>Pm.Tliq|theta0(dikepts(end)+1:end)<Pm.Tsol) = 0;
 end
 if sum(isfinite(meltdk))>0
     dfdT(dikepts(meltdk)) = Pm.bd*((theta0(meltdk)-Pm.Tsold)./(Pm.Tliqd-Pm.Tsold)).^(Pm.bd-1); 
     %= interp1q(Pm.Teval',Pm.dfbdT',theta0(meltdk)); 
 end
 %zero out any points that are larger than basalt liquidus
 dfdT(isnan(dfdT)) = 0;

    cpEff=cp; 
     %effective heat capacity to account for Latent heat
    cpEff(dikepts(meltdk))=cpEff(dikepts(meltdk)) + Pm.Lf.*dfdT(dikepts(meltdk))/(Pm.Tliqd-Pm.Tsold);
    cpEff(nondikepts(meltcr))=cpEff(nondikepts(meltcr)) + Pm.Lf.*dfdT(nondikepts(meltcr))/(Pm.Tliq-Pm.Tsol) ;

     if sum(imag(cpEff))>0 || sum(isnan(cpEff))>0
         disp('something wrong with cpEff')
         disp(cpEff)
         disp(isnan(cpEff))
         error('ending simulation')
         
     end
%      disp(sum(meltcr))
%      disp(sum(meltdk))
%      if t/Pm.year>Pm.tcenter
%          keyboard
%      end
     
    k = Pm.k*ones(size(kappa)); %thermal conductivity, W/m K
    dkdx = zeros(size(cp));%gradient(k,Pm.h);

if Pm.dononuniform
    alpha = dxdZinv./(Pm.rho.*cpEff); %coefficients in transformed coordinates
    beta = dxdZinv.*k; %coefficients in transformed coordinates
end
    Flux = (1-tanh(10*(t/Pm.year - Pm.tcenter)/Pm.twidth))/2;
    %finite difference operators
    theta1=zeros(N+1,1);
    
if Flux>.01 %t<Pm.tdike
%time varying temperature within dike (assume uniform for now)
    if theta0(1)>theta0(length(dikepts)+1)
        theta1(dikepts)=dT0dt(t,theta0,Pm);
    else
        theta1(dikepts)=0;
    end

    %0; %constant temperature LHS
%interior nodes
    for i=length(dikepts)+1:N
        if Pm.dononuniform
            theta1(i) = alpha(i)/(2*dZ^2) *((beta(i+1)+beta(i))*theta0(i+1) - ...
                (beta(i+1)+2*beta(i)+beta(i-1))*theta0(i) + (beta(i)+beta(i-1))*theta0(i-1));
        else
            theta1(i)= dkdx(i)/(Pm.rho* cpEff(i)) * ((theta0(i+1)-theta0(i-1))/(2*Pm.h)) + ...
               k(i)/(Pm.rho *cpEff(i)) * (theta0(i+1)-2*theta0(i)+theta0(i-1))/Pm.h^2 ;
        end
    end
        if Pm.RightBC==1%far-field BC is insulating
            theta1(N+1) = 2*beta(N+1)*alpha(N+1) * (theta0(N)-theta0(N+1))/dZ^2 ; 
            %theta1(N+1) = 2*k(N+1)/(Pm.rho*cpEff(N+1)) * (theta0(N)-theta0(N+1))/Pm.h^2 ;
        elseif Pm.RightBC==2 %far-field BC is constant temp
            theta1(N+1) = 0;
        end
        
else
%if dike is no longer providing constant temperature, need to change the LHS boundary condition    
    if Pm.dononuniform
        theta1(1) = 2*beta(1)*alpha(1) * (theta0(2)-theta0(1))/dZ^2 ;    
    else
        theta1(1) = 2*k(1)/(Pm.rho*cpEff(1)) * (theta0(2)-theta0(1))/Pm.h^2 ;
    end
%interior nodes
    for i=2:N
        if Pm.dononuniform
            theta1(i) = alpha(i)/(2*dZ^2) *((beta(i+1)+beta(i))*theta0(i+1) - ...
                (beta(i+1)+2*beta(i)+beta(i-1))*theta0(i) + (beta(i)+beta(i-1))*theta0(i-1));
        else
            theta1(i)= dkdx(i)/(Pm.rho*cpEff(i)) * ((theta0(i+1)-theta0(i-1))/(2*Pm.h)) + ...
                k(i)/(Pm.rho*cpEff(i)) * (theta0(i+1)-2*theta0(i)+theta0(i-1))/Pm.h^2;
        end
    end
        if Pm.RightBC==1 %far-field BC is insulating
            theta1(N+1) = 2*beta(N+1)*alpha(N+1) * (theta0(N)-theta0(N+1))/dZ^2 ;    
        %theta1(N+1) = 2*k(N+1)/(Pm.rho*cpEff(N+1)) * (theta0(N)-theta0(N+1))/Pm.h^2 ;
        elseif Pm.RightBC==2 %far-field BC is constant temp
            theta1(N+1) = 0;
        end
        
    
end


end

function rhsDT0Dt = dT0dt(t,~,Pm)
%calculate time derivative of boundary temperature, in this case with form
%T = (TDike-Tbkd)*cos(2*pi*t/T)*(1-tanh((10*(t-b)/c))/2 + Tbkd
%only keeping positive part
%         if cos(2*pi*t/Pm.tperiod)>0
%     rhsDT0Dt=Pm.tperiod.^(-1).*Pm.twidth.^(-1).*(Pm.Tbackground+(-1).*Pm.TDike).*(5.*a.*cos(2.*a.^(-1).*pi.*t).* ...
%   sech(10.*Pm.twidth.^(-1).*((-1).*Pm.tcenter+t)).^2+(-1).*Pm.twidth.*pi.*sin(2.*Pm.tperiod.^(-1).* ...
%   pi.*t).*((-1)+tanh(10.*Pm.twidth.^(-1).*((-1).*Pm.tcenter+t))));
%         else
%     rhsDT0Dt=0;
%         end
        
%just tanh part
    rhsDT0Dt= -5.*Pm.twidth.^(-1).*Pm.year.^(-1).*(Pm.Tliqd-Pm.Tbackground).*sech(10.*Pm.twidth.^(-1).*((-1).*Pm.tcenter+t/Pm.year)).^2;  
end