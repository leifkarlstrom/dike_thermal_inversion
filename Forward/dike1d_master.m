%     This script implements solves 1D unsteady heat equation implicitly, 
%     using stretched coordinates with finite differences.
%     Second part of the script computes chemical diffusion at series of query points to predict     
%     one or more of: (U-Th)/He diffusion, Bt/Ar diffusion, Ap Fission Track, 
%     Magnetic geothermometers
%
%
%     Copyright (C) July 2022  Leif Karlstrom leif@uoregon.edu
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%add data path
%%%% CHANGE THIS PATH TO WHATEVER IT IS LOCALLY
addpath('../../Data')

%% model parameters

%control parameters varied in inversions
Pm.Tbackground = 46;%82.6;  %background temperature in C
Pm.twidth = 2.9; %vary this parameter to change the rate at which dike temperatures decrease
Pm.tcenter = 3.3; %scale for dike shut-off time (years)
Pm.k = 5;%8.1;%5.5;%9.2074; %thermal conductivity (assumed same between materials)

%kinetic parameters
Tm.EaZr=167.48; %167.2;  %activation energy Zr
Tm.EaAp=132.41;%130.6; %activation energy Ap
Tm.EaBt=222.23;%197; %activation energy Bi

%%%%%%%%%%%%%

Pm.DikeThick=4.5; %dike half-thickness in m; 
Pm.L = 595; %total length of domain in meters (this matters when RightBC=2)

%boundary condition in far field away from dike. 
%RightBC=1 is insulating (Neumann condition), RightBC=2 is Dirichlet (fixed temperature at initial conditions)
Pm.RightBC=2; 
Pm.h = 7;%7;%Pm.L/100; %grid spacing, m
Pm.dip = sind(90); % correct for dip of dike

Pm.year = 3600*24*365;  %one year in seconds

Pm.rho=2600 ;%host rock density, kg/m3
Pm.rhom=2800 ;%magma density, kg/m3
Pm.Lf=310000; %latent heat fusion J/Kg
Pm.Tsol=725; %tonalite solidus in C
Pm.Tliq=1100; %tonalite liquidus in C
Pm.Tsold=1015; %dike material solidus in C
Pm.Tliqd=1165; %dike material liquidus in C
Pm.R=8.314*1e-3; %gas constant kJ/mol K

Pm.b=0.5; %exponent for melt fraction temperature curve. b-->1 is "more mafic"
Pm.bd=1.7; %exponent for basaltic dike
Pm.cp = 1100; %specific heat capacity (assumed same between materials)

Pm.kappa = Pm.k/(Pm.cp*Pm.rho); %thermal diffusivity

%%%%%%%%%%%
%% Total time of model run
%choose simulation time based on diffusion distance of farthest sample(approximately 100 m from dike)
diffdist = 100^2/Pm.kappa /Pm.year; 

Pm.TotT = 5.4*diffdist;  %total time of simulation in years

%%
%coordinate transform for increased spatial resolution near dike
Pm.dononuniform=1; %flag for uniform or non-uniform grid spacing
Pm.ell=Pm.L/10; %scale parameter for non-uniform grid (arctangent function) - smaller values are more uniform grid spacing
Pm.N = Pm.L/Pm.h ; %number of nodes = domain length in km/spacing

if Pm.dononuniform
    Pm.dZ=1/(Pm.N);
    Pm.dikeThickZ=atan(Pm.DikeThick/Pm.ell)/atan(Pm.L/Pm.ell);
    Pm.dikepts = 1:round(Pm.dikeThickZ/Pm.dZ);
    Pm.Z = 0:Pm.dZ:1; %distance vector in scaled coordinates
    Pm.Z=Pm.Z';
    Pm.x = Pm.ell*tan(atan(Pm.L/Pm.ell)*Pm.Z);
    
    Pm.dxdZinv = 1./(Pm.ell*atan(Pm.L/Pm.ell).*sec(Pm.Z.*atan(Pm.L/Pm.ell)).^2);    

    Pm.beta = Pm.dxdZinv.*Pm.k.*ones(length(Pm.x),1); %coefficients in transformed coordinates

    e2=(Pm.beta(2:Pm.N)+Pm.beta(3:Pm.N+1));
    e1=(Pm.beta(2:Pm.N)+Pm.beta(1:Pm.N-1));
    e0=-(2*Pm.beta(2:Pm.N)+Pm.beta(3:Pm.N+1)+Pm.beta(1:Pm.N-1));
    Dint=1./(2.*Pm.dZ^2)*(diag(e1(2:Pm.N-1),-1)+diag(e2(1:Pm.N-2),1)+diag(e0,0));
    
    D=zeros(Pm.N+1);
%sparse representation - doesnt speed things up unless N is large    
    e11=[e1(2:Pm.N-1);0];e22=[0;e1(2:Pm.N-1)];
    Dint2 = 1./(2.*Pm.dZ^2) .*spdiags([e11 e0 e22],-1:1,Pm.N-1,Pm.N-1);
%     
    D(2:Pm.N,2:Pm.N)=Dint2;
    D(2,1)=1./(2.*Pm.dZ^2)*(Pm.beta(2)+Pm.beta(1));
    D(Pm.N,Pm.N+1)=1./(2.*Pm.dZ^2)*(Pm.beta(Pm.N+1)+Pm.beta(Pm.N));
    
    Pm.D=D;
    clear D Dint e0 e1 e2 
else
    Pm.dikepts = 1:round(Pm.DikeThick/Pm.h);
    Pm.x = 0:Pm.h:Pm.L; %distance vector
end

%% run diffusion model
tic
[x,time,theta] = diffusion_1d_dike(Pm);

Pts = linspace(2,100,100);%[2 4 10 40];

%calculate time-temperature histories at query points
TimeTemp = timetemphistory(x,theta,Pts);
toc
disp('Done with thermal model')
disp(['Total simulation time: ' num2str(Pm.TotT)])

%% plot summary from thermal model
%uncomment to display data associated with the thermal model only

%summaryplots_dike1d(time,theta,x,Pm)

%% thermochron forward calculations
%compare with data, assuming all available dataset here

load('JacksonA_meanages_3_2022.mat');

ZPts = 1:1:220;%Pm.DikeThick + zhe_dist*Pm.dip; %evaluation points for Zr
APts = 1:1:220;%Pm.DikeThick + ahe_dist*Pm.dip; %evaluation points for Ap
BPts = 1:1:220;%Pm.DikeThick + ahe_dist*Pm.dip; %evaluation points for Ap
AFTPts = linspace(1,220,40);

%define data structure associated with thermochron model(s)
Tm.trZ = 15;%age of resetting in Ma
Tm.tiZ= 112;%age of background in Ma

Tm.trA = 15;%age of resetting in Ma
Tm.tiA= 86;%age of background in Ma

Tm.trB = 15;%age of resetting in Ma
Tm.tiB= 125;%age of background in Ma

Tm.trAFT = 15;
Tm.tiAFT = 97;

%more kinetic parameters for Bt
D0Bt=7.5e-2; % Bi reference diffusivity, cm^2/s
A2_Bt=(750e-4)^2; %radius squared, cm
Tm.D0a2_Bt = D0Bt/A2_Bt;

%parameters for AFT in Ketcham 1999 model (from Reiners 2008) lc,mod
Tm.c0 = -19.844;%-26.039;
Tm.c1 = 0.38951;%0.53168;
Tm.c2 = -51.253;%-62.319;
Tm.c3 = -7.6423;%-7.8935;
Tm.alpha = -0.12327;%-0.20196;
Tm.beta = -11.988;%-7.4224;

%parameters for ZFT from Ch. 3 of 'Fission Track Thermochronology and it's
%application to Geology' (Ketcham 2018)
% Tm.c0 = -91.659;
% Tm.c1 = 2.09366;
% Tm.c2 = -314.94;
% Tm.c3 = -14.287;
% Tm.alpha = -0.05721;
% Tm.beta = -1;

rmr0 =  0.869;%0.8204;

%parameterizes D0^2 in terms of Ea, returns larger Tm
Tm = fracreset_paramfit(Tm);

%Zircon
zhe_dist=zhe_jacka_agemean(:,1);
zhe_age=zhe_jacka_agemean(:,2);
zhe_err=zhe_jacka_agemean(:,3);
Zr_agemodel = thermochron_calc(zhe_dist,theta,x,time,Tm.D0a2_Zr,Tm.EaZr,Tm.trZ,Tm.tiZ,Tm,Pm,1);

%Apatite
ahe_dist=ahe_jacka_agemean(:,1);
ahe_age=ahe_jacka_agemean(:,2);
ahe_err=ahe_jacka_agemean(:,3);
Ap_agemodel = thermochron_calc(ahe_dist,theta,x,time,Tm.D0a2_Ap,Tm.EaAp,Tm.trA,Tm.tiA,Tm,Pm,1);

%Biotite
bt_dist=btar_jacka_agemean(:,1);
bt_age=btar_jacka_agemean(:,2);
bt_err=btar_jacka_agemean(:,3);
Bt_agemodel = thermochron_calc(bt_dist,theta,x,time,Tm.D0a2_Bt,Tm.EaBt,Tm.trB,Tm.tiB,Tm,Pm,2);

%AFT
aft_dist=aft_jacka_agemean(:,1);
aft_age=aft_jacka_agemean(:,2);
aft_err=aft_jacka_agemean(:,3);
AFT_agemodel = thermochron_calc(aft_dist,theta,x,time,rmr0,Tm.EaBt,Tm.trAFT,Tm.tiAFT,Tm,Pm,3);

%ZFT
% zft_dist=zft_jacka_agemean(:,1);
% zft_age=zft_jacka_agemean(:,2);
% zft_err=zft_jacka_agemean(:,3);
% ZFT_agemodel = thermochron_calc(zft_dist,theta,x,time,Tm.D0a2_Bt,Tm.EaBt,Tm.trB,Tm.tiB,Tm,Pm,4);

%PaleoMagnetic Resetting
%find the 580 deg isotherm
CurieT = 580;
MaxT=max(theta,[],2);

%assuming hanging wall Dc range from Biasi and Karlstrom and taking average...

mgt_dist=Pm.DikeThick+(5.4+3.7)/2; 
mgt_err=(5.4-3.7)/2;%assuming error corresponding to range 

%finds the location in wall rocks at which the maximum temperature reaches
%Curie point during run
X_MGT = interp1(MaxT(x>Pm.DikeThick),x(x>Pm.DikeThick),CurieT);

disp('Done with chemical model')

%% Calculate residuals and make more plots
%calculate residuals
Zresidual=zhe_age - Zr_agemodel;
Aresidual=ahe_age - Ap_agemodel;
Bresidual=bt_age - Bt_agemodel;
AFTresidual=aft_age - AFT_agemodel;
% ZFTresidual=aft_age - ZFT_agemodel;
MGTresidual=mgt_dist - X_MGT;

%errors, covariance matrix
Ap_covinv = diag(1./ahe_err.^2);
Zr_covinv = diag(1./zhe_err.^2);
Bt_covinv = diag(1./bt_err.^2);
AFT_covinv = diag(1./aft_err.^2);
%ZFT_covinv = diag(1./zft_err.^2);
MGT_covinv = diag(1./mgt_err.^2);

disp('norm of residuals for Zr and Ap mean ages')

ZResTot = (Zresidual' * Zr_covinv * Zresidual)/2;
AResTot = (Aresidual' * Ap_covinv * Aresidual)/2;
BResTot = (Bresidual' * Bt_covinv * Bresidual)/2;
AFTResTot = (AFTresidual' * AFT_covinv * AFTresidual)/2;
%ZFTResTot = (ZFTresidual' * ZFT_covinv * ZFTresidual)/2;
MGTResTot = (MGTresidual' * MGT_covinv * MGTresidual)/2;

disp((Zresidual' * Zr_covinv * Zresidual)/2)
disp((Aresidual' * Ap_covinv * Aresidual)/2)

Zr_agemodel_pl = thermochron_calc(ZPts,theta,x,time,Tm.D0a2_Zr,Tm.EaZr,Tm.trZ,Tm.tiZ,Tm,Pm,1);
Ap_agemodel_pl = thermochron_calc(APts,theta,x,time,Tm.D0a2_Ap,Tm.EaAp,Tm.trA,Tm.tiA,Tm,Pm,1);
Bt_agemodel_pl = thermochron_calc(BPts,theta,x,time,Tm.D0a2_Bt,Tm.EaBt,Tm.trB,Tm.tiB,Tm,Pm,2);
AFT_agemodel_pl = thermochron_calc(AFTPts,theta,x,time,rmr0,Tm.EaBt,Tm.trAFT,Tm.tiAFT,Tm,Pm,3);
%ZFT_agemodel_pl = thermochron_calc(APts,theta,x,time,Tm.D0a2_Bt,Tm.EaBt,Tm.trB,Tm.tiB,Tm,Pm,4);

%plot up data vs model
figure(1)
plot(ZPts,Zr_agemodel_pl,'r',APts,Ap_agemodel_pl,'b')
hold on
plot(BPts,Bt_agemodel_pl,'k')
plot(AFTPts,AFT_agemodel_pl,'g')
%plot(APts,ZFT_agemodel_pl,'m')
errorbar(zhe_dist,zhe_age,zhe_err,'ro')
errorbar(ahe_dist,ahe_age,ahe_err,'bx')
errorbar(bt_dist,bt_age,bt_err,'ks')
errorbar(aft_dist,aft_age,aft_err,'g*')
%errorbar(zft_dist,zft_age,zft_err,'m*')


labelT{1}=['Zr model, residual: ' num2str(ZResTot)];labelT{2}=['Ap model, residual: ' num2str(AResTot)];labelT{3}=['Bt model, residual: ' num2str(BResTot)];
labelT{5}=['Zr data'];labelT{6}=['Ap data'];labelT{7}=['Bt data'];labelT{4}=['AFT model, residual:' num2str(AFTResTot)];labelT{8}=['AFT data'];%label{9}=['ZFT data'];label{10}=[ZFT model, residual'];
legend(labelT)
xlabel('Distance from contact (m)')
ylabel('Age (Ma)')
hold off

% %% plotting for comparison
% 
% XPts = [20,50,100]; %distances at which you want to compare T-t histories
% TtXPts = rescompare(x,theta,XPts);
% 
% Ap_age_XPts = thermochron_calc(XPts,theta,x,time,Tm.D0a2_Ap,Tm.EaAp,Tm.trA,Tm.tiA,Tm,Pm,3);
% 
% figure(5)
% for i=1:length(XPts)
% plot(time/Pm.year,TtXPts(:,i));
% hold on
% disp(['At distance  = ' num2str(XPts(i)) ' m, predicted age = ' num2str(Ap_age_XPts(i)) ' Ma'])
% end
% disp(['Total simulation time: ' num2str(Pm.TotT)])

figure(2);
plot(x,MaxT,'o-')
hold on
plot(X_MGT,580,'v-')
errorbar(mgt_dist,580,0,0,mgt_err,mgt_err,'*')
labelMGT{1}=['Max Temperature reached during model'];
labelMGT{2}=['Predicted magnetic reset distance'];
labelMGT{3}=['Observed magnetic reset distance, residual: ' num2str(MGTResTot)];
legend(labelMGT)
xlabel('Distance from dike center (m)')
ylabel('Temperature (C)')
xlim([0,20])
hold off

%% save stuff
%try to minimize saving large arrays to keep file size down
%save('Forwardmodel.mat','x','theta','time','Pm','zhe_dist','zhe_age','zhe_err');%,'hammer','params')
