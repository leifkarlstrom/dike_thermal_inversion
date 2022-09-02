%     This script implements a Markov Chain Monte Carlo inversion of one or more
%     thermochronometric datasets assuming a 1D heating event associated with dike intrusion
%     Solves 1D unsteady heat equation implicitly, using stretched coordinates with finite differences
%     
%     Currently set up to deal with one or more of: (U-Th)/He diffusion, Bt/Ar diffusion, Ap Fission Track, 
%     Magnetic geothermometers
%
%     Uses the MCMC Hammer inversion code written by Kyle Anderson
%     which is a required dependency (Matlab toolbox addon) 
%     https://gitlab.com/kander/mcmc_hammer
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
clear 
close all

tic

%add folders containing data and forward models
%modify these to whatever they are locally!!
addpath('../Forward')%('../../Data','../../Forward')

%specify output path
savedir = '../Runs/';

%% Set up the inversion
% Use parallel computing toolbox for inversion? 1=yes, 0=no
parpool(4) %uncomment if you want to change the number of cores used
runparallel = 1;                

%shuffle random number generator or use "default"
rng("shuffle")

% Which data are we interested in inverting?
%1 = apatite/zircon (U-Th)/He thermochron data
%2 = magnetic geothermometry data
%3 = full suite of multiple chronometers including fission track
%4 = full suite of thermochronometers + magnatic geothermometry data 
Pm.DataFlag = 2;

%using most recent version of the Jackson A dataset here
load('JacksonA_meanages_3_2022.mat');

% MCMC Hammer parameters
nIts = 1e3; %total number of model iterations to be completed in a MCMC chain
nBurn = 1e2; %number of iterations to discard (as "Burn In"). Typically ~10% of nIts 
nThin = 1; %only keep every nThin interation, values greater than one represent "thinned" models, helps reduce autocorrelation in chains (sortof)
updateInt = 5; %how often (after how many iterations) to print out status report

Pm.nChains = 1; %number of MCMC chains

%% Define parameters 
%parameters fixed during MCMC inversion

%geometric parameters and simulation parameters
Pm.L = 595; %total length of domain in meters (this matters when RightBC=2)
%boundary condition in far field away from dike. 
%RightBC=1 is insulating (Neumann condition), RightBC=2 is Dirichlet (fixed temperature at initial conditions)
Pm.RightBC=2; 
Pm.h = 7;%7;%Pm.L/100; %grid spacing, m
Pm.dip = sind(90); % to correct for dike dip
Pm.DikeThick=4.5; %dike half-thickness in m; Biasi says 9m

%coordinate transform for increased spatial resolution near dike
Pm.dononuniform=1; %flag for uniform or non-uniform grid spacing
Pm.ell=Pm.L/10; %scale parameter for non-uniform grid (arctangent function) - smaller values are more uniform grid spacing

%these will get reset during inversions 
Pm.Tbackground = 50;  %background temperature in C
Pm.twidth = 0.1; %vary this parameter to change the rate at which dike temperatures decrease
Pm.tcenter = 1; %scale for dike shut-off time (years)
Pm.k = 3; %thermal conductivity (assumed same between materials)

%physical constants
Pm.R=8.314*1e-3; %gas constant kJ/mol K
Pm.rho=2600 ;%host rock density, kg/m3
Pm.rhom=2800 ;%magma density, kg/m3
Pm.Lf=310000; %latent heat fusion J/Kg
Pm.Tsol=725; %tonalite solidus in C
Pm.Tliq=1100; %tonalite liquidus in C
Pm.Tsold=1015; %dike material solidus in C
Pm.Tliqd=1165; %dike material liquidus in C

Pm.b=0.5; %exponent for melt fraction temperature curve. b-->1 is "more mafic"
Pm.bd=1.7; %exponent for basaltic dike
Pm.cp = 1100; %specific heat capacity (assumed same between materials)
Pm.kappa = Pm.k/(Pm.cp*Pm.rho); %thermal diffusivity

%% set total time of each conduction simulation
Pm.year = 3600*24*365;  %one year in seconds
%choose simulation time based on diffusion distance of farthest sample (approximately 100 m from dike)
diffdist = 100^2/Pm.kappa /Pm.year; %diffusion distance depends on thermal diffusivity (kappa)
Pm.TotT = 5.4*diffdist;  %total time of simulation in years

%% 
%precompute grid transform
Pm.N = Pm.L/Pm.h ; %number of nodes = domain length in km/spacing

if Pm.dononuniform
    Pm.dZ=1/(Pm.N);
    Pm.dikeThickZ=atan(Pm.DikeThick/Pm.ell)/atan(Pm.L/Pm.ell);
    Pm.dikepts = 1:round(Pm.dikeThickZ/Pm.dZ);
    Pm.Z = 0:Pm.dZ:1; %distance vector in scaled coordinates
    Pm.Z=Pm.Z';
    Pm.x = Pm.ell*tan(atan(Pm.L/Pm.ell)*Pm.Z);    
    Pm.dxdZinv = 1./(Pm.ell*atan(Pm.L/Pm.ell).*sec(Pm.Z.*atan(Pm.L/Pm.ell)).^2);    
else
    Pm.dikepts = 1:round(Pm.DikeThick/Pm.h);
    Pm.x = 0:Pm.h:Pm.L; %distance vector
end

%% load data, define MCMC parameters

if Pm.DataFlag == 1 || Pm.DataFlag == 3 || Pm.DataFlag == 4
    %zircon data
    zhe_dist=zhe_jacka_agemean(:,1);
    zhe_age=zhe_jacka_agemean(:,2);
    zhe_err=zhe_jacka_agemean(:,3);
    %apatite data
    ahe_dist=ahe_jacka_agemean(:,1);
    ahe_age=ahe_jacka_agemean(:,2);
    ahe_err=ahe_jacka_agemean(:,3);
    
    %get the data into datastructure
    params.data.ZPts = zhe_dist; %evaluation points for Zr 
    params.data.APts = ahe_dist; %evaluation points for Ap 
    
    params.data.ApAge = ahe_age;
    params.data.Ap_covinv = diag(1./ahe_err.^2);
    params.data.ZrAge = zhe_age;
    params.data.Zr_covinv = diag(1./zhe_err.^2);
    
    % compute the chemical parameters: we use different ages for each
    % system to focus on the common RESETTING part of the thermochron history
    Tm.trZ = 15;%age of resetting in Ma
    Tm.tiZ= 112;%age of background in Ma

    Tm.trA = 15;%age of resetting in Ma
    Tm.tiA= 86;%age of background in Ma
    
    %define data structure associated with thermochron model(s)
    Tm.EaZr=167.2;  %activation energy Zr
    Tm.EaAp=130.6; %activation energy Ap
    
    %parameterizes D0^2 in terms of Ea, returns larger Tm
    Tm = fracreset_paramfit(Tm); %for zircon and apatite
end
if Pm.DataFlag == 3 || Pm.DataFlag == 4
    %biotite data
    btar_dist=btar_jacka_agemean(:,1);
    btar_age=btar_jacka_agemean(:,2);
    btar_err=btar_jacka_agemean(:,3);
    %AFT
    aft_dist=aft_jacka_agemean(:,1);
    aft_age=aft_jacka_agemean(:,2);
    aft_err=aft_jacka_agemean(:,3);
    %ZFT
    aft_dist=aft_jacka_agemean(:,1);
    aft_age=aft_jacka_agemean(:,2);
    aft_err=aft_jacka_agemean(:,3);
    
    params.data.BPts = btar_dist; %evaluation points for Bt 
    params.data.AFTPts = aft_dist; %evaluation points for Bt 

    params.data.BtAge = btar_age;
    params.data.Bt_covinv = diag(1./btar_err.^2);

    params.data.AFTAge = aft_age;
    params.data.AFT_covinv = diag(1./aft_err.^2);
    
    %ages for Bt and AFT systems
    Tm.trB = 15;%age of resetting in Ma
    Tm.tiB= 125;%age of background in Ma

    Tm.trAFT = 15;
    Tm.tiAFT = 97;
    
    Tm.EaBt=230;%197; %activation energy Bi

    D0Bt=7.5e-2; %biotite reference diffusivity, cm^2/s
    A2_Bt=(750e-4)^2; %radius squared, cm
    Tm.D0a2_Bt = D0Bt/A2_Bt;

    %parameters for AFT in Ketcham 1999 model (from Reiners 2008) lc,mod
    Tm.c0 = -19.844;
    Tm.c1 = 0.38951;
    Tm.c2 = -51.253;
    Tm.c3 = -7.6423;
    Tm.alpha = -0.12327;
    Tm.beta = -11.988;
end

%Magnetic resetting data

%assuming hanging wall Dc range from Biasi and Karlstrom and taking average...
%note that we probably should view these constraints loosely and use the
%whole range implied between both hanging and footwall measurements - this
%will increase both the predicted distance and the errors.
if Pm.DataFlag == 2 || Pm.DataFlag == 4
    mgt_dist=Pm.DikeThick+(5.4+3.7)/2;  %from Biasi and Karlstrom 2021
    mgt_err=(5.4-3.7)/2;%assuming error corresponding to range 

    Pm.CurieT = 580; %Curie temperature (C)
    
    params.data.MGTdist = mgt_dist;
    params.data.MGT_covinv = diag(1./mgt_err.^2);
end


%get parameter ranges assuming uniform distributions. 
%Note that if you change the variables in here you have to also modify the
%hammer.setparam below
PmRg = ParameterRanges(Pm); 

%thermal model parameters
params.Pm=Pm;

if Pm.DataFlag~=2
%thermochron parameters
    params.Tm=Tm;
end

%% Run the inversion

hammer = mcmc_hammer('NumIterations',nIts,'NumBurn',nBurn,'NumThin',nThin,'UpdateInterval',updateInt,'PerturbModelFunc','gaussian','Cascade',1,...
    'RandomStart',1,'RandomStartCheck',1,...
    'StepWidthGlobalDiv',   110,...
    'UserParams',           params,...
    'LikelihoodFunc',       @dike_likelihood);

%set the parameters and ranges, assuming uniform distribution. Note which
%ones are in logspace and which are not!!!
hammer = hammer.setparam('log10BgT','lims',[log10(PmRg(1,1)) log10(PmRg(1,2))],'title','Background T (C)');
hammer = hammer.setparam('log10Tw','lims',[log10(PmRg(2,1)) log10(PmRg(2,2))],'title','Dike unsteadiness');
hammer = hammer.setparam('log10Tc','lims',[log10(PmRg(3,1)) log10(PmRg(3,2))],'title','Dike longevity');
hammer = hammer.setparam('k','lims',[(PmRg(4,1)) (PmRg(4,2))],'title','thermal conductivity k');


if Pm.DataFlag == 1 || Pm.DataFlag == 3 || Pm.DataFlag == 4
hammer = hammer.setparam('Ea_Z','lims',[PmRg(6,1) PmRg(6,2)],'title','Ea Zr');
hammer = hammer.setparam('Ea_A','lims',[PmRg(5,1) PmRg(5,2)],'title','Ea Ap');
end

if Pm.DataFlag == 3 || Pm.DataFlag == 4
hammer = hammer.setparam('Ea_B','lims',[PmRg(7,1) PmRg(7,2)],'title','Ea Bt');
hammer = hammer.setparam('rmr0_AFT','lims',[PmRg(8,1) PmRg(8,2)],'title','AFT rmr0');
end

% %post processing
% hammer = hammer.addpostprocessedvar('Tf',0.1*(10.*10.^hammer.getsamples('log10Tc') + 10.^hammer.getsamples('log10Tw') .*atanh(1-0.02)));
% hammer = hammer.addpostprocessedvar('Tbg',10.^hammer.getsamples('log10BgT'));
% %hammer = hammer.addpostprocessedvar('k',10.^hammer.getsamples('log10k'));
% hammer = hammer.addpostprocessedvar('Tc',10.^hammer.getsamples('log10Tc'));
% hammer = hammer.addpostprocessedvar('Tw',10.^hammer.getsamples('log10Tw'));
% %hammer = hammer.addpostprocessedvar('Tp',10.^hammer.getsamples('log10Tp'));

%run MCMC hammer
if runparallel
    hammer = hammer.runparallel;
else
    hammer = hammer.run;
end

savefilename = [savedir 'JackAllThermo_' date '.mat'];

save(savefilename);
toc

