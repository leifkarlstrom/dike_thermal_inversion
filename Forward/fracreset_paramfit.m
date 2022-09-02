function Tm = fracreset_paramfit(Tm)
%Oct 2021

%parameter fit of ln(D0/a^2) as a function of activation energy Ea and Ez
%from experimental data

logD0a2_Ap = 0.00045*Tm.EaAp.^2.083;

logD0a2_Zr = 2e-11*Tm.EaZr.^5.074;

Tm.D0a2_Zr=exp(logD0a2_Zr); %Diffusivity Zr
Tm.D0a2_Ap=exp(logD0a2_Ap);  %Diffusivity Ap





