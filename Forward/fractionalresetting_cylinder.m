function f = fractionalresetting_cylinder(T,t,D0a2,Ea,R)
%predicts the fractional resetting, assuming spherical mineral grain (using
%series solution of Wolff et al., 1998)

%inputs:
%T time-Temperature history
%t time
%D diffusivity
%a spherical grain radius
%Ea activation energy 
%R %gas constant

%ratio of diffusion time to duration of heating, temperature in Kelvin!
tau = D0a2 * trapz(t,exp(-Ea/R *1./(T+273)));

%fractional resetting
f=0;   

%number of terms in series to evaluate. Is it enough??
num=250;

%calculate Bessel zeros
J0s = besselzero(0, num);

for n=1:num %calculate approx based on num terms in infinite series
    alphan = J0s(n);
    f = (1/alphan^2)*exp(-alphan^2*tau) + f;
end

f=1-4*f;

%approximation in Reiners 2009
% if f>=0.6
%     fapprox = 1 - 9/13*exp(-5.78*tau);
% else
%     fapprox = 4/sqrt(pi) * sqrt(tau) - tau;
% end
%keyboard
end

