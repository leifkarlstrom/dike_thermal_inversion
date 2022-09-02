function f = fractionalresetting_sphere(T,t,D0a2,Ea,R)
%predicts the fractional resetting, assuming spherical mineral grain (using
%series solution of Wolff et al., 1998)

%inputs:
%T time-Temperature history
%t time
%D diffusivity
%a spherical grain radius
%Ea activation energy 
%R %gas constant

%ratio of diffusion time to duration of heating
tau = D0a2 * trapz(t,exp(-Ea/R *1./(T+273)));

%fractional resetting
f=0;   
for n=1:1000 %calculate approx based on 1000 terms in infinite series
    f = (1/n^2)*exp(-n^2*pi^2*tau) + f;
end

f=1-(6/pi^2)*f;

%approximation in Reiners 2009
% if f>=0.85
%     fapprox = 1 - (6/pi^2)*exp(-pi^2*tau);
% else
%     fapprox = 6*sqrt(tau/pi) - 3*tau;
% end
%keyboard
end

