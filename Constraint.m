%Constraint Analysis
clc;clear all;close all


Wing_Loading = 500:250:3000;

%Landing Approach Speed
Clmax = 2.7;
Vapp = 100*0.514444; %Max Approach Speed
rhoSL = 1.225; 
Wing_Loading_App = Clmax*(rhoSL/2)*(Vapp/1.3)^2;

%Take off field length
ClmaxTO = 2.7;
TOFL = 1000; %Take off field length (runway length)
a = Wing_Loading.^2 * (0.225/ClmaxTO^2);
b = (11.8/ClmaxTO)*Wing_Loading;
c = -TOFL;

Power_Loading1 =1./( (-b+sqrt(b.^2 - 4.*a.*c))./(2.*a) );

%Cruise speed Constraint
V = 0.25*328.387;
rho = 0.904637;
q = 0.5*rho*V^2;
Clcruise = (1/q).*Wing_Loading;
prop_eff = 0.77;
Cd0 = 0.0191;
A = 20;
e = 0.827;
K = 1/(3.1415*e*A);
%Ps = 0;
%n = 1;

Power_Loading2 = (V/prop_eff).*((Cd0./Clcruise)+(K.*Clcruise));

%Service ceiling
Vv_max = (100*0.3048)/60; %Rate of Climb at ceiling (min)
rhomax = 0.849137;
L_D = 1/ sqrt(4*Cd0*K);

Power_Loading3 = (Vv_max/prop_eff) + (2/(prop_eff*rhomax)) * sqrt(K/(3*Cd0)) .* (Wing_Loading.^0.5) * (1.155/L_D);

%Climb Gradient
sec_seg = 0.015; %Minimum 2nd-segment climb gradient w/ 1 engine
c = sqrt(1.2^2*2/(rhoSL*ClmaxTO));
Cl = ClmaxTO/(1.2^2);
Cd0G = 0.028;
Cd = Cd0G+(Cl^2/(3.1415*A*e));
Power_Loading4 = (2* (sec_seg + (1.2^2)*Cd / ClmaxTO) / (prop_eff/c)).*Wing_Loading.^0.5;


%ROC @ Cruise
Vv = 1.524;
Power_Loading5 = (Vv/prop_eff) + (2/(prop_eff*rho)) * sqrt(K/(3*Cd0)) .* (Wing_Loading.^0.5) * (1.155/L_D);

%CL*
Clstar = sqrt(Cd0/K);
Wing_Loading_Star = q*Clstar;

%n-limit
n=3.5;
Clmaxclean = 1.9;
Wing_Loading_n = (q/n)*Clmaxclean;

hold on;
plot(Wing_Loading,Power_Loading1,Wing_Loading,Power_Loading2,Wing_Loading,Power_Loading3,Wing_Loading,Power_Loading4,Wing_Loading,Power_Loading5);

SP = Wing_Loading_App;
y1 = get(gca,'ylim');
hold on
plot([SP SP],y1,[Wing_Loading_Star Wing_Loading_Star],y1,[Wing_Loading_n Wing_Loading_n],y1)
legend('TOFL','Cruise speed','Ceiling','Climb grad','ROC @ Cruise','landing','Cl*','nlimit');
xlabel('W0/S')
ylabel('P0/W0')