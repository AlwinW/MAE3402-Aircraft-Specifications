%%
clear all
close all
clc
%code written by AJ Kusangaya 26071371
%Uses design for proposal parameters to find the acceptable design range
%Uses aerodynamic inputs from the energy requirement spreadsheet
%% Aerodynamic and Flight parameters
K = 0.0208;
CD0 = 0.01648;
CL_star = sqrt((CD0/K));
LD_star = 1/(2*sqrt((CD0*K)));
pcruise = 0.9091;
W0_S = 500: 250 : 3000;
P0_S = 0:0.5:25;
CLto = 2.3;
CLmax = 2.8;
Vapp = 51.444;
SSCtheta = 8.53;
npr = 0.85;
Vcruise = 82.1;
W0_Sclstar = 2656.4.*ones(1,length(P0_S));
%%Landing Wing Loading
W0_Slanding = CLmax*(1.225/2)*((Vapp/1.3)^2);
W0_Sland = W0_Slanding.*ones(1,length(P0_S));
%% Cruise at any Cl (L/D = 27)
CL = (2/pcruise).*W0_S.*(1/(Vcruise^2));
Kalt = (CL - (27*CD0))./(27.*(CL.^2));
P_W_cruise = Vcruise.*((CD0./CL) + (Kalt.*CL));
A = 1./(pi.*Kalt.*0.85);
%% Second Segment Climb One engine inoperative
CD = CD0 + (Kalt.*CL.*CL);
P_W_ssc = 2.*((0.015 + ((1.44/CLmax).*CD))./ ((npr/1.138).*((W0_S).^-0.5)));
%% TOFL
P_Wtofl = zeros(1,length(W0_S));
for i = 1:length(W0_S);
    TOFL = @(x) 1000 - (11.8*(x*W0_S(i)*(1/CLto))) - (0.0255*((x*W0_S(i)*(1/CLto))^2));
    root = false_pos(TOFL,0,50,0.001);
    P_Wtofl(i) = 1/root;
end
plot(W0_S,P_W_cruise,'r', ...
    W0_Sclstar,P0_S,'b--', ...
    W0_Sland,P0_S,'c', ...
    W0_S,P_W_ssc,'y', ...
    W0_S,P_Wtofl,'g')
title('Design Constraints Wing and Power Loading L/D 27 For Cruise (Aspect Ratio)')
xlabel('Wing Loading Pa')
ylabel('Power Loading W/N (Aspect Ratio)')
legend('Cruise ','Wing Loading Energy Analysis, min AR','Landing Wing Loading','Second Segment Climb OEI','TOFL 1km','100ft/min Cruise Climb','cruise at current cl')
% figure
% P_A = 10:0.5:23;
% W0_Sclstar = 2656.4.*ones(1,length(P_A));
% plot(W0_S,A,'g',W0_Sclstar,P_A,'b--')
% title('Aspect Ratio Vs Wing Loading')
% xlabel('Wing Loading Pa')
% ylabel('Aspect Ratio')