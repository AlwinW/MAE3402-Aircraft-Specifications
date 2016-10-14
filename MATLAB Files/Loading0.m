%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loading/Bending Moment Analysis %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% n = 0 %%
clc; clear; clear all;

% Common Parametres %
g = 9.8065;
n = 3.5;
s_f = 1.5;
b = 22.135;
a = 11.0675;
Stress_Allowable = 950*10^6;

% Chord length/Thickness %
for z = 1:1107;
        c(z) = 1.58114.*[1 - 0.0542.*abs(z./1107.*a)];
        t(z) = 0.18.*c (z);
end

% Location of Spar 1/2 %
s1 = 0.15.*c;
s2 = 0.65.*c;

%%%% Spar 1 %%%%
for v = 1:1107;
        h(v) = 230/130*0.1.*t(v)./0.18 - 0.004; %% Height of the spar, since we need to figure out the thickness
        h2(v) = 0.1.*t(v)./0.18;
end

%% Loading Parametres %%
L_Aero = zeros(1,1107);

L_Weight = ones(1,1107);
L_Weight = 2800.*g.*L_Weight;

L_Batt = ones(1,1107);
L_Batt = c./mean(c).*(0.8.*3500 + 300).*g./2./1107.*L_Batt;

L_Eng = zeros(1,1107);
for x1 = 2.75:0.01:3.06;
    L_Eng(x1*100) = 150*g/301;
end
M_Eng = zeros(1,1107);

L_Wing = ones(1,1107);
L_Wing = 186.3.*g./2./1107.*L_Wing;

%% Moment Derivation %%
for t = 1:1107;
    for t1 = t:1107;
        M_Aero_Neg (t1) = (1107 - t1)./1107.*a.*L_Aero(1108 - t1);
    end 
    for t2 = t:1107;
        M_Batt_Neg (t2) = (1107 - t2)./1107.*a.*L_Batt(1108 - t2);
    end 
    for t3 = t:1107;
        M_Wing_Neg (t3) = (1107 - t3)./1107.*a.*L_Wing(1108 - t3);
    end
    for t4 = t:305;
        M_Eng_Neg (t4) = (305 - t4)./1107.*a.*L_Eng(306 - t4);
    end
    
    M_Aero(t) = sum(M_Aero_Neg(1:1108-t));
    M_Batt(t) = sum(M_Batt_Neg(1:1108-t));
    M_Wing(t) = sum(M_Wing_Neg(1:1108-t));
    if t < 306;
          M_Eng(t) = sum(M_Eng_Neg(1:291-t));
    end
    
    Moment (t) = 0 - M_Aero(t) - M_Batt(t) - M_Wing(t) - M_Eng(t);
end

%% Determining Bending Load %%

for q = 1:1107;
    B(q) = 6.*Moment(q)./Stress_Allowable./[h(q)].^2;
end
M_Wing = sum(B.*h).*a.*2.*1590./1106;

%% Iterations to refine the wing weight %%
n = 0;
while n < 10;
    n = n+1;
    
    L_Wing = ones(1,1107);
    L_Wing = L_Wing.*B.*h.*g.*1590./1107;
 
    for t = 1:1107;

    for t1 = t:1107;
        M_Aero_Neg (t1) = (1107 - t1)./1107.*a.*L_Aero(1108 - t1);
    end 
    for t2 = t:1107;
        M_Batt_Neg (t2) = (1107 - t2)./1107.*a.*L_Batt(1108 - t2);
    end 
    for t3 = t:1107;
        M_Wing_Neg (t3) = (1107 - t3)./1107.*a.*L_Wing(1108 - t3);
    end
    for t4 = t:305;
        M_Eng_Neg (t4) = (305 - t4)./1107.*a.*L_Eng(306 - t4);
    end
    
    M_Aero(t) = sum(M_Aero_Neg(1:1108-t));
    M_Batt(t) = sum(M_Batt_Neg(1:1108-t));
    M_Wing(t) = sum(M_Wing_Neg(1:1108-t));
    if t < 306;
          M_Eng(t) = sum(M_Eng_Neg(1:291-t));
    end
    
    Moment (t) = 0 - M_Aero(t) - M_Batt(t) - M_Wing(t) - M_Eng(t);
    end  
end

hold on
plot(0:0.01:a, Moment, 'y')
title('Moment distribution over the wing, n = 0')
xlabel('Root chord to tip chord')
ylabel('Bending Moment')
%% Net Force Determination %%
Force = L_Aero - L_Batt - L_Wing - L_Eng;
figure
plot(0:0.01:a, 100*Force, 'k')
title('Force distribution over the wing, n = 0')
xlabel('Root chord to tip chord')
ylabel('Net Force')




