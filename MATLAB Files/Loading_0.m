%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loading/Bending Moment Analysis %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; clear all;

% Common Parametres %
g = 9.8065;
n = 3.5;
s_f = 1.5;
b = 22.135;
a = 11.0675;
S_Shear_Max = 550 * 10^6; % This results in 12.5%% strain (8.5 not including the safety factor) which is not ideal
S_Shear_Acceptable = S_Shear_Max/2; %4.25 strain here%
Stress_Allowable = 1793*10^6;

% Chord length/Thickness %
for z = 1:1107;
        c(z) = 1.58114.*[1 - 0.0542.*abs(z./1107.*a)];
        t(z) = 0.18.*c (z);
end

% Location of Spar 1/2 %
s1 = 0.20.*c;
s2 = 0.70.*c;

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

L_Eng = 150.*g; %% Located at around 2.9 meters (subject to change)
M_Eng = zeros(1,1107);

L_Wing = ones(1,1107);
L_Wing = 186.3.*g./2./1107.*L_Wing;

%% Moment Derivation %%
for t = 1:1107;
    for t1 = t:1107;
        M_Aero_Pos (t1) = (1107 - t1)./1107.*a.*L_Aero(1108 - t1);
    end 
    for t2 = t:1107;
        M_Batt_Neg (t2) = (1107 - t2)./1107.*a.*L_Batt(1108 - t2);
    end 
    for t3 = t:1107;
        M_Wing_Neg (t3) = (1107 - t3)./1107.*a.*L_Wing(1108 - t3);
    end
    for t4 = t:290;
        M_Eng_Neg (t4) = (290 - t4)./1107.*a.*L_Eng;
    end
    
    M_Aero(t) = sum(M_Aero_Pos(1:1108-t));
    M_Batt(t) = sum(M_Batt_Neg(1:1108-t));
    M_Wing(t) = sum(M_Wing_Neg(1:1108-t));
    if t < 291;
          M_Eng(t) = sum(M_Eng_Neg(1:291-t));
    end
    
    Moment (t) = M_Aero(t) - M_Batt(t) - M_Wing(t) - M_Eng(t);
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
        M_Aero_Pos (t1) = (1107 - t1)./1107.*a.*L_Aero(1108 - t1);
    end 
    for t2 = t:1107;
        M_Batt_Neg (t2) = (1107 - t2)./1107.*a.*L_Batt(1108 - t2);
    end 
    for t3 = t:1107;
        M_Wing_Neg (t3) = (1107 - t3)./1107.*a.*L_Wing(1108 - t3);
    end
    for t4 = t:290;
        M_Eng_Neg (t4) = (290 - t4)./1107.*a.*L_Eng;
    end
    
    M_Aero(t) = sum(M_Aero_Pos(1:1108-t));
    M_Batt(t) = sum(M_Batt_Neg(1:1108-t));
    M_Wing(t) = sum(M_Wing_Neg(1:1108-t));
    if t < 291;
          M_Eng(t) = sum(M_Eng_Neg(1:291-t));
    end
    
    Moment (t) = M_Aero(t) - M_Batt(t) - M_Wing(t) - M_Eng(t);
    end  
end

plot(0:0.01:a, Moment+M_Eng)


