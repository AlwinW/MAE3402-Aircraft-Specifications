%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loading/Bending Moment Analysis %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% n = 3.5 %%
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
y = 0;

for x = 0:0.01:a;
    y = y+1;
    L_Aero (y) = Aero_function (x, a);
end

L_Aero = real(L_Aero); 

L_Weight = ones(1,1107);
L_Weight = 2800.*g.*L_Weight;

L_Batt = ones(1,1107);
L_Batt = c./mean(c).*(0.8.*3500 + 300).*g./2./1107.*L_Batt;

L_Eng = zeros(1,1107);
for x1 = 2.75:0.01:3.06;
    L_Eng(x1*100) = 162.4*g/301;
end
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
    for t4 = t:305;
        M_Eng_Neg (t4) = (305 - t4)./1107.*a.*L_Eng(306 - t4);
    end
    
    M_Aero(t) = sum(M_Aero_Pos(1:1108-t));
    M_Batt(t) = sum(M_Batt_Neg(1:1108-t));
    M_Wing(t) = sum(M_Wing_Neg(1:1108-t));
    if t < 306;
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
    for t4 = t:305;
        M_Eng_Neg (t4) = (305 - t4)./1107.*a.*L_Eng(306 - t4);
    end
    
    M_Aero(t) = sum(M_Aero_Pos(1:1108-t));
    M_Batt(t) = sum(M_Batt_Neg(1:1108-t));
    M_Wing(t) = sum(M_Wing_Neg(1:1108-t));
    if t < 306;
          M_Eng(t) = sum(M_Eng_Neg(1:291-t));
    end
    
    Moment (t) = M_Aero(t) - M_Batt(t) - M_Wing(t) - M_Eng(t);
    end  
end

hold on
plot(0:0.01:a, Moment, 'b')
title('Moment distribution over the wing, n = 3.5')
xlabel('Root chord to tip chord')
ylabel('Bending Moment')
%% Net Force Determination %%
Force = L_Aero - L_Batt - L_Wing - L_Eng;
figure
plot(0:0.01:a, 100*Force, 'g')
title('Force distribution over the wing, n = 3.5')
xlabel('Root chord to tip chord')
ylabel('Net Force')

%% Finding Spar Weight Optimisation %%

r1 = 0.0035; r1_2 = 0.0025;
r2 = 0.020; r2_2 = 0.01;
z1 = 0.16; z1_2 = 0.1;
z2 = 0.08; z2_2 = 0.053;

I1 = z1 * (r1.^3)./6 + 2.*z1*r1*[r2+0.5*r1].^2;
I2 = z2*(r2.^3)./12;
I_1 = I1 + I2;

I1_2 = z1_2 * (r1_2.^3)./6 + 2.*z1_2*r1_2*[r2_2+0.5*r1_2].^2;
I2_2 = z2_2*(r2_2.^3)./12;
I_2 = I1_2 + I2_2;
I = (I_1 + (2*r1*z1 + r2*z2)*(0.02603)^2 + I_2 + (2*r1_2.*z1_2 + r2_2*z2_2)*(0.01258)^2);

q = r2 + r2_2 +(r1+r1_2)/2;
    for t = 1:1107;
        Stress(t) = Moment(t).*q./I;
    end

Stress = 1.5.*Stress;
max_s = max(Stress);
calibration = Stress./max(Stress);

r1 = r1*ones(1,1107);
r2 = r2*ones(1,1107);
z1 = z1*ones(1,1107);
z2 = z2*ones(1,1107);

r1_2 = r1_2*ones(1,1107);
r2_2 = r2_2*ones(1,1107);
z1_2 = z1_2*ones(1,1107);
z2_2 = z2_2*ones(1,1107);

r1 = r1.*calibration;
r2 = r2.*calibration;
z1 = z1.*calibration;
z2 = z2.*calibration;

r1_2 = r1_2.*calibration;
r2_2 = r2_2.*calibration;
z1_2 = z1_2.*calibration;
z2_2 = z2_2.*calibration;

for t = 1:1107;
    Area(t) = 2*r1(t).*z1(t) + r2(t).*z2(t) + 2*r1_2(t).*z1_2(t) + r2_2(t).*z2_2(t);
end

M_Wing = mean(Area)*2*a*4420;

M_Wing_Spars_Ti = M_Wing
%% Wing Skin Determination %%

for t = 1:1106;
    Shear_Force(t) = abs(Moment(t+1)-Moment(t))./0.01;
end

M_Wing_Skin = 24.5*2*0.003*1900

%% Ribs Sizing/Weight %%
% Ribs are taken to occur every 0.25 meters %
Rib_Num = round(2*a/0.25);

% Area taken to be approxiamately rectangular, originating from the nose of
% airfoil back to the second spar %
h_Ribs = 0.8.*h;

l_Ribs = 0.65.*c;
Area_Ribs = h_Ribs .* l_Ribs;

Vol_Total_Ribs = mean(Area_Ribs).*0.0024225.*Rib_Num.*2700.*0.85 % 2.4225mm thickness, with a typical/conservative 85%% filled volume %

%% Total Wing Mass %%

Mass_Wing = M_Wing_Skin + Vol_Total_Ribs + M_Wing_Spars_Ti


