function [ L_Aero ] = Untitled2( x, a );
% Function: Loading for the Aerodynamic Equation

fun = @(y) 7686.632.*sqrt(1 - y.^2./a.^2) + 8338.302.* [1 + (1.58114.*[1 - 0.0542 .* abs(y)] -0.7906)./1.58114];
L_Aero = integral (fun, x , x + 0.01);
L_Aero = L_Aero/1.5;
end

