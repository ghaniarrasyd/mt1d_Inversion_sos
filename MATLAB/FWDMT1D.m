function [apparentResistivity, phase] = FWDMT1D(resistivities, thicknesses, frequency)

for i=1:length(frequency)
mu = 4*pi*1E-7; %Magnetic Permeability (H/m)  
w(i) = 2 * pi * frequency(i); %Angular Frequency (Radians);
n=length(resistivities); %Number of Layers

impedances = zeros(n,1);
%Layering in this format
%  Layer     j
% Layer 1    1
% Layer 2    2
% Layer 3    3
% Layer 4    4
% Basement   5
%

% Steps for modelling (for each geoelectric model and frequency)
% 1. Compute basement impedance Zn using sqrt((w * mu * resistivity))
% 2. Iterate from bottom layer to top(not the basement)
    % 2.1. Calculate induction parameters
    % 2.2. Calculate Exponential factor from intrinsic impedance
    % 2.3 Calculate reflection coeficient using current layer
    %          intrinsic impedance and the below layer impedance
    
% 3. Compute apparent resistivity from top layer impedance
        %   apparent resistivity = (Zn^2)/(mu * w)

%Symbols
% Zn - Basement Impedance
% Zi - Layer Impedance
% wi - Intrinsic Impedance
% di - Induction parameter
% ei - Exponential Factor
% ri - Reflection coeficient
% re - Earth R.C.
        
%Step 1 : Calculate basement impedance  
Zn(i) = sqrt(sqrt(-1)*w(i)*mu*resistivities(n)); 
impedances(n) = Zn(i); 

%Iterate through layers starting from layer j=n-1 (i.e. the layer above the basement)        
for j = n-1:-1:1
    resistivity = resistivities(j);
    thickness = thicknesses(j);
                
    % 3. Compute apparent resistivity from top layer impedance
    %Step 2. Iterate from bottom layer to top(not the basement) 
    % Step 2.1 Calculate the intrinsic impedance of current layer
    dj(i) = sqrt(sqrt(-1)* (w(i) * mu * (1/resistivity)));
    wj(i) = dj(i) * resistivity;
        
    % Step 2.2 Calculate Exponential factor from intrinsic impedance
    ej(i) = exp(-2*thickness*dj(i));                     

    % Step 2.3 Calculate reflection coeficient using current layer
    %          intrinsic impedance and the below layer impedance
    belowImpedance = impedances(j + 1);
    rj(i) = (wj(i) - belowImpedance)/(wj(i) + belowImpedance); 
    re(i) = rj(i)*ej(i); 
    Zj(i) = wj(i) * ((1 - re(i))/(1 + re(i)));
    impedances(j) = Zj(i);               
end
% Step 3. Compute apparent resistivity from top layer impedance
Z(i) = impedances(1);
absZ(i) = abs(Z(i)); 
apparentResistivity(i) = (absZ(i) * absZ(i))/(mu * w(i));
phase(i) = atan2(imag(Z(i)),real(Z(i)))*180/pi;
end
    