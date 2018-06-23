classdef UnitConvert
    
%UnitConvert
%Author       : Stepen (stepen3112@gmail.com)
%Version      : 1.00
%Dependencies : N/A
%UnitConvert class provides constants to convert various unit of dimension.

%% ConstantPropertiesDeclaration------------------------------------------%
%Declaring constants for time unit conversion
properties(Constant = true)
    
    %Conversion factor from second to hour
    S2H = 1/3600;
    
end

%Declaring constants for length unit conversion
properties(Constant = true)
    
    %Conversion factor from meter to kilometer
    M2KM = 1/1000;
    %Conversion factor from meter to feet
    M2FT = 3.28084;
    %Conversion factor from kilometer to mile
    KM2MI = 0.621371;
    %Conversion factor from kilometer to nautical mile
    KM2NMI = 0.539957;
    
end

%Declaring constants for velocity unit conversion
properties(Constant = true)
    
    %Conversion factor from meter/second to knots
    MPS2KTS = UnitConvert.M2KM*UnitConvert.KM2NMI/(UnitConvert.S2H);
    
end
% EndOfConstantPropertiesDeclaration---------------------------------------%

end