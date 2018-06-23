classdef Atmosphere

%Atmosphere
%Author       : Stepen (stepen3112@gmail.com)
%Version      : 1.00
%Dependencies : UnitConvert.m
%Atmosphere class provides functions to calculate atmosphere properties
%based on ISA model (US 1976) and other atmospheric related operations such
%as airspeed conversions.

%% ConstantPropertiesDeclaration------------------------------------------%
%Declaring constants for ISA atmosphere properties
properties(Constant = true)
    
    %Specific Heat Capacity Ratio
    GAMMA = 1.4;
    %Ideal Gas Constant in N.m/kk.K
    R     = 287.05287;
    %Gravity Acceleration in m/sec^2
    G0    = 9.80665;
    %Standard Atmosphere Property Table (US 1976)
    PROP  = table([0;11000;20000;32000;47000;51000;71000],... %m
                  [11000;20000;32000;47000;51000;71000;84000],... %m
                  [288.15;216.65;216.65;228.65;270.65;270.65;214.65],... %K
                  [-0.0065;0;0.001;0.0028;0;-0.0028;-0.002],... %K/m
                  [101325;22632.0400950078;5474.87742428105;...
                   868.015776620216;110.90577336731;66.9385281211797;...
                   3.9563921603966],... %Pa
                  'VariableNames',{'LowerAltitude';...
                                   'UpperAltitude';...
                                   'BaseTemperature';...
                                   'TemperatureLapseRate';
                                   'BasePressure'},...
                  'RowNames',{'Troposphere';'Tropopause';...
                              'Stratosphere1';'Stratosphere2';...
                              'Stratopause';'Mesosphere1';
                              'Mesosphere2'});
    
end

%Declaring constants for user interface settings
properties(Constant=true,Hidden=true)
    
    %RESERVED
    
end
% EndOfConstantPropertiesDeclaration--------------------------------------%

%% ConstructorMethodDeclaration-------------------------------------------%
methods
    
    function this = Atmosphere()
    %Draws a plot of basic atmosphere properties over altitude changes.
        
        %Executing default methods to generate atmosphere chart
        Atmosphere.drawAtmosphere();
        %Suppress output argument
        clear('this');
        
    end
    
end
% EndOfConstructorMethodDeclaration---------------------------------------%

%% StaticMethodsDeclaration-----------------------------------------------%
methods(Static = true)

    function layerId = layer(h)
    % Finds atmospheric layers for given geopotential altitudes.
    % layerId = layer(h) returns atmospheric layer index for given
    % geopotential altitude h (in meters).
    
        %Checking for input argument with row array shape
        rowIn = isrow(h);
        %Transposing input argument with row array shape
        if rowIn
            h = h';
        end
        %Determining atmospheric layer for each geopotential altitudes
        layerId = zeros(size(h));
        for id_layer = 1:height(Atmosphere.PROP)
            layerId((h>=Atmosphere.PROP.LowerAltitude(id_layer))&...
                    (h<Atmosphere.PROP.UpperAltitude(id_layer))) = ...
                id_layer;
        end
        %Transposing output argument to match input argument's shape
        if rowIn
            layerId = layerId';
        end
    
    end
    
    function T = T(h)
    % Calculates ISA temperatures for given geopotential altitudes.
    % T = T(h) returns ISA temperatures (in Kelvin) for given geopotential
    % altitudes h (in meters).

        %Checking for input argument with row array shape
        rowIn = isrow(h);
        %Transposing input argument with row array shape
        if rowIn
            h = h';
        end
        %Determining atmosphere layer
        layer = Atmosphere.layer(h);
        %Calculating temperature
        T = Atmosphere.PROP.BaseTemperature(layer) + ...
              (Atmosphere.PROP.TemperatureLapseRate(layer).*...
               (h-Atmosphere.PROP.LowerAltitude(layer)));
        %Transposing output argument to match input argument's shape
        if rowIn
            T = T';
        end
        
    end

    function p = p(h)
    % Calculates ISA ambient pressures for given geopotential altitudes.
    % p = p(h) returns ISA ambient pressures (in Pascal) for given
    % geopotential altitudes h (in meters).
        
        %Checking for input argument with row array shape
        rowIn = isrow(h);
        %Transposing input argument with row array shape
        if rowIn
            h = h';
        end
        %Determining atmosphere layer
        layer = Atmosphere.layer(h);
        %Finding layer with no temperature lapse rate
        nolapse = (Atmosphere.PROP.TemperatureLapseRate(layer)==0);
        %Preallocating array for output
        p = zeros(size(h));
        %Calculating pressure
        p(nolapse)  = Atmosphere.PROP.BasePressure(layer(nolapse)).*...
                      exp(-Atmosphere.G0.*...
                          (h(nolapse)-...
                           Atmosphere.PROP.LowerAltitude(...
                               layer(nolapse)))./...
                          (Atmosphere.PROP.BaseTemperature(...
                               layer(nolapse)).*...
                           Atmosphere.R));
        p(~nolapse) = Atmosphere.PROP.BasePressure(layer(~nolapse)).*...
                      ((1+(Atmosphere.PROP.TemperatureLapseRate(...
                               layer(~nolapse)).*...
                           (h(~nolapse)-...
                            Atmosphere.PROP.LowerAltitude(...
                                layer(~nolapse)))./...
                           Atmosphere.PROP.BaseTemperature(...
                               layer(~nolapse)))).^...
                       (-Atmosphere.G0./...
                        (Atmosphere.PROP.TemperatureLapseRate(...
                             layer(~nolapse)).*...
                         Atmosphere.R)));
        %Transposing output argument to match input argument's shape
        if rowIn
            p = p';
        end
        
    end

    function delta = delta(h)
    % Calculates ISA pressure ratios for given geopotential altitudes.
    % delta = delta(h) returns ISA pressure ratios for given geopotential
    % altitudes h (in meters).
        
        delta = Atmosphere.p(h)./Atmosphere.p(0);
        
    end

    function rho = rho(h)
    % Calculates ISA air densities for given geopotential altitudes.
    % rho = rho(h) returns ISA air densities (in kg/m^3) for given
    % geopotential altitudes h (in meters).
        
        rho = Atmosphere.p(h)./(Atmosphere.R*Atmosphere.T(h));
        
    end

    function sigma = sigma(h)
    % Calculates ISA density ratios for given geopotential altitudes.
    % delta = delta(h) returns ISA density ratios for given geopotential
    % altitudes h (in meters).
        
        sigma = Atmosphere.rho(h)./Atmosphere.rho(0);
        
    end

    function a = a(h)
    % Calculates speeds of sound at given geopotential altitudes.
    % a = a(h) returns speed of sounds (in m/s) at given geopotential
    % altitudes h (in meters).
    
        a = sqrt(Atmosphere.GAMMA.*Atmosphere.R.*Atmosphere.T(h));
        
    end

    function tas = mach2tas(mach,h)
    % Converts Mach number to TAS at given geopotential altitudes.
    % tas = mach2tas(mach,h) returns true airspeed conversion (in m/s) of
    % given Mach number, mach, at given geopotential altitudes h (in
    % meters).
        
        %TODO: Input size check
        tas = mach.*Atmosphere.a(h);
        
    end

    function mach = tas2mach(tas,h)
    % Converts TAS to Mach number at given geopotential altitudes.
    % mach = tas2mach(tas,h) returns Mach number conversion of given true
    % airspeed tas (in m/s) at given geopotential altitudes h (in
    % meters).
    
        %TODO: Input size check
        mach = tas./Atmosphere.a(h);
        
    end

    function eas = tas2eas(tas,h)
    % Converts TAS to EAS at given geopotential altitudes.
    % eas = tas2eas(tas,h) returns equivalent airspeed conversion of given
    % true airspeed tas at given geopotential altitudes h (in meters). Note
    % that eas returned will be in the same unit as the given tas unit.
        
        %TODO: Input size check
        eas = tas.*sqrt(Atmosphere.sigma(h));
        
    end

    function tas = eas2tas(eas,h)
    % Converts EAS to TAS at given geopotential altitudes.
    % tas = eas2tas(eas,h) returns true airspeed conversion of given
    % equivalent airspeed eas at given geopotential altitudes h (in
    % meters). Note that tas returned will be in the same unit as the given
    % eas unit.
        
        %TODO: Input size check
        tas = eas./sqrt(Atmosphere.sigma(h));
        
    end

    function [cas,ccv] = eas2cas(eas,h)
    % Converts EAS to CAS at given geopotential altitudes.
    % [cas,ccf] = eas2cas(eas,h) returns calibrated airspeed conversion of
    % given equivalent airspeed eas at geopotential altitudes h (in
    % meters) and its compressibility correction value ccv. Note that cas
    % and ccv returned will be in the same unit as the given eas unit.
        
        %TODO: Input size check
        M   = Atmosphere.tas2mach(Atmosphere.eas2tas(eas,h),h);
        d   = Atmosphere.delta(h);
        cas = eas.*(1+...
                    ((1-d).*(M.^2)/8)+...
                    (3*(1-(10.*d)+(9.*d.^2)).*(M.^4)/640));
        ccv = eas - cas;
        
    end

    function out = cas2eas(cas,h)
    % Converts CAS to EAS at given geopotential altitudes.
        
        %TODO: To be added
        out = [];
        
    end

    function DrawAtmosphere
    % Draws a plot of basic atmosphere properties over altitude changes.
        
        %Declaring constants for chart generation
        WINDOW_AR     = 0.5;
        WINDOW_MARGIN = 0.2;
        UI_MARGIN_PX  = 40;
        UI_AXES_LOWERPOS_NORM = [0,0.1,0.2];
        UI_AXES_LIMIT_T       = [200,300];
        UI_AXES_LIMIT_P       = [0,110000];
        UI_AXES_LIMIT_RHO     = [0,1.25];
        %Generating chart figure
        screenSize   = get(0,'ScreenSize');
        windowWidth  = WINDOW_AR*(1-WINDOW_MARGIN)*screenSize(4);
        windowHeight = (1-WINDOW_MARGIN)*screenSize(4)*0.75;
        uiWindow = figure('Menubar','none',...
                          'WindowStyle','normal',...
                          'DockControls','off',...
                          'Units','Pixels',...
                          'Position',[0.5*(screenSize(3)-windowWidth),...
                                      0.5*(screenSize(4)-windowHeight),...
                                      windowWidth,windowHeight]);
        %Generating chart axes
        drawableArea = get(uiWindow,'InnerPosition');
        usableWidth  = 0.7*drawableArea(3) - (2*UI_MARGIN_PX);
        usableHeight = 0.7*drawableArea(4) - (2*UI_MARGIN_PX);
        uiAxesMain = axes('Parent',uiWindow,...
                          'XColor','r',...
                          'XLim',UI_AXES_LIMIT_T,...
                          'Units','pixels',...
                          'Position',[UI_MARGIN_PX,...
                                      UI_MARGIN_PX+...
                                      (UI_AXES_LOWERPOS_NORM(3)*...
                                       usableHeight),...
                                      usableWidth,...
                                      (1-UI_AXES_LOWERPOS_NORM(3))*...
                                       usableHeight]);
        xlabel('Outside Air Temperature (K)');
        ylabel('Geopotential Altitude (m)');
        uiAxesSec1 = axes('Parent',uiWindow,...
                          'XColor','g',...
                          'XLim',UI_AXES_LIMIT_P,...
                          'Units','pixels',...
                          'Position',[UI_MARGIN_PX,...
                                      UI_MARGIN_PX+...
                                      (UI_AXES_LOWERPOS_NORM(2)*...
                                       usableHeight),...
                                      usableWidth,...
                                      1]);
        xlabel('Ambient Pressure (Pa)');
        uiAxesSec2 = axes('Parent',uiWindow,...
                          'XColor','b',...
                          'XLim',UI_AXES_LIMIT_RHO,...
                          'Units','pixels',...
                          'Position',[UI_MARGIN_PX,...
                                      UI_MARGIN_PX+...
                                      (UI_AXES_LOWERPOS_NORM(1)*...
                                       usableHeight),...
                                      usableWidth,...
                                      1]);
        xlabel('Air Density (kg/m3)');
        %Generating basic atmosphere properties data for plotting
        h     = 0:100:80000;
        T     = Atmosphere.T(h);
        p     = Atmosphere.p(h);
        p_s   = UI_AXES_LIMIT_T(1) + ...
                (diff(UI_AXES_LIMIT_T)/diff(UI_AXES_LIMIT_P)*...
                 (p-UI_AXES_LIMIT_P(1)));
        rho   = Atmosphere.rho(h);
        rho_s = UI_AXES_LIMIT_T(1) + ...
                (diff(UI_AXES_LIMIT_T)/diff(UI_AXES_LIMIT_RHO)*...
                 (rho-UI_AXES_LIMIT_RHO(1)));
        %Drawing atmosphere properties trend lines
        line(uiAxesMain,T,h,'Color','r');
        line(uiAxesMain,p_s,h,'Color','g');
        line(uiAxesMain,rho_s,h,'Color','b');
        
    end

    function DrawCompressibilityCorrection(varargin)
    % Draws compressibility correction chart for EAS determination.
        
        %Parsing input arguments
        switch nargin
            case 0
                h_ft = [0:1000:4000,5000:5000:60000]; % ft
                mach = 0:0.025:0.8;
        end
        %Generating meshgrid
        [grid_h,grid_M] = meshgrid((1/UnitConvert.M2FT)*h_ft,mach);
        %Calculating equivalent and calibrated airspeed
        eas_mps = Atmosphere.tas2eas(...
                      Atmosphere.mach2tas(grid_M,grid_h),...
                      grid_h);
        cas_mps = Atmosphere.eas2cas(eas_mps,grid_h);
        %Converting airspeed and calculating its difference
        eas_kts = eas_mps*UnitConvert.MPS2KTS;
        cas_kts = cas_mps*UnitConvert.MPS2KTS;
        dV_kts  = cas_kts - eas_kts;
        %Generating chart
        figure();
        hold on;
        line(cas_kts,dV_kts,'Color','k');
        hold off;
        
    end

end
% EndOfStaticMethodsDeclaration-------------------------------------------%

end