classdef Roots < handle
    %Root parameters for water uptake
    
    properties
        rootWaterUptake     % do root water uptake?
        rootDistribution    % file: load from file; exp: exponential function
        rootDensity         % root density 
        rootLength          % depth of roots in soil
        feddes              % values of the feddes function
    end
    
    methods
        %% constructor
        function obj = Roots(doUptake,dist,charDistValue,z,H,delzCS,feddes)
            obj.rootWaterUptake = doUptake;
            if doUptake
                obj.rootDistribution = dist;
                if strcmp(dist,'file')
                    density = charDistValue;
                elseif strcmp(dist,'exp')
                    obj.rootLength = charDistValue;
                    density = exp(-(H-z)/charDistValue);
                end
                % normalize, sum has to be equal to 1
                obj.rootDensity = density./(sum(density)*delzCS);

                obj.feddes = feddes;
            end
        end
    end
    
end

