classdef UZParameterization < handle
    %Parameterization Soil hydraulic parameters
    
    properties
        model               % retention model (VG: van Genuchten, RG: Russo Gardner)
        K                   % saturated hydraulic conductivity
        Alpha               % model parameter
        N                   % model parameter
        aKR                 % model parameter
        Ssat                % maximal saturation
        Sres                % residual saturation
        poro                % porosity
        storage             % storage capacity
        sDiff               % effective saturation range
        M                   % model parameter
        KrelWeight          % How are the relative permeabilities weighted? 
                            % (1: upstream; 2: arithmetic; 3: harmonic; 4: geometric)
        KrelMin = 1e-10;    % Minimum value of relative permeability 
        boundaryKR          % Krel at the lower boundary
        Ksat                % harmonic averaged saturated conductivity
        useK_GW             % use saturated hydraulic conductivity of groundwater model
    end
    
    methods
        % constructor
        function obj = UZParameterization(model,K,Alpha,N,aKR,Ssat,Sres,poro,...
                storage,useK_GW)
            obj.KrelWeight = 1;
            obj.model = model;
            obj.K = K;
            obj.Alpha = Alpha;
            obj.N = N;
            obj.aKR = aKR;
            obj.Ssat = Ssat;
            obj.Sres = Sres;
            obj.poro = poro;
            obj.storage = storage;
            obj.sDiff = Ssat-Sres;
            obj.M = 1-1./N;
            obj.useK_GW = useK_GW;
            getHarmonicAverages(obj);
        end
        %% calculateBoundaryKR
        function obj = calculateBoundaryKR(obj,bHead)
            if strcmp(obj.model,'VG')
                % van Genuchten model
                aB=obj.Alpha(1)*abs(bHead);
                bB=(1+aB^obj.N(1))^obj.M(1);
                obj.boundaryKR=((1 -((aB)^(obj.N(1)-1))/bB)^2)*(bB^obj.aKR(1));
            elseif strcmp(obj.model,'RG')
                % Russo-Gardener model
                obj.boundaryKR=exp(-obj.Alpha(1)*abs(bHead));
            end
        end
        %% getHarmonicAverages
        function obj = getHarmonicAverages(obj)
            obj.Ksat=zeros(length(obj.K)+1,1);
            obj.Ksat(2:end-1)=2*obj.K(2:end).*obj.K(1:end-1)./(obj.K(2:end)+obj.K(1:end-1));
            % bottom and top are the originals
            obj.Ksat(1)=obj.K(1);
            obj.Ksat(end)=obj.K(end);
        end
        %% updateK
        function obj = updateK(obj,iEnd,K_GW)
            % update saturated K in saturated zone
            obj.K(:) = obj.K(end);
            obj.K(1:iEnd) = K_GW;
        end
    end
    
end

