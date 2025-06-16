classdef UZModel < handle
    %UZModel Calculates 1d unsaturated flow solving the Richards equation
    
    properties
        % geometry
        H               % height of the soil column
        % discretization
        nZ              % number of cells
        delzCS          % height of the cells
        delzNN          % distances between the cell centered nodes
        z               % level of the cell centered nodes
        nrStates        % number of model states
        % parameter
        parameter       % soil parameter
        % roots
        roots           % parameter for root water uptake
        % sources and sinks
        ss              % sources and sinks   
        % boundary conditions
        boundaries      % contains all boundary conditions and timing
        % initial condition
        IC              % initial condition
        IC_GW           % get full initial condition from groundwater data
        % states
        h               % states (pressure heads)
        % water content
        WC              % water content
        % solver
        solver          % numerical solver
        % timing
        dt              % time step size
        doAdaptiveTime  % adaptive time step size  
        dtMax = 3600*24;% max time step size
        dtMin = 1;      % min time step size
    end
    
    methods
        %% constructor
        function obj = UZModel(H, delzCS, parameter, boundaries, ss, IC, IC_GW, ...
                solver, dt, doAdaptiveTime)
            obj.H = H;
            obj.delzCS = delzCS;
            obj.nZ = length(delzCS);
            [obj.delzNN, obj.z] = obj.initDiscretization();
            obj.parameter = parameter;
            obj.boundaries = boundaries;
            obj.ss = ss;
            obj.IC = IC;
            obj.IC_GW = IC_GW;
            obj.solver = solver;
            obj.dt = dt;
            obj.doAdaptiveTime = doAdaptiveTime;
            obj.h = IC;
            obj.nrStates = length(IC);
            if ~isempty(IC)
                obj.calculateWC();
            end
        end
        %% initDiscretization
        function [nn, z] = initDiscretization(obj)
            nn=zeros(obj.nZ+1,1);
            nn(1)=obj.delzCS(1)/2;
            nn(2:end-1) = (obj.delzCS(1:end-1)+obj.delzCS(2:end))/2;
            nn(end)=obj.delzCS(end)/2;
            z = cumsum(nn(1:end-1));
        end
        %% setIC
        function obj = setIC(obj, h, minValue)
            % groundwater table elevation given
            % assume steady state above
            if obj.IC_GW
                obj.IC = h-obj.z;
                obj.h = h-obj.z;
            % pressure heads are given directly
            else
                obj.IC = h;
                obj.h = h;
            end
            % apply minimum value
            obj.IC(obj.h<minValue) = minValue;
            obj.h(obj.h<minValue) = minValue;
            
            % update depending variables
            obj.nrStates = length(obj.IC);
            obj.calculateWC();
        end
        %% setRoots
        function setRoots(obj,doUptake,dist,charDistValue,feddes)
            obj.roots = Roots(doUptake,dist,charDistValue,obj.z,obj.H,...
                obj.delzCS,feddes);
        end
        %% calculateWC
        function calculateWC(obj)
            [sat,~] = calcSat(obj.h,obj.parameter);
            obj.WC = sat.*obj.parameter.poro;
            obj.WC(isinf(obj.h)|isnan(obj.h)) = NaN;
        end
        %% advance
        function [h_new,R,hu] = advance(obj,tStart,tStop,source,hGW)
            % run model forward
            [h_new,R,hu] = obj.solver.solveUZ(tStart,tStop,source,obj.dt,obj,hGW);
        end
        %% getH
        function H = getH(obj)
            % calculate position of groundwater table (linear interpolation)
            full = false;
            if obj.h(1)>0
                hi = find(obj.h<0,1,'first');
                if isempty(hi)
                    full = true;
                else
                    hs = obj.h(hi-1);
                    zs = obj.z(hi-1);
                    hu = obj.h(hi);
                    zu = obj.z(hi);
                end
            else
                hs = obj.boundaries.currHeadB;
                zs = obj.z(1)-obj.delzCS(1)/2;
                hu = obj.h(1);
                zu = obj.z(1);
            end
            if ~full
                H = zs+(zu-zs)*(0-hs)/(hu-hs);
            else
                H = obj.z(end)+obj.h(end);
            end
        end
        %% updateK
        function obj = updateK(obj,K_GW)
            % update saturated K in saturated zone
            if obj.parameter.useK_GW
                iEnd = find(obj.z<obj.getH(),1,'last');
                if (isempty(iEnd) || iEnd == obj.nZ)
                    iEnd = obj.nZ-1;
                end
                obj.parameter.updateK(iEnd,K_GW);
            end
        end
        %% getPressureHeads
        function P = getPressureHeads(obj)
            P = obj.h;
        end
        %% getWaterContent
        function WC = getWaterContent(obj)
            obj.calculateWC();
            WC = obj.WC;
        end
    end
    
end

