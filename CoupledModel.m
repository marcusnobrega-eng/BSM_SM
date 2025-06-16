classdef CoupledModel < handle
    % CoupledModel 2.5-D
    % Calculates 1d unsaturated flow solving the Richards equation for
    % unsaturated soil columns and the depth averaged 2d groundwater
    % equation for the saturated part

    properties
        %         uzmodel         % unsaturated soil solumns (EDIT)
        gwmodel         % groundwater model
        xi              % cells of the groundwater model with soil column
        yi              % cells of the groundwater model with soil column
        nrStates        % number of model states
        eps = 1e-6;     % convergence criterion for coupling
        dH_R            % groundwater table fluctuation due to recharge
        zoneMap         % maps groundwater zones and soil columns
        maxIter = 41;   % maximum number of iterations for coupling
        discharge = 0.0;
        storage = 0.0;
    end

    methods
        %% constructor
        function obj = CoupledModel(gwmodel)
            %             obj.uzmodel = uzmodel; %EDIT
            obj.gwmodel = gwmodel;
            %             obj.xi = xi;
            %             obj.yi = yi;
            obj.dH_R = zeros(gwmodel.ny,gwmodel.nx);
            %             obj.zoneMap = zoneMap;
            obj.nrStates = getNumberOfStates(obj);
            obj.discharge = 0.0;
            obj.storage = 0.0;
        end
        %% getNumberOfStates
        function nr = getNumberOfStates(obj)
            nr = obj.gwmodel.nrStates;
        end
        %% advance
        function obj = advance(obj,tStart,tStop)

            %% initialize iteration variables
            lost_mass_mm = 0;
            % differences of groundwater table elevation
            dh = 1;
            % iteration counter
            it = 0;

            % groundwater table elevation of groundwater model
            h_GW_last = obj.gwmodel.h; % Groundwater Elevation from Datum Reference at the beginning
            
            % specific yield
            Sy = obj.gwmodel.parameter.Sy;

            % Implementing Depth-dependent Specific Yield
            % Hillberts (2005) formulation
            % -- VG Parameters
            porosity = obj.gwmodel.parameter.Poro; 
            alpha = obj.gwmodel.parameter.Alpha;
            n = obj.gwmodel.parameter.N;            
            % alpha = -0.9;
            % n = 1.7608;
            % f'(h) = (theta_s - theta_r)(1  - (1 + (alpha(h -
            % z))^n)^((n+1)/n)
            Sy = porosity.*(1 - (1 + (alpha.*(max(h_GW_last - obj.gwmodel.z0,0) - obj.gwmodel.parameter.H)).^n).^(-((n+1)./n)));
            % Sy = porosity; % DELETE
            % Sy_last = Sy;
            % Sy(logical(obj.gwmodel.boundaries.SeepageFace)) = porosity(1,1);
            % Sy = max(Sy,0.05);
            % time step size
            dt = tStop-tStart;

            %% Run groundwater model forward

            % update specific yield
            obj.gwmodel.parameter.updateSy(Sy);

            % calculate groundwater recharge (no flooding)
            Rflooded = (obj.gwmodel.h+obj.dH_R)>obj.gwmodel.lsurf;
            obj.dH_R(Rflooded) = obj.gwmodel.lsurf(Rflooded) - ...
                obj.gwmodel.h(Rflooded);
            R = obj.dH_R(:).*0; % EDIT. Recharge is null and Q is actually the recharge value
            dirchlet_cells = obj.gwmodel.boundaries.perimI;
            dirchlet_cells(isnan(dirchlet_cells)) = 0;
            Perim_Dir = logical(dirchlet_cells);

            % run groundwater model with Q
            % Correct Q to consider delay in recharge (Method 1)
            idx = find(obj.gwmodel.boundaries.changeTQ<= tStart,1,'last');
            % Delay in Recharge
            % R = P*exp(-beta_r*(d_r - h)), h <= d_r
            % R = P, h > d_r
            % Acharya et al., (2012)
            % Loamy Sand
            % beta_r = 0.08 cm, d_r = 45 cm
            % h: water depth [cm]; P: rainfall rate [cm/h]
            % beta_r = 0.16; 
            % d_r = 40/100; % m
            % depth = max(h_GW_last - obj.gwmodel.z0,0); % m
            % factor = exp(-beta_r.*(max(d_r - depth,0)));   
            % Recharge = P.*factor; % Recharge            
            % --- Changing P based on the observed runoff coefficient
            damping_factor = obj.gwmodel.parameter.Irr_d; % Based on observations
            % damping_factor = 1;
            P = damping_factor*obj.gwmodel.boundaries.Q(:,:,idx); % Precipitation in m/s

            % --- Running a linear reservoir model
            k_UZ_reservoir = obj.gwmodel.parameter.k;
            % k_UZ_reservoir = 1e-4*ones(obj.gwmodel.ny,obj.gwmodel.nx); % You can calibrate
            
            % Delay Function
            [Recharge,obj.gwmodel.parameter.S_UZ] = irrigation_delay(P,k_UZ_reservoir,dt,obj.gwmodel.parameter.S_UZ);

            % Refreshing Recharge
            obj.gwmodel.boundaries.Q(:,:,idx) = Recharge; % New Recharge Rate

            % Run Model Forward
            h_GW = obj.gwmodel.advance(tStart,tStop,-R);
            zzz = obj.gwmodel.z0; zzz = zzz(:); zzz(~isnan(zzz)) = h_GW;
            h_GW = zzz;
            h_GW = reshape(h_GW,obj.gwmodel.ny,obj.gwmodel.nx);

            [X,Y] = meshgrid(1:1:obj.gwmodel.nx,1:1:obj.gwmodel.ny);
            X = X*obj.gwmodel.dx;
            Y = Y*obj.gwmodel.dy;
            % for i = 1:length(obj.gwmodel.boundaries.perimIdir)
            %     dirchlet_mask = (obj.gwmodel.boundaries.perimI == obj.gwmodel.boundaries.perimIdir(i));
            %     Perim_Dir = logical(Perim_Dir + dirchlet_mask);
            % end

            % Check Flooding
            % eliminate surface water
            flooding = h_GW>obj.gwmodel.lsurf;
            h_GW(flooding) = obj.gwmodel.lsurf(flooding);
            % if sum(sum(flooding)) > 0
            %     error('Flooding. Flooded Volume (mm/h)');
            % end
            % 
            % [~] = surf_plot(max(max(obj.gwmodel.parameter.H)),tStop/86400,'h_g','m',h_GW - obj.gwmodel.z0,0,1,256,0.9,0,[0,90],X,Y);
            % pause(0.0000001);

            % stop iteration if groundwater model crashed
            % EDIT: throw an error instead
            if sum(sum(isnan(h_GW))) == numel(h_GW) % Edit Marcus
                error('GW model crashed!')
                %break;
            end


            % sprintf('Mass lost (mm) in a GW time-step: %f', lost_mass_mm)
            % pause(0.15)
            %% end of iteration

            % Discharge Calculation
            mask = ~isnan(obj.gwmodel.z0);
            Acell = obj.gwmodel.dx*obj.gwmodel.dy;
            idx = find(obj.gwmodel.boundaries.changeTQ<= tStart,1,'last');

            Sy_actual = porosity.*(1 - (1 + (alpha.*(max(h_GW - obj.gwmodel.z0,0) - obj.gwmodel.parameter.H)).^n).^(-((n+1)./n)));
            h_GW_half = 1/2*(h_GW + h_GW_last); % Head at 1/2 time-step
            Sy_ = porosity.*(1 - (1 + (alpha.*(max(h_GW_half - obj.gwmodel.z0,0) - obj.gwmodel.parameter.H)).^n).^(-((n+1)./n)));            
            % Sy_ = porosity; % DELETE
            % surf(Sy_);
            % pause(0.00001)
            obj.discharge = mean(mean(-obj.gwmodel.boundaries.Q(:,:,idx) - (h_GW - h_GW_last).*Sy_/dt)); % m/s
            % obj.discharge = mean(R(Perim_Dir(:)) - (h_GW(Perim_Dir(:)) - obj.gwmodel.parameter.Sy(Perim_Dir(:)).*obj.gwmodel.h(Perim_Dir(:)))/dt)
            discharge_mmh = obj.discharge*1000*3600; % mm/h
            obj.discharge = max(0,obj.discharge); % In cases discharge becomes positive, we assume 0

            if discharge_mmh > 100
                error('Discharge too high');
            end

            % Vol_last = Acell*sum(sum((h_GW_last(mask) -obj.gwmodel.z0(mask)).*Sy_(mask))); % m3
            % Recharge_flux = (-obj.gwmodel.boundaries.Q(:,:,idx))*Acell*dt; % m3 per cell
            % Recharge = sum(sum(Recharge_flux(mask))); % m3 per domain
            % Vol_t = Acell*sum(sum((h_GW(mask) -obj.gwmodel.z0(mask)).*Sy_(mask))); % m3 per domain
            % Outflow = (Recharge - (Vol_t - Vol_last))/dt; % m3/s
            % obj.discharge = Outflow/(Acell*numel(h_GW)); % m/sl
            % discharge_mmh_2 = (1000*3600)*Outflow/(Acell*numel(h_GW)); % m/s

            % update groundwater table elevation in groundwater model
            obj.gwmodel.h = h_GW;
            % update groundwater model parameters (if needed)
            obj.gwmodel.parameter.updateParam(h_GW-obj.gwmodel.z0);

            % crashed
            if any(isnan(h_GW))
                % stop simulation
                % EDIT: throw error instead
                error('--------------Coupling did not converge, breaking out!')
                % too many iteration steps
            end

        end
        %% getPressureHeads
        function P = getPressureHeads(obj)
            P = zeros(obj.nrStates,1);
            counter = 1;
            P(counter:end) = obj.gwmodel.getPressureHeads();
        end

        %% getDischarge
        function discharge = getDischarge(obj)
            discharge=obj.discharge;
        end

        %% getStorage
        function storage = getStorage(obj)
            storage=obj.storage;
        end
    end

end

