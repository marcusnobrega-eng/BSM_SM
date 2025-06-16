classdef GWParameterization < handle
    %Parameterization Soil hydraulic parameters
    
    properties
        K       % saturated hydraulic conductivity
        Sy      % specific yield
        H       % Groundwater Tickness (EDIT)
        K_3D    % 3D satuared hydraulic conductivity field
        dz_3D   % 3D grid cell sizes 
        z_3D    % 3D grid coordinates of upper cell boundaries
        N       % Van-Genutchen n parameter
        Poro    % UZ Porosity [-]
        Alpha   % UZ Alpha [1/m]
        S_UZ    % Current UZ storage [mm]
        k       % UZ linear reservoir damping factor [1/sec]
        Irr_d   % Irrigation damping: Irrigation = Irrigation * Irr_d
    end
    
    methods
        %% constructor 
        function obj = GWParameterization(K,Sy,H,N,Poro,Alpha,S_UZ,S_k,Irr_d) % EDIT
            obj.K = K;
            obj.Sy = Sy;
            obj.H = H; % EDIT
            obj.N = N; % EDIT
            obj.Poro = Poro; % EDIT
            obj.Alpha = Alpha; % EDIT
            obj.S_UZ = S_UZ; % EDIT
            obj.k = S_k; % EDIT
            obj.Irr_d = Irr_d; % EDTI
        end
        %% set3D_Data
        function set3D_data(obj, K, dz)
            obj.K_3D = K;
            obj.dz_3D = dz;
            obj.z_3D = cumsum(dz,3)-dz/2;
        end
        %% update specific yield
        function updateSy(obj,Sy)
            obj.Sy = Sy;
            % apply minimum value
            obj.Sy(obj.Sy<5e-5) = 5e-5;
        end
        %% updateParam
        function updateParam(obj,h)
            % if 3D fields are given
            if ~isempty(obj.dz_3D)
                % get saturated cells and cell sizes
                z_sat = repmat(h,1,1,size(obj.dz_3D,3))-obj.z_3D;
                dz_sat = obj.dz_3D;
                dz_sat(z_sat<0) = NaN;
                % average ksat in the saturated zone according to
                % cell sizes
                obj.K = obj.K_3D.*dz_sat;
                obj.K = (nansum(obj.K,3)./nansum(dz_sat,3)); 
                K_bottom = obj.K_3D(:,:,1);
                obj.K(isnan(obj.K)) = K_bottom(isnan(obj.K));
            end
        end
    end
    
end

