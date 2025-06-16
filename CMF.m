classdef CMF < handle
    %CMF Coupled Model Framework

    properties
        % versioning
        author = 'Brandhorst';  % original author
        lastEdit = '11/20';     % month of last update
        version = '2.0';        % actual version
        % for saving
        name                    % run name
        saveOption              % specifies how and which variables are saved
        tSave                   % times at which the results are saved
        % timing
        tStart                  % start time of the forward models
        tEnd                    % end time of the forward models
        t                       % current time
        % model
        coupledModel            % coupled model
        % matrices for saving
        pressureHeadStore       % pressure heads
        SyStore                 % specific yield
        dischargeStore          % Discharge Store
        storageStore            % Storage Store
        parameterStore          % Ksat (for particle filter only)
        dischargeStore_PF       % Discharge stored (for particle filter only)       
        GW_Depth_Store_PF       % GW Depth (for particle filter only)
        Sy_Store_PF             % Particle Filter Sy (for particle filter only)
    end

    methods
        %% constructor
        function obj = CMF(name, saveOption, tSave, tStart, tEnd, coupledModel)
            obj.name = name;
            obj.saveOption = saveOption;
            obj.tSave = tSave;
            obj.tStart = tStart;
            obj.t = obj.tStart;
            obj.tEnd = tEnd;
            obj.coupledModel = coupledModel;
            obj.pressureHeadStore = zeros(coupledModel.nrStates,length(tSave));
            obj.SyStore = zeros(coupledModel.gwmodel.nrStates,length(tSave));
            obj.dischargeStore = zeros(1,length(tSave));
            obj.storageStore = zeros(1,length(tSave));
            obj.parameterStore = [];          
            obj.dischargeStore_PF = [];  
            obj.GW_Depth_Store_PF = [];
            obj.Sy_Store_PF = [];
        end
        %% run and save hydrograph (Edit Marcus)
        function [Qmod,GWdepth,dx,dy,DX,DY,nx,ny] = run_boussinesqModel(obj)
            % for all time steps
            %             try

            for i=1:length(obj.tSave)
                tic
                % start and end time of current time step
                time = obj.t;
                ti = obj.tSave(i);
                disp([obj.name ' - ' num2str(ti/(3600*24)) 'd/' num2str(obj.tSave(end)/(3600*24)) 'd'])

                % advance coupled model
                obj.coupledModel.advance(time,ti);

                % store model states
                obj.pressureHeadStore(:,i) = obj.coupledModel.getPressureHeads();
                obj.SyStore(:,i) = obj.coupledModel.gwmodel.getSy();
                obj.dischargeStore(i) = obj.coupledModel.getDischarge();
                obj.storageStore(i) = obj.coupledModel.getStorage();
                % update current time
                obj.t = ti;

                time_sim(i) = toc;

                if i > 5
                    if sum(time_sim((i-3):i)) > 240 % Seconds
                        error('Simulation found an error plato.')
                    end
                end


                if sum(time_sim) > 30*60
                    error('Total simulation time exceeded 30-min')
                end

                dx = obj.coupledModel.gwmodel.dx;
                dy = obj.coupledModel.gwmodel.dy;
                DY = obj.coupledModel.gwmodel.ny*obj.coupledModel.gwmodel.dy;
                DX = obj.coupledModel.gwmodel.nx*obj.coupledModel.gwmodel.dx;
                ny = obj.coupledModel.gwmodel.ny;
                nx = obj.coupledModel.gwmodel.nx;


                % if abs(obj.dischargeStore(i))*1000*3600 > 1.5
                %     error('Outflow too large')
                % end

            end
            %             catch ME
            %                 warning('Saving results up to this time')
            %             end
            GWdepth =  obj.pressureHeadStore;
            Qmod = obj.dischargeStore; % Saving Discharge
        end
        %% run
        function [] = run(obj)

            % for all time steps
            for i=1:length(obj.tSave)
                if i == 14*24*60/30
                    ttt = 1;
                end

                % start and end time of current time step
                time = obj.t;
                ti = obj.tSave(i);
                disp([obj.name ' - ' num2str(ti/(3600*24)) 'd/' num2str(obj.tSave(end)/(3600*24)) 'd'])

                % UZ Soil Moisture (initial)
                S_UZ_original = obj.coupledModel.gwmodel.parameter.S_UZ;

                % Original Recharge
                if i == 1
                    original_irrigation = obj.coupledModel.gwmodel.boundaries.Q; % Original values for irrigation
                end
                
                % Particle Filter Setup
                flag_PF = obj.coupledModel.gwmodel.flag_PF; % 1 runs the particle filter model
                PF_Setup();

                % advance coupled model
                if flag_PF == 1
                    % Correct UZ soil moisture
                    obj.coupledModel.gwmodel.parameter.S_UZ = S_UZ_original;   
                    try
                        obj.coupledModel.advance(time,ti); % Sometimes we have errors in the original model
                    end
                else
                    obj.coupledModel.advance(time,ti);
                end

                if flag_PF == 1
                    for j = 1:size(parameter_PF,2)
                        obj.parameterStore(i,:,j) = parameter_PF(:,j);
                    end
                    obj.dischargeStore_PF(i,:) = discharge_PF;
                    obj.GW_Depth_Store_PF(:,:,i) = GW_Depth_PF;
                    obj.Sy_Store_PF(:,i) = Sy_PF;
                end                

                % store model states
                obj.pressureHeadStore(:,i) = obj.coupledModel.getPressureHeads();
                obj.SyStore(:,i) = obj.coupledModel.gwmodel.getSy();
                obj.dischargeStore(i) = obj.coupledModel.getDischarge();
                if obj.dischargeStore(i)*1000*3600 > 10
                    ttt = 1;
                end
                obj.storageStore(i) = obj.coupledModel.getStorage();
                % update current time
                obj.t = ti;
            end
        end
        %% saveResults
        function saveResults(obj)
            runName = obj.name;
            runName(runName=='.') = [];
            % cells of groundwater model with soil column
            if obj.saveOption == 1
                XY = [obj.coupledModel.xi;obj.coupledModel.yi];
            else
                XY = NaN;
            end
            % save only groundwater levels
            GW = [];
            if obj.saveOption == 2 || obj.saveOption == 3
                GW = obj.pressureHeadStore(obj.coupledModel.nrStates- ...
                    obj.coupledModel.gwmodel.nrStates+1:end,:);
                GW = reshape(GW,obj.coupledModel.gwmodel.ny, ...
                    obj.coupledModel.gwmodel.nx,length(obj.tSave));
            end
            % save

            % Adding groundwate to the data
            GW = obj.pressureHeadStore(obj.coupledModel.nrStates- ...
                obj.coupledModel.gwmodel.nrStates+1:end,:);
            GW = reshape(GW,obj.coupledModel.gwmodel.ny, ...
                obj.coupledModel.gwmodel.nx,length(obj.tSave));

            Save.doSave(obj.saveOption,obj.tSave,runName,XY, ...
                obj.pressureHeadStore, ...
                obj.dischargeStore, GW,obj.SyStore,...
                obj.coupledModel.gwmodel.z0,obj.coupledModel.gwmodel.lsurf,...
                obj.coupledModel.gwmodel.dx,obj.coupledModel.gwmodel.dy,obj.coupledModel.gwmodel.nx,obj.coupledModel.gwmodel.ny,obj.coupledModel.gwmodel.parameter.H, ...
                obj.coupledModel.gwmodel.boundaries.Q,...
                obj.coupledModel.gwmodel.boundaries.changeTQ,obj.storageStore,obj.parameterStore,obj.dischargeStore_PF,obj.GW_Depth_Store_PF,obj.Sy_Store_PF);            
               
        end
    end

    methods (Static)
        % main (to be executed)
        function [] = main(setup)

            addpath Setups/
            addpath UZ/
            addpath GW/
            addpath Output/
            addpath Extra/
            addpath 'Particle Filter'\

            % initalize framework
            eval(['cmf = ' setup '();'])
            clc;

            % start
            cmf.run();

            % save
            cmf.saveResults();
        end
    end

end

