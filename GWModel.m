classdef GWModel < handle
    %GWModel calculates the depth averaged 2d groundwater equation
    
    properties
        % geometry
        z0              % bottom elevation
        lsurf           % land surface
        % discretization
        nx              % number of cells in x-direction
        ny              % number of cells in y-direction
        dx              % size of cells in x-direction
        dy              % size of cells in y-direction
        nrStates        % number of model states
        % parameter
        parameter       % soil parameter
        zones           % zones
        % boundary conditions
        boundaries      % contains all boundary conditions and timing
        % wells
        wells           % wells in the domain
        % rivers
        rivers          % rivers in the domain
        % states
        h               % states (pressure heads)
        % initial condition
        IC              % initial condition
        % solver
        solver          % numerical solver
        % timing
        dt              % time step size
        doAdaptiveTime  % adaptive time step size  
        flag_PF         % Activate particle filter
        model_dir      % Excel spreadsheet directory
    end
    
    methods
        %% constructor
        function obj = GWModel(z0,lsurf,nx,ny,dx,dy,parameter,boundaries, ...
                wells,rivers,IC,solver,dt,doAdaptiveTime,flag_PF,dir)
            obj.z0 = z0;
            obj.lsurf = lsurf;
            obj.nx = nx;
            obj.ny = ny;
            obj.dx = dx;
            obj.dy = dy;
            obj.parameter = parameter;
            obj.boundaries = boundaries;
            obj.wells = wells;
            obj.rivers = rivers;
            obj.IC = IC;
            obj.h = IC;
            obj.nrStates = nx*ny;
            obj.solver = solver;
            obj.dt = dt;
            obj.doAdaptiveTime = doAdaptiveTime;
            obj.flag_PF = flag_PF;
            obj.model_dir = dir;
        end
        %% advance
        function h_new = advance(obj,tStart,tStop,R)
            % run model forward
            h_new = obj.solver.solveGW(tStart,tStop,obj.dt,obj,R);
        end
        %% getPressureHeads
        function P = getPressureHeads(obj)
            P = obj.h(:)-obj.z0(:);
        end
        %% getSy
        function Sy = getSy(obj)
            Sy = obj.parameter.Sy(:);
        end
    end
    
end

