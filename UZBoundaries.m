classdef UZBoundaries < handle
    %Boundary conditions for the unsaturated zone model
    
    properties
        % top boundary
        fluxT                % boundary condition at the top (Neumann)
        currFluxT            % current boundary condition at the top (Neumann)
        changeT              % times at which boundary condition at the top is changed
        % bottom boundary
        headB                % boundary condition at the bottom (Dirichlet)
        noFluxB              % no flow boundary at the bottom
        currHeadB            % current boundary condition at the bottom (Dirichlet)
        changeB              % times at which boundary condition at the bottom is changed
        % transpiration
        trans                % transpiration time series
        currTrans            % current transpiration rate
        changeTr             % times at which transpiration in changed
        % general settings
        ETBreak = -1000;     % cut outflow
        doFreeDrain = false; % free drainage at the bottom
        freeDrainOption = 1; % 1: gravity drain; 2: equal flux; 3: on/off boundary condition
        noBackFlow = false;  % should water be stopped from flowing back into the system?
    end
    
    methods
        %% constructor
        function obj = UZBoundaries(fluxT,changeT,headB,noFluxB,changeB,trans,changeTr)
            obj.fluxT = fluxT;
            obj.currFluxT = fluxT;
            obj.changeT = changeT;
            obj.headB = headB;
            obj.noFluxB = noFluxB;
            obj.currHeadB = headB;
            obj.changeB = changeB;
            obj.trans = trans;
            obj.currTrans = trans;
            obj.changeTr = changeTr;
        end
        %% update boundary conditions
        function [obj] = updateBoundaries(obj,t)
            if ~obj.noFluxB
                tmp = obj.headB(obj.changeB<=t);
                obj.currHeadB = tmp(end);
            end
            tmp = obj.fluxT(obj.changeT<=t);
            obj.currFluxT = tmp(end);
            tmp = obj.trans(obj.changeTr<=t);
            obj.currTrans = tmp(end);
        end
    end
    
end

