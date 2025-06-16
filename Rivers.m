classdef Rivers < handle
    %Rivers in the groundwater model
    
    properties
        doRivers            % rivers in the domain?
        riverMask           % position of the rivers
        riverMaskDomain     % position of the rivers (logical)
        dx                  % discretization in x-direction
        dy                  % discretization in y-direction
        changeT             % changing times of the river water level
        K                   % hydraulic conductivity of the river bed
        d                   % thickness of the river bed
        h                   % water level
    end
    
    methods
        %% constructor
        function obj = Rivers(doRivers,riverMask,dx,dy,changeT,K,d,h)
            obj.doRivers = doRivers;
            obj.riverMask = riverMask;
            obj.dx = dx;
            obj.dy = dy;
            obj.changeT = changeT;
            obj.K = K;
            obj.d = d;
            obj.h = h;
        end
        %% update rivers
        function hR = updateRivers(obj,t)
            tmp = obj.h(obj.changeT<=t);
            hR = tmp(end);
        end
    end
    
end

