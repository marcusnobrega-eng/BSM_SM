classdef GWBoundaries < handle
    %Boundary conditions for the groundwater model
    
    properties
        mask            % assign active and non-active cells
        Q               % recharge (not the one calculated internally)
        changeTQ        % times at which recharge is changed
        perimI          % perimeter indices (rest no flow)
        perimIdir       % indices of Dirichlet boundaries
        perimIneu       % incides of Neumann boundaries
        doInternal      % internal dirichlet boundaries
        changeTDir      % times at which Dirichlet boundary is changed
        changeTNeu      % times at which Neumann boundary is changed
        boundDirValues  % values of the Dirichlet boundaries
        boundNeuValues  % values of the Neumann boundaries
        SeepageFace     % values equals 1 in seepage face
    end
    
    methods
        %% constructor
        function obj = GWBoundaries(mask,Q,changeTQ,perimI,perimIdir,perimIneu,...
                changeTDir,changeTNeu,boundDirValues,boundNeuValues)
            obj.mask = mask;
            obj.Q = Q;
            obj.changeTQ = changeTQ;
            obj.perimI = perimI;
            obj.perimIdir = perimIdir;
            obj.perimIneu = perimIneu;
            obj.doInternal = any(any(perimI(2:end-1,2:end-1)>0));
            obj.changeTDir = changeTDir;
            obj.changeTNeu = changeTNeu;
            obj.boundDirValues = boundDirValues;
            obj.boundNeuValues = boundNeuValues;
        end
        %% update boundary conditions
        function [q,dir,neu] = updateBoundaries(obj,t,iPerim,indexTot)
            dir = nan(length(indexTot),1);
            neu = zeros(length(indexTot),1);
            if ~isempty(obj.boundDirValues)
                tmp = obj.boundDirValues(:,obj.changeTDir<=t);
                for j=1:length(obj.perimIdir)
                    dir(iPerim==obj.perimIdir(j)) = tmp(j,end);
                end
            end
            if ~isempty(obj.boundNeuValues)
                tmp = obj.boundNeuValues(:,obj.changeTNeu<=t);
                for j=1:length(obj.perimIneu)
                    neu(iPerim==obj.perimIneu(j)) = tmp(j,end);
                end
            end
            tmp = obj.Q(:,:,obj.changeTQ<=t);
            tmp = tmp(:,:,end);
            q = tmp(indexTot);
        end
    end
    
end

