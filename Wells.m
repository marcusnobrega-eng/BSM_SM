classdef Wells < handle
    %Wells in the groundwater model
    
    properties
        doWells     % wells in the domain?
        nrWells     % number of wells
        x           % x-positions of wells
        y           % y-positions of wells
        tStart      % start time of wells
        tStop       % stop time of wells
        Q           % pumpring rate (negative is extraction)     
    end
    
    methods
        %% constructor
        function obj = Wells(doWells,nrWells,x,y,tStart,tStop,Q)
            obj.doWells = doWells;
            obj.nrWells = nrWells;
            obj.x = x;
            obj.y = y;
            obj.tStart = tStart;
            obj.tStop = tStop;
            obj.Q = Q;
        end
        %% update wells
        function W = updateWells(obj,t,ny,nx,wXY)
            % if not in the start-file, use old approach
            W=zeros(ny,nx);
            % Check if any wells are active
            activeW=t-obj.tStart >= 0 & t-obj.tStop < 0;
            if any(activeW)
                for iw=1:obj.nrWells
                    if activeW(iw)
                        W(wXY==iw)=obj.Q(iw);
                    end
                end
            end
        end
    end
    
end

