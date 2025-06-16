classdef Save
    % Save contains all save options
    
    properties
    end
    
    methods (Static)
        function [] = doSave(saveOption,tSave,runName,XY,pressureHead, ...
                discharge,GW,Sy,z0,lsurf,dx,dy,nx,ny,H,TopBC,ChangeT,storage,parameterStore,dischargeStore_PF,GW_Depth_Store_PF,Sy_Store_PF)
            switch saveOption
                case 1
                    % save all data
                    save(['Output/' runName '.mat'], 'tSave','XY','pressureHead', ...
                        'waterContent')
                case 2
                    % save only groundwater data
                    save(['Output/' runName '.mat'], 'GW','Sy')
                case 3
                    % save all data
                    save(['Output/' runName '.mat'], 'tSave','XY','pressureHead', ...
                        'waterContent','discharge','GW','Sy','z0','lsurf','dx','dy','nx','ny','parameterStore','dischargeStore_PF','GW_Depth_Store_PF','Sy_Store_PF')
                case 4
                    % save more data
                    save(['Output/' runName '.mat'], 'tSave','XY','pressureHead', ...
                        'discharge','GW','Sy','z0','lsurf','dx','dy','H', ...
                        'TopBC','ChangeT','storage','parameterStore','dischargeStore_PF','GW_Depth_Store_PF','Sy_Store_PF')                    
            end
        end
    end
    
end

