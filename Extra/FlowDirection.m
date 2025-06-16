function [f_dir,idx_fdir] = FlowDirection(DEM,Delta_x,Delta_y,coord_outlet)
% Developer: Marcus Nobrega Gomes Junior and Marcio H. Giacomoni
% 11/02/2021
% Input Data
% This function develops the information regarding the V-Tilted Catchment
% Goal - Assign parameters for cells according to the land use and land
% cover given by the imperviousness map
% Last Update - 11/2/2021

dimensions = size(DEM);
row_out = coord_outlet(1); col_out = coord_outlet(2);
n_Col = dimensions(2);
n_Lin = dimensions(1);
% Grid Dimension
Delta_xY = sqrt((Delta_x) ^ 2 + (Delta_y) ^ 2);
% Preallocating Array
f_dir = zeros(n_Lin,n_Col);
% Initializing Variables
a = 0 ; b = 0 ; c = 0 ; d = 0 ; e = 0 ; f = 0 ; g = 0 ; h = 0 ; i = 0;
Slope1 = 0 ; Slope2 = 0 ; Slope3 = 0 ; Slope4 = 0 ; Slope5 = 0 ; Slope6 = 0 ; Slope7 = 0 ; Slope8 = 0;
Steepest_Slope = 0;

for lin = 1 : n_Lin
    for col = 1 : n_Col
        e = DEM(lin, col);
        Steepest_Slope = 1000000000;
        if ~isnan(e)
            if ((lin == 1) && (col == 1))
                a = e;
                b = e;
                c = e;
                d = e;
                f = DEM(lin, col + 1);
                g = e;
                h = DEM(lin + 1, col);
                i = DEM(lin + 1, col + 1);
            elseif ((lin == 1) && (col > 1) && (col < n_Col))
                a = e;
                b = e;
                c = e;
                d = DEM(lin, col - 1);
                f = DEM(lin, col + 1);
                g = DEM(lin + 1, col - 1);
                h = DEM(lin + 1, col);
                i = DEM(lin + 1, col + 1);
            elseif ((lin == 1) && (col == n_Col))
                a = e;
                b = e;
                c = e;
                d = DEM(lin, col - 1);
                f = e;
                g = DEM(lin + 1, col - 1);
                h = DEM(lin + 1, col);
                i = e;
            elseif ((lin > 1) && (lin < n_Lin) && (col == 1))
                a = e;
                b = DEM(lin - 1, col);
                c = DEM(lin - 1, col + 1);
                d = e;
                f = DEM(lin, col + 1);
                g = e;
                h = DEM(lin + 1, col);
                i = DEM(lin + 1, col + 1);
            elseif ((lin > 1) && (lin < n_Lin) && (col == n_Col))
                a = DEM(lin - 1, col - 1);
                b = DEM(lin - 1, col);
                c = e;
                d = DEM(lin, col - 1);
                f = e;
                g = DEM(lin + 1, col - 1);
                h = DEM(lin + 1, col);
                i = e;
            elseif ((lin == n_Lin) && (col == 1))
                a = e;
                b = DEM(lin - 1, col);
                c = DEM(lin - 1, col + 1);
                d = e;
                f = DEM(lin, col + 1);
                g = e;
                h = e;
                i = e;
            elseif ((col > 1) && (col < n_Col) && (lin == n_Lin))
                a = DEM(lin - 1, col - 1);
                b = DEM(lin - 1, col);
                c = DEM(lin - 1, col + 1);
                d = DEM(lin, col - 1);
                f = DEM(lin, col + 1);
                g = e;
                h = e;
                i = e;
            elseif ((col == n_Col) && (lin == n_Lin))
                a = DEM(lin - 1, col - 1);
                b = DEM(lin - 1, col);
                c = e;
                d = DEM(lin, col - 1);
                f = e;
                g = e;
                h = e;
                i = e;
            else
                a = DEM(lin - 1, col - 1);
                b = DEM(lin - 1, col);
                c = DEM(lin - 1, col + 1);
                d = DEM(lin, col - 1);
                f = DEM(lin, col + 1);
                g = DEM(lin + 1, col - 1);
                h = DEM(lin + 1, col);
                i = DEM(lin + 1, col + 1);
            end

            Slope1 = (f - e) / Delta_x;
            Slope2 = (i - e) / Delta_xY;
            Slope3 = (h - e) / Delta_y;
            Slope4 = (g - e) / Delta_xY;
            Slope5 = (d - e) / Delta_x;
            Slope6 = (a - e) / Delta_xY;
            Slope7 = (b - e) / Delta_y;
            Slope8 = (c - e) / Delta_xY;

            if (Slope1 >= 0 || isnan(Slope1)) && (Slope2 >= 0 || isnan(Slope2)) && (Slope3 >= 0 || isnan(Slope3)) && (Slope4 >= 0 || isnan(Slope4)) && (Slope5 >= 0 || isnan(Slope5)) && (Slope6 >= 0 || isnan(Slope6)) && (Slope7 >= 0 || isnan(Slope7)) && (Slope8 >= 0 || isnan(Slope8)) 
                f_dir(lin,col) = 0;
            else
                if (Slope1 < Steepest_Slope)
                    f_dir(lin, col) = 1;
                    Steepest_Slope = Slope1;
                end
                if (Slope2 < Steepest_Slope)
                    f_dir(lin, col) = 2;
                    Steepest_Slope = Slope2;
                end
                if (Slope3 < Steepest_Slope)
                    f_dir(lin, col) = 4;
                    Steepest_Slope = Slope3;
                end
                if (Slope4 < Steepest_Slope)
                    f_dir(lin, col) = 8;
                    Steepest_Slope = Slope4;
                end
                if (Slope5 < Steepest_Slope)
                    f_dir(lin, col) = 16;
                    Steepest_Slope = Slope5;
                end
                if (Slope6 < Steepest_Slope)
                    f_dir(lin, col) = 32;
                    Steepest_Slope = Slope6;
                end
                if (Slope7 < Steepest_Slope)
                    f_dir(lin, col) = 64;
                    Steepest_Slope = Slope7;
                end
                if (Slope8 < Steepest_Slope)
                    f_dir(lin, col) = 128;
                    Steepest_Slope = Slope8;
                    if isnan(DEM(lin-1,col+1))
                        error('Please make sure your DEM is hydrologic corrected or if the watershed resolution is properly defined.')
                    end
                    % Acrescentar mensagem de erro quando Bh não está bem definida.
                end
            end
        else
            f_dir(lin,col) = nan;
        end
    end
    f_dir(row_out,col_out) = 0; % Outlet consideration
    idx_fdir = f_dir>0; % Only cells with flow direction
    idx_fdir(row_out,col_out) = 1; % Considering the outlet in the mass balance
    idx_fdir = idx_fdir(:);
end