function correct_dem = extract_upstream_cells(dem, outlet_i, outlet_j,moore,X,Y)
% Initialize a matrix to mark visited cells

% Get the dimensions of the DEM
[rows, cols] = size(dem);
idx_ups = nan(rows,cols);
for i = 1:rows
    for j = 1:cols
        dem_cell = dem(outlet_i,outlet_j);
        % Left
        try
            if dem(i,j-1) > dem_cell
                idx_ups(i,j-1) = 1;
            end
        end
        % Right
        try
            if dem(i,j+1) > dem_cell
                idx_ups(i,j+1) = 1;
            end
        end
        % Up
        try
            if dem(i-1,j) > dem_cell
                idx_ups(i-1,j) = 1;
            end
        end
        % Down
        try
            if dem(i+1,j) > dem_cell
                idx_ups(i+1,j) = 1;
            end
        end
        if moore == 1
            % NE
            try
                if dem(i-1,j+1) > dem_cell
                    idx_ups(i-1,j+1) = 1;
                end
            end
            % NW
            try
                if dem(i-1,j-1) > dem_cell
                    idx_ups(i-1,j-1) = 1;
                end
            end
            % SE
            try
                if dem(i+1,j+1) > dem_cell
                    idx_ups(i+1,j+1) = 1;
                end
            end
            % SW
            try
                if dem(i+1,j-1) > dem_cell
                    idx_ups(i+1,j-1) = 1;
                end
            end
        end
    end
end
correct_dem = dem.*idx_ups;
surf_plot(max(max(correct_dem)),0,'z','m',correct_dem,1,'DEM',256,0.9,1,[0,90],X,Y)

end

