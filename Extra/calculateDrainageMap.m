function drainageMap = calculateDrainageMap(DEM, i, j)
    [rows, cols] = size(DEM);
    drainageMap = false(rows, cols); % Initialize drainage map
    
    % Calculate flow direction
    flowDir = calculateFlowDirection(DEM);

    % Define a visited matrix to keep track of visited cells
    visited = false(rows, cols);

    % Recursively mark cells that drain to the outlet
    markDrainingCells(i, j);

    % Define recursive function to mark draining cells
    function markDrainingCells(x, y)
        % Mark current cell as visited
        visited(y, x) = true;

        % Mark current cell as draining
        drainageMap(y, x) = true;

        % Get downstream neighbor
        [nx, ny] = getDownstreamCell(x, y, flowDir);

        % Check if neighbor is within bounds and not already visited
        if nx >= 1 && nx <= cols && ny >= 1 && ny <= rows && ~visited(ny, nx)
            % Recursively mark downstream cells
            markDrainingCells(nx, ny);
        end
    end

    % Inner function to calculate downstream cell
    function [nx, ny] = getDownstreamCell(x, y, flowDir)
        dir = flowDir(y, x);
        switch dir
            case 1 % Northwest
                nx = x - 1;
                ny = y - 1;
            case 2 % North
                nx = x;
                ny = y - 1;
            case 3 % Northeast
                nx = x + 1;
                ny = y - 1;
            case 4 % East
                nx = x + 1;
                ny = y;
            case 5 % Southeast
                nx = x + 1;
                ny = y + 1;
            case 6 % South
                nx = x;
                ny = y + 1;
            case 7 % Southwest
                nx = x - 1;
                ny = y + 1;
            case 8 % West
                nx = x - 1;
                ny = y;
            otherwise
                nx = x;
                ny = y;
        end
    end
end

function flowDir = calculateFlowDirection(DEM)
    % D8 algorithm to calculate flow direction
    [rows, cols] = size(DEM);
    flowDir = zeros(rows, cols);

    for y = 1:rows
        for x = 1:cols
            [dx, dy] = getNeighborDeltas(x, y, rows, cols);
            [row, col] = getDownstreamCell(x, y, DEM, dx, dy);
            flowDir(y, x) = getDirection(row, col, x, y);
        end
    end
end

function [dx, dy] = getNeighborDeltas(x, y, rows, cols)
    % Get indices of neighboring cells
    dx = [x-1, x-1, x-1, x, x+1, x+1, x+1, x];
    dy = [y-1, y, y+1, y+1, y+1, y, y-1, y-1];
    % Ensure indices are within bounds
    dx = max(1, min(dx, cols));
    dy = max(1, min(dy, rows));
end

function [row, col] = getDownstreamCell(x, y, DEM, dx, dy)
    % Find the downstream cell based on the steepest slope
    [slopes, indices] = sort((DEM(y, x) - DEM(dy, dx)) ./ sqrt((x - dx).^2 + (y - dy).^2), 'descend');
    row = dx(indices(1));
    col = dy(indices(1));
end

function dir = getDirection(row, col, x, y)
    % Calculate flow direction based on downstream cell
    if row == x - 1 && col == y - 1
        dir = 1; % Northwest
    elseif row == x && col == y - 1
        dir = 2; % North
    elseif row == x + 1 && col == y - 1
        dir = 3; % Northeast
    elseif row == x + 1 && col == y
        dir = 4; % East
    elseif row == x + 1 && col == y + 1
        dir = 5; % Southeast
    elseif row == x && col == y + 1
        dir = 6; % South
    elseif row == x - 1 && col == y + 1
        dir = 7; % Southwest
    elseif row == x - 1 && col == y
        dir = 8; % West
    else
        dir = 0; % No flow
    end
end
