function [h_t, u_x, u_y, q_exf, q_river] = Boussinesq_2D_explicit(dt,dx,dy,h,z0,S_s,R,K,river_mask,K_river,h_river,z_river,Courant,h_soil,catchment_mask,dirichlet_mask,h_dirichlet,free_flow_mask)
% Solves the 2D Boussinesq equation using an explicit finite difference method
% with adaptive time-stepping based on the Courant stability condition.
%
% Inputs:
% - dt: Initial model time-step [sec]
% - dx: Grid spacing in x-direction [m]
% - dy: Grid spacing in y-direction [m]
% - h: Initial hydraulic head [m] (Nx, Ny)
% - z0: Ground elevation [m] (Nx, Ny)
% - S_s: Specific storage coefficient [1/m]
% - R: Recharge rate [m/s] (Nx, Ny)
% - K: Hydraulic conductivity [m/s] (Nx, Ny)
% - river_mask: logical mask indicating where rivers exist
% - K_river: hydraulic conductivity of the river [m/s]
% - h_river: hydraulic head at the rivers [matrix in DEM size]
% - z_river: topographic river elevtion [mm]
% - Courant: Stability factor (0 < Courant <= 1)
% - h_soil: Soil depth [m] (Nx, Ny)
% - catchment_mask: Binary mask indicating computational domain (Nx, Ny)
% - dirichlet_mask: Binary mask for Dirichlet boundary conditions (Nx, Ny)
% - h_dirichlet: Prescribed hydraulic head at Dirichlet boundaries [m] (Nx, Ny)
% - free_flow_mask: Binary mask for free-flow boundaries (Nx, Ny)
%
% Outputs:
% - h_t: Final hydraulic head [m] (Nx, Ny)
% - u_x: Groundwater velocity in x-direction [m/s] (Nx-1, Ny)
% - u_y: Groundwater velocity in y-direction [m/s] (Nx, Ny-1)
% - q_exf: Exfiltration rate at cells [m/s] (Nx, Ny)
% - q_river: River exfiltrtion [m/s]. Positive leaves the aquifer. Negative,
% enters the aquifer
%

%% Surface Water Level
h_surf = z0 + h_soil;  % Compute the surface elevation from topography

%% INITIAL CONDITION
q_exf = zeros(size(h)); % Initialize exfiltration rate matrix

%% River Flux Length 
L_river = 1/2*(dx+dy);

%% TIME STEPPING LOOP
t = 0; % Initialize time
T = dt; % Set total simulation time to initial dt

while t < T
    % Compute new stable dt at each step
    dt_GW = compute_stable_dt(K, S_s, dx, dy, Courant); % Adaptive time step
    
    % Ensure dt does not exceed remaining simulation time
    dt = min(dt_GW,dt);

    h_t = h; % Copy current head for updating

    % Compute groundwater discharge to rivers
    q_river = zeros(Nx, Ny);  % Initialize discharge array
    river_cells = (river_mask == 1);  % Identify river locations
    q_river(river_cells) = K_river .* ((h(river_cells) - h_river(river_cells)) ./ L_river); % m/s

    % River water depth constraint
    river_depth = h_river(river_cells) - z_river(river_cells);
    q_river(q_river<0 & river_depth <= 0) = 0; % No influx in the aquifer from zero-depth rivers
    
    if max(abs(q_river)) > 0 % At least one river cell has fluxes
        % Updated Groundwater Head
        h(river_cells) = h(river_cells) - dt*q_river(river_cells);
    end
    
    % Apply Dirichlet boundary conditions
    h_t(dirichlet_mask) = h_dirichlet(dirichlet_mask); 
    
    % Apply free-flow boundary conditions
    h_t = apply_free_flow_bc(h_t, free_flow_mask);

    % Compute second derivatives using central differencing
    [d2h_dx2, d2h_dy2] = compute_gradients(h, K, dx, dy);

    % Update hydraulic head using explicit finite difference
    h_t(2:end-1,:) = h(2:end-1,:) + dt/S_s .* d2h_dx2; % x contribution
    h_t(:,2:end-1) = h_t(:,2:end-1) + dt/S_s .* d2h_dy2; % y contribution
    h_t = h_t + dt/S_s .* R; % Recharge contribution

    % Ensure groundwater head does not fall below ground level
    h_t = max(h_t, z0);

    % Compute exfiltration where groundwater head exceeds surface elevation
    q_exf = max(0, (h - h_surf) / dt);

    % Groundwater Head Considering Exfiltration
    h_t = min(h_t, h_surf);

    % Compute groundwater velocities
    [u_x, u_y] = compute_boussinesq_velocities(h_t, K, dx, dy);

    % Apply catchment mask (ensuring calculations remain within domain)
    h_t(~catchment_mask) = 0;
    u_x(~catchment_mask) = 0;
    u_y(~catchment_mask) = 0;
    q_exf(~catchment_mask) = 0;

    % Visualization (Remove for efficiency in large simulations)
    surf(h_t);
    title(sprintf('2D Groundwater Flow (Time = %.2f s)', t));
    view(0,90);
    colormap('turbo');
    colorbar;
    pause(0.00001);

    % Update variables for next time step
    h = h_t;
    t = t + dt;
end


end

%% LOCAL FUNCTIONS

% Compute harmonic mean of two values
function h_mean = mean_function(K1, K2,flag_mean)
if K1 + K2 == 0
    h_mean = 0; % Avoid division by zero
else
    if flag_mean == 1
        h_mean = 2 * (K1 .* K2) ./ (K1 + K2);
    else
        h_mean = 1/2*(K1 + K2);
    end
end
end


% Compute the second derivatives (gradient) using matrix operations
function [d2h_dx2, d2h_dy2] = compute_gradients(h, K, dx, dy)
% Use finite difference with convolution for second derivatives

% Compute harmonic mean conductivities at interfaces
Kx_p = mean_function(K(2:end-1, :), K(1:end-2, :),1); % Right side conductivity
Kx_m = mean_function(K(1:end-2, :), K(2:end-1, :),1); % Left side conductivity
hx_p = mean_function(h(2:end-1, :), h(1:end-2, :),2); % Right side conductivity
hx_m = mean_function(h(1:end-2, :), h(2:end-1, :),2); % Left side conductivity

% Central difference for y-direction
Ky_p = mean_function(K(:, 2:end-1), K(:, 1:end-2),1); % Upward conductivity
Ky_m = mean_function(K(:, 1:end-2), K(:, 2:end-1),1); % Downward conductivity
hy_p = mean_function(h(:, 2:end-1), h(:, 1:end-2),2); % Upward conductivity
hy_m = mean_function(h(:, 1:end-2), h(:, 2:end-1),2); % Downward conductivity

% Compute second derivatives using correct finite differences
d2h_dx2 = 1/dx*(Kx_p .* hx_p.*(h(3:end,:) - h(2:end-1,:))/dx - Kx_m .* hx_m .* (h(2:end-1,:) - h(1:end-2,:))/dx);
d2h_dy2 = 1/dy*(Ky_p .* hy_p.*(h(:, 3:end) - h(:,2:end-1))/dy - Ky_m .* hy_m .* (h(:,2:end-1) - h(:,1:end-2))/dy);

end


% Compute the maximum stable time step based on the Courant condition
function dt_max = compute_stable_dt(K, S_s, dx, dy, Courant)
K_min = min(K(:)); % Find the smallest K in the domain
dt_max = (S_s / K_min) * ((dx^2 * dy^2) / (2 * (dx^2 + dy^2)));
dt_max = Courant * dt_max; % Apply a safety factor
end

function h = apply_free_flow_bc(h, mask)
[Nx, Ny] = size(h);

for i = 2:Nx-1
    for j = 2:Ny-1
        if ~mask(i, j) % If it's outside the catchment
            % free_flow BC: Set to nearest interior neighbor
            if mask(i-1, j), h(i, j) = h(i-1, j); end
            if mask(i+1, j), h(i, j) = h(i+1, j); end
            if mask(i, j-1), h(i, j) = h(i, j-1); end
            if mask(i, j+1), h(i, j) = h(i, j+1); end
        end
    end
end
end



function [Sx, Sy] = compute_slopes_from_dem(DEM, dx, dy)
% COMPUTE_SLOPES_FROM_DEM - Computes the slopes in the x and y directions from a DEM.
%
% Inputs:
% - DEM: Digital Elevation Model matrix [m] (nx, ny)
% - dx: Grid spacing in the x-direction [m]
% - dy: Grid spacing in the y-direction [m]
%
% Outputs:
% - Sx: Slope in the x-direction (nx, ny)
% - Sy: Slope in the y-direction (nx, ny)

% Get size of DEM
[nx, ny] = size(DEM);

% Initialize slope matrices
Sx = zeros(nx, ny);
Sy = zeros(nx, ny);

% Compute slope in x-direction (central difference)
% Sx(2:nx-1, :) = (DEM(3:nx, :) - DEM(1:nx-2, :)) / (2 * dx);
Sx(2:nx-1, :) = (DEM(2:nx-1, :) - DEM(1:nx-2, :)) / (dx);

% Compute slope in y-direction (central difference)
% Sy(:, 2:ny-1) = (DEM(:, 3:ny) - DEM(:, 1:ny-2)) / (2 * dy);
Sy(:, 2:ny-1) = (DEM(:, 2:ny-1) - DEM(:, 1:ny-2)) / (dy);

% Boundary conditions (one-sided difference at the boundaries)
Sx(1, :) = (DEM(2, :) - DEM(1, :)) / dx;  % Forward difference for left boundary
Sx(nx, :) = (DEM(nx, :) - DEM(nx-1, :)) / dx;  % Backward difference for right boundary
Sy(:, 1) = (DEM(:, 2) - DEM(:, 1)) / dy;  % Forward difference for bottom boundary
Sy(:, ny) = (DEM(:, ny) - DEM(:, ny-1)) / dy;  % Backward difference for top boundary

end

function [vx, vy] = compute_boussinesq_velocities(h, K, dx, dy)
% COMPUTE_BOUSSINESQ_VELOCITIES - Computes groundwater velocities from the Boussinesq equation
% using Darcy's law.
%
% Inputs:
% - h: Hydraulic head matrix [m] (nx, ny)
% - K: Hydraulic conductivity matrix [m/s] (nx, ny)
% - dx: Grid spacing in x-direction [m]
% - dy: Grid spacing in y-direction [m]
%
% Outputs:
% - vx: Groundwater velocity in x-direction [m/s] (nx-1, ny)
% - vy: Groundwater velocity in y-direction [m/s] (nx, ny-1)

[nx, ny] = size(h);

% Compute harmonic mean of hydraulic conductivity at cell interfaces
Kx = 2 ./ (1./K(1:nx-1, :) + 1./K(2:nx, :)); % Interface values in x-direction
Ky = 2 ./ (1./K(:, 1:ny-1) + 1./K(:, 2:ny)); % Interface values in y-direction

% Compute hydraulic head gradients using central differences
dhdx = (h(2:nx, :) - h(1:nx-1, :)) / dx;  % Gradient in x-direction
dhdy = (h(:, 2:ny) - h(:, 1:ny-1)) / dy;  % Gradient in y-direction

% Compute velocities using Darcy's law
vx = -Kx .* dhdx; % x-component of groundwater velocity
vy = -Ky .* dhdy; % y-component of groundwater velocity

end