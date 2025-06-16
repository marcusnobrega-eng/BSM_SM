function groundwater_flow_explicit()
% Solves the 2D Boussinesq equation using an explicit finite difference method
% with adaptive time-stepping based on the Courant stability condition.

%% PARAMETERS
Nx = 50;         % Number of grid points in x-direction
Ny = 50;         % Number of grid points in y-direction
Lx = 100;        % Domain length in x (meters)
Ly = 100;        % Domain length in y (meters)
dx = Lx / Nx;    % Grid spacing in x
dy = Ly / Ny;    % Grid spacing in y
S_s = 0.1;       % Specific yield coefficient
P = 1e-4*ones(Nx, Ny);        % Rainfall rate (m/s)
k = 1e-5*ones(Nx, Ny);        % Reservoir Constant
I = 1e-7 * ones(Nx, Ny); % Infiltration rate [m/s]
manning_n = 0.03*ones(Nx, Ny); % Manning’s roughness coefficient

% Soil K
K = ones(Nx, Ny)*1/3600; % m/s

% Total simulation time
T = 60000;          % Simulation time (seconds)
t = 0;           % Start time
Courant = 0.25;

%% Surface
z0 = zeros(Nx,Ny); % Bedrock elevation [m]
z0(:, Ny/3:2*Ny/3) = -0.25;     % Initial perturbation at the center
soil_depth = ones(Nx,Ny); % Soil depth [m]
h_surf = z0 + soil_depth;

%% INITIAL CONDITION
h = zeros(Nx, Ny);  % Initializing the hydraulic head field
h(:, Ny/3:2*Ny/3) = 0.75;     % Initial perturbation at the center
h_s = zeros(Nx, Ny);  % Initial surface water depth [m]
exfiltration_rate = zeros(Nx, Ny); % Exfiltration rate
initial_soil_moisture = 10/1000*ones(Nx,Ny); % Initil soil moisture [mm]
u = zeros(Nx, Ny); % Velocity in x [m/s]
v = zeros(Nx, Ny); % Velocity in y [m/s]

%% Slopes
[Sx, Sy] = compute_slopes_from_dem(h_surf + h_s, dx, dy);

%% BOUNDARY CONDITIONS
% Define Dirichlet boundary condition mask
dirichlet_mask = false(Nx, Ny);
dirichlet_mask(1, :) = true;  % Fixed head on the left
dirichlet_mask(Nx, :) = true; % Fixed head on the right

% Define fixed head values at Dirichlet boundaries
% mask = create_catchment_mask(Nx, Ny);
mask = zeros(Nx,Ny); mask(:,1) = 1; mask(:,end) = 1; mask(1,:) = 1; mask(end,:) = 1;
h_dirichlet = zeros(Nx, Ny);
h_dirichlet(1, :) = 0;   % Fixed head on the left boundary
h_dirichlet(Nx, :) = 0;  % Fixed head on the right boundary


%% TIME STEPPING LOOP
while t < T

    % Compute new stable dt at each step
    dt_GW = compute_stable_dt(K, S_s, dx, dy, Courant);

    [Sx, Sy] = compute_slopes_from_dem(h_surf + h_s, dx, dy);
    dt_OLF = compute_overland_timestep(h_s, u, v, Sx, Sy, dx, dy, manning_n, 9.81, Courant);

    % New timestep
    dt = min(min(dt_GW,dt_OLF));

    % Surface Water Balance
    I = 0.5*P;

    % Compute overland flow routing
    [h_s_new, u, v] = solve_overland_flow(h_s, z0,  u, v, P, I, manning_n, Sx, Sy, dx, dy, dt, exfiltration_rate);

    % Compute recharge
    [R, updated_soil_moisture] = simulate_groundwater_recharge(I, initial_soil_moisture, k, dt);

    h_new = h; % Copy current head for updating

    % Apply boundary conditions using the mask
    h_new(dirichlet_mask) = h_dirichlet(dirichlet_mask); % Dirichlet BC
    h_new = apply_neumann_bc(h_new, mask);               % Neumann BC

    % Compute gradients using vectorized operations
    [d2h_dx2, d2h_dy2] = compute_gradients(h, K, dx, dy);


    h_new(2:end-1,:) = h(2:end-1,:) + dt/S_s*d2h_dx2;     % x contribution
    h_new(:,2:end-1) = h_new(:,2:end-1) + dt/S_s*d2h_dy2; % y contribution
    h_new = h_new + dt/S_s.*R;                            % Recharge Contribution

    h_new = max(h_new, z0); % Avoiding negative depths

    % Exfiltration occurs if groundwater head exceeds surface elevation
    exfiltration_rate = max(0, (h - h_surf) / dt);        % m/s

    % Groundwater Head Considering Exfiltration
    h_new = min(h_new,h_surf);

    subplot(2,1,1)
    surf(h);
    title(sprintf('2D Groundwater Flow (Time = %.2f s)', t)); % Display time in title
    view(0,90); colormap('turbo'); colorbar;
    subplot(2,1,2)
    surf(h_s_new);
    title(sprintf('2D Overland Flow (Time = %.2f s)', t)); % Display time in title
    view(0,90); colormap('turbo'); colorbar;
    pause(0.00001)
    h = h_new; % Update the field for the next step
    h_s = h_s_new;
    initial_soil_moisture = updated_soil_moisture;
    t = t + dt; % Advance time
end


end

%% LOCAL FUNCTIONS

% Compute second derivative in x-direction
function d2h_dx2 = second_derivative_x(h, i, j, dx)
d2h_dx2 = (h(i+1, j) - 2*h(i, j) + h(i-1, j)) / dx^2;
end

% Compute second derivative in y-direction
function d2h_dy2 = second_derivative_y(h, i, j, dy)
d2h_dy2 = (h(i, j+1) - 2*h(i, j) + h(i, j-1)) / dy^2;
end

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

% Create a catchment mask (1 = active, 0 = inactive)
function mask = create_catchment_mask(Nx, Ny)
mask = true(Nx, Ny);
for i = 1:Nx
    for j = 1:Ny
        if (i < Nx/4 || i > 3*Nx/4) && (j < Ny/4 || j > 3*Ny/4)
            mask(i, j) = false; % Remove corners to make an irregular shape
        end
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

function h = apply_neumann_bc(h, mask)
[Nx, Ny] = size(h);

for i = 2:Nx-1
    for j = 2:Ny-1
        if ~mask(i, j) % If it's outside the catchment
            % Neumann BC: Set to nearest interior neighbor
            if mask(i-1, j), h(i, j) = h(i-1, j); end
            if mask(i+1, j), h(i, j) = h(i+1, j); end
            if mask(i, j-1), h(i, j) = h(i, j-1); end
            if mask(i, j+1), h(i, j) = h(i, j+1); end
        end
    end
end
end

function [h_s_new, u, v] = solve_overland_flow(h_s, z0, u, v,rainfall, infiltration, manning_n, Sx, Sy, dx, dy, dt, exfiltration_rate)
% SOLVE_OVERLAND_FLOW Computes the new water depths and velocities using the diffusive wave equation
%
% Inputs:
% - h_s: Surface water depth matrix [m]
% - v: Velocity matrix [m/s]
% - rainfall: Rainfall rate matrix [m/s]
% - infiltration: Infiltration rate matrix [m/s]
% - manning_n: Manning's roughness coefficient
% - slope: Surface slope matrix (assumed constant or spatially varying)
% - dx, dy: Grid spacing in x and y directions [m]
% - dt: Time step [s]
%
% Outputs:
% - h_s_new: Updated water depth [m]
% - v_new: Updated velocity [m/s]

nx = size(h_s,1);
ny = size(h_s,2);

tol = 1e-6;

h_s = max(h_s,tol);

h_s_x = max(h_s(2:end-1, :) + z0(2:end-1, :), h_s(1:end-2, :) + z0(1:end-2, :)) - max(z0(2:end-1, :),z0(1:end-2, :)); % Right side conductivity

h_s_x = [h_s(1,:); h_s_x;  h_s(end,:)];

h_s_y = max(h_s(:, 2:end-1) + z0(:, 2:end-1), h_s(:,1:end-2) + z0(:,1:end-2)) - max(z0(:, 2:end-1), z0(:,1:end-2)); % Right side conductivity

h_s_y = [h_s(:,1), h_s_y, h_s(:,end)];

% Compute fluxes in x and y directions using upwind differences
q_x = h_s_x .* u;  % Discharge per unit width in x direction
q_y = h_s_y .* v;  % Discharge per unit width in y direction

% Using Local-Inertial Formulation
dq_x = (q_x - 9.81*h_s_x.*dt.*Sx)./(1 + 9.81*dt.*manning_n.^2*abs(q_x)./(h_s_x.^(7/3))); % meter per unit width
dq_y = (q_y - 9.81.*h_s_y.*dt.*Sy)./(1 + 9.81*dt.*manning_n.^2*abs(q_y)./(h_s_y.^(7/3))); % meter per unit width

dq_x(h_s_x <= 2*tol) = 0;
dq_y(h_s_y <= 2*tol) = 0;


% Update water depth using mass conservation
Vol_Flux = ([zeros(ny,1) , 1/dy*dq_x(:,1:(nx-1))] - 1/dy*dq_x ...
                         + 1/dx*dq_y - [1/dx*dq_y(2:end,:); zeros(1,nx)]);

Vol_Flux(2:end-1,2:end-1) = 0;

% Outlet Flux
out_flux = 1./manning_n.*h_s.^(2/3).*sqrt(0.02).*(dx + dy)/2*h_s./(dx*dy); % m/s

out_flux(2:end-1,2:end-1) = 0;

% out_flux(h_s<=2*tol) = 0;

h_s_new = h_s + dt * (rainfall + exfiltration_rate - infiltration + Vol_Flux - out_flux);

% Ensure no negative water depth
h_s_new(h_s_new < 0) = 0;

% Compute refreshed velocities using Manning's equation (fully vectorized)
u = dq_x*dx;
v = dq_y*dy;

end


function dt_new = compute_overland_timestep(h, u, v, Sx, Sy, dx, dy, n, g, C)
% COMPUTE_OVERLAND_TIMESTEP - Computes stable time step for 2D overland flow using CFL and diffusion conditions.
%
% Inputs:
% - h: Water depth matrix [m] (nx, ny)
% - u: Velocity in x-direction [m/s] (nx, ny)
% - v: Velocity in y-direction [m/s] (nx, ny)
% - Sx: Slope in x-direction (nx, ny)
% - Sy: Slope in y-direction (nx, ny)
% - dx: Grid spacing in x-direction [m]
% - dy: Grid spacing in y-direction [m]
% - n: Manning's roughness coefficient
% - g: Gravitational acceleration [m/s²]
% - C: Safety factor (0.5 - 0.9)
%
% Output:
% - dt_new: Stable time step [s]

% Compute wave celerity (vectorized)
c = sqrt(g * h); % Wave speed [m/s]

% Compute CFL condition (vectorized)
dt_CFL_x = dx ./ (u + c + 1e-9); % Avoid division by zero
dt_CFL_y = dy ./ (v + c + 1e-9);

% Compute diffusion stability condition (vectorized)
K_diff_x = (h.^(5/3)) ./ n .* abs(Sx);
K_diff_y = (h.^(5/3)) ./ n .* abs(Sy);
dt_diff = min(dx^2 ./ (K_diff_x + 1e-9), dy^2 ./ (K_diff_y + 1e-9));

% Compute final time step (vectorized)
dt_new = C * min([dt_CFL_x(:); dt_CFL_y(:); dt_diff(:)]);

% Ensure time step is positive
dt_new = max(dt_new, 1e-4);
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
