function [recharge_rate, updated_soil_moisture] = simulate_groundwater_recharge(infiltration_rate, initial_soil_moisture, alpha, dt)
    % Simulate groundwater recharge using a linear reservoir approach.
    %
    % INPUTS:
    % infiltration_rate      - The infiltration rate at the surface (mm/h)
    % initial_soil_moisture  - The initial soil moisture in the unsaturated zone (mm)
    % alpha                  - The linear coefficient representing the recharge rate (1/h)
    % dt                     - The time step (hours)
    %
    % OUTPUTS:
    % recharge_rate          - The computed recharge rate (mm/h)
    % updated_soil_moisture  - The updated soil moisture (mm)
    
    % Calculate the recharge rate based on the linear reservoir approach
    recharge_rate = alpha .* initial_soil_moisture;
    
    % Update soil moisture state
    updated_soil_moisture = initial_soil_moisture + dt .* (infiltration_rate - recharge_rate);
    
    % Ensure that soil moisture does not go below zero
    updated_soil_moisture = max(updated_soil_moisture, 0);
    
end
