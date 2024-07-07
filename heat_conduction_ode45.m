function heat_conduction_ode45

    % Input Parameters
    nx = 100; % Number of grid points in x-direction
    ny = nx; % Number of grid points in y-direction
    Lx = 1; % Length of the domain in x-direction
    Ly = 1; % Length of the domain in y-direction
    alpha = 1.1e-3; % Thermal diffusivity
    nt = 1400; % Total number of time steps
    T_L = 400; % Left boundary temperature
    T_R = 600; % Right boundary temperature
    T_T = 800; % Top boundary temperature

    T_B = 900; % Bottom boundary temperature
    T_initial = 300; % Initial temperature for interior points
    
    % Create spatial grid
    x = linspace(0, Lx, nx);
    y = linspace(0, Ly, ny);
    [X, Y] = meshgrid(x, y);
    
    % Define initial condition
    T0 = T_initial * ones(nx, ny); % Initial temperature at all interior points
    
    % Set boundary conditions
    T0(:, 1) = T_L; % Left boundary
    T0(:, end) = T_R; % Right boundary
    T0(1, :) = T_T; % Top boundary
    T0(end, :) = T_B; % Bottom boundary
    
    % Reshape the temperature matrix into a column vector for ODE solver
    T0_vec = reshape(T0, [], 1);
    
    % Time span for ODE solver
    tspan = linspace(0, nt, nt+1); % Time steps
    
    % Solve the ODE using ode45
    [~, T_sol] = ode45(@heat_conduction_rhs, tspan, T0_vec);
    
    % Reshape the solution back to temperature matrix
    T_solution = reshape(T_sol', nx, ny, []);
    
    % Plotting the temperature contour over time
    figure;
    for k = 1:nt+1
        contourf(X, Y, T_solution(:, :, k));
        colorbar;
        colormap(jet);
        xlabel('X-Axis');
        ylabel('Y-Axis');
        title(sprintf('Temperature Contour at Time Step %d', k-1));
        pause(0.1); % Pause to display each time step
    end

end

function dTdt = heat_conduction_rhs(~, T_vec)
    % Function to compute the right-hand side of the ODE system
    
    % Reshape the temperature vector back to matrix
    nx = 100; % Number of grid points in x-direction
    ny = nx; % Number of grid points in y-direction
    T_matrix = reshape(T_vec, nx, ny);
    
    % Calculate spatial derivatives (using central differences)
    dx = 1 / (nx - 1); % Grid spacing in x-direction
    dy = 1 / (ny - 1); % Grid spacing in y-direction
    alpha = 1.1e-3; % Thermal diffusivity
    
    dTdt = zeros(nx, ny);
    
    for i = 2:nx-1
        for j = 2:ny-1
            d2T_dx2 = (T_matrix(i+1,j) - 2*T_matrix(i,j) + T_matrix(i-1,j)) / dx^2;
            d2T_dy2 = (T_matrix(i,j+1) - 2*T_matrix(i,j) + T_matrix(i,j-1)) / dy^2;
            dTdt(i,j) = alpha * (d2T_dx2 + d2T_dy2);
        end
    end
    
    % Flatten the matrix into a vector for ODE solver
    dTdt = reshape(dTdt, [], 1);

end
