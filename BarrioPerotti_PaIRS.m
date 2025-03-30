%% Water propulsion function
function dydt_water = waterODEs(t, params, variables)
    % Unpacking variables and constants
    m_w = variables(1);
    v_b = variables(2);
    v_wn = variables(3);
    p_a = variables(4);

    m_a = params.m_a;
    m_b = params.m_b;
    S = params.S;
    S_b = params.S_b;
    p_atm = params.p_atm;
    rho_w = params.rho_w;
    rho_atm = params.rho_atm;
    V = params.V;
    V_w0 = params.V_w0;
    p_a0 = params.p_a0;
    mu = params.mu;
    gamma = params.gamma;
    theta = params.theta;
    g = params.g;
    CD = params.CD;
    
    % Calculating forces
    F_D = 0.5 * rho_atm * CD * S_b * v_b^2;
    F_M = 0; % (Uncomment and adjust if friction is considered)
    W_x = (m_b + m_a + m_w) * g * sin(theta);
    
    % Equation (18): Mass conservation
    ode_m_w = -rho_w * v_wn * S;
    
    % Equation (19): Rocket velocity v_b
    ode_v_b = -(F_D + F_M + W_x) / (m_b + m_a) ...
               + ( S_b * (p_a - p_atm) - ((S_b - S)^2 / (2 * S_b)) * rho_w * v_wn^2 ) / (m_b + m_a);
    
    % Equation (20): Water exit velocity v_wn
    term1 = (S_b^2 * (p_a - p_atm)) / S * (1 / m_w + 1 / (m_b + m_a));
    term2 = (v_wn^2 * rho_w) / (2 * S) * (S_b^2 / m_w + (S_b - S)^2 / (m_b + m_a));
    term3 = S_b / (S * (m_b + m_a)) * (F_D + F_M + W_x);
    ode_v_wn = term1 - term2 - term3;
    
    % Equation (21): Air pressure p_a evolution
    ode_p_a = -gamma * (p_a^((gamma+1)/gamma)) / ((V - V_w0) * p_a0^(1/gamma)) * v_wn * S;
    
    % Packing into vector
    dydt_water = [ode_m_w; ode_v_b; ode_v_wn; ode_p_a];
end

%% Air propulsion function
function dydt_air = airODEs(t, params, variables)
    % Unpacking variables and state variables
    v_b = variables(1);
    m_a = variables(2);  % Air mass (state variable)
    p_a = variables(3);

    m_b = params.m_b;
    S = params.S;
    S_b = params.S_b;
    p_atm = params.p_atm;
    rho_w = params.rho_w;
    rho_atm = params.rho_atm;
    V = params.V;
    V_w0 = params.V_w0;
    p_a0 = params.p_a0;
    mu = params.mu;
    gamma = params.gamma;
    theta = params.theta;
    g = params.g;
    CD = params.CD;
    rho_a0 = params.rho_a0;
    m_w = 0; % Assume that all water is depleted

    % Calculating forces
    F_D = 0.5 * rho_atm * CD * S_b * v_b^2;
    F_M = 0; % (Uncomment and adjust if friction is considered)
    W_x = (m_b + m_a + m_w) * g * sin(theta);

    % Compute nozzle pressure p_an (this is an assumed relationship)
    p_an = p_a * (2 / (gamma+1))^(gamma/(gamma-1));
    
    % Equation (29): Rocket acceleration during air propulsion
    ode_v_b = S / (m_b + m_a) * ( (2*gamma/(gamma-1)) * (p_a^((gamma-1)/gamma) - p_an^((gamma-1)/gamma)) * p_an^(1/gamma) ) ...
              - (F_D + F_M + W_x) / (m_b + m_a);
    
    % Equation (30): Air mass evolution
    ode_m_a = - S * p_an^(1/gamma) * sqrt( (2*gamma)/(gamma-1) * rho_a0 / (p_a0^(1/gamma)) * (p_a^((gamma-1)/gamma) - p_an^((gamma-1)/gamma) ) );
    
    % Equation (31): Air pressure evolution    
    ode_p_a = - (S / V) * sqrt( (2*gamma)/(gamma-1) * (p_a0^(1/gamma) / rho_a0) ) * ...
              ( p_a^((gamma-1)/gamma) - p_an^((gamma-1)/gamma) )^(1/2) * p_an^(1/gamma) * p_a^((gamma-1)/gamma);
    
    % Packing into vector
    dydt_air = [ode_v_b; ode_m_a; ode_p_a];
end

%% Ballistic flight function
function dydt_ballistic = ballisticODEs(t, params, variables)
    % Unpacking variables and constants
    position_x = variables(1);
    position_y = variables(2);
    vx = variables(3);
    vy = variables(4);

    m_b = params.m_b;
    % Other parameters (m_a, etc.) may not be relevant in the ballistic phase
    rho_atm = params.rho_atm;
    g = params.g;
    CD = params.CD;
    S_b = params.S_b;
    
    % ODEs for positions
    ode_x = vx;
    ode_y = vy;

    % Equation (32): Ballistic flight in x axis
    ode_vx = -((rho_atm * CD * S_b) / (2 * m_b)) * vx * sqrt(vx^2 + vy^2);

    % Equation (33): Ballistic flight in y axis
    ode_vy = -g - ((rho_atm * CD * S_b) / (2 * m_b)) * vy * sqrt(vx^2 + vy^2);

    dydt_ballistic = [ode_x; ode_y; ode_vx; ode_vy];
end

%% Main Script

% Define constants & initial conditions
V = 0.0020; % !
S_b = (((0.10)/2)^2 * pi); % !

V_w0 = V*(1/3);
rho_w = 998;
rho_atm = 1.225;
m_payload = V_w0*rho_w;
m_w = V_w0*rho_w;
m_b = m_payload + 1.00; % !
S = (((20.5e-3)/2)^2) * pi;
p_atm = 101325;

p_a0 = 10 * p_atm;

mu = 0.4; % Not needed.
gamma = 1.4; % Better guess

theta = 82 * (pi/180);
g = 9.81;
CD = 0.4; % Guess

% Compute ambient air density (rename to match parameter usage)
rho_a0 = p_a0 / (287 * (10 + 273.15));
m_a_initial = rho_a0 * (V - V_w0);

% Assemble parameter structure
params.m_a = m_a_initial; % initial air mass (may change during air phase)
params.m_b = m_b;
params.S = S;
params.S_b = S_b;
params.p_atm = p_atm;
params.rho_w = rho_w;
params.rho_atm = rho_atm;
params.V = V;
params.V_w0 = V_w0;
params.p_a0 = p_a0;
params.mu = mu;
params.gamma = gamma;
params.theta = theta;
params.g = g;
params.CD = CD;
params.rho_a0 = rho_a0;

%% Phase 1: Water propulsion
% Initial water mass is computed as water volume * density (kg)
variables_water0 = [V_w0 * rho_w, 0, 0, p_a0]';
% Use a finite tspan, e.g. [0, 5] seconds, and rely on event functions to stop early
tspan_water = [0 5];

% Event function: stop when water mass reaches 0
options_water = odeset('Events',@mass_is_zero);

[t_water, variables_water] = ode45(@(t,y) waterODEs(t, params, y), tspan_water, variables_water0, options_water);
t_water = t_water;  % Save time vector for water phase

%% Phase 2: Air propulsion
% Use the final water phase values to initialize the air phase.
% Use the final rocket velocity and air pressure from water phase.
variables_air0 = [variables_water(end,2), m_a_initial, variables_water(end,4)]';
tspan_air = [0 5]; % finite tspan

options_air = odeset('Events',@pressure_is_atm);

[t_air, variables_air] = ode45(@(t,y) airODEs(t, params, y), tspan_air, variables_air0, options_air);
t_air = t_air;  % Save time vector for air phase

%% Phase 3: Ballistic flight
v_net = variables_air(end,1);
% Initial positions at (0,0) and velocity components based on rocket angle
variables_ballistic0 = [0, 0, v_net * cos(theta), v_net * sin(theta)]';
tspan_ballistic = [0 20]; % Choose a tspan long enough for flight

options_ballistic = odeset('Events', @y_is_zero);

[t_ballistic, variables_ballistic] = ode45(@(t,y) ballisticODEs(t, params, y), tspan_ballistic, variables_ballistic0, options_ballistic);
t_ballistic = t_ballistic;  % Save time vector for ballistic phase

disp('Total horizontal distance travelled:')
disp(variables_ballistic(end, 1))

%% Event Functions (Nested at the end of the file)
function [value, isterminal, direction] = mass_is_zero(t, variables)
    % Stop integration when water mass is nearly zero.
    value = variables(1);  % When m_w is zero.
    isterminal = 1;
    direction = -1;
end

function [value, isterminal, direction] = pressure_is_atm(t, variables)
    % Stop integration when air pressure reaches atmospheric pressure.
    value = variables(3) - 101325;
    isterminal = 1;
    direction = -1;
end

function [value, isterminal, direction] = y_is_zero(t, variables)
    % Stop integration when vertical position returns to ground.
    value = variables(2);  % when y is zero.
    isterminal = 1;
    direction = -1;
end


%% Plotting - Water Propulsion Phase
figure;
subplot(2,2,1);
plot(t_water, variables_water(:,1), 'LineWidth',1.5);
title('Water Phase: Water Mass (m_w)');
xlabel('Time (s)'); ylabel('Mass (kg)');
grid on;

subplot(2,2,2);
plot(t_water, variables_water(:,2), 'LineWidth',1.5);
title('Water Phase: Rocket Velocity (v_b)');
xlabel('Time (s)'); ylabel('Velocity (m/s)');
grid on;

subplot(2,2,3);
plot(t_water, variables_water(:,3), 'LineWidth',1.5);
title('Water Phase: Water Exit Velocity (v_{wn})');
xlabel('Time (s)'); ylabel('Velocity (m/s)');
grid on;

subplot(2,2,4);
plot(t_water, variables_water(:,4), 'LineWidth',1.5);
title('Water Phase: Air Pressure (p_a)');
xlabel('Time (s)'); ylabel('Pressure (Pa)');
grid on;

%% Plotting - Air Propulsion Phase
figure;
subplot(3,1,1);
plot(t_air, variables_air(:,1), 'LineWidth',1.5);
title('Air Phase: Rocket Velocity (v_b)');
xlabel('Time (s)'); ylabel('Velocity (m/s)');
grid on;

subplot(3,1,2);
plot(t_air, variables_air(:,2), 'LineWidth',1.5);
title('Air Phase: Air Mass (m_a)');
xlabel('Time (s)'); ylabel('Mass (kg)');
grid on;

subplot(3,1,3);
plot(t_air, variables_air(:,3), 'LineWidth',1.5);
title('Air Phase: Air Pressure (p_a)');
xlabel('Time (s)'); ylabel('Pressure (Pa)');
grid on;

%% Plotting - Ballistic Flight Phase
figure;
% Trajectory (x vs. y)
subplot(2,1,1);
plot(variables_ballistic(:,1), variables_ballistic(:,2), 'LineWidth',1.5);
title('Ballistic Phase: Trajectory');
xlabel('Horizontal Position (m)'); ylabel('Vertical Position (m)');
grid on;
axis equal;  

% Velocities over time (v_x and v_y)
subplot(2,1,2);
plot(t_ballistic, variables_ballistic(:,3), 'b', t_ballistic, variables_ballistic(:,4), 'r', 'LineWidth',1.5);
title('Ballistic Phase: Velocity vs. Time');
xlabel('Time (s)'); ylabel('Velocity (m/s)');
legend('v_x', 'v_y');
grid on;
axis equal; 