function four_magnet_rigid_body
% 4-Magnet Rigid Body — 3 DOF (x, theta, phi) Control
% 13 states: [x, th, ph, dx, dth, dph, eix, eith, eiph, I1, I2, I3, I4]

clear physics_loop;

% --- setup ---
P.g    = 9.81;
t_end  = 10.0;
P.dt   = 0.01; %time between updates
opts   = odeset('RelTol',1e-6,'AbsTol',1e-8);

% --- mechanical stuff ---
P.m      = 4;       % kg
P.Jy     = 0.1;     % pitch inertia
P.Jx     = 0.1;     % roll inertia
P.Lx     = 0.1;     % m (front/rear dist)
P.Ly     = 0.08;    % m (left/right dist)

% --- magnet model ---
P.Imax  = 15;
P.nI    = 2.0264;
P.nD    = 0.5416;
P.Cmag  = 0.1061;

% --- electrical lag (L/R) ---
P.L     = 0.1;   % H
P.R     = 2.0;   % Ohm
P.tau   = P.L / P.R; % 50ms lag

% --- reference ---
deq    = 0.008;     % m, goal hover gap
d1     = deq+0.01;
d2     = deq;
t_step = 0.5;
P.d_nom  = deq;
P.x_ref     = @(t) ((t < t_step).*(P.d_nom - d1) + (t >= t_step).*(P.d_nom - d2));
P.x_dot_ref = @(t) 0*t;
% safety
P.d_min  = 0.001;
P.d_max  = 0.050;

% --- controller gains (3 PIDs) ---
P.Kpx = 1e3; P.Kdx = 4e2; P.Kix = 3e2; % x
P.Kpth = 75; P.Kdth = 20; P.Kith = 10; % theta
P.Kpph = 75; P.Kdph = 20; P.Kiph = 10; % phi

% --- initial conditions (13 states) ---
x0      = P.d_nom - d1;   % start at d1
th0     = 0; ph0 = 0;
dx0     = 0; dth0 = 0; dph0 = 0;
eix0    = 0; eith0 = 0; eiph0 = 0;
% initial current to hold weight
F_hold  = P.m * P.g / 4; % per magnet
I0 = ((F_hold * d1^P.nD) / P.Cmag)^(1/P.nI);
X0 = [x0; th0; ph0; dx0; dth0; dph0; eix0; eith0; eiph0; I0; I0; I0; I0];

% --- simulate ---
[t, X] = ode45(@(t,X) physics_loop(t,X,P), [0 t_end], X0, opts);

% --- unpack, post-process ---
th    = X(:,2);
ph      = X(:,3);
I_act   = X(:, 10:13);
[F_mag, I_cmd, gaps] = post_proc(t, X, P);
dref_gap  = arrayfun(@(t) (t < t_step).*d1 + (t >= t_step).*d2, t);

% --- plot it ---
figure(1); clf;
subplot(4,1,1)
plot(t, 1e3*gaps,'LineWidth',1.8); grid on
hold on; plot(t, 1e3*dref_gap,'--r','LineWidth',1.4)
ylabel('Gaps d_i [mm]'); title('4 mag sim')
legend('d1','d2','d3','d4','ref','Location','best')

subplot(4,1,2)
plot(t, th * 180/pi,'LineWidth',1.8); grid on; hold on
plot(t, ph * 180/pi,'LineWidth',1.8);
ylabel('Angle [deg]');
legend('Pitch','Roll','Location','best')

subplot(4,1,3)
plot(t, I_act,'LineWidth',1.8); grid on;
hold on; 
% NEW: use 'stairs' plot for I_cmd to show the sample-and-hold
stairs(t, I_cmd, '--', 'LineWidth', 1.0); 
yline(P.Imax,'--k');
ylabel('Current I [A]')
legend('I1 act','I2 act','I3 act','I4 act','Location','best')

subplot(4,1,4)
plot(t, F_mag, 'LineWidth',1.8); grid on
yline(P.m*P.g/4,'--'); ylabel('F_{mag} [N]')
xlabel('Time [s]');
legend('F1','F2','F3','F4','Location','best')
end

% main loop
function dX = physics_loop(t, X, P)
    % persistent variables, hold state between controller updates
    persistent t_last_update I_cmd_held
    
    % reset persistent vars if t=0 or t is less than last update (a new run)
    if isempty(t_last_update) || (t < t_last_update)
        t_last_update = -P.dt; % force update on first step
        I_cmd_held = zeros(4,1);
    end

    % unpack states
    x=X(1); th=X(2); ph=X(3);
    dx=X(4); dth=X(5); dph=X(6);
    eix=X(7); eith=X(8); eiph=X(9);
    I_act = X(10:13);
    
    % sample+hold logic
    if (t - t_last_update) >= P.dt
        t_last_update = t;
        
        % kinematics (get gaps)
        x1 = x + P.Lx*th + P.Ly*ph; % FL
        x2 = x + P.Lx*th - P.Ly*ph; % FR
        x3 = x - P.Lx*th + P.Ly*ph; % RL
        x4 = x - P.Lx*th - P.Ly*ph; % RR
        x_vec = [x1; x2; x3; x4];
        d_vec = P.d_nom - x_vec;
        
        % 3 PID controllers
        x_ref = P.x_ref(t);
        e_x  = x_ref - x;     edot_x = P.x_dot_ref(t) - dx;
        e_th = 0 - th;        edot_th = 0 - dth;
        e_ph = 0 - ph;        edot_ph = 0 - dph;
        
        Fz_cmd  = P.m*P.g - (P.Kpx*e_x + P.Kdx*edot_x + P.Kix*eix);
        Tau_th_cmd = 0    + (P.Kpth*e_th + P.Kdth*edot_th + P.Kith*eith); % sign fix
        Tau_ph_cmd = 0    + (P.Kpph*e_ph + P.Kdph*edot_ph + P.Kiph*eiph); % sign fix
        
        % allocation (virtual -> 4 forces)
        F_base   = Fz_cmd / 4;
        F_pitch  = Tau_th_cmd / (4 * P.Lx);
        F_roll   = Tau_ph_cmd / (4 * P.Ly);
        
        F_cmd_vec = [F_base - F_pitch - F_roll;  % F1
                     F_base - F_pitch + F_roll;  % F2
                     F_base + F_pitch - F_roll;  % F3
                     F_base + F_pitch + F_roll]; % F4
    
        % actuator (force -> current)
        I_cmd_target_vec = zeros(4,1);
        for i = 1:4
            F_cmd_clamped = max(0, F_cmd_vec(i));
            d_eff = min(max(d_vec(i), P.d_min), P.d_max);
            I_cmd_target = ((F_cmd_clamped * d_eff^P.nD) / P.Cmag)^(1/P.nI);
            I_cmd_target_vec(i) = max(0, min(P.Imax, I_cmd_target));
        end
        
        % hold this command
        I_cmd_held = I_cmd_target_vec;
    end

    % runs every ode45 step
    I_cmd_target_vec = I_cmd_held;

    % electrical lag
    dI_dt_vec = (I_cmd_target_vec - I_act) / P.tau;
    
    % force model (actual current -> actual force)
    % needs the gaps, must recalculate them for physics
    x1_phys = x + P.Lx*th + P.Ly*ph;
    x2_phys = x + P.Lx*th - P.Ly*ph;
    x3_phys = x - P.Lx*th + P.Ly*ph;
    x4_phys = x - P.Lx*th - P.Ly*ph;
    d_vec_phys = P.d_nom - [x1_phys; x2_phys; x3_phys; x4_phys];

    F_act_vec = zeros(4,1);
    for i = 1:4
        d_eff = min(max(d_vec_phys(i), P.d_min), P.d_max);
        F_act_vec(i) = P.Cmag * (I_act(i)^P.nI) / (d_eff^P.nD);
    end
    F1=F_act_vec(1); F2=F_act_vec(2); F3=F_act_vec(3); F4=F_act_vec(4);
    
    % disturbances
    F_disturb = 0; T_pitch_disturb = 0; T_roll_disturb = 0;
    if (t > 2.0 && t < 2.2)
        F_disturb = 15;     % vertical hit
        T_pitch_disturb = 0.2;  % pitch hit
    end
    if (t > 3.0 && t < 3.2)
        T_roll_disturb = 2;   % roll hit
    end
    if (t > 4.0 && t < 4.2)
        T_pitch_disturb = -2; % combo hit
        T_roll_disturb = 2;
    end
    
    % EOMs (physics)
    xdd   = P.g - (F1 + F2 + F3 + F4) / P.m + F_disturb / P.m;
    thdd = (P.Lx * (F3 + F4 - F1 - F2)) / P.Jy + T_pitch_disturb / P.Jy;
    phdd = (P.Ly * (F2 + F4 - F1 - F3)) / P.Jx + T_roll_disturb / P.Jx;

    % integrator derivatives (calculated inside  x ms loop)
    % need to update even if the controller doesn't run
    x_ref = P.x_ref(t);
    e_x  = x_ref - x;
    e_th = 0 - th;
    e_ph = 0 - ph;
    eix_dot  = e_x;
    eith_dot = e_th;
    eiph_dot = e_ph;

    % pack em up
    dX = [dx; dth; dph; ...
          xdd; thdd; phdd; ...
          eix_dot; eith_dot; eiph_dot; ...
          dI_dt_vec];
end

%
% plots processing
%
function [F_mag, I_cmd, gaps] = post_proc(t, X, P)
persistent t_last_update_plot I_cmd_held_plot
I_cmd_held_plot = zeros(4,1); % reset for plots
t_last_update_plot = -P.dt;

n = numel(t);
F_mag  = zeros(n, 4);
I_cmd = zeros(n, 4);
gaps  = zeros(n, 4);

for k = 1:n
    % get states
    x=X(k,1); th=X(k,2); ph=X(k,3);
    dx=X(k,4); dth=X(k,5); dph=X(k,6);
    eix=X(k,7); eith=X(k,8); eiph=X(k,9);
    I_act_vec = X(k, 10:13)';
    
    % get gaps
    x1 = x + P.Lx*th + P.Ly*ph;
    x2 = x + P.Lx*th - P.Ly*ph;
    x3 = x - P.Lx*th + P.Ly*ph;
    x4 = x - P.Lx*th - P.Ly*ph;
    x_vec = [x1; x2; x3; x4];
    d_vec = P.d_nom - x_vec;
    gaps(k,:) = d_vec;
   
    if (t(k) - t_last_update_plot) >= P.dt
        t_last_update_plot = t(k);
        
        % get PID commands
        x_ref = P.x_ref(t(k));
        e_x  = x_ref - x;     edot_x = P.x_dot_ref(t(k)) - dx;
        e_th = 0 - th;        edot_th = 0 - dth;
        e_ph = 0 - ph;        edot_ph = 0 - dph;
        Fz_cmd  = P.m*P.g - (P.Kpx*e_x + P.Kdx*edot_x + P.Kix*eix);
        Tau_th_cmd = 0    + (P.Kpth*e_th + P.Kdth*edot_th + P.Kith*eith);
        Tau_ph_cmd = 0    + (P.Kpph*e_ph + P.Kdph*edot_ph + P.Kiph*eiph);
        
        % get allocated forces
        F_base   = Fz_cmd / 4;
        F_pitch  = Tau_th_cmd / (4 * P.Lx);
        F_roll   = Tau_ph_cmd / (4 * P.Ly); % <-- TYPO FIX
        F_cmd_vec = [F_base - F_pitch - F_roll;
                     F_base - F_pitch + F_roll;
                     F_base + F_pitch - F_roll;
                     F_base + F_pitch + F_roll];
        
        % get target currents
        I_cmd_target_vec = zeros(4,1);
        for i = 1:4
            F_cmd_clamped = max(0, F_cmd_vec(i));
            d_eff = min(max(d_vec(i), P.d_min), P.d_max);
            I_cmd_target = ((F_cmd_clamped * d_eff^P.nD) / P.Cmag)^(1/P.nI);
            I_cmd_target_vec(i) = max(0, min(P.Imax, I_cmd_target));
        end
        I_cmd_held_plot = I_cmd_target_vec;
    end
    
    I_cmd(k, :) = I_cmd_held_plot; % store held command
    
    % get actual forces
    for i = 1:4
        d_eff = min(max(d_vec(i), P.d_min), P.d_max);
        F_mag(k, i) = P.Cmag * (I_act_vec(i)^P.nI) / (d_eff^P.nD);
    end
end
end