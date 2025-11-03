function one_magnet_height_control
% single electromagnet — 1 dof gap control
% states: x (downward displacement), dx, ei (integral error), i_actual (current)
% gap:    d = d_nom - x   [m]   (keep d positive)
% model:  fmag = cmag * i^ni / d^nd
% tau * di/dt = i_cmd - i_actual

%% --- sim ---
g      = 9.81;        % [m/s^2]
t_end  = 5.0;         % [s]
opts   = odeset('RelTol',1e-6,'AbsTol',1e-8);

%% --- mechanical ---
m_eff    = 4;           % [kg] effective levitated mass
k_struct = 0;           % [n/m] stiffness
c_struct = 0;           % [n*s/m] damping

%% --- magnet model ---
% fmag = cmag * i^ni / d^nd
imax   = 15;          % [a] current supply limit
ni     = 2.0264;      % exponent on current (fit)
nd     = 0.5416;      % exponent on gap distance (fit)
cmag   = 0.1061;      % [n*m^nd / a^ni] fitted constant
tau_l  = 0.01;        % [s] l/r time constant for current lag

%% --- gap reference ---
deq    = 0.004;       % [m] nominal hover gap
d1     = deq + 0.01;  % [m] initial gap (0..t_step)
d2     = deq;         % [m] target after t_step
t_step = 0.5;         % [s]
d_ref_fun    = @(t) (t < t_step).*d1 + (t >= t_step).*d2;
d_refdot_fun = @(t) 0*t;   % hold position
% safety range (physics/sensor limits)
d_min  = 0.001;       % [m]
d_max  = 0.050;       % [m]

%% --- controller (gap pid on d) ---
kp  = 2e4;    
kd  = 2e3;    
ki  = 1e3;    
kaw = 200;    % [1/s] anti-windup back-calc gain

%% --- pack parameters ---
p.g = g; p.m = m_eff; p.k = k_struct; p.c = c_struct;
p.cmag = cmag; p.ni = ni; p.nd = nd; p.imax = imax;
p.d_min = d_min; p.d_max = d_max;
p.kp = kp; p.kd = kd; p.ki = ki; p.kaw = kaw;
p.d_ref_fun = d_ref_fun; p.d_refdot_fun = d_refdot_fun;
p.d_nom = deq;        % x=0 ↔ d=deq
p.tau = tau_l;        % [s] l/r time constant

%% --- initial conditions ---
x0  = p.d_nom - d1;   % so d(0) = d1
dx0 = 0;
ei0 = 0;
% initial current to hold weight at d1
f_eq = p.m * p.g;
i0 = ((f_eq * d1^p.nd) / p.cmag)^(1/p.ni);
i0 = min(i0, p.imax); % saturate if needed

x0  = [x0; dx0; ei0; i0]; % add i_actual as 4th state

%% --- simulate ---
[t, x] = ode45(@(t,x) dyn_gap_pid(t,x,p), [0 t_end], x0, opts);

% unpack
x_pos    = x(:,1);
x_vel    = x(:,2);
ei       = x(:,3);
i_actual = x(:,4); % i_actual
d        = p.d_nom - x_pos;
dref     = arrayfun(p.d_ref_fun, t);

% postcompute commanded signals & actual force
[~, i_cmd, f_cmd] = post_signals(t, x_pos, x_vel, ei, p);
% calculate actual force from simulated actual current
d_eff_vec = min(max(d, p.d_min), p.d_max);
% clamp i_actual for post-processing too, just in case
fmag_actual = p.cmag * (max(0,i_actual).^p.ni) ./ (d_eff_vec.^p.nd);

%% --- plots ---
figure(1); clf
subplot(4,1,1)
plot(t, 1e3*d,'LineWidth',1.8); grid on
hold on; plot(t, 1e3*dref,'--','LineWidth',1.4)
ylabel('gap d [mm]'); title('single magnet — gap control with current lag')
legend('d','d_{ref}','Location','best')

subplot(4,1,2)
plot(t, i_cmd,'--','LineWidth',1.4); grid on; hold on
plot(t, i_actual, 'LineWidth', 1.8);                 
yline(p.imax,'--r'); ylabel('current i [a]')
legend('i_{cmd}', 'i_{actual}', 'Location','best')    

subplot(4,1,3)
plot(t, fmag_actual,'LineWidth',1.8); grid on         
yline(p.m*p.g,'--'); ylabel('f_{mag} (actual) [n]') 
legend('f_{actual}', 'm*g', 'Location','best')

subplot(4,1,4)
plot(t, f_cmd,'LineWidth',1.8); grid on
xlabel('time [s]'); ylabel('f_{cmd} [n]')
end

%% ===== dynamics =====
function dx = dyn_gap_pid(t, x_state, p)
    % unpack 4 states
    x        = x_state(1);
    dx       = x_state(2);
    ei       = x_state(3);
    i_actual = max(0, x_state(4)); % prevent negative current -> complex num -> nan

    % gap and rate
    d  = p.d_nom - x;
    dd = -dx;

    % reference and errors
    d_ref   = p.d_ref_fun(t);
    d_ref_d = p.d_refdot_fun(t);
    e       = d_ref   - d;
    edot    = d_ref_d - dd;

    % unsaturated force command
    f_cmd_unsat = p.m*p.g + p.kp*e + p.kd*edot + p.ki*ei;

    % saturate force, not just current
    d_eff   = min(max(d, p.d_min), p.d_max);
    f_max   = p.cmag * (p.imax^p.ni) / (p.d_min^p.nd); % physical max at imax, smallest gap
    f_cmd   = min(max(f_cmd_unsat, 0), f_max);         % clamp force command

    % map clamped force to coimmanded current
    if f_cmd <= 0
        i_cmd = 0;
    else
        i_cmd = ((f_cmd * d_eff^p.nd) / p.cmag)^(1/p.ni);
        i_cmd = max(0, min(p.imax, i_cmd));
    end

    % realized magnetic force from acctual current
    fmag = p.cmag * (i_actual.^p.ni) / (d_eff.^p.nd);

    % back calc antiwindup using force clamp
    ei_dot = e + p.kaw * (f_cmd - f_cmd_unsat);

    % plant dynamics
    f_disturb = 0;
    if (t > 2.0 && t < 2.5)
        f_disturb = 15; % 15 n, pushing the mass away from mag
    end
    
    % +g (down), -fmag (up)
    xdd = (p.m*p.g - fmag - p.k*x - p.c*dx + f_disturb)/p.m;

    % 1st order lag
    di_dt = (1/p.tau) * (i_cmd - i_actual);

    % pack 4 d/dts
    dx = [dx; xdd; ei_dot; di_dt];
end

%% postcomutee
function [fmag, i_cmd, f_cmd] = post_signals(t, x, dx, ei, p)
% calculates commanded signals for plotting
n = numel(t);
fmag = zeros(n,1); i_cmd = zeros(n,1); f_cmd = zeros(n,1);
for k = 1:n
    d  = p.d_nom - x(k);
    dd = -dx(k);
    e    = p.d_ref_fun(t(k)) - d;
    edot = p.d_refdot_fun(t(k)) - dd;
    
    % unsaturated force command
    f_cmd_unsat = p.m*p.g + p.kp*e + p.kd*edot + p.ki*ei(k);
    
    % saturate force
    f_max   = p.cmag * (p.imax^p.ni) / (p.d_min^p.nd);
    f_cmd(k,1) = min(max(f_cmd_unsat, 0), f_max);
    
    % map clamped force to current
    d_eff = min(max(d, p.d_min), p.d_max);
    if f_cmd(k,1) <= 0
        i = 0;
    else
        i = ((f_cmd(k,1) * d_eff^p.nd) / p.cmag)^(1/p.ni);
        i = max(0, min(p.imax, i));
    end
    i_cmd(k,1) = i;
    
    %  commanded force
    fmag(k,1)  = p.cmag * (i.^p.ni) / (d_eff.^p.nd);
end
end