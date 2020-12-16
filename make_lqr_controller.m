% Good Q values:
%    10.00000           0           0           0           0           0           0           0           0           0           0           0
%           0    10.00000           0           0           0           0           0           0           0           0           0           0
%           0           0   100.00000           0           0           0           0           0           0           0           0           0
%           0           0           0   100.00000           0           0           0           0           0           0           0           0
%           0           0           0           0   100.00000           0           0           0           0           0           0           0
%           0           0           0           0           0   100.00000           0           0           0           0           0           0
%           0           0           0           0           0           0     1.00000           0           0           0           0           0
%           0           0           0           0           0           0           0     1.00000           0           0           0           0
%           0           0           0           0           0           0           0           0     1.00000           0           0           0
%           0           0           0           0           0           0           0           0           0     0.10000           0           0
%           0           0           0           0           0           0           0           0           0           0     0.10000           0
%           0           0           0           0           0           0           0           0           0           0           0     0.10000

% Good R values:
%   0.010000          0          0          0
%          0   0.010000          0          0
%          0          0   0.010000          0
%          0          0          0   0.010000

function [quad, quad_controlled, K_prime, L] = make_lqr_controller(x_f, y_f, z_f, Q, R)
  arm_mass = 0.2; %kg
  arm_length = 0.6; %meter
  arm_lengths = [arm_length, arm_length, arm_length, arm_length];
  motor_mass = 0.065; %kg
  body_mass = arm_mass*2 + motor_mass*4; %kg
  arms_inertia = zeros(3, 3);
  arms_inertia(1, 1) = (1/12)*arm_mass*arm_length^2;
  arms_inertia(2, 2) = (1/12)*arm_mass*arm_length^2;
  arms_inertia(3, 3) = (1/12)*arm_mass*arm_length^2;
  motors_inertia = zeros(3, 3);
  motors_inertia_tensor(1, 1) = 2*motor_mass*arm_length**2;
  motors_inertia_tensor(2, 2) = 2*motor_mass*arm_length**2;
  motors_inertia_tensor(3, 3) = 4*motor_mass*arm_length**2;
  
  J_inert = arms_inertia + motors_inertia;
  
  g = 9.81;
  neutral_buoyancy_pulse = 750;
  k_force = body_mass*g/(4*neutral_buoyancy_pulse^2);
  
  prop_radius = 0.1016; % meters
  prop_pitch = 0.1143; % meters
  c_d = 0.09; % airfoil
  % Max body angular velocity about z axis
  body_max_omega_z = 3600*(pi/180); % radians/second;
  
  k_torque = 0.5*(prop_radius^4)*prop_pitch*c_d*(body_max_omega_z^2)/(10^2);
  
  usca = 9.8/0.3;
  u = [usca/4, usca/4, usca/4, usca/4];
  [A, B, C, D] = linear_dynamics_euler(body_mass, J_inert, arm_length*k_force, k_torque, 0, 0, 0, k_force, [0 0 0 0 0 0 1 0 0 0 0 0 0.0], u);
  quad = ss(A, B, C);
  rank(ctrb(quad)) == 12;
##  Q = 10*eye(12);
##  Q(1, 1) = 10;
##  Q(2, 2) = 10;
##  Q(3, 3) = 10;
##  Q(4, 4) = 100;
##  Q(5, 5) = 100;
##  Q(6, 6) = 100;
##  Q(7, 7) = 10;
##  Q(8, 8) = 10;
##  Q(9, 9) = 10;
##  R = 5*eye(4);
  [K, P, L] = lqr(quad, Q, R);
  % Get rid of extremely small values in the feed back gain matrix.
  mask = K > 1e-9;
  K_prime = K.*mask;
  
  quad_controlled = ss((A - B*K_prime), -A + B*K_prime, C);
  
  t = 0:0.1:500;
  u = max(0, min(t-1, 1));
  u_vec = [x_f*u; y_f*u; z_f*u; 0*u; 0*u; 0*u; 0*u; 0*u; 0*u; 0*u; 0*u; 0*u;]';
  [Y, t_out, X] = lsim(quad_controlled, u_vec, t);
  figure; plot(t_out, Y(:, 1)); figure; plot(t_out, Y(:, 2)); figure; plot(t_out, Y(:, 3));
  
  figure; plot(Y(:, 1), Y(:, 3));
  figure; plot(Y(:, 1), Y(:, 2));
  
endfunction
