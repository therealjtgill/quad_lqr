function [state_history, euler_angles, stimuli] = runner
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
  
  state_history = [];
  dt = 0.001/2;
  old_states = [0, 0, 10, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0];
  stimuli = [];
  run_length = 20000;
%  state_history = [state_history; old_states];
  state_history = zeros(length(old_states), run_length);
  euler_angles = zeros(3, run_length);
  for i = 1:20000
    T = neutral_buoyancy_pulse;
    lambda = 0;
    zeta = 0;
    q = old_states(7:10);
    phi = atan2(2*(q(1)*q(2) + q(3)*q(4)), 1 - 2*(q(2)^2 + q(3)^2));
    theta = asin(2*(q(1)*q(3) - q(2)*q(4)));
    psi = atan2(2*(q(1)*q(4) + q(2)*q(3)), 1 - 2*(q(3)^2 + q(4)^2));
    if i < 2000
      lambda = 10;
      zeta = (T - sqrt(T^2 - 4*lambda^2 - 4*T*lambda))/2;
    end
    u_motor_speeds = [...
      T + lambda,...
      T + lambda,...
      T - zeta,...
      T - zeta ...
    ];
    [ext_force_b, ext_torque_b] = external_stimuli(...
      u_motor_speeds,...
      k_force,...
      k_torque,...
      arm_lengths...
    );
    stimuli = [stimuli; transpose(ext_force_b), transpose(ext_torque_b)];
    new_states = forward_dynamics(...
      old_states,...
      ext_force_b,...
      ext_torque_b,...
      body_mass,...
      J_inert,...
      dt...
    );
%    state_history = [state_history; new_states];
    euler_angles(:, i) = [phi, theta, psi];
    state_history(:, i) = new_states;
    old_states = new_states;
  end  
endfunction  