function [A, B, C, D] = linear_dynamics(...
  mass,...
  inertia_tensor,...
  nu,...    % thrust coefficient for motors
  eta,...   % body yaw coefficient
  rho,...   % air density
  c_d,...   % damping coefficient for air frame
  a_ref,... % reference area for air frame
  k,...     % thrust ratio
  current_state,...
  control_input...
  )
  % q_dot = f(q_b, w_b)
  % w_b_dot = g(u, w_b)
  % v_dot = (1/m)*h(q_b, v_g, u)
  % Force exerted in body frame is a sum of individual control inputs, u.
  xg = transpose(current_state(1:3));
  vg = transpose(current_state(4:6));
  qb2g = transpose(current_state(7:10));
  wb = transpose(current_state(11:13));
  
  % [x y z vz vy vz phi theta psi wbx wby wbz]
  
  % Force in body frame is in z_hat for a quadrotor.
  force_body = [0; 0; k*sum(control_input)];
  f_t = force_body(3);
  
  f_w = zeros(4, 3);
  f_w(1, :) = -0.5*qb2g(2:end);
  f_w(2:end, :) = 0.5*(qb2g(1)*eye(3) + cpm3(qb2g(2:end)));
  
  f_q_vec = zeros(4, 3);
  f_q_vec(1, :) = -0.5*wb;
  f_q_vec(2:end, :) = -0.5*cpm3(wb);

  f_q_sca = zeros(4, 1);
  f_q_sca(2:end, 1) = 0.5*wb;

  % Hack: delete effects of air resistance on body rates.
  lambda = 0.0;
  g_u = inv(inertia_tensor)*[[0 nu 0 -nu]; [nu 0 -nu 0]; [-eta eta -eta eta]];
  g_w = inv(inertia_tensor)*(-lambda*eye(3) - cpm3(wb)*(inertia_tensor - eye(3)));
  
  q_sca = qb2g(1);
  q_vec = qb2g(2:end);
  h_q_sca = 2*qb2g(1)*force_body - 2*cross(q_vec, force_body);
  
  h_q_vec = -2*q_vec*force_body' + 2*q_sca*cpm3(force_body) + 2*force_body*q_vec' + 2*dot(q_vec, force_body)*eye(3);
  
  h_v = zeros(3);
  if norm(vg) > 1e-6
    h_v = -rho*norm(vg)*c_d*a_ref*eye(3) - rho*c_d*a_ref*(norm(vg)^-1)*vg*vg';
  else
    h_v = zeros(3);
  end  
  
  
  k
  % Use product rule to get partial of h with respect to u.
  h_fb = (q_sca^2 - norm(q_vec)^2)*eye(3) - 2*q_sca*cpm3(q_vec) + 2*q_vec*q_vec';
  fb_u = zeros(3, 4);
  fb_u(3, :) = [k k k k];
  h_u = h_fb*fb_u;
  
  A = [ [0 0 0 1 0 0 0 0 0 0 0 0];...
        [0 0 0 0 1 0 0 0 0 0 0 0];...
        [0 0 0 0 0 1 0 0 0 0 0 0];...
        [0 0 0 0 0 0 0 -f_t 0 0 0 0];...
        [0 0 0 0 0 0 f_t 0 0 0 0 0];...
        [0 0 0 0 0 0 0 0 0 0 0 0];...
        [0 0 0 0 0 0 0 0 0 1 0 0];...
        [0 0 0 0 0 0 0 0 0 0 1 0];...
        [0 0 0 0 0 0 0 0 0 0 0 1];...
        [0 0 0 0 0 0 0 0 0 0 0 0];...
        [0 0 0 0 0 0 0 0 0 0 0 0];...
        [0 0 0 0 0 0 0 0 0 0 0 0];...
      ];

  B = [ [0 0 0 0];...
        [0 0 0 0];...
        [0 0 0 0];...
        [0 0 0 0];...
        [0 0 0 0];...
        [1 1 1 1];...
        [0 0 0 0];...
        [0 0 0 0];...
        [0 0 0 0];...
        [g_u(1, :)];...
        [g_u(2, :)];...
        [g_u(3, :)];...
      ];
      
  C = zeros(12);
  % Want to control quad position, so we include this term.
  C(1, 1) = 1;
  C(2, 2) = 1;
  C(3, 3) = 1;
  
  D = zeros(4);
  
##  xg = transpose(current_state(1:3));
##  vg = transpose(current_state(4:6));
##  qb2g = transpose(current_state(7:10));
##  wb = transpose(current_state(11:13));
endfunction