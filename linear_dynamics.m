function A, B, C, D = lineardynamics(current_state)
  % q_dot = f(q_b, w_b)
  % w_b_dot = g(tau_b, w_b)
  % v_dot = (1/m)*h(q_b, v_g, u)
  % Force exerted in body frame is a sum of individual control inputs, u.
  xg = transpose(current_state(1:3));
  vg = transpose(current_state(4:6));
  qb2g = transpose(current_state(7:10));
  wb = transpose(current_state(11:13));
  
  F_w = zeros(4, 3);
  F_w(1, :) = -0.5*qb2g;
  F_w(2:, :) = 0.5*(qb2g(1)*eye(3) - cpm3(q(2:)));
  
  F_q_vec = zeros(4, 3);
  F_q_vec(1, :) = -0.5*wb;
  F_q_vec(2:, :) = 0.5*cpm(q(2:));
  
  F_q_sca = zeros(4, 1);
  F_q_sca(2:, 1) = 0.5*wb;
  
  
endfunction