function [quad_controlled, K_prime] = make_lqr_controller(x_f, y_f, z_f, Q, R)
  J = 2*eye(3);
  J(3, 3) = 1.01;
  usca = 9.8/0.3;
  u = [usca/4, usca/4, usca/4, usca/4];
  [A, B, C, D] = linear_dynamics_euler(1, J, 0.11, 0.12, 0.13, 0.2, 0.5, 0.3, [0 0 0 0 0 0 1 0 0 0 0 0 0.0], u);
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
  
  t = 0:0.1:20;
  u = max(0, min(t-1, 1));
  u_vec = [x_f*u; y_f*u; z_f*u; 0*u; 0*u; 0*u; 0*u; 0*u; 0*u; 0*u; 0*u; 0*u;]';
  [Y, t_out, X] = lsim(quad_controlled, u_vec, t);
  figure; plot(t_out, Y(:, 1)); figure; plot(t_out, Y(:, 2)); figure; plot(t_out, Y(:, 3));
  
endfunction
