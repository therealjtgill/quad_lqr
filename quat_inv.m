function qinv = quat_inv(q)
  qinv = [q(1); -q(2:4)];
endfunction  