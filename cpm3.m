function a_cpm = cpm3(a)
  a_cpm = zeros(3, 3);
  a_cpm(1, 2) = -a(3);
  a_cpm(1, 3) = a(2);
  a_cpm(2, 3) = -a(1);
  
  a_cpm += -1*transpose(a_cpm);
endfunction  