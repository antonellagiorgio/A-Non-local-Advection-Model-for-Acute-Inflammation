function [alphaT_local, k_loc] = critvaluelocal3(L, D1, D2, k1, k2, k3, k4, n)

% Note L is the domain total length 

k = 2*pi.*n./L; 
q6 = (D1*D2).*k.^6;
q4 = (((k1*k4 + k2*k3)*D1*D2)/(k2*(k3+k4)) + k2*D2 + k4*D1).*k.^4;
q2 = ( (k4*(k1-k2)*k2*D2)/(k2*(k3+k4))+ (k4*(k1*k4+k2*k3)*D1)/(k2*(k3+k4)) + k2*k4).*k.^2;
q0 = (k4*(k1*k4 + k2*k3)*(k1-k2))/(k1*(k3+k4));
num = q6 +q4 +q2 +q0;
den = (D2*k4*(k1-k2))/(k3+k4).*k.^4 + (k4^2*(k1-k2)*(k1*k4 + k2*k3))/(k1*(k3+k4)^2).*k.^2;
alphaloc_k = num ./ den;
[alphaT_local, c_loc] = min(alphaloc_k); 
k_loc = k(c_loc); 

end


