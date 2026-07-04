function [alphaT_nonlocal, k_nonloc] = critvaluenonlocal3(L, R, D1, D2, k1, k2, k3, k4, n)

% Note L is the domain total length 

k = 2*pi.*n./L;
K = sin(k.*R)./(k.*R);
knl_validi = k(K>0);
K_validi = K(K>0);
q6 = (D1*D2).*knl_validi.^6;
q4 = (((k1*k4 + k2*k3)*D1*D2)/(k2*(k3+k4)) + k2*D2 + k4*D1).*knl_validi.^4;
q2 = ( (k4*(k1-k2)*k2*D2)/(k2*(k3+k4))+ (k4*(k1*k4+k2*k3)*D1)/(k2*(k3+k4)) + k2*k4).*knl_validi.^2;
q0 = (k4*(k1*k4 + k2*k3)*(k1-k2))/(k1*(k3+k4));
num = q6 +q4 +q2 +q0;
den = K_validi.*((D2*k4*(k1-k2))/(k3+k4).*knl_validi.^4 + (k4^2*(k1-k2)*(k1*k4 + k2*k3))/(k1*(k3+k4)^2).*knl_validi.^2);
alphanonloc_k = num ./ den;
[alphaT_nonlocal, c_nl] = min(alphanonloc_k); 
k_nonloc = knl_validi(c_nl); 

end