function [alphaT_nonlocal, k_nonloc] = critvaluenonlocal2(L, R, D, k2, ni, n)

% Note L is the domain total length 

k = 2*pi.*n./L;
K = sin(k.*R)./(k.*R);
knl_validi = k(K>0);
K_validi = K(K>0);
num_nl = D*knl_validi.^4 + (D*ni./k2 + k2).*knl_validi.^2 + ni - k2;
den_nl = K_validi.*(ni - k2).*knl_validi.^2;
alphanonloc_k = num_nl./den_nl;
[alphaT_nonlocal, c_nl] = min(alphanonloc_k); 
k_nonloc = knl_validi(c_nl); 

end