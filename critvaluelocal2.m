function [alphaT_local, k_loc] = critvaluelocal2(L, D, k2, ni, n)

% Note L is the domain total length 

k = 2*pi.*n./L; 
num_loc = D*k.^4 + (D*ni./k2 + k2).*k.^2 + ni - k2;
den_local = (ni - k2).*k.^2;
alphaloc_k = num_loc./den_local;
[alphaT_local, c_loc] = min(alphaloc_k); 
k_loc = k(c_loc); 

end




