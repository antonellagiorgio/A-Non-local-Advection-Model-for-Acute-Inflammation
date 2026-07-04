clear
clc


%% Difference between alpha_T non local and alpha_T local varying R (F-S)
Rval = [ 0.5, 1, 1.5, 2];
L = 15;
D = 1; 
ni = 2; 
k2 = 1;
n = 1:20;                      
[alphaT_local, k_loc] = critvaluelocal2(2*L, D, k2, ni, n);
fprintf('alphaT local = %.4f \n', alphaT_local)
fprintf('wave number local = %.4f \n', k_loc)

for R = Rval
    [alphaT_nonlocal, k_nonloc] = critvaluenonlocal2(2*L, R, D, k2, ni, n);
    diff = alphaT_nonlocal - alphaT_local;
    sprintf(('alphaT non local for R = %.4f : %.4f \n'), R, alphaT_nonlocal)
    sprintf(('difference for R = %.4f: %.4f \n'), R, diff)
    sprintf(('wave number critico for R = %.4f: %.4f \n'), R, k_nonloc)
end


%% Difference between alpha_T non local and alpha_T local varying R (Full model)
Rval = [ 0.5, 1, 1.5, 2];
L = 15;
D1 = 1;
D2 = 1;
k1 = 2;
k2 = 1;
k3 = 1;
k4 = 2; 
n = 1:20;
[alphaT_local, k_loc] = critvaluelocal3(2*L, D1, D2, k1, k2, k3, k4, n);
fprintf('alphaT local = %.4f \n', alphaT_local)
fprintf('wave number local = %.4f \n', k_loc)
for R = Rval
    [alphaT_nonlocal, k_nonloc] = critvaluenonlocal3(2*L, R, D1, D2, k1, k2, k3, k4, n);
    diff = alphaT_nonlocal - alphaT_local;
    sprintf(('alphaT non local for R = %.4f : %.4f \n'), R, alphaT_nonlocal)
    sprintf(('difference for R = %.4f: %.4f \n'), R, diff)
    sprintf(('wave number critico for R = %.4f: %.4f \n'), R, k_nonloc)
end






