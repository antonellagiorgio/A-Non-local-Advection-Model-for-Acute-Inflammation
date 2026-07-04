function kernel = kernel_generation(R,L, S, dx)
    kernel  = zeros(S,S);
    for j = 1:S
        for i = 1:S-1
            dist = min(abs(i-j)*dx, 2*L-abs(i-j)*dx);
            if dist <= R 
                kernel(j,i) = 1/(2*R);
            end
        end
    end
end 