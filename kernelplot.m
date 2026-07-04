% Kernel examples

sigma_val = [1, 1.5, 2, 2.5];
for sigma = sigma_val 
    laplace = @(x) sqrt(2)/(2*sigma).*exp(-sqrt(2).*abs(x) ./ sigma);
    
    gaussian = @(x) 1/(sqrt(2*pi)*sigma) .* exp(- x.^2 ./ (2*sigma^2));
    
    tophat = @(x) (abs(x) > sqrt(3)*sigma).*(0) + (abs(x) <= sqrt(3)*sigma).*(1/(2*sqrt(3)*sigma));
    
    x = -10:0.00001:10;
    figure()
    hold on
    grid on
    box on
    set(gcf, 'Color', 'w');
    set(gca, 'Color', 'w');
    set(gca, 'FontSize', 14);
    set(gca, 'XColor', 'k', 'YColor', 'k');
    set(gca, 'GridColor', 'k');
    set(gca, 'MinorGridColor', 'k');
    lap = plot(x,laplace(x), 'g', 'LineWidth', 1.8);
    gau = plot(x,gaussian(x),'LineWidth', 1.8);
    top = plot(x, tophat(x), 'b', 'LineWidth', 1.8);
    legend([lap, gau, top], 'Laplace', 'Gaussian', 'Top-hat');
    axis([-7 7 0 0.75])
    xlabel('x')
    ylabel('K(x)')
    
end

