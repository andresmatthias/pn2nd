function henyeyGreensteinCheck()
    clear
    close all
    addpath('../PN2nd/basisFunctions/')
    g = 0.3;
    HG = @(theta) 1 / 4 / pi * (1 - g^2) ./ ((1 + g^2 - 2 * g * cos(theta)).^(3/2)); 
    f = @(theta) HG(theta);
    theta = linspace(0, pi, 500);
    
    plot([0, 0], [0, 0.1], 'k', 'LineWidth', 3)
    hold on
    plot([0, -0.01], [0.1, 0.09], 'k', 'LineWidth', 3)
    plot([0, 0.01], [0.1, 0.09], 'k', 'LineWidth', 3)
    axis equal

    
    plotPolar(theta, f(theta), 'k');
    for n = 5
        gain = HGApprox(n, g, theta);
        plotPolar(theta, gain, 'g');
    end
    figure
    plot(theta, f(theta), 'k')
    hold on
    plot(theta, gain)
end

function plotPolar(theta, gain, varargin)
    x = sin(theta) .* gain;
    y = cos(theta) .* gain;
    plot(x,y, varargin{:})
end

function y = HGApprox(n, g, theta)
% https://www.sciencedirect.com/science/article/pii/0022407388900787
    y = theta * 0;
    for k = 0 : n
       y = y + (2 * k + 1) * g^k * legendreP(k, cos(theta));  
    end
    y = y / 4 / pi;
end