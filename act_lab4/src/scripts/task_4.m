% plant parameters
t = linspace(-2*pi, 2*pi, 1000);
ampl = 1;
m = 3;

g_approx = zeros(size(t));

for k = 1:m
    n = 2*k - 1;
    bn = (4*ampl)/(pi*n);
    g_approx = g_approx + bn * sin(n * t);
end

% true square wave
g_true = ampl * sign(sin(t));

% plot
plot(t, g_true, 'k--', 'LineWidth', 1.2); hold on;
plot(t, g_approx, 'r', 'LineWidth', 2);
legend('g', 'fourier g (m=3)');
xlabel('t'); ylabel('g(t)');
title('Fourier approximation');
grid on;
