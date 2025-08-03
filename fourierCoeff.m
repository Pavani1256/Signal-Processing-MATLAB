
T = 1;
N = 5;
t1 = 0;
t2 = T;
syms t;
xt = 2*cos(2*pi*t) + cos(6*pi*t);
time_grid = -0.5:0.01:0.5;


function F = compfourierCoeff(t, xt, T, t1, t2, N)
    F = sym(zeros(2*N+1, 1));  % Preallocate symbolic array
    k_vals = -N:N;
    for idx = 1:length(k_vals)
        k = k_vals(idx);
        F(idx) = (1/T) * int(xt * exp(-1i * 2 * pi * k * t / T), t, t1, t2);
    end
end

% Call the coefficient function
F = compfourierCoeff(t, xt, T, t1, t2, N);

% Function to compute partial sum
function y = partialsum(F, T, time_grid, N)
    y = zeros(size(time_grid));  % Fix: use zeros instead of zero
    k_vals = -N:N;
    for idx = 1:length(k_vals)
        k = k_vals(idx);
        y = y + double(F(idx)) * exp(1i * 2 * pi * k * time_grid / T);
    end
end


y_approx = partialsum(F, T, time_grid, N);


FS_idx = -N:N;
figure;
stem(FS_idx, double(F)); grid on;
title('Fourier Coefficients');
xlabel('Harmonic Index'); ylabel('Coefficient');

% Plot partial sum approximation
figure;
plot(time_grid, real(y_approx), 'LineWidth', 1.5); grid on;
title('Partial Sum Approximation of xt');
xlabel('Time'); ylabel('xt reconstructed');