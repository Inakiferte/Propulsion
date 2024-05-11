function plotPerformanceV(X,F_s_m,I_sp_m,eta_p_m,eta_th_m,eta_o_m)
% Plotting parameters
fs = 25;

% Plotting the values
figure;

subplot(1,2,1)
plot(X, F_s_m * 1e-3, 'Color', 'r', 'LineWidth', 2);

% Adding labels and title
xlabel('$T_{Fuel}^{In}$ [C]', 'Interpreter', 'latex', 'FontSize', fs)
ylabel('$F_{s}$ [kN]', 'Interpreter', 'latex', 'FontSize', fs);
xlim([X(1) X(end)]);
set(gca, 'ticklabelinterpreter', 'latex','FontSize', fs);

text(0.1, 0.98, 'a)', 'Interpreter', 'latex', 'Units', 'normalized', 'FontSize', fs, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');

% Display the grid
grid on;

% Plotting the values
subplot(1,2,2)
plot(X, I_sp_m, 'r', 'Color', 'b', 'LineWidth', 2);

% Adding labels and title
xlabel('$T_{Fuel}^{In}$ [C]', 'Interpreter', 'latex', 'FontSize', fs);
ylabel('$I_{sp}$ [s]', 'Interpreter', 'latex', 'FontSize', fs);
xlim([X(1) X(end)]);
set(gca, 'ticklabelinterpreter', 'latex','FontSize', fs);
text(0.1, 0.98, 'b)', 'Interpreter', 'latex', 'Units', 'normalized', 'FontSize', fs, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');

% Display the grid
grid on;

figure;

% Plotting the values
subplot(1,3,1)

plot(X, eta_p_m * 100, 'r', 'Color', 'r', 'LineWidth', 2);


% Adding labels and title
xlabel('$T_{Fuel}^{In}$ [C]', 'Interpreter', 'latex', 'FontSize', fs);
ylabel('$\eta_{p}$ [$\%$]', 'Interpreter', 'latex', 'FontSize', fs);
xlim([X(1) X(end)]);
set(gca, 'ticklabelinterpreter', 'latex','FontSize', fs);
text(0.2, 0.98, 'a)', 'Interpreter', 'latex', 'Units', 'normalized', 'FontSize', fs, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');

% Display the grid
grid on;

% Plotting the values
subplot(1,3,2)

plot(X, eta_th_m * 100, 'r', 'Color', 'b', 'LineWidth', 2);


% Adding labels and title
xlabel('$T_{Fuel}^{In}$ [C]', 'Interpreter', 'latex', 'FontSize', fs);
ylabel('$\eta_{th}$ [$\%$]', 'Interpreter', 'latex', 'FontSize', fs);
xlim([X(1) X(end)]);
set(gca, 'ticklabelinterpreter', 'latex','FontSize', fs);
text(0.2, 0.98, 'b)', 'Interpreter', 'latex', 'Units', 'normalized', 'FontSize', fs, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');

% Display the grid
grid on;

% Plotting the values
subplot(1,3,3)

plot(X, eta_o_m * 100, 'r', 'Color', 'g', 'LineWidth', 2);

% Adding labels and title
xlabel('$T_{Fuel}^{In}$ [C]', 'Interpreter', 'latex', 'FontSize', fs);
ylabel('$\eta_{o}$ [$\%$]', 'Interpreter', 'latex', 'FontSize', fs);
xlim([X(1) X(end)]);
set(gca, 'ticklabelinterpreter', 'latex','FontSize', fs);
text(0.2, 0.98, 'c)', 'Interpreter', 'latex', 'Units', 'normalized', 'FontSize', fs, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');

% Display the grid
grid on;
end