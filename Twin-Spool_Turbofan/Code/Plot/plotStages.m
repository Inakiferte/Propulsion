function [] = plotStages(Stages,po,To)
% Plotting parameters
fs = 20;

% Plotting the values
figure;

plot(Stages, po * 1e-3, 's--', 'MarkerSize', 8, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'b', 'Color', 'b', 'LineWidth', 2)

% Adding labels and title
xlabel('Turbofan Stages', 'Interpreter', 'latex', 'FontSize', fs);
ylabel('Po [kPa]', 'Interpreter', 'latex', 'FontSize', fs);
title('Stagnation Pressure Through Stages', 'Interpreter', 'latex', 'FontSize', fs);
set(gca, 'ticklabelinterpreter', 'latex');

% Display the grid
grid on;

% Plotting the values
figure;

plot(Stages, To, 's--', 'MarkerSize', 8, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'Color', 'r', 'LineWidth', 2)

% Adding labels and title
xlabel('Turbofan Stages', 'Interpreter', 'latex', 'FontSize', fs);
ylabel('To [K]', 'Interpreter', 'latex', 'FontSize', fs);
title('Stagnation Temperature Through Stages', 'Interpreter', 'latex', 'FontSize', fs);
set(gca, 'ticklabelinterpreter', 'latex');

% Display the grid
grid on;

end