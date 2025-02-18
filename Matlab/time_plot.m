% Load the data from the saved .mat file (if not already in the workspace)
% Uncomment the next line if you are loading the saved data
load('imu_data.mat', 'data');

% Plot gravity (g) vs time
figure;
subplot(2,1,1);
plot(data.time, data.g(:,1), 'r', 'LineWidth', 1.5); hold on;
plot(data.time, data.g(:,2), 'g', 'LineWidth', 1.5);
plot(data.time, data.g(:,3), 'b', 'LineWidth', 1.5);
title('Gravity Vector Components (g) vs Time');
xlabel('Time (s)');
ylabel('g (m/s^2)');
legend('g[0]', 'g[1]', 'g[2]');
grid on;

% Plot angular velocity (dps) vs time
subplot(2,1,2);
plot(data.time, data.dps(:,1), 'r', 'LineWidth', 1.5); hold on;
plot(data.time, data.dps(:,2), 'g', 'LineWidth', 1.5);
plot(data.time, data.dps(:,3), 'b', 'LineWidth', 1.5);
title('Gyroscope Angular Velocity Components (dps) vs Time');
xlabel('Time (s)');
ylabel('dps (deg/s)');
legend('dps[0]', 'dps[1]', 'dps[2]');
grid on;
