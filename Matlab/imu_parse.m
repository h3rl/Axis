% This file is used to parse serial data from the micro controller into a
% format that is suitable for matlab
com_port = "COM3";
baud_rate = 115200;
tag_name = "IMU_Serial";

% Check if the serial port is already open
s = serialportfind("Tag", tag_name);
if isempty(s)
    s = serialport(com_port, baud_rate);
    s.Tag = tag_name;
    configureTerminator(s, "LF");
    fprintf("Opened %s successfully with tag '%s'.\n", com_port, tag_name);
else
    fprintf("Serial port %s is already open with tag '%s'.\n", com_port, tag_name);
end

s.flush()

% Initialize the structure to store the data
data = struct('time', [], 'g', [], 'dps', []);


% Read data in an infinite loop
while true
    if s.NumBytesAvailable > 0
        % Read the line from the serial port
        ln = readline(s);
        % Display the raw data (for debugging)
        fprintf("%s", ln);
        
        % Parse the line into individual components
        values = sscanf(ln, '%f %f %f %f %f %f %f');
        
        % Ensure that the line contains exactly 7 values: [time, g[0], g[1], g[2], dps[0], dps[1], dps[2]]
        if numel(values) == 7
            % Create a new entry in the structure for the current data
            data.time(end+1) = values(1);
            data.g(end+1, :) = values(2:4);  % g[0], g[1], g[2]
            data.dps(end+1, :) = values(5:7);  % dps[0], dps[1], dps[2]
        end
    end
end
% For saving the data
% shift so time starts at 0
data.time = data.time - data.time(1);
% convert ms to seconds
data.time = data.time / 1000;
% save the data
save('imu_data.mat', 'data');