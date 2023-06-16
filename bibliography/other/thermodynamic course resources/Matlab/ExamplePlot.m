% Script that plots data with squares

% using a matrix makes it easy to enter tabular data
% first column is x, second y, third Pressure
% using tabs can make it easy to enter data
% set tab length in 'File>Preferences'
xyPData = [
    0   0    0.2939;
    0.1 0.3 0.3616;
    0.5 0.75 0.475;
    1   1    0.55;];

% enter 'help plot' in the command window for a quick overview of codes for lines and symbols.
% common marker options, 'o', 's', '^', 'v'; see 'help plot'
plot(xyPData(:,1),xyPData(:,3),'-s')
hold on % hold to add a second data set
plot(xyPData(:,2),xyPData(:,3),'--ro')
hold off
xlabel('x1 values')
ylabel('y1 values')
% put legend in the best location
legend('first set','second set', 'Location', 'Best')
