%% Import Data from xlsx table
filename = 'C:\Users\JimmyZhang\OneDrive - personalmicrosoftsoftware.uci.edu\School Work\PhD\PhD Research\Methamphetamine Study\Large Exp\Graphs Large Exp.xlsx';
T = xlsread(filename,'Matlab');
ControlT = T(1:7,:);
TwentyFiveT = T(11:17,:);
FourtyT = T(21:27,:);
FiftyT = T(31:37,:);
% Columns = Day, HR (BPM),QT (ms) QTc (ms), Ratio, RR (s), RT (s), RTc (s), ST (ms)

% Percentage    
for j = 2:8
    for i = 1:7
PerControlT(i,j) = ControlT(i,j)/ControlT(1,j);
PerTwentyFiveT(i,j) = TwentyFiveT(i,j)/TwentyFiveT(1,j);
PerFourtyT(i,j) = FourtyT(i,j)/FourtyT(1,j);
PerFiftyT(i,j) = FiftyT(i,j)/FiftyT(1,j);
    end
end

% HR plot
plot(ControlT(:,1),ControlT(:,2))
hold on
plot(ControlT(:,1),TwentyFiveT(:,2))
plot(ControlT(:,1),FourtyT(:,2))
plot(ControlT(:,1),FiftyT(:,2))
hold off
legend('Control','25 \muM', '40 \muM', '50 \muM', 'location', 'southeast');
title('Heart Rate vs Day of Methamphetamine Treatment');
xlabel('Day After First Methamphetamine Treatment');
ylabel('Heart Rate (Beats per Minute)');
figure(2)

% Percent HR Plot
plot(ControlT(:,1),PerControlT(:,2))
hold on
plot(ControlT(:,1),PerTwentyFiveT(:,2))
plot(ControlT(:,1),PerFourtyT(:,2))
plot(ControlT(:,1),PerFiftyT(:,2))
hold off
legend('Control','25 \muM', '40 \muM', '50 \muM', 'location', 'southeast');
title('Percentage Heart Rate vs Day of Methamphetamine Treatment');
xlabel('Day After First Methamphetamine Treatment');
ylabel('Percent Heart Rate (% of Baseline)');
figure(3)

% QTc plot
plot(ControlT(:,1),ControlT(:,4))
hold on
plot(ControlT(:,1),TwentyFiveT(:,4))
plot(ControlT(:,1),FourtyT(:,4))
plot(ControlT(:,1),FiftyT(:,4))
hold off
legend('Control','25 \muM', '40 \muM', '50 \muM', 'location', 'southeast');
title('QTc vs Day of Methamphetamine Treatment');
xlabel('Day After First Methamphetamine Treatment');
ylabel('QTc (milliseconds)');
figure(4)

% Percent QTc Plot
plot(ControlT(:,1),PerControlT(:,4))
hold on
plot(ControlT(:,1),PerTwentyFiveT(:,4))
plot(ControlT(:,1),PerFourtyT(:,4))
plot(ControlT(:,1),PerFiftyT(:,4))
hold off
legend('Control','25 \muM', '40 \muM', '50 \muM', 'location', 'southeast');
title('Percent QTc vs Day of Methamphetamine Treatment');
xlabel('Day After First Methamphetamine Treatment');
ylabel('Percent QTc (% of Baseline)');