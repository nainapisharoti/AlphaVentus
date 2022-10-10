clear
clc
close all

turbine = 'M10';

filepath = fullfile('C:\Users\pisharoti1\Downloads\alpha_ventus_data-20220719T191212Z-001\alpha_ventus_data\',turbine)
cd(filepath)

filename = sprintf('%s_dataset_compiled.csv', turbine);
Data = readtable(filename);

%% Plotting Direction
id_angle_wind = Data.Wind_Direction>180;
id_angle_nacelle = Data.Nacelle_Position>180;
id_angle_rotor = Data.Rotor_Position>180;

for i = 1:size(Data.Wind_Direction,1)
if id_angle_wind(i) == 1
Wind_Direction(i,:) = Data.Wind_Direction(i)-360;
else
Wind_Direction(i,:) = Data.Wind_Direction(i);
end
end

for i = 1:size(Data.Nacelle_Position,1)
if id_angle_nacelle(i) == 1
Nacelle_Position(i,:) = Data.Nacelle_Position(i)-360;
else
Nacelle_Position(i,:) = Data.Nacelle_Position(i);
end
end

for i = 1:size(Data.Rotor_Position,1)
if id_angle_rotor(i) == 1
Rotor_Position(i,:) = Data.Rotor_Position(i)-360;
else
Rotor_Position(i,:) = Data.Rotor_Position(i);
end
end

%% Plotting for the existing window

%Desired Time Window
Time = Data.Time;
tstart = datetime(2010,05,16,00,00,00, 'Format','yyyy-MMM-dd HH:mm:ss Z');
tend = datetime(2010,05,16,06,00,00, 'Format','yyyy-MMM-dd HH:mm:ss Z');

id = tstart<=Time & tend>=Time;

%Plotting only "Quality Data"
id_quality = Data.Data_Quality==1 | Data.Data_Quality ==2;

%% Plots for the whole time frame
% Create figure Power-Vs-Time
figure1 = figure;
axes1 = axes('Parent',figure1);
plottools('on')
plot(Data.Time, Data.Active_Power)
hold on
plot(Data.Time, Data.Reactive_Power)
% Create ylabel
ylabel('Power (kW)');
% Create xlabel
xlabel('Time');
set(axes1,'FontSize',14,'FontWeight','bold');

% Create figure Property-Vs-Time
figure1 = figure;
axes1 = axes('Parent',figure1);
plottools('on')
plot(Data.Time, Data.Wind_Speed)
% Create ylabel
ylabel('Wind Speed (m/s)');
% Create xlabel
xlabel('Time');
set(axes1,'FontSize',14,'FontWeight','bold');

% Create figure Property-Vs-Time
figure1 = figure;
axes1 = axes('Parent',figure1);
plottools('on')
plot(Data.Time, Data.Generator_Speed)
hold on
plot(Data.Time, Data.Rotor_Speed)
% Create ylabel
ylabel('Speed (m/s)');
% Create xlabel
xlabel('Time');
set(axes1,'FontSize',14,'FontWeight','bold');

% Create figure Property-Vs-Time
figure1 = figure;
axes1 = axes('Parent',figure1);
plottools('on')
plot(Data.Time, Wind_Direction)
hold on
plot(Data.Time, Nacelle_Position)
plot(Data.Time, Rotor_Position)
% Create ylabel
ylabel('Angle (^o)');
% Create xlabel
xlabel('Time');
set(axes1,'FontSize',14,'FontWeight','bold')

%Create figure Property-Vs-Property
figure3 = figure;
axes3 = axes('Parent',figure3);
plottools('on')
plot(Data.Wind_Speed, Data.Active_Power)
% Create ylabel
ylabel('Power (kW)');
% Create xlabel
xlabel('Wind Speed (m/s)');
set(axes3,'FontSize',14,'FontWeight','bold');

%Create figure Property-Vs-Property
figure3 = figure;
axes3 = axes('Parent',figure3);
plottools('on')
scatter(Data.Generator_Speed, Data.Active_Power)
hold on
scatter(Data.Rotor_Speed, Data.Active_Power)
% Create ylabel
ylabel('Power (kW)');
% Create xlabel
xlabel('Speed (m/s)');
set(axes3,'FontSize',14,'FontWeight','bold');

%% Plots for Wills's simulated time frame

% Create figure Power-Vs-Time
figure1 = figure;
axes1 = axes('Parent',figure1);
plottools('on')
plot(Data.Time, Data.Active_Power)
hold on
plot(Data.Time, Data.Reactive_Power)
xlim([tstart tend])
% Create ylabel
ylabel('Power (kW)');
% Create xlabel
xlabel('Time');
set(axes1,'FontSize',14,'FontWeight','bold');

% Create figure Property-Vs-Time
figure1 = figure;
axes1 = axes('Parent',figure1);
plottools('on')
plot(Data.Time, Data.Wind_Speed)
xlim([tstart tend])
% Create ylabel
ylabel('Wind Speed (m/s)');
% Create xlabel
xlabel('Time');
set(axes1,'FontSize',14,'FontWeight','bold');

% Create figure Property-Vs-Time
figure1 = figure;
axes1 = axes('Parent',figure1);
plottools('on')
plot(Data.Time, Wind_Direction)
hold on
plot(Data.Time, Nacelle_Position)
plot(Data.Time, Rotor_Position)
xlim([tstart tend])
% Create ylabel
ylabel('Angle (^o)');
% Create xlabel
xlabel('Time');
set(axes1,'FontSize',14,'FontWeight','bold');

% Create figure Property-Vs-Time
figure1 = figure;
axes1 = axes('Parent',figure1);
plottools('on')
plot(Data.Time, Data.Generator_Speed)
hold on
plot(Data.Time, Data.Rotor_Speed)
% Create ylabel
ylabel('Speed (m/s)');
% Create xlabel
xlabel('Time');
set(axes1,'FontSize',14,'FontWeight','bold');

%Create figure Property-Vs-Property
figure3 = figure;
axes3 = axes('Parent',figure3);
plottools('on')
plot(Data.Wind_Speed(id), Data.Active_Power(id))
% Create ylabel
ylabel('Power (kW)');
% Create xlabel
xlabel('Wind Speed (m/s)');
set(axes3,'FontSize',14,'FontWeight','bold');

%Create figure Property-Vs-Property
figure3 = figure;
axes3 = axes('Parent',figure3);
plottools('on')
scatter(Data.Generator_Speed(id), Data.Active_Power(id))
hold on
scatter(Data.Rotor_Speed(id), Data.Active_Power(id))
% Create ylabel
ylabel('Power (kW)');
% Create xlabel
xlabel('Speed (m/s)');
set(axes3,'FontSize',14,'FontWeight','bold');
%% Plots based on Quality flags

clean_quality = Data.Data_Quality<10;

% Create figure Power-Vs-Time
figure1 = figure;
axes1 = axes('Parent',figure1);
plottools('on')
plot(Data.Time(id_quality), Data.Active_Power(id_quality))
hold on
plot(Data.Time(id_quality), Data.Reactive_Power(id_quality))
scatter(Data.Time(clean_quality), 100.*Data.Data_Quality(clean_quality))
% Create ylabel
ylabel('Power (kW)');
% Create xlabel
xlabel('Time');
set(axes1,'FontSize',14,'FontWeight','bold');

% Create figure Property-Vs-Time
% figure1 = figure;
% axes1 = axes('Parent',figure1);
figure(4)
hold on
plottools('on')
scatter(Time(id_quality), Wind_Direction(id_quality))
hold on
scatter(Time(id_quality), Nacelle_Position(id_quality))
scatter(Time(id_quality), Rotor_Position(id_quality))
% scatter(Data.Time(clean_quality), Data.Data_Quality(clean_quality))
% Create ylabel
ylabel('Angle (^o)');
% Create xlabel
xlabel('Time');
set(axes1,'FontSize',14,'FontWeight','bold');





