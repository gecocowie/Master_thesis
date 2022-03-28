%% Clear and Close Figure
clear all
close all
file_list = [
"B080817185719.txt"
"B100817185216.txt"
"B170817185244.txt"
"B330817185204.txt"
"B680817185148.txt"
"B830817185154.txt"
"B080817194711.txt"
"B100817194656.txt"
"B170817194728.txt"
"B330817194646.txt"
"B680817194630.txt"
"B830817194634.txt"
"B840817204053.txt"
"B080817203909.txt"
"B100817203625.txt"
"B170817203929.txt"
"B330817203612.txt"
"B680817203610.txt"
"B830817203600.txt"
];

for i=1:length(file_list)

%% Clear and Close Figure
%clear all
%close all

%addpath('C:\Users\Tuhtan\OneDrive - TTU\MATLAB\dependencies'); % dependency: sort_nat.m creates list of files in a specified folder.
%filePath = '/home/laura/Documents/TalTechRobotics/TestIMUGlacier'; % FOLDER WHERE B-SERIES FILES ARE LOCATED
%filePath = '/home/laura/Documents/TalTechRobotics/Drifter_Tube_SmallSensor12032020/Tube/B85/'; % FOLDER WHERE B-SERIES FILES ARE LOCATED
%filePath = '/home/laura/Documents/TalTechRobotics/Glaciers_Tracking/Data_Sets/Glaciers/Data/Supraglacial/Supraglacial_All0702/txt_files_run_new/'
filePath = "/Users/georgecowie/Documents/Master/Masteroppgave/data/2020/all_deployments_17082020/";
fileNameTxt = file_list{i}%'B800818151707.txt'; %r_B630630201247 r_B630630203630
fileFull = fullfile(filePath,fileNameTxt);

% contents = cellstr(ls(filePath));
% contents = sort_nat(contents(3:end,1)); % uses FEX SORT_NAT to get correct sorting otherwise 1,10,2,24,etc. will be used by MATLAB

% for itFile = 1:size(contents,1) % This iterates over all open B-series files, ending with *.txt
%     fileNameTxt = char(contents(itFile));
%     fileFull = strcat(filePath,fileNameTxt);
%     disp(['Importing and transforming file: ',fileNameTxt,' ...']);

fileNameNoExt = fileNameTxt(1:end-4); % remove classification ending

%% Mairo's code inserted here:
tubeVer = '1B'; % Options: 1A or 1B
sampleRate = 100; % NDOF mode:250 or ACC+GYRO mode:100

% Definitions
timeDiffMax = 15;
Pmin = 950;         % Pressure MIN limit
Pmax = 1400;        % Pressure MAX limit
Tmin = 18;          % Temperature MIN limit
Tmax = 35;          % Temperature MAX limit
Calmin = 0;
Calmax = 3;

%% Set packet sieze 
if sampleRate == 100
    packetSize = 100;
else
    packetSize = 60;
end

% Load file
fid = fopen(fileFull,'r','b');
% Calculate file size
fstat1 = dir(fileFull);
flen = floor((fstat1.bytes/packetSize))-1;

len = flen; % Number of lines to extract from binary file
if strcmp(tubeVer, '1B')
    if sampleRate == 100
        A = [len, 25];
    else
        A = [len, 18];
    end
else
    A = [len, 27];
end
   
for i=1:len
    j = 0;  
    %percent = floor((i/flen)*100);
    
    if mod(i,1000) == 0
        fprintf('%d of %d\n', i, flen);
    end
    
    fread(fid, 1, '*uint8'); % Check digits, not to be included in data matrix!
    fread(fid, 1, '*uint8');
    %fread(fid, 1, '*uint8')
    
    if strcmp(tubeVer, '1B')
        j=j+1;
        A(i,j) = fread(fid, 1, '*uint16', 'ieee-le');  % sample rate
        % disp(num2str(A(i,j))); % display sample rate
    end
    
    j=j+1;    A(i,j) = fread(fid, 1, '*uint32', 'ieee-le');   % timestamp
    j=j+1;    A(i,j) = fread(fid, 1, '*float', 'ieee-le');   % P1
    j=j+1;    A(i,j) = fread(fid, 1, '*float', 'ieee-le');   % T1
    j=j+1;    A(i,j) = fread(fid, 1, '*float', 'ieee-le');   % P2
    j=j+1;    A(i,j) = fread(fid, 1, '*float', 'ieee-le');   % T2
    j=j+1;    A(i,j) = fread(fid, 1, '*float', 'ieee-le');   % P3
    j=j+1;    A(i,j) = fread(fid, 1, '*float', 'ieee-le');   % T3
    if sampleRate == 100 
        j=j+1;    A(i,j) = fread(fid, 1, '*float', 'ieee-le');   % eul head
        j=j+1;    A(i,j) = fread(fid, 1, '*float', 'ieee-le');   % eul roll
        j=j+1;    A(i,j) = fread(fid, 1, '*float', 'ieee-le');   % eul pitch
        j=j+1;    A(i,j) = fread(fid, 1, '*float', 'ieee-le');   % quat w
        j=j+1;    A(i,j) = fread(fid, 1, '*float', 'ieee-le');   % quat x
        j=j+1;    A(i,j) = fread(fid, 1, '*float', 'ieee-le');   % quat y
        j=j+1;    A(i,j) = fread(fid, 1, '*float', 'ieee-le');   % quat z
        j=j+1;    A(i,j) = fread(fid, 1, '*float', 'ieee-le');   % mag x
        j=j+1;    A(i,j) = fread(fid, 1, '*float', 'ieee-le');   % mag y
        j=j+1;    A(i,j) = fread(fid, 1, '*float', 'ieee-le');   % mag z
    end
    j=j+1;    A(i,j) = fread(fid, 1, '*float', 'ieee-le');   % acc x
    j=j+1;    A(i,j) = fread(fid, 1, '*float', 'ieee-le');   % acc y
    j=j+1;    A(i,j) = fread(fid, 1, '*float', 'ieee-le');   % acc z
    j=j+1;    A(i,j) = fread(fid, 1, '*float', 'ieee-le');   % gyro x
    j=j+1;    A(i,j) = fread(fid, 1, '*float', 'ieee-le');   % gyro y
    j=j+1;    A(i,j) = fread(fid, 1, '*float', 'ieee-le');   % gyro z
    j=j+1;    A(i,j) = fread(fid, 1, '*uint8');   % cal mag
    j=j+1;    A(i,j) = fread(fid, 1, '*uint8');   % cal acc
    j=j+1;    A(i,j) = fread(fid, 1, '*uint8');   % cal gyro
    j=j+1;    A(i,j) = fread(fid, 1, '*uint8');   % cal imu
end
fclose(fid);
j = 1; % restart column position counter

fs = A(:,j); j=j+1;         % 1
timestamp = A(:,j); j=j+1;  % 2
P1 = A(:,j); j=j+1;         % 3
T1 = A(:,j); j=j+1;         % 4
P2 = A(:,j); j=j+1;         % 5
T2 = A(:,j); j=j+1;         % 6
P3 = A(:,j); j=j+1;         % 7
T3 = A(:,j); j=j+1;         % 8
if sampleRate == 100 
    eul_head = A(:,j); j=j+1;   % 9
    eul_roll = A(:,j); j=j+1;   % 10
    eul_pitch = A(:,j); j=j+1;  % 11
    quat_w = A(:,j); j=j+1;     % 12
    quat_x = A(:,j); j=j+1;     % 13
    quat_y = A(:,j); j=j+1;     % 14
    quat_z = A(:,j); j=j+1;     % 15
    mag_x = A(:,j); j=j+1;      % 16
    mag_y = A(:,j); j=j+1;      % 17
    mag_z = A(:,j); j=j+1;      % 18
end
acc_x = A(:,j); j=j+1;      % 9/22
acc_y = A(:,j); j=j+1;      % 10/23
acc_z = A(:,j); j=j+1;      % 11/24
gyro_x = A(:,j); j=j+1;     % 12/19
gyro_y = A(:,j); j=j+1;     % 13/20
gyro_z = A(:,j); j=j+1;     % 14/21
cal_mag = A(:,j); j=j+1;    % 15/25
cal_acc = A(:,j); j=j+1;    % 16/26
cal_gyro = A(:,j); j=j+1;   % 17/27
cal_imu = A(:,j); j=j+1;    % 18/28

% Analyze data
for i=1:len
    if (i>1) % Skip first line
       timediff(i) = timestamp(i)-timestamp(i-1);
       if timediff(i) > timeDiffMax
           fprintf('Time diff too large at line %d: %d \n', i, timediff(i));
       end
    end
   if (i==2)
       timediff(i-1) = timediff(i);
   end
end
%%

tm = (timestamp - timestamp(1)) ./ (60*1000);
ts = tm .* 60; % time vector in seconds

%% Export data as csv including header

if sampleRate == 100
    dataExport(:,1) = ts; % set first time stamp to zero
    dataExport(:,2) = P1; % pressure data including atmoshperic correction
    dataExport(:,3) = T1;
    dataExport(:,4) = P2;
    dataExport(:,5) = T2;
    dataExport(:,6) = P3;
    dataExport(:,7) = T3;
    dataExport(:,8) = eul_head;
    dataExport(:,9) = eul_roll;
    dataExport(:,10) = eul_pitch;
    dataExport(:,11) = quat_w;
    dataExport(:,12) = quat_x;
    dataExport(:,13) = quat_y;
    dataExport(:,14) = quat_z;
    dataExport(:,15) = mag_x;
    dataExport(:,16) = mag_y;
    dataExport(:,17) = mag_z;
    dataExport(:,18) = acc_x;
    dataExport(:,19) = acc_y;
    dataExport(:,20) = acc_z;
    dataExport(:,21) = gyro_x;
    dataExport(:,22) = gyro_y;
    dataExport(:,23) = gyro_z;
    dataExport(:,24) = cal_mag;
    dataExport(:,25) = cal_acc;
    dataExport(:,26) = cal_gyro;
    dataExport(:,27) = cal_imu;

else if sampleRate == 250
    dataExport(:,1) = ts; % set first time stamp to zero
    dataExport(:,2) = P1; % pressure data including atmoshperic correction
    dataExport(:,3) = T1;
    dataExport(:,4) = P2;
    dataExport(:,5) = T2;
    dataExport(:,6) = P3;
    dataExport(:,7) = T3;
    dataExport(:,8) = acc_x;
    dataExport(:,9) = acc_y;
    dataExport(:,10) = acc_z;
    dataExport(:,11) = gyro_x;
    dataExport(:,12) = gyro_y;
    dataExport(:,13) = gyro_z;
    dataExport(:,14) = cal_mag;
    dataExport(:,15) = cal_acc;
    dataExport(:,16) = cal_gyro;
    dataExport(:,17) = cal_imu;
    end
end
        
dataExport = dataExport(1:end-201,:);

if sampleRate == 100
acc_mag = sqrt(acc_x.^2 + acc_y.^2 + acc_z.^2);
mag_mag = sqrt(mag_x.^2 + mag_y.^2 + mag_z.^2);
elseif sampleRate == 250 % if 250 Hz  
else
end

% calculated normalized magnitude data
acc_magMax = max(acc_mag);
acc_magMin = min(acc_mag);
acc_magNorm = (acc_mag - acc_magMin)./(acc_magMax - acc_magMin);

%% Convert body to Earth reference frame
quat = quaternion([dataExport(:,11),dataExport(:,12),dataExport(:,13),dataExport(:,14)]);
accBody = [dataExport(:,18),dataExport(:,19),dataExport(:,20)]; 
%accEarth = rotatepoint(quat,accBody);

%dataExport(:,18) = accEarth(:,1);
%dataExport(:,19) = accEarth(:,2);
%dataExport(:,20) = accEarth(:,3);
%%

%% STEP 4: Export as csv text file for reporting / sharing raw data
if sampleRate == 100
cHeader = {'Time [s]','PL [hPa]','TL [C]','PC [hPa]','TC [C]','PR [hPa]','TR [C]','EX [deg]','EY [deg]','EZ [deg]','QW [-]','QX [-]','QY [-]','QZ [-]','MX [microT]','MY [microT]','MZ [microT]','AXEarth [m/s2]','AYEarth [m/s2]','AZEarth [m/s2]','RX [rad/s]','RY [rad/s]','RZ [rad/s]','CSM','CSA','CSR','CSTOT'}; %dummy header
elseif sampleRate == 250 % if 250 Hz   
cHeader = {'Time [s]','PL [hPa]','TL [C]','PC [hPa]','TC [C]','PR [hPa]','TR [C]','AXEarth [m/s2]','AYEarth [m/s2]','AZEarth [m/s2]','RX [rad/s]','RY [rad/s]','RZ [rad/s]','CSM','CSA','CSR','CSTOT'}; %dummy header    
else
end

commaHeader = [cHeader;repmat({','},1,numel(cHeader))]; %insert commaas
commaHeader = commaHeader(:)';
textHeader = cell2mat(commaHeader); %cHeader in text with commas
textHeader = textHeader(1:end-1);

%write header to file
exportFile = strcat(filePath,'\',fileNameNoExt,'.csv');
fid = fopen(exportFile,'wt'); 
fprintf(fid,'%s\n',textHeader)
fclose(fid)

%write data to end of file
dlmwrite(exportFile,dataExport,'-append');
%%

%% STEP 5: Export .mat file
%exportMatFile = strcat(filePath,'\',fileNameNoExt,'.mat');
%if sampleRate == 100 % if 100 Hz
%save(exportMatFile,'ts','P1','T1','P2','T2','P3','T3','eul_head','eul_roll','eul_pitch','quat_x','quat_y','quat_z','quat_w','mag_x','mag_y','mag_z','acc_x','acc_y','acc_z','gyro_x','gyro_y','gyro_z','cal_mag','cal_acc','cal_gyro','cal_imu');
%elseif sampleRate == 250 % if 250 Hz
%save(exportMatFile,'ts','P1','T1','P2','T2','P3','T3','acc_x','acc_y','acc_z','gyro_x','gyro_y','gyro_z','cal_mag','cal_acc','cal_gyro','cal_imu');    
%else
%end
%%

%% Figures
gridStatus = 'on';
figFontSize = 10; % axis text font size
figLineWidth = 1;
fig1 = figure('units','normalized','outerposition',[0 0 1 1]);
left_color = [0 0 0];
right_color = [1 0 1];
set(fig1,'defaultAxesColorOrder',[left_color; right_color]);

yyaxis right
%subplot(3,1,1)
% plot(ts,sqrt(acc_x.^2 + acc_y.^2 + acc_z.^2) ,'m','LineWidth',1);
% ylabel('Acceleration Magnitude (m/s�)');
% yyaxis left

%plot(ts,P1,'-r','LineWidth',figLineWidth);
%ax = gca;
%hold on
%plot(ts,P2,'-b','LineWidth',figLineWidth);
%plot(ts,P3,'-k','LineWidth',figLineWidth);
%legend('Left','Center','Right');
%ylabel('Total Pressure (mbar)')
%xlabel('Time (s)');
% xlim([15 60]);
% ylim([900 1100]);
%ax.XGrid = gridStatus;
%ax.YGrid = gridStatus;
%set(gca,'fontsize',figFontSize)
%title(fileNameTxt);

% subplot(3,1,2)
% plot(ts,eul_head,'k','LineWidth',figLineWidth);
% ax = gca;
% hold on
% plot(ts,eul_roll,'r','LineWidth',figLineWidth);
% plot(ts,eul_pitch,'b','LineWidth',figLineWidth);
% legend('Heading','Roll','Pitch');
% ylabel('Euler angle (deg)')
% xlabel('Time (s)');
% ax.XGrid = gridStatus;
% ax.YGrid = gridStatus;
% % xlim([15 60]);
% set(gca,'fontsize',figFontSize)

% subplot(3,1,2)
% plot(ts,(eul_pitch),'b','LineWidth',figLineWidth);
% ax = gca;
% legend('Pitch');
% ylabel('Euler angle (deg)')
% xlabel('Time (s)');
% ax.XGrid = gridStatus;
% ax.YGrid = gridStatus;
% % xlim([15 60]);
% set(gca,'fontsize',figFontSize)

%subplot(3,1,2)
%plot(ts,(mag_x),'k','LineWidth',figLineWidth);
%ax = gca;
%hold on
%plot(ts,(mag_y),'r','LineWidth',figLineWidth);
%plot(ts,(mag_z),'b','LineWidth',figLineWidth);
% legend('Pitch');
%ylabel('Magnetometer (micro T)')
%xlabel('Time (s)');
%ax.XGrid = gridStatus;
%ax.YGrid = gridStatus;
% ylim([-80 80]);
% xlim([15 60]);
%set(gca,'fontsize',figFontSize)

% subplot(3,1,3)
% plot(ts,gyro_x .* 180/pi,'k','LineWidth',3);
% hold on
% plot(ts,gyro_y .* 180/pi,'r','LineWidth',2);
% plot(ts,gyro_z .* 180/pi,'b','LineWidth',1);
% legend('Heading','Roll','Pitch');
% ylabel('Rate gyro (deg/s)')
% xlabel('Time (s)');
% % xlim([15 60]);
% set(gca,'fontsize',figFontSize)

% Plot linear accel, minus gravity vector (100 Hz fusion mode only!)
% subplot(3,1,3)
% plot(ts,acc_x,'k','LineWidth',figLineWidth);
% ax = gca;
% hold on
% plot(ts,acc_y,'r','LineWidth',figLineWidth);
% plot(ts,acc_z,'b','LineWidth',figLineWidth);
% legend('Acc X','Acc Y','Acc Z');
% ylabel('Linear acceleration (m/s^2)')
% xlabel('Time (s)');
% ax.XGrid = gridStatus;
% ax.YGrid = gridStatus;
% % xlim([15 60]);
% set(gca,'fontsize',figFontSize)

% Plot linear accel MAGNITUDE
%subplot(3,1,3)
%plot(ts,acc_mag,'r','LineWidth',figLineWidth);
%ax = gca;
%hold on
% plot(ts,abs(eul_pitch),'b','LineWidth',figLineWidth);
%legend('Acc XYZ');
%ylabel('Accel Mag (m/s^2)')
%xlabel('Time (s)');
%ax.XGrid = gridStatus;
%ax.YGrid = gridStatus;
% xlim([15 60]);
%set(gca,'fontsize',figFontSize)

% subplot(3,1,3)
% tubespectro % separate script for plotting spectrogram

%% Export current figure as image
%figName = erase(fileNameTxt,'.txt');
%figureExport = strcat(filePath,figName,'.pdf');
% saveas(gcf,figureExport);
% print(figureExport,'-dpdf','-fillpage')
% orient(fig1,'landscape')
%print(fig1,figureExport,'-dpdf','-bestfit')

%pause(1)

%% Create output for dead reckoning
deadReckCSV = 0;
if deadReckCSV == 1
g = 9.818559; % https://www.ptb.de/cartoweb3/SISproject.php

out = ((1:size(acc_x)))';
out(:,2) = gyro_x;
out(:,3) = gyro_y;
out(:,4) = gyro_z;
out(:,5) = acc_x ./ g;
out(:,6) = acc_y ./ g;
out(:,7) = acc_z ./ g;
out(:,8) = mag_x .* 0.01; % convert from microTesla to Gauss
out(:,9) = mag_x .* 0.01; % convert from microTesla to Gauss
out(:,10) = mag_x .* 0.01; % convert from microTesla to Gauss
csvwrite('E:\MATLAB\Dead_Reckoning\Gait-Tracking-With-x-IMU-master\Gait-Tracking-With-x-IMU-master\Gait Tracking With x-IMU\Datasets\test.csv',out);
else 
end
%%

close all
clear dataExport
%% Clear and Close Figure
close all

% end % loop over all files in current folder

fclose('all');
end