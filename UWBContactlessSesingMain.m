clear all; 
close all; 
clc

%% Query user for logfile
[fnm,dnm] = uigetfile('*.csv');
current_dir = '';
C = strsplit(dnm, '/')
C{1,2}
for i=1:length(C)-2;
    current_dir = strcat(current_dir, C{1,i})
    current_dir = strcat(current_dir, '/')
end
fprintf('Reading logfile %s\n',fullfile(dnm,fnm));
[cfg,req,scn,det] = readMrmRetLog(fullfile(dnm,fnm));

 %% Pull out the raw scans (if saved)
scansI = find([scn.Nfilt] == 1);
scansV = reshape([scn(scansI).scn],[],length(scansI))';

%% Create the waterfall horizontal (also know as range or fast time) and vertical axes (also known as scan number or slow time)
Tbin = 32/(512*1.024);  % ns   (T2-T1)/N=Tbin  from the formulas in the biosensors paper
T0 = 0; % ns
c = 0.29979;  % m/ns
Rbin = c*(Tbin*(0:size(scansV,2)-1) - T0)/2;  % Range Bins in meters
Idata = 1:size(scansV,1);

%% Plot raw data as a waterfall
figure;
imagesc(Rbin,Idata,scansV);
xlabel('Range (m)')
ylabel('Scan Number')
title('Waterfall plot of raw scans')
colorbar
drawnow

%% bandpass filter
[b,a] = filt_coefs;
bpscans = filter(b,a,scansV,[],2);

figure; % Visualizing the raw bandpass data
imagesc(Rbin,Idata,bpscans);
xlabel('Range (m)')
ylabel('Scan Number')
title('Waterfall plot of bandpass filtered scans')
colorbar
drawnow

%% filtering to remove clutter
% We apply the following differece equation to remove static clutter
b = [1 -0.6 -0.3 -0.1]; 
a = 1;

figure; freqz(b,a); %FIR4 highpass 0.245Hz:cutoff freq, 8Hz=Fs 

filterscan = filter(b,a,bpscans,[],1); %filter each column % This we we got rid of static clutter and also preserve motion information (note we might have lost vital sign information though)

figure; % Visualizing the clutter removed data
imagesc(Rbin,Idata,filterscan);
xlabel('Range (m)')
ylabel('Scan Number')
title('Waterfall plot of bandpass filtered scans')
colorbar
drawnow

%% find envelope of the nonclutter rawscans

[b,a] = butter(6,0.4);  
%freqz(b,a) % This is a low pass filter.
filteredScansV = filter(b,a,abs(filterscan),[],2);% envelope of each scan or each row of the noclutterscansV
filteredScansV = max(filteredScansV,0);% now just ignoring anything less than zero

figure;
imagesc(Rbin,Idata,filteredScansV);
xlabel('Range (m)')
ylabel('Scan Number')
title('Waterfall plot of envelope of no clutter scan')
colorbar
drawnow

[maxes, max_indexes] = max(filteredScansV,[],2);
range = Rbin(max_indexes);
range = smoothdata(range,'movmedian',5);

figure;
plot(Idata, range);
ylabel('Range in m')
xlabel('Scan')
title('Estimated Range/Time')
drawnow

files = ["participant1_drink.csv", 1, "participant1_opendoor.csv", "participant1_sitdown.csv", "participant1_standup.csv", "participant1_walkslow.csv", ...
         "participant2_drink.csv", 1, "participant2_opendoor.csv", "participant2_sitdown.csv", "participant2_standup.csv", "participant2_walkslow.csv", ...
         "participant3_drink.csv", 1, "participant3_opendoor.csv", "participant3_sitdown.csv", "participant3_standup.csv", "participant3_walkslow.csv", ...
         "participant4_drink.csv", 1, "participant4_opendoor.csv", "participant4_sitdown.csv", "participant4_standup.csv", "participant4_walkslow.csv", ...
         "participant5_drink.csv", 1, "participant5_opendoor.csv", "participant5_sitdown.csv", "participant5_standup.csv", "participant5_walkslow.csv", ...
         "participant6_drink.csv", 1, "participant6_opendoor.csv", "participant6_sitdown.csv", "participant6_standup.csv", "participant6_walkslow.csv", ...
         "participant7_drink.csv", 1, "participant7_opendoor.csv", "participant7_sitdown.csv", "participant7_standup.csv", "participant7_walkslow.csv", ...
         "participant8_drink.csv", 1, "participant8_opendoor.csv", "participant8_sitdown.csv", "participant8_standup.csv", "participant8_walkslow.csv", ...
         "participant9_drink.csv", 1, "participant9_opendoor.csv", "participant9_sitdown.csv", "participant9_standup.csv", "participant9_walkslow.csv", ...
         "participant10_drink.csv", 1, "participant10_opendoor.csv", "participant10_sitdown.csv", "participant10_standup.csv", "participant10_walkslow.csv"];  
x_axis = [];
y_axis = [];

for file = files.'
    fileName = char(file(1))
    temp = strsplit(fileName, '_')
    class(fullfile(current_dir, temp{1,1} ,fileName))
    class(dnm)
    [cfg,req,scn,det] = readMrmRetLog(fullfile(current_dir, temp{1,1} ,fileName));
    scansI = find([scn.Nfilt] == 1);
    scansV = reshape([scn(scansI).scn],[],length(scansI))';
    Tbin = 32/(512*1.024);
    T0 = 0;
    c = 0.29979;
    Rbin = c*(Tbin*(0:size(scansV,2)-1) - T0)/2;
    Idata = 1:size(scansV,1);
    [b,a] = filt_coefs;
    bpscans = filter(b,a,scansV,[],2);
    b = [1 -0.6 -0.3 -0.1]; 
    a = 1;
    filterscan = filter(b,a,bpscans,[],1);
    [b,a] = butter(6,0.4);  
    filteredScansV = filter(b,a,abs(filterscan),[],2);
    filteredScansV = max(filteredScansV,0);
    [maxes, max_indexes] = max(filteredScansV,[],2);
    range = Rbin(max_indexes);
    features = [];
    labels = [];
    label = file(2);
    
    for i = 0:((length(range) - 25)/5)
        section_max = maxes((i * 5) + 1 : (i * 5) + 25);
        section_range = range((i * 5) + 1 : (i * 5) + 25);
        if max(section_max) > 1000
            features = [features; var(section_range), std(section_range), mean(gradient(section_range)), max(section_max), min(section_max), mean(section_max)];
            labels = [labels; label];
        end     
    end
    x_axis = [x_axis; features];
    y_axis = [y_axis; labels]; 
end

final = fitctree(x_axis, y_axis);

prediction = predict(final, x_axis);
cmmatrix = confusionmat(cellstr(y_axis),cellstr(prediction));