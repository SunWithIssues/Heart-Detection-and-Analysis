% Title: HeartRateAnalysis2.m
% Undergraduates: Alvin, Anya
% Graduates: Akwasi, Pascal
% Professor guidance by: Dr. Deo, Dr. Albin
% Purpose:
%    This should be able to recieve an ECG data, filter it, and 
%      determine the QRS Peaks. It returns the estimated beats per minute
%      and if the heart rate stays the same throughout the sample (i.e 
%      if the heart has an arrhythmia issue).
% Notes/Problems to Fix:
%    1. Try to implement a Bandpass without the explicit butterworth
%       parameter. This means use a high pass after the low pass data.
%    2. Need to program an adaptive threshold, for detection purposes.
%         Plan: FIXED detec -> Weighted (top/bottom) detec -> Adaptive
% Method:
%    The program goes through all CSV data inside the directory. 
% Outputs:
%    Length of file. Beats per Minutes.
%    Possibility of arrhythmia

clear all
close all
cd '/Users/Looks/Documents/MATLAB/New_Data/Jailyn/' % The directory of CSV data
dr = datastore(pwd);
for i = 1:length(dr.Files) % Goes through all files inside the directory

    dloc = sprintf('d%d', i); 
    Data.(dloc) = importdata(char(dr.Files(i)));
    t = Data.(dloc).data(:,1);    % time values
    
    % We assume that there are two channels, and the amplitude of ECG is in
    % one of them. Which one, is not defined, before processing.
    ch1 = Data.(dloc).data(:,2);  % channel 1 values
    ch2 = Data.(dloc).data(:,3);  % channel 2 values
    
    % =========== Checks for channel of interest ================ %
    % The channel used has both negative and positive amplitudes,
    % therefore, the range is greater for the channel which the data
    % was obtained
    % =========================================================== %
    if range(ch2) > range(ch1) %channel 2 should be observed
        Amp = ch2(:);
    else                       %channel 1 should be observed
        Amp = ch1(:);
    end
    
    % =================== Plot, Raw Data vs Time ================ %
    figure
    plot(t,Amp);
    xlabel('time (s)')
    ylabel('Voltage(V)')
    title('sensor voltage vrs time')
    
    
    % ============== Single-Sided Spectrum, Raw Data ============ %
    % Uses Fourier Transform
    %    Things to Remember:
    %       * fourier transform -> changes time domain to frequency domain
    %       * inverse fourier transform -> changes freq to time
    % Honestly, I have very little understanding of these next lines of
    %    code (that is, lines 64-71). These 8 lines were written by Akwasi, 
    %    they can be seen elsewhere in the project, with not much change.
    % =========================================================== %
    Fs = 20;          % sampling frequency
    T = 1/Fs;         % sampling period
    L = length(t);    % length of data
    Y = fft(Amp);     % returns complex form
    P2 = abs(Y/L);    % This is the magnitude
    P1 = P2(1:L/2+1); % rescale axis due to double sideband
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(L/2))/L;
    
    % Plots the Raw, Single-Sided Amplitude Spectrum
    figure
   
    subplot(4,2,[5,6]);
    plot(f,P1);
    title('Single-Sided Amplitude Spectrum of X(t)')
    xlabel('f (Hz)');
    ylabel('|P1(f)|');

    % ================= Testing Low Pass filter ================= %
    % A low pass filter gets rid of higher frequencies. This uses 
    % the butterworth function as well as the zerophase filter.
    %    Things to look into: 
    %       1. Cutoff frequency -> Wn
    % =========================================================== %
    Wn = 0.05;                    % cutoff frequency
    N = 3;                        % 3rd order
    [a,b] = butter(N, Wn, 'Low'); % a,b represent the coefficients
    low_f = filtfilt(a,b,Amp);    % this is the result for the low pass filter
    
    % Plot, Low Pass Filter
    subplot(4,2,[1,2]);
    plot(low_f);
    title('Low Pass Filter');
    
    % Plot, Raw Data
    subplot(4,2,[3,4]);
    plot(Amp);
    title('Raw');
    
    
    % ================ Transfering filter data ==================== %
    % This stub comments were used for transfering graphs and analyis
    % to subgroup 3 (i.e. android application portion).
    % ============================================================= %
    % Creates a table of the data processed.
    % singleSet = table(t, low_f);
    % name = "data" + dloc + ".csv";
    % writetable(singleSet, name);
    
    
    % ============ Testing Low Pass, Single-Sided Spectrum ========= %
    Y = fft(low_f);     %returns complex form
    P2 = abs(Y/L);      %this is the magnitude
    P1_L = P2(1:L/2+1); %rescale axis due to double sideband
    P1_L(2:end-1) = 2*P1_L(2:end-1);
    f = Fs*(0:(L/2))/L;
    
    % Plot, Low Pass Filter's Single-Sided Amplitude Spectrum
    subplot(4,2,[7,8]);
    plot(f,P1_L);
    title('Single-Sided Amplitude Spectrum of X(t) after Low')
    xlabel('f (Hz)');
    ylabel('|P1_L(f)|');
    
    % ======= Testing HIGH pass Filter, Bandpass Version ======== %
    % This should use the low pass infomation, i.e. the low_f, to 
    %    create a bandpass with the butterworth function instead 
    %    of outright using the BANDPASS filter.
    % =========================================================== %
    
    Wn = 0.05;                    % cutoff frequency
    N = 3;                         % 3rd order
    [a,b] = butter(N, Wn, 'High'); % a,b represent the coefficients
    hl_f = filtfilt(a,b,low_f);    % this is the result for the high pass filter
    
    
    % ============ Testing BANDPASS pass filter ================= %
    % This uses the butterworth function. A BANDPASS pass filter gets
    % rid of lower frequencies.
    %    Things to look into: 
    %       1. Cutoff frequency -> Wn
    %       2. Do I need a zero phase filter -> filtfilt()
    %       3. Affects on Amplitude Spectrum
    % =========================================================== %
    Wn = [0.02 0.05];              % cutoff frequency
    N = 3;                         % 3rd order
    [a,b] = butter(N, Wn);         % a,b represent the coefficients
    band_f = filtfilt(a,b,Amp);    % this is the result for the bandpass filter
    
    % Plot, Bandpass Filter 
    figure
    subplot(4,2,[1,2]);
    plot(band_f);
    title('Bandpass Filter');
    
    % Plot, Raw Data
    subplot(4,2,[3,4]);
    plot(Amp);
    title('Raw');
    
    % ============ Testing BANDPASS, Single-Sided Spectrum ========= %
    Y = fft(band_f);    % returns complex form
    P2 = abs(Y/L);      % this is the magnitude
    P1_B = P2(1:L/2+1); % rescale axis due to double sideband
    P1_B(2:end-1) = 2*P1_B(2:end-1);
    f = Fs*(0:(L/2))/L;
    
    
    % Plot, Bandpass Single-Sided Amplitude Spectrum
    subplot(4,2,[7,8]);
    plot(f,P1_B);
    title('Single-Sided Amplitude Spectrum of X(t) after Bandpass')
    xlabel('f (Hz)');
    ylabel('|P1_B(f)|');

    % Plot, Raw Data Single-Sided Amplitdue Spectrum
    subplot(4,2,[5,6]);
    plot(f,P1);
    title('Single-Sided Amplitude Spectrum of X(t)')
    xlabel('f (Hz)');
    ylabel('|P1(f)|');
     
    % ======= Detection - Heartbeat per Minute Algorithm ========== %
    % This will hopefully be using an adaptive threshold for QRS
    %     detection. I just got to figure that part out. hahaha.
    % Facts to Remember:
    %     * A heart CANNOT have two heartbeats in 200 ms.
    % ============================================================= %
    [pks, locs] = findpeaks(band_f, 'MINPEAKDISTANCE', 250); % Is this 200 ms
    lenPks = length(pks);
    
    devi_pks = std(pks);      % The standard deviation of the peaks
    devi_pks = devi_pks / 2; 
    
    mean_pks = mean(pks);     % The mean value of the peaks
    
    realPks  = zeros(lenPks); % Preallocating a size of lenPks
    realLocs = zeros(lenPks); % Preallocating a size of lenPks
    count = 1;                % Number of heartbeats
    
    for pk = 1:lenPks
        % If the peak in pks is w/i the standard deviation range then ...
        if pks(pk) <= (mean_pks + devi_pks) && pks(pk) >=  (mean_pks - devi_pks)
            realPks(count)  = pks(pk);   % add that pk into realPks,
            realLocs(count) = locs(pk);  % add that location into realLocs,
            count = count + 1;           % increase count
        end
    end
    
    realPks = realPks(1:count-1);
    realLocs = realLocs(1:count-1);
    
    % Plots the locations of each peak of interest
    figure
    plot(t,band_f);
    hold on
    scatter(t(realLocs), band_f(realLocs),20,'filled')
    hold off
    
    time  = range(t);
    denom = 60 / time;
    
    fprintf("The time in %f secs. Heartbeat per minute: %d \n", single(time), int16(count*denom));
    
    % ========= Detection - Regular Heartbeat Continuance ========= %
    % Detects if three pairs of heatbeats have distances, vaguely 
    %   similar, through the euclidean distances. Therefore detecting
    %   arrhythmia.
    % ============================================================= %
    
    % Area 1
    p1  = [realLocs(1) realPks(1)];
    p12 = [realLocs(2) realPks(2)];
    
    d1 = norm(p1 - p12); % euclideans distance of area 1
    
    % Area 2
    k = int32(2/3 * count);
    p2  = [realLocs(k)   realPks(k)];
    p22 = [realLocs(k+1) realPks(k+1)];
    
    d2 = norm(p2 - p22); % euclideans distance of area 2
    
    % Area 3
    p3  = [realLocs(count-2)  realPks(count-2)];
    p32 = [realLocs(count-1) realPks(count-1)];
    
    d3 = norm(p3 - p32); % euclideans distance of area 3
    
    % Statistical Information
    dist = [d1 d2 d3];
    dist_mean = mean(dist); 
    dist_std  = std(dist);
    
    % Checks for arrhythmia possibility
    for n = 1:length(dist)
        if (dist(n) <= dist_mean + dist_std) && (dist(n) >= dist_mean - dist_std)
           continue; % Within range 
        else
            fprintf("\n Possible arrhythmia. Dist %d, mean %d, std %d \n", int16(dist(n)), int16(dist_mean), int16(dist_std));
            break;
        end
    end
    
end