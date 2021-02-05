%% This code shows the work flow for a receiver side dsp: data generation,
%% synchronization, signal extraction, PAM4 decision and BER calculation.
% Rx_DSP for PAM signal after direct detection
% Lei Xue
% July 2020
% V1.1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
clc;

%% change the current directory to the folder which contains this m file
cd(fileparts(which(mfilename)));
% Signal parameters
M = 4; % M-PAM 
bits_per_sym = log2(M);
Fd = 40e9;    % Sample rate of DSO
Fb = 25e9;  % Baud rate of the aimed signal
Fs = 4*Fb; % Desired Sample rate
real_time=0;

% Label of algorithms
FFE_label = 1;
Volterra_label = 0;
MLSD_label = 0;

%% Resample to desired sample rate
if real_time==1
    % Directly get data from DSO
    filename =['-9'];
    Rx_wfm=osc_test(filename);
    Rx_wfm=Rx_wfm';
    Rx_wfm =Rx_wfm-mean(Rx_wfm);
else
    % load the csv file;
    %Rx_wfm = csvread('.\Sampled data\dml with different output power\C1dml6dbm00002.csv',3,0); 
    %Rx_wfm = csvread('.\Sampled data\27_July_SOA2\C1-1dbm00000.csv'); 
    load('.\Sampled data\Receive sensitivity with SOA\250mA\rx-0dBm.mat') 
    Rx_wfm=data_temp1';%%11
    Rx_wfm=Rx_wfm;
    Rx_wfm =Rx_wfm-mean(Rx_wfm);
end 
Rx_wfm = resample(Rx_wfm,Fs,Fd); % downsample只是整数倍的减少;resample to interger times of signal data rate
% eyediagram of sampled signal
eyediagram(Rx_wfm(1:2^15),3*4);
title('eyediagram of sampled signal');
% H(Rx_wfm)
%% spectrum of sampled signal
% Ex_spec3 = Rx_wfm;
% Npoints3 = length(Ex_spec3);
% FFT_Ex_3 = fftshift(fft(Ex_spec3));
% FFT_Ex3 = abs(FFT_Ex_3)./(length(Ex_spec3));
% Frek3 = (Fs*(-(Npoints3)/2:((Npoints3/2)-1)))/Npoints3;
% figure
% plot(Frek3./1e9,10*log10(FFT_Ex3.^2));
% title('Spectrum of sampled signal');

%% Match filtering
Nsym = 10;           % Filter span in symbols
sps = 4;             % Samples per symbol
txrolloff = 0.15;    % The compressed bandwidth is (1+)
Hd = rcosdesign(txrolloff,Nsym,sps,'sqrt');% Square-root rasise cosin filter
rxFilt = upfirdn(Rx_wfm,Hd,1,4);

%% Synchronization
% Load original data
  load 'original_data1.mat'
  OriginalData = tx_sig*(-1);
% OriginalData = tx_sig;
  SynSeed = upsample(OriginalData,1);
  CorrelationResult = conv(rxFilt(1:end), conj(SynSeed(end:-1:1)));
  [Max,Index] = max(CorrelationResult);
  ExtractedSignal = rxFilt(Index-length(OriginalData)+1:Index); 
  figure
  plot(CorrelationResult);
 
%% Downsampling to 1sps
% save the data for Neural network training
 ExtractedSignal =downsample(ExtractedSignal,1);
%  csvwrite('.\ml\Rxsignal2.csv',ExtractedSignal); 
  
%% Pre_equalization with 2 sps
% DFE 
if FFE_label ==1
   Rx_eq = FFE_Equalizer(ExtractedSignal, OriginalData, 'lms',117, 0.001,5);
end
csvwrite('.\ml\Rxsignal1.csv',Rx_eq); 
%% Eyediargam after FFE
H = comm.EyeDiagram('SampleRate',25E9,'SamplesPerSymbol',1,...
    'DisplayMode','2D color histogram',...
    'YLimits',[-4,4],...
    'OversamplingMethod' ,'Input interpolation','ShowGrid',0);
eyeObj.ColorScale = 'Logarithmic';
H(Rx_eq*(-1))
% Volterra
if Volterra_label ==1
    Rx_eq = Volterra_Equalize(ExtractedSignal, OriginalSignal, 'lms', 10, 47, 0.1, 11, [], 5, [], false);
end
% MLSD
if MLSD_label ==1
    AlphaTap = 0.001;%0.4
    BitPerSym = 3;
    RxSym = Rx_eq;
    RxSym = conv(RxSym,[1 AlphaTap]);
    RxSym = RxSym(1:end-1);
    RxSym = MLSD(RxSym,BitPerSym,AlphaTap);
    RxSym=RxSym';
end

%% Downsample to 1sps and BER caculation
[BitErrorRate, SymErrorRate, BitErrorNum,SymErrorNum] = Decision_Cal_Ber(Rx_eq,OriginalData,M);
fprintf('Sample error number : %d \n', SymErrorNum);
fprintf('Bit error number : %d \n', BitErrorNum);
fprintf('SER: %e\n', SymErrorRate);
fprintf('BER: %e\n', BitErrorRate);

%% save the original data as target for NN training
% OriginalData(find(OriginalData==-3))=0;
% OriginalData(find(OriginalData==1))=2;
% OriginalData(find(OriginalData==-1))=1;

% OriginalData = repmat(OriginalData,2,1); 
% csvwrite('.\ml\Txsignal2.csv',OriginalData);

Rxsignal0=csvread('.\ml\Rxsignal0.csv'); 
Rxsignal1=csvread('.\ml\Rxsignal1.csv'); 
Rxsignal=[Rxsignal0;Rxsignal1];
csvwrite('.\ml\Rxsignal.csv',Rxsignal);

% OriginalData(find(OriginalData==-7))=0;
% OriginalData(find(OriginalData==5))=6;
% OriginalData(find(OriginalData==3))=5;
% OriginalData(find(OriginalData==1))=4;
% OriginalData(find(OriginalData==-1))=3;
% OriginalData(find(OriginalData==-3))=2;
% OriginalData(find(OriginalData==-5))=1;
% csvwrite('.\ml\Txsignal.csv',OriginalData);














