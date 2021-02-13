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

%% Resample to desired sample rate
% load the csv file;
Rx_wfm1 = csvread('.\Sampled data\C1_1dbm00000.csv'); 
Rx_wfm1 =Rx_wfm1-mean(Rx_wfm1);
Rx_wfm1 = resample(Rx_wfm1,Fs,Fd); % downsample只是整数倍的减少
Rx_wfm2 = csvread('.\Sampled data\C1-10dbm00001.csv'); 
Rx_wfm2 =Rx_wfm2-mean(Rx_wfm2);
Rx_wfm2 = resample(Rx_wfm2,Fs,Fd); % downsample只是整数倍的减少

%% Match filtering
Nsym = 10;           % Filter span in symbols
sps = 4;            % Samples per symbol
txrolloff = 0.15;
Hd = rcosdesign(txrolloff,Nsym,sps,'sqrt');% Square-root rasise cosin filter
rxFilt1 = upfirdn(Rx_wfm1,Hd,1,4);
rxFilt2 = upfirdn(Rx_wfm2,Hd,1,4);
%% Synchronization
% Load original data
  load 'original_data.mat'
  OriginalData = tx_sig*(-1);
  SynSeed = upsample(OriginalData,1);
% Extract two segment of data
  CorrelationResult = conv(rxFilt1(1:end), conj(SynSeed(end:-1:1)));
  [Max,Index] = max(CorrelationResult);
  ExtractedSignal1 = rxFilt1(Index-length(OriginalData)+1:Index); 
  CorrelationResult = conv(rxFilt2(1:end), conj(SynSeed(end:-1:1)));
  [Max,Index2] = max(CorrelationResult);
  figure
  plot(CorrelationResult);
  ExtractedSignal2 = rxFilt2(Index2-length(OriginalData)+1:Index2);
% Combine the data and save as csv file
  ExtractedSignal = [ExtractedSignal1;ExtractedSignal2];
  csvwrite('.\ml\Rxsignal.csv',ExtractedSignal);
  OriginalData = repmat(OriginalData,2,1);
  OriginalData(find(OriginalData==-3))=0;
  OriginalData(find(OriginalData==1))=2;
  OriginalData(find(OriginalData==-1))=1;
  csvwrite('.\ml\Txsignal.csv',OriginalData);