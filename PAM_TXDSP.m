% M-PAM waveform generator for the AWG with Pulse shaping and pre-equalization
% Keysight 8195A
% Lei Xue
% June 2020
% V1.2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mkdir('AWG_DSO')
% addpath(genpath('AWG_DSO'))
%% Initialization
clear, close all
clc
modulation='PAM';%'Duobinary','External_load'
% System parameters
Fb = 50e9;  % Baud rate of the aimed signal
Tb = 1/Fb;  % Symbol time interval
Fs = 4*Fb;  % Original Sampling frequency 16
Ts = 1/Fs;
Nss = Fs/Fb;% Oversimpling rate
num_symb = 2^16; %Number symbols
Fs_awg = 90e9; % AWG sampling rate for resampling

% signal parameters
M = 4;  % M-PAM %Duobinary it is still 2!!!!!!!!!!
bits_per_sym = log2(M);
PRBS_order = 'random';

% Pulse shaping initializaiton
PulseShaping = 'RRC'; % 'Nyquist'|'Rcos' |'No_Shaping'|'RRC'
txrolloff = 0.15;     % Rolloff factor 0.15
% Pre-compensation initializaiton
pre_comp = 'Nolinear_THP'; % 'FFE','Linear_THP','NN','Nolinear_THP'
%% Signal generation
switch PRBS_order
    case 7
        load PRBS\PRBS_7
    case 9
        load PRBS\PRBS_9
    case 15
        load PRBS\PRBS_15.mat
    case 19
        load PRBS\PRBS_19
    case 'random'
        tx_data = randi([0 1],num_symb*bits_per_sym,1);
        %save('tx_data.mat','tx_data'); % the orignal bits should be saved for the fisrt time
end
%  rep_PRBS = floor(num_symb*bits_per_sym/(length(PRBS)));
%  tx_data = repmat(PRBS,[rep_PRBS 1]);
L = length(tx_data)/bits_per_sym;

switch modulation
    
    case 'Duobinary'
        
        tx_data = reshape(tx_data,bits_per_sym,L);
        tx_data2(1)=0;
        for i=2:length(tx_data)
            tx_data2(i)=double(xor(tx_data2(i-1),tx_data(i-1)));
        end
        tx_data2=tx_data2';
        tx_data3=tx_data2+[tx_data2(2:end);tx_data2(1)];
        mods = modem.pammod('M',3,'SymbolOrder','bin');
        tx_sig = modulate(mods,tx_data3');
    case 'External_load'
        
        load('original_Data.mat');
%         tx_sig=tx_sig-mean(tx_sig);
%         tx_sig=(tx_sig-min(tx_sig))/(max(tx_sig)-min(tx_sig));
        %csvwrite('.\50G pam4.vtmu_pack\Inputs\Original_Data.csv',tx_sig);
    otherwise
        tx_data = reshape(tx_data,bits_per_sym,L);
        mods = modem.pammod('M',M,'SymbolOrder','bin');
        tx_sig = modulate(mods, bi2de(tx_data','left-msb'));
        %save('original_data.mat','tx_sig');
        
end
tx_up = upsample(tx_sig,Nss); %  数据间不具备相关性所以spectrum是白噪声
num = ones(Nss,1);
den = 1;
tx_wfm=filter(num,den,tx_up); %Idea rectangular pulse shaping
%% Pre-compensation
% Downsample the data rate to AWG sampling rate, then use
% FFE to obtain and save the ffetaps
switch pre_comp
    case 'FFE'
        load('ffetaps.mat');
        h_AWG_DACRate = reshape(w,1,[]);
        length_prefix = floor(length(h_AWG_DACRate)/2);
        tempX = [tx_sig(end -length_prefix+1:end );tx_sig;tx_sig(1:length_prefix)];
        tx_wfm_precom = conv(tempX, fliplr(h_AWG_DACRate));
        tx_wfm_precom(end-2*length_prefix+1:end) = [];
        tx_wfm_precom(1:2*length_prefix) = [];
        tx_wfm_precom = upsample(tx_wfm_precom,Nss);
    case 'NN'
        load('./ml/nn_pre.csv');
        tx_wfm_precom=nn_pre;
        %tx_wfm_precom =tx_wfm_precom-mean(tx_wfm_precom);
        tx_wfm_precom = upsample(tx_wfm_precom,Nss);
    case 'Linear_THP'
        %% THP precoding
        load fb_filter.mat % load the taps from DFE
        Thp_filter = fb_filter;
        tx_wfm_precom=zeros(1,length(tx_sig));
        Thp_filter_in =zeros(1,length(Thp_filter));
        Bk=zeros(1,length(tx_sig));% The values of modulo operation
        error=0; % init value of feedback output
        for i=1:length(tx_sig)
            Pre_mod = tx_sig(i)-error;
            % modulo operation
            Bk(i)=round(Pre_mod/8);
            tx_wfm_precom(i) = Pre_mod-Bk(i)*8;
            
            Thp_filter_in(2:end)=Thp_filter_in(1:end-1);
            Thp_filter_in(1)=tx_wfm_precom(i);
            error=Thp_filter*Thp_filter_in.';
        end
        label_ffe =tx_sig-8*Bk;
        %save('THP_out.mat','label_ffe') % save pre-coded sequence
        tx_wfm_precom = upsample(tx_wfm_precom,Nss);
    case 'Nolinear_THP'
        load fb_volfilter.mat % load the taps from DFE in the receiver
        Thp_volfilter = fb_volfilter;
        tx_wfm_precom=zeros(1,length(tx_sig));
        fb_filter_ip = zeros(1,ch1a);%the input of feedback volterra before caculation
        Thp_volfilter_in =zeros(1,length(Thp_volfilter));
        Bk=zeros(1,length(tx_sig));% The values of modulo operation
        error=0; % init value of feedback output
        for i=1:length(tx_sig)
            Pre_mod = tx_sig(i)-error;
            % modulo 2M operation
            Bk(i)=round(Pre_mod/8);
            tx_wfm_precom(i) = Pre_mod-Bk(i)*8;
            % feedback update
            fb_filter_ip(2:end)=fb_filter_ip(1:end-1);
            fb_filter_ip(1)=tx_wfm_precom(i);
            Thp_volfilter_in = Input_gen(fb_filter_ip,ch1a,ch2a);
            error=Thp_volfilter*Thp_volfilter_in.';
        end
        label_ffe =tx_sig-8*Bk';
        %save('THP_out.mat','label_ffe') % save pre-coded sequence
        tx_wfm_precom = upsample(tx_wfm_precom',Nss);
    otherwise
        tx_wfm_precom = tx_up;
end

%% Nyquist/Raise cosine pulse shaping
switch PulseShaping
    case 'RRC'
        Nsym = 10;           % Filter span in symbols
        sps = Nss;            % Samples per symbol/upsampling rate
        Hd = rcosdesign(txrolloff,Nsym,sps,'sqrt');%square-root rasise cosin filter
        %tx_wfm_precom_ps =conv(tx_wfm_precom,Hd);
        tx_wfm_precom_ps = upfirdn(tx_wfm_precom, Hd, 1);% the data length is len(Hd)+len(tx_sig)*NSs-1
    case 'No_Shaping'
        tx_wfm_precom_ps=tx_wfm;
end

%% Eyediagram before pulse shaping
eyediagram(tx_wfm(1:2^15),3*Nss);
title('Eyediagram before pulse shaping');

%% Eyediagram after pulse shaping
% H = comm.EyeDiagram('SampleRate',100E9,'SamplesPerSymbol',4,...
%     'DisplayMode','2D color histogram',...
%     'YLimits',[-4,4],...
%     'OversamplingMethod' ,'Input interpolation','ShowGrid',0);
% eyeObj.ColorScale = 'Logarithmic';
% H(tx_wfm_precom_ps)
eyediagram(tx_wfm_precom_ps(3*Nss:2^15),5*Nss);
title('Eyediagram after pulse shaping');

%% Resampling the signal with AWG sampling rate
if Fs ~= Fs_awg
    tx_wfm_precom_ps_re = resample(tx_wfm_precom_ps,Fs_awg,Fs);
else
    tx_wfm_precom_ps_re = tx_wfm_precom_ps;
end

%% Plot and compare signal before and after resampling
% Plotting symbols
figure(8);
subplot(4,1,1);
plot(tx_wfm(1:300),'LineWidth',2);
title('No pulse shaping','FontWeight','bold');
% xlabel('first 200 samples');
ylabel('Amplitude');

subplot(4,1,2);
plot(tx_wfm_precom (1:300),'k','LineWidth',2);
title('After pre-compensation','FontWeight','bold');
%xlabel('first 300 samples');
ylabel('Amplitude');

subplot(4,1,3);
plot(tx_wfm_precom_ps(1:300),'r','LineWidth',2);
title('After pulse shaping','FontWeight','bold');
% xlabel('first 200 samples');
ylabel('Amplitude');

subplot(4,1,4);
plot(tx_wfm_precom_ps_re(1:300),'g','LineWidth',2);
title('After re-sampling','FontWeight','bold');
xlabel('first 200 samples');
ylabel('Amplitude');

% Plot spectrum
figure(9);
data1 = double(tx_wfm);
data2 = double(tx_wfm_precom);
data3 = double(tx_wfm_precom_ps);
data4 = double(tx_wfm_precom_ps_re);

Ex_spec1 = data1;
Npoints1 = length(Ex_spec1);
FFT_Ex_1 = fftshift(fft(Ex_spec1));
FFT_Ex1 = abs(FFT_Ex_1)./(length(Ex_spec1));
Frek1 = (Fs*(-(Npoints1)/2:((Npoints1/2)-1)))/Npoints1;

Ex_spec2 = data2;
Npoints2 = length(Ex_spec2);
FFT_Ex_2 = fftshift(fft(Ex_spec2));
FFT_Ex2 = abs(FFT_Ex_2)./(length(Ex_spec2));
Frek2 = (Fs*(-(Npoints2)/2:((Npoints2/2)-1)))/Npoints2;

Ex_spec3 = data3;
Npoints3 = length(Ex_spec3);
FFT_Ex_3 = fftshift(fft(Ex_spec3));
FFT_Ex3 = abs(FFT_Ex_3)./(length(Ex_spec3));
Frek3 = (Fs*(-(Npoints3)/2:((Npoints3/2)-1)))/Npoints3;

Ex_spec4 = data4;
Npoints4 = length(Ex_spec4);
FFT_Ex_4 = fftshift(fft(Ex_spec4));
FFT_Ex4 = abs(FFT_Ex_4)./(length(Ex_spec4));
Frek4 = (Fs_awg*(-(Npoints4)/2:((Npoints4/2)-1)))/Npoints4;

subplot(2,2,1);
plot(Frek1./1e9,10*log10(FFT_Ex1.^2),'.');
title('Spectrum before pulse shaping');

subplot(2,2,2);
plot(Frek2./1e9,10*log10(FFT_Ex2.^2),'.');
title('Spectrum after Pre-equalization');

subplot(2,2,3);
plot(Frek3./1e9,10*log10(FFT_Ex3.^2));
title('Spectrum after Pulse shaping');

subplot(2,2,4);
plot(Frek4./1e9,10*log10(FFT_Ex4.^2));
title('Spectrum after Resampling');

%% Save the waveform in .MAT format
file_prefix = 'PAM';
BaudRT_int = floor(Fb/1e9);
BaudRT_dec = (Fb/1e9-BaudRT_int)*10;
awg_file = [num2str(M) file_prefix '_' num2str(BaudRT_int) 'p' num2str(BaudRT_dec) 'Gbaud' '_PRBS' num2str(PRBS_order) '_' PulseShaping '_' pre_comp]; % Baseline name for Tek AWG files
Waveform_Name = awg_file;
Waveform_Data =  tx_wfm_precom_ps_re; %already a double array
Waveform_Sampling_Rate = Fs_awg;      % Saving sampling rate
%Waveform_Amplitude_1 = 0.300;          %and amplitude in V
%save(['awg/' awg_file '.mat'], 'Waveform_*', '-v6');
Waveform_Data=Waveform_Data-mean(Waveform_Data);
Waveform_Data=(Waveform_Data-min(Waveform_Data))/(max(Waveform_Data)-min(Waveform_Data));
a=floor(length(Waveform_Data)/256); % 数据前加一段0,使长度为128的整数倍
b=(a+3)*256-length(Waveform_Data);
c=zeros(b,1);
Y=[c;Waveform_Data];
Y =single(Y);
csvwrite('.\Txsignal.csv',Y);%有效信号作为同步
XDelta= 1/Fs_awg;
%save(['.\AWG\' awg_file '.mat'], 'Y','XDelta', '-v6');

%% load data to AWG
ArbConfig = loadArbConfig_M8195A();
Setting = loadSetting_M8195A(length(Y));
IQdownload_M8195A(ArbConfig, Setting.fs,Y,  Setting.marker1, ...
    Setting.marker2, Setting.segmNum, Setting.keepOpen, Setting.chMap,...
    Setting.sequence, Setting.run,Setting.segmentLength,Setting.segmentOffset);

