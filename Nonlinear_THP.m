%% Nonlinear Tomlinson-Harashima precoding simultion in a AWGN channel
%  Base on Decision feedback equalizer
%  Lei Xue
%  2021-02-11
% Designed for PAM4 signal

clear all
close all
clc

%% Parameter initialization
load original_data.mat % load pam4 sequence
training_seq = tx_sig';
training_len = length(training_seq);%length of the training sequence

snr_dB = 30; % snr in dB
snr = 10^(0.1*snr_dB); % power 'w'
noise_var = 1/(2*snr); % noise variance

ch1 = 7;ch2 =5;ch3=0; %Volterra ffe kernel size
ch1a = 9;ch2a = 5;ch3a=0;% Volterra dfe kernel size
THP=0; %THP control label
figure(1)
hist(training_seq);
title('Histogram before precoder');
%% Volterra THP precoding 
if THP==1
    load fb_volfilter.mat % load the taps from DFE in the receiver
    Thp_volfilter = fb_volfilter;
    Input_precode=zeros(1,length(training_seq));
    fb_filter_ip = zeros(1,ch1a);%the input of feedback volterra before caculation
    Thp_volfilter_in =zeros(1,length(Thp_volfilter));
    Bk=zeros(1,length(training_seq));% The values of modulo operation
    error=0; % init value of feedback output
    for i=1:training_len
        Pre_mod = training_seq(i)-error;
        % modulo 2M operation
        Bk(i)=round(Pre_mod/8);
        Input_precode(i) = Pre_mod-Bk(i)*8;
        % feedback update
        fb_filter_ip(2:end)=fb_filter_ip(1:end-1);
        fb_filter_ip(1)=Input_precode(i);
        Thp_volfilter_in = Input_gen(fb_filter_ip,ch1a,ch2a);
        error=Thp_volfilter*Thp_volfilter_in.';
    end
    label_ffe =training_seq-8*Bk; %Generate the ideal output of volterra THP coder, ref.Po-digital communication
    
    figure(2)
    histogram(Input_precode);
    title('Histogram after precoder');
end
%% Channel transmission
% impulse response of the channel 
fade_chan =  [1 0.234 0.407 0.815 0.407];%(PROAKIS B CHANNEL)[1 0.234 0.407 0.815 0.407]
fade_chan = fade_chan/norm(fade_chan);  %sqrt(fade_chan^2)
chan_len = length(fade_chan);
% awgn
noise = normrnd(0,sqrt(noise_var),1,training_len+chan_len-1);
% channel output
if THP==1
    chan_op= conv(fade_chan,Input_precode);
else
    chan_op= conv(fade_chan,training_seq);
end
% chan_op= conv(fade_chan,training_seq)+noise;
for j=1:length(chan_op)
        chan_op(j)=chan_op(j)+0.06*chan_op(j)^2;
end
chan_op=chan_op+noise;
% chan_op=chan_op-mean(chan_op);
% chan_op=(chan_op-min(chan_op))/(max(chan_op)-min(chan_op))*6-3;

%% ------------ Volterra_DFE filter------------------------------------------
%taps update
[Kernelsize1,max1] = Kernel_cal(ch1,ch2,ch3); % calculate the kernel size
[Kernelsize2,max2] = Kernel_cal(ch1a,ch2a,ch3a);

ff_volfilter =  zeros(1,Kernelsize1);% Volterra filter initialization
ff_filter_ip = zeros(1,ch1);
ff_volfilter_ip = zeros(1,Kernelsize1); % feedforward volterra filter input after calculation

fb_volfilter =zeros(1,Kernelsize2);
fb_filter_ip = zeros(1,ch1a);
fb_volfilter_ip = zeros(1,Kernelsize2); % feedback volterra filter input
fb_volfilter_op = 0; % feedback filter output symbol
% chan_op=(chan_op-min(chan_op))/(max(chan_op)-min(chan_op));
InputSignalZP = [zeros(1,floor(ch1/2)) chan_op(1:end-4) zeros(1,floor(ch1/2))];% add zero for entire output
dec_seq = zeros(1,length(InputSignalZP)-ch1+1);% Output  from Volterra dfe

epoch=10;
costs = zeros(epoch, 1);
step_size = 0.00001;
% estimating the autocorrelation of received sequence at zero lag
% Rvv0 = (chan_op*chan_op')/(training_len+chan_len-1);
% maximum step size
% max_step_size = 2/(ff_filter_len*Rvv0+fb_filter_len*(1));
% step_size = 0.125*max_step_size; % step size 
%% loop for Volterra DFE
% % DFE training for taps update
if THP==0
    for n = 1 : epoch
        for i = 1: length(InputSignalZP) - ch1 + 1
            
            ff_filter_ip =InputSignalZP(i : i + ch1 - 1); % input sequence
            ff_volfilter_ip = Input_gen(ff_filter_ip,ch1,ch2);% the first and second order input for Volterra
            ff_volfilter_op = ff_volfilter * ff_volfilter_ip.';
            
            ff_and_fb = ff_volfilter_op-fb_volfilter_op; % difference between ffe and dfe output
            error = ff_and_fb-training_seq(i);
            % hard decision
            dec_seq(i)=decision(ff_and_fb);
            
            %taps updates
            ff_volfilter = ff_volfilter - step_size * error * ff_volfilter_ip;
            fb_volfilter = fb_volfilter+step_size*error*fb_volfilter_ip;
            %feedback output
            fb_filter_ip(2:end) = fb_filter_ip(1:end-1);
            fb_filter_ip(1)=dec_seq(i);
            fb_volfilter_ip=Input_gen(fb_filter_ip,ch1a,ch2a);
            fb_volfilter_op = fb_volfilter*fb_volfilter_ip.';
            
            costs(n) = costs(n) + 0.5 * (error^ 2); %Accumulated error
        end
        %Record the cost/error of each epoch
        costs(n) = costs(n) / (length(InputSignalZP) - ch1 + 1);% mean error per bit
    end
    ber = nnz(dec_seq-training_seq)/training_len   
    %     Save Volterra DFE taps
    save('fb_volfilter.mat','fb_volfilter','ch1a','ch2a');
    %% loop for Volterra FFE
else
    for n = 1 : epoch
        for i = 1: length(InputSignalZP) - ch1 + 1
            
            ff_filter_ip =InputSignalZP(i : i + ch1 - 1); % input sequence
            ff_volfilter_ip = Input_gen(ff_filter_ip,ch1,ch2);% the first and second order input for Volterra
            dec_seq(i) = ff_volfilter * ff_volfilter_ip.';
            error = dec_seq(i)-label_ffe(i);
            %taps updates
            ff_volfilter = ff_volfilter - step_size * error * ff_volfilter_ip;
            
            costs(n) = costs(n) + 0.5 * (error^ 2); %Accumulated error
        end
        %Record the cost/error of each epoch
        costs(n) = costs(n) / (length(InputSignalZP) - ch1 + 1);% mean error per bit
    end
    figure(3)
    histogram(dec_seq);
    title('the histogram of Volterra FFE output')
    %% THP decode
    dec_seq_demod = dec_seq-round(dec_seq/8)*8;
    % decision
    for i=1:length(dec_seq_demod)
        a(i)=decision(dec_seq_demod(i));
    end
    %% Bit error rate
    % ber = nnz(dec_seq(1:end-ff_filter_len+1)-training_seq(1:end-ff_filter_len+1))/(training_len-ff_filter_len+1)
    ber = nnz(a-training_seq)/training_len
end

