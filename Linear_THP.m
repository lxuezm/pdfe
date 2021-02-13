%% Linear Tomlinson-Harashima precoding simultion in a AWGN channel
%  Base on Decision feedback equalizer
%  Lei Xue
%  2021-02-10
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
ff_filter_len = 21; % feedforward filter length, should be a
fb_filter_len =9; % feedback filter length

figure(1)
hist(training_seq);
title('Histogram before precoder');
%% THP precoding 
load fb_filter.mat % load the taps from DFE in the receiver
Thp_filter = fb_filter;
Input_precode=zeros(1,length(training_seq));
Thp_filter_in =zeros(1,length(Thp_filter));
Bk=zeros(1,length(training_seq));% The values of modulo operation
error=0; % init value of feedback output
for i=1:training_len
    Pre_mod = training_seq(i)-error;
    % modulo operation
    Bk(i)=round(Pre_mod/8);
    Input_precode(i) = Pre_mod-Bk(i)*8;
    
    Thp_filter_in(2:end)=Thp_filter_in(1:end-1);
    Thp_filter_in(1)=Input_precode(i);
    error=Thp_filter*Thp_filter_in.';
end
label_ffe =training_seq-8*Bk;%Generate the ideal output of THP coder, ref.Po

figure(2)
histogram(Input_precode);
title('Histogram after precoder');
%% Channel transmission
% impulse response of the channel 
fade_chan =  [1 0.234 0.407 0.815 0.407];%(PROAKIS B CHANNEL)
fade_chan = fade_chan/norm(fade_chan);  %sqrt(fade_chan^2)
chan_len = length(fade_chan);
% awgn
noise = normrnd(0,sqrt(noise_var),1,training_len+chan_len-1);
% channel output
chan_op= conv(fade_chan,Input_precode)+noise;
%% ------------ DFE filter------------------------------------------
%tap update
InputSignalZP = [zeros(1,floor(ff_filter_len/2)) chan_op(1:end-4) zeros(1,floor(ff_filter_len/2))];% add zero for entire output
ff_filter = zeros(1,ff_filter_len); % feedforward filter initialization
fb_filter =zeros(1,fb_filter_len);
ff_filter_ip = zeros(1,ff_filter_len); % feedforward filter input
fb_filter_ip = zeros(1,fb_filter_len); % feedback filter input
fb_filter_op = 0; % feedback filter output symbol
dec_seq = zeros(1,length(InputSignalZP)-ff_filter_len+1);% output sequence from dfe
epoch=10;
costs = zeros(epoch, 1);
step_size = 0.001;
% estimating the autocorrelation of received sequence at zero lag
% Rvv0 = (chan_op*chan_op')/(training_len+chan_len-1);
% maximum step size
% max_step_size = 2/(ff_filter_len*Rvv0+fb_filter_len*(1));
% step_size = 0.125*max_step_size; % step size 
%% loop for DFE
% DFE training for taps update
% for n = 1 : epoch
% 	for i = 1 : length(InputSignalZP) - ff_filter_len + 1
%           ff_filter_ip=InputSignalZP(i : i + ff_filter_len - 1);
% 		    ff_filter_op = ff_filter * ff_filter_ip.';
%           
%           ff_and_fb = ff_filter_op-fb_filter_op; % difference between ffe and dfe output
%           error = ff_and_fb-training_seq(i);
%           hard decision
%           dec_seq(i)=decision(ff_and_fb);
%           
%           taps updates
%           ff_filter = ff_filter - step_size * error * ff_filter_ip;
%           fb_filter = fb_filter+step_size*error*fb_filter_ip;
%           feedback output
%           fb_filter_ip(2:end) = fb_filter_ip(1:end-1);
%           fb_filter_ip(1)=dec_seq(i);
%           fb_filter_op = fb_filter*fb_filter_ip.';
%           
% 		  costs(n) = costs(n) + 0.5 * (error^ 2); %Accumulated error
%     end
% 		  Record the cost/error of each epoch
% 		  costs(n) = costs(n) / (length(InputSignalZP) - ff_filter_len + 1);% mean error per bit
% end
%% loop for FFE
% load ff_filter
% tx_wfm_precom = conv(chan_op, fliplr(ff_filter));
% FFE training for taps update
for n = 1 : epoch
	for i = 1 : length(InputSignalZP) - ff_filter_len + 1
          ff_filter_ip=InputSignalZP(i : i + ff_filter_len - 1);
		  dec_seq(i) = ff_filter * ff_filter_ip.';
          error = dec_seq(i)-label_ffe(i); % Here the label is generated after precoder
          % taps updates
          ff_filter = ff_filter - step_size * error * ff_filter_ip; 
		  costs(n) = costs(n) + 0.5 * (error^ 2); %Accumulated error
    end
		  % Record the cost/error of each epoch
		  costs(n) = costs(n) / (length(InputSignalZP) - ff_filter_len + 1);% mean error per bit
end
figure(3)
histogram(dec_seq);
title('the histogram of FFE output')
%% Save DFE taps
% save('fb_filter.mat','fb_filter')
% save('ff_filter.mat','ff_filter')
%% Demapping symbols back to bits
% dec_a = dec_seq<0;
%% THP decode
dec_seq_demod = dec_seq-round(dec_seq/8)*8;
% decision
for i=1:length(dec_seq_demod)
    a(i)=decision(dec_seq_demod(i));
end
%% Bit error rate
% ber = nnz(dec_seq(1:end-ff_filter_len+1)-training_seq(1:end-ff_filter_len+1))/(training_len-ff_filter_len+1)
ber = nnz(a-training_seq)/training_len
