function output = DFE_Equalizer(InputSignal, TrainingSignal,ff_filter_len,fb_filter_len,step_size,epoch)																							 FFETaps, alpha, epoch)
% This function performs the feed decision equalization with LMS algorithm.

%% Signal Normalization and Duplication
%   InputSignal = InputSignal -mean(InputSignal);
InputSignal = 6/(max(InputSignal)-min(InputSignal))*(InputSignal-min(InputSignal))-3;
% 	plot(InputSignal*(-1),'.')
% Both signal is duplicated for better performance
% 	InputSignalDup = repmat(InputSignal, 2, 1);
%   TrainingSignalDup = repmat(TrainingSignal, 2, 1);
TrainingSignalDup = upsample(TrainingSignal, 1);
InputSignalZP = [zeros(floor(FFETaps/2), 1); InputSignal; zeros(floor(FFETaps/2), 1)];% Zero Padding for input signal for stable output
%% Weights Initializing
% TODO choose on : randomly init weights or init to 0
% w = zeros(FFETaps, 1);
% w(floor(length(w)/2) + 1) = 1;
ff_filter = zeros(1,ff_filter_len); % feedforward filter initialization
fb_filter =zeros(1,fb_filter_len);
ff_filter_ip = zeros(1,ff_filter_len); % feedforward filter input
fb_filter_ip = zeros(1,fb_filter_len); % feedback filter input
fb_filter_op = 0; % feedback filter output symbol
dec_seq = zeros(1,length(InputSignalZP)-ff_filter_len+1);% output sequence from dfe
% epoch=10;
costs = zeros(epoch, 1);
% step_size = 0.001;
%% TrainingLoop for DFE
for n = 1 : epoch
	for i = 1 : length(InputSignalZP) - ff_filter_len + 1
          ff_filter_ip=InputSignalZP(i : i + ff_filter_len - 1);
		  ff_filter_op = ff_filter * ff_filter_ip.';
          
          ff_and_fb = ff_filter_op-fb_filter_op; % difference between ffe and dfe output
          error = ff_and_fb-training_seq(i);
          % hard decision
          dec_seq(i)=decision(ff_and_fb);
          
          % taps updates
          ff_filter = ff_filter - step_size * error * ff_filter_ip;
          fb_filter = fb_filter+step_size*error*fb_filter_ip;
          % feedback output
          fb_filter_ip(2:end) = fb_filter_ip(1:end-1);
          fb_filter_ip(1)=dec_seq(i);
          fb_filter_op = fb_filter*fb_filter_ip.';
          
		  costs(n) = costs(n) + 0.5 * (error^ 2); %Accumulated error
    end
% 		  Record the cost/error of each epoch
		  costs(n) = costs(n) / (length(InputSignalZP) - ff_filter_len + 1);% mean error per bit
end

%% Save the taps for pre-eq
%save('fb_filter.mat','fb_filter')
%save('ff_filter.mat','ff_filter')
%% Convergence trends figure
figure;
plot(costs);
title('Curve of Convergence');
xlabel('Epoch'); ylabel('Cost');
%%  Choose a half as the output
output = downsample(dec_seq,1);
% %  output = y(1:length(y)/2);
% figure;
% %plot(output*(-1),'.');
% a=[1:1:65536];
% binscatter(a',output*(-1),[250,250]);
% colormap(gca,'hot');
% plot(output*(-1),'.');
end
