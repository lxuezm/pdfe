function [output, TrainingSignal] = FFE_Equalizer(InputSignal, TrainingSignal, AlgType, ...
																							 FFETaps, alpha, epoch)
	% This function performs the feed forward equalization with LMS or RLS algorithm.
	% First, the InputSignal and TrainingSignal will be normalized to 0-1.
	% The training will use all the InputSignal and will be performed epoch times.
	% The alpha is the learning rate of LMS or the forgetting factor of RLS, 
	% which should be chosen carefully withthe help of curve of convergence. 
	% After training, the equalization will be performed and then the result will 
	% be returned.
	%
	% input: 
	%     InputSignal
	%       The input signal to be equalized.
	%     TrainingSignal
	%       The actual signal to be equalized to.
	%     AlgType
	%       'lms' for LMS or 'rls' for RLS.
	%     FFETaps (optional)
	%       The numbers of FFE taps which must be odd.
	%       Default: 5
	%     alpha (optional)
	%       The learning rate of LMS algorithm or the forgetting factor of RLS.
	%       Default: 0.01 for AlgType = 'lms', 0.99 for AlgType = 'rls'
	%     epoch (optional)
	%       The epoch of the learning of LMS through all the input signal.
	%       Default: 1
	% output:
	%     output
	%       The equalized signal with the same length of InputSignal
	%       Size: length(InputSignal), 1
	%     w
	%       Weights of FFE
	%       Size: FFETaps, 1
	%     costs
	%       The costs after each training epoch, which is used to draw a 
	%       curve of convergence and thus determine the best learning rate.
	
	%% Parameter Checking
	narginchk(3, 6);
	
	if ~exist('FFETaps','var') || isempty(FFETaps)
		FFETaps = 5;
	end
	
	if ~exist('alpha','var') || isempty(alpha)
		if AlgType == 'lms'
			alpha = 0.01;
		elseif AlgType == 'rls'
			alpha = 0.99;
		end
	end
% 	
% 	if ~exist('epoch','var') || isempty(epoch)
% 		epoch = 1;
% 	end
% 	
% 	% TODO add some parameter checking
% 	if mod(FFETaps, 2) == 0
% 		error('linearFeedForwardEqualize:argChk', 'FFE taps must be odd');
% 	end
% 	
% 	if ~strcmp(AlgType, 'lms') && ~strcmp(AlgType, 'rls')
% 		error('linearFeedForwardEqualize:argChk', 'AlgType must be lms or rls');
% 	end
% 
% 	if FFETaps <= 0
% 		error('linearFeedForwardEqualize:argChk', 'FFE taps must be bigger than 0');
% 	end
% 	
 %% Signal Normalization and Duplication
% 	% InputSignal and TrainingSignal is normalized to the range between -3 to 3.
%     InputSignal = InputSignal -mean(InputSignal);
    InputSignal = 6/(max(InputSignal)-min(InputSignal))*(InputSignal-min(InputSignal))-3;
% 	plot(InputSignal*(-1),'.')
    % Both signal is duplicated for better performance
% 	InputSignalDup = repmat(InputSignal, 2, 1);
%   TrainingSignalDup = repmat(TrainingSignal, 2, 1);
    TrainingSignalDup = upsample(TrainingSignal, 1);
	% Zero Padding for input signal
	InputSignalZP = [zeros(floor(FFETaps/2), 1); InputSignal; zeros(floor(FFETaps/2), 1)];
%   InputSignalZP = [zeros(floor(FFETaps/2), 1); InputSignal; zeros(floor(FFETaps/2), 1)];
	%% Weights Initializing
	% TODO choose on : randomly init weights or init to 0
	% w = zeros(FFETaps, 1);
	% w(floor(length(w)/2) + 1) = 1;
	rng('shuffle');
	w = rand(FFETaps, 1);
	w = 2 * w - 1;
	
	%% Training
	costs = zeros(epoch, 1);
	y = zeros(size(TrainingSignalDup));
	if AlgType == 'lms'
		% The LMS learning algorithm
		for n = 1 : epoch
			for i = 1 : length(InputSignalZP) - FFETaps + 1
				y(i) = w' * InputSignalZP(i : i + FFETaps - 1);
				w = w - alpha * (y(i) - TrainingSignalDup(i)) * InputSignalZP(i : i + FFETaps - 1);
				costs(n) = costs(n) + 0.5 * ((y(i) - TrainingSignalDup(i)) ^ 2);%Accumulated error
			end
			% Record the cost/error of each epoch
			costs(n) = costs(n) / (length(InputSignalZP) - FFETaps + 1);% mean error per bit
		end
	elseif AlgType == 'rls'
		% The RLS learning algorithm
		Sd = eye(FFETaps);
		for n = 1 : epoch
			for i = 1 : length(InputSignalZP) - FFETaps + 1
				x = InputSignalZP(i : i + FFETaps - 1);
				e = TrainingSignalDup(i) - w' * x;
				phi = Sd * x;
				Sd = (1 / alpha) * (Sd - (phi * phi') / (alpha + phi' * x));
				w = w + e * Sd * x;
				costs(n) = costs(n) + 0.5 * (e ^ 2);
			end
			% Record the cost/error of each epoch
			costs(n) = costs(n) / (length(InputSignalZP) - FFETaps + 1);
		end
    end 
%% Use DFE module in matlab
% dfeObj=dfe(10,5,lms(0.0001));
% dfeObj.nSampPerSym = 1;
% dfeObj.SigConst=[-3,-1,1,3];
% dfeObj.ResetBeforeFiltering = 0;
% dfeObj.Weights=[1,zeros(1,14)];
% for n=1:epoch
%     [eq_sig,~,err] = equalize(dfeObj,InputSignal,TrainingSignalDup);
% end    
% output = equalize(dfeObj,InputSignal);
% plot(output,'.')
    
%% DFE self-writing
% dec_seq = zeros(1,data_len-ff_filter_len+1);% output from dfe
% ff_filter_ip = zeros(1,ff_filter_len); % feedforward filter input
% fb_filter_ip = zeros(1,fb_filter_len); % feedback filter input
% fb_filter_op = 0; % feedback filter output symbol
% for i1=1:data_len-ff_filter_len+1 % steady state part
%          ff_filter_ip(2:end)=ff_filter_ip(1:end-1);
%          ff_filter_ip(1) = chan_op(i1);
%          ff_filter_op = ff_filter*ff_filter_ip.';
%          
%          ff_and_fb = ff_filter_op-fb_filter_op;
%         
%          % hard decision
%          temp = ff_and_fb<0;
%          dec_seq(i1) = 1-2*temp;
%          
%          fb_filter_ip(2:end)=fb_filter(1:end-1);
%          fb_filter_ip(1) = dec_seq(i1);
%          
%          fb_filter_op = fb_filter*fb_filter_ip.'; % feedback filter output
% end
% % demapping symbols back to bits
% dec_a = dec_seq<0;

%% Save the taps for pre-eq
save('ffetaps.mat','w')

%% Convergence trends figure
figure;
plot(costs);
title('Curve of Convergence');
xlabel('Epoch'); ylabel('Cost');

%% Using Trained Weights to Equalize Data
for i = 1 : length(InputSignalZP) - FFETaps + 1
    y(i) = w' * InputSignalZP(i : i + FFETaps - 1);
end
%%  Choose a half as the output
output = downsample(y,1);
% %  output = y(1:length(y)/2);
% figure;
% %plot(output*(-1),'.');
% a=[1:1:65536];
% binscatter(a',output*(-1),[250,250]);
% colormap(gca,'hot');
% plot(output*(-1),'.');

