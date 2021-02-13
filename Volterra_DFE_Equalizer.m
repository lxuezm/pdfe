function output = Volterra_DFE_Equalizer(InputSignal, TrainingSignal)
% This function performs the volterra feedforward equalization with LMS or RLS algorithm.
% For now only 1st-3rd order are supported. 1st order must be included, while 2nd
% and 3rd orders are optional and controlled by the ChanLen2nd and ChanLen3rd
% flags. The LMS learning rate of different order can be different by adjusting the
% Alpha1st, Alpha2nd and Alpha3rd. The equalizer will be trained on the
% InputSignal epoch times and will then perform a equalization.
%% Signal Normalization 
InputSignal=InputSignal-mean(InputSignal);
InputSignal=(InputSignal-min(InputSignal))/(max(InputSignal)-min(InputSignal))*6-3;
%% Volterra weights initialization
ch1 = 21;ch2 =5;ch3=0; %Volterra ffe kernel size
ch1a = 15;ch2a = 0;ch3a=0;% Volterra dfe kernel size
step_size = 0.00001; % learning rate

[Kernelsize1,max1] = Kernel_cal(ch1,ch2,ch3); % calculate the feedforward kernel size
[Kernelsize2,max2] = Kernel_cal(ch1a,ch2a,ch3a);%calculate the feedback kernel size

ff_volfilter =  zeros(1,Kernelsize1);% Volterra filter initialization
ff_filter_ip = zeros(1,max1);
ff_volfilter_ip = zeros(1,Kernelsize1); % feedforward volterra filter input after calculation

fb_volfilter =zeros(1,Kernelsize2);
fb_filter_ip = zeros(1,max2);
fb_volfilter_ip = zeros(1,Kernelsize2); % feedback volterra filter input
fb_volfilter_op = 0; % feedback filter output symbol

InputSignalZP = [zeros(1,floor(ch1/2)) InputSignal zeros(1,floor(ch1/2))];% add zero for entire output
dec_seq = zeros(1,length(InputSignalZP)-ch1+1);% Output  from Volterra dfe
training_seq = TrainingSignal;
epoch=10;
costs = zeros(epoch, 1);
%% Loop of Volterra-DFE
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
    costs(n) = costs(n) / (length(InputSignalZP) - ff_filter_len + 1);% mean error per bit
end
ber = nnz(dec_seq-training_seq)/training_len
%     Save Volterra DFE taps
save('fb_volfilter.mat','fb_volfilter');
%% Convergence trends figure
figure;
plot(costs);
title('Curve of Convergence');
xlabel('Epoch'); ylabel('Cost');
% TODO choose a half of the output
output = dec_seq;
figure;
plot(output,'.');
end
