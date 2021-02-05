function [BitErrorRate, SymErrorRate, BitErrorNum, SymErrorNum] = Decision_Cal_Ber(InputSignal, OriginalData,M, threshold)
	% This function performs the PAM4 decision and bit error ratio calculation for PAM4.
	% First the PAM4 decision will be performed and then the error counting and error ratio
	%	calculation. The %threshold% for PAM4 decision is a optional parameter, which have the
	% default value [0.25; 0.5; 0.75], allowing outside to perform some optimization on
	% the decision threshold of PAM4.
	%
	% input:
	%     InputSignal
	%       The input signal to be decided and error-calculated.
	%     OriginalData
	%       The origin data to be compared to the decided input signal to get the errors.
	%     threshold (optional)
	%       The PAM4 decision threshold, which can be adjusted outside this function.
	%       Default: [0.25; 0.5; 0.75]
	% output:
	%     BitErrorRate
	%       The bit error ratio of the %InputSignal% comparing to %OriginalData%.
	%     SymErrorRate
	%       The symbol error ratio of the %InputSignal% comparing to %OriginalData%.
	%     BitErrorNum
	%       The number of bit error, which is 1/sym when only 1 bit changes and 2/sym
	%       when both bits change.
	narginchk(2, 3);
if M==4
	% Input Signal Normalization
    InputSignal = InputSignal -mean(InputSignal);
    InputSignal = 6/(max(InputSignal)-min(InputSignal))*(InputSignal-min(InputSignal))-3;
%   OriginalData = OriginalData - mean(OriginalData);

	if ~exist('threshold','var') || isempty(threshold)
		threshold = zeros(3, 1);
		threshold(2) = mean(InputSignal);
		threshold(1) = mean(InputSignal(find(InputSignal <= threshold(2))));
		threshold(3) = mean(InputSignal(find(InputSignal > threshold(2))));
		% threshold = [0.25; 0.5; 0.75];
	end

	%% Input Signal Decision
	InputSignal(find(InputSignal > threshold(3))) = 3;
	InputSignal(find(InputSignal > threshold(2) & InputSignal <= threshold(3))) = 1;
	InputSignal(find(InputSignal > threshold(1) & InputSignal <= threshold(2))) = -1;
	InputSignal(find(InputSignal <= threshold(1))) = -3;

	%% Error Counting
	SymErrorNum = length(find(InputSignal ~= OriginalData));
	% The bit error number is 1/sym when only 1 bit changes and 2/sym when both bits change.
% 	BitErrorNum = SymErrorNum + length(find((InputSignal == 3) & (OriginalData == 0))) ...
% 														+ length(find((InputSignal == 2) & (OriginalData == 1))) ...
% 														+ length(find((InputSignal == 1) & (OriginalData == 2))) ...
% 														+ length(find((InputSignal == 0) & (OriginalData == 3)));
    BitErrorNum = SymErrorNum + length(find((InputSignal == 3) & (OriginalData == -3))) ...
														+ length(find((InputSignal == 1) & (OriginalData == -1))) ...
														+ length(find((InputSignal == -1) & (OriginalData == 1))) ...
														+ length(find((InputSignal == -3) & (OriginalData == 3)));
    SymErrorRate = SymErrorNum / length(OriginalData);
	BitErrorRate = BitErrorNum / (2 * length(OriginalData));
elseif M==8
    % Input Signal Normalization
    InputSignal = InputSignal -mean(InputSignal);
    InputSignal = 14/(max(InputSignal)-min(InputSignal))*(InputSignal-min(InputSignal))-7;
%   OriginalData = OriginalData - mean(OriginalData);

	if ~exist('threshold','var') || isempty(threshold)
		threshold = zeros(7, 1);
		threshold(4) = mean(InputSignal);
        threshold(2) = mean(InputSignal(find(InputSignal > threshold(4))));
        threshold(1) = mean(InputSignal(find(InputSignal > threshold(2))));
		threshold(3) = mean(InputSignal(find(threshold(4)<=InputSignal& InputSignal<= threshold(2))));
        
        threshold(6) = mean(InputSignal(find(InputSignal < threshold(4))));
        threshold(7) = mean(InputSignal(find(InputSignal < threshold(6))));
		threshold(5) = mean(InputSignal(find(threshold(6)<=InputSignal& InputSignal <= threshold(4))));

		% threshold = [6; 4; 2; 0; -2; -4; -6];
	end

	%% Input Signal Decision
	InputSignal(find(InputSignal > threshold(1))) = 7;
	InputSignal(find(InputSignal > threshold(2) & InputSignal <= threshold(1))) = 5;
	InputSignal(find(InputSignal > threshold(3) & InputSignal <= threshold(2))) = 3;
    InputSignal(find(InputSignal > threshold(4) & InputSignal <= threshold(3))) = 1;
    
    InputSignal(find(InputSignal > threshold(5) & InputSignal <= threshold(4))) = -1;
	InputSignal(find(InputSignal > threshold(6) & InputSignal <= threshold(5))) = -3;
    InputSignal(find(InputSignal > threshold(7) & InputSignal <= threshold(6))) = -5;
	InputSignal(find(InputSignal <= threshold(7))) = -7;

	%% Error Counting
    SymErrorNum = length(find(InputSignal ~= OriginalData));
    deco_input=pamdemod(InputSignal,M);
    deco_input=de2bi(deco_input,'left-msb');
    deco_input=deco_input';
    deco_input1=reshape(deco_input,1,3*length(InputSignal));
    deco_origidata =pamdemod(OriginalData,M);
    deco_origidata=de2bi(deco_origidata,'left-msb');
    deco_origidata=deco_origidata';
    deco_origidata1=reshape(deco_origidata,1,3*length(InputSignal));
	% The bit error number is 1/sym when only 1 bit changes and 2/sym when both bits change.
    BitErrorNum = length(find(deco_input1 ~= deco_origidata1));
    SymErrorRate = SymErrorNum / length(OriginalData);
	BitErrorRate = BitErrorNum / length(deco_origidata1);
end
        