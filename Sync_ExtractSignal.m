function [ExtractedSignalUS, OriginalSignalUS] = Sync_ExtractSignal(SampledData, OriginalData, OverSamplingRatio, UpSamplingRatio)
	% This function performs the synchronization and signal extraction for the 
	% sampled signal from DSO. The %SampledData% will be down sampled by the 
	% %OverSamplingRatio%/%UpSamplingRatio% and correlated with the %OriginalData%
	% upsampled by %UpSamplingRatio%. Then the index of max correlating coeffients 
	% will be found and the synchronized data will be extracted.
	%
	% input:
	% 		SampledData
	%				The data sampled by DSO.
	%		  OriginalData
	%				The data which is transmitted with 1 symbol per second.
	%		  OverSamplingRatio
	%				The ratio of down sampling progress, which is usually the DSO sample 
	% 			rate dividing signal symbol rate
	%		  UpSamplingRatio (optional)
	%				The ratio of upsampling progress.
	%				Default value: 1
	% output:
	%     ExtractedSignalUS
	%       The synchronized and extracted signal upsampled by %UpSamplingRatio%.
	%       Size: length(OriginalData)*UpSamplingRatio, 1
	%     OriginalSignalUS
	%       The original signal upsampled by %UpSamplingRatio%.
	%       Size: length(OriginalData)*UpSamplingRatio, 1
	
	% Parameters checking
	narginchk(3,4);

	if ~exist('UpSamplingRatio','var') || isempty(UpSamplingRatio)
% 		UpSamplingRatio = 1;
      UpSamplingRatio = 4;
	end
	
	% Downsampling
% 	DownSampledData = SampledData(11:OverSamplingRatio/UpSamplingRatio:end, 1);
	DownSampledData = SampledData(1:end);
	% Preparing transmitted data
% 	OriginalDataRemapped = (OriginalData - 1.5) * 2;
    OriginalDataRemapped = OriginalData;
    OriginalSignalUS  = upsample(OriginalDataRemapped,4);
% 	OriginalSignalUS = reshape(repmat(OriginalDataRemapped, 1, UpSamplingRatio)', UpSamplingRatio * numel(OriginalDataRemapped), 1);
	% Correlation
	CorrelationResult = conv(DownSampledData(1:end), conj(OriginalSignalUS(end:-1:1)));
	[MaxCorr, index] = max(CorrelationResult);
    figure
    plot(CorrelationResult);
	ExtractedSignalUS = DownSampledData(index-length(OriginalData)*UpSamplingRatio+1 : index);

    Ex_spec1 = ExtractedSignalUS ;
    Npoints1 = length(Ex_spec1);
    FFT_Ex_1 = fftshift(fft(Ex_spec1));
    FFT_Ex1 = abs(FFT_Ex_1)./(length(Ex_spec1));
    Frek1 = (80e9*(-(Npoints1)/2:((Npoints1/2)-1)))/Npoints1;
    plot(Frek1./1e9,10*log10(FFT_Ex1.^2));
    title('Spectrum of original signal');
%     ExtractedSignalUS = DownSampledData(index+1 : index+length(OriginalData)*UpSamplingRatio);
    ExtractedSignalUS=ExtractedSignalUS(1:UpSamplingRatio:end);
    OriginalSignalUS=OriginalSignalUS(1:UpSamplingRatio:end);
% 	ExtractedSignalUS = SampledData(index*16-length(OriginalData)*OverSamplingRatio+1 : index);
%     nbits=floor(length(SampledData)/OverSamplingRatio);
%     for i=1:OverSamplingRatio
%         Re_SampledData(:,i)=SampledData(i:OverSamplingRatio:OverSamplingRatio*nbits);
%     end
%     ExtractedSignal= Re_SampledData(index-length(OriginalData)+1:index,:)';
%     ExtractedSignalUS=reshape(ExtractedSignalUS ,524272,1);
%     ExtractedSignalUS= ExtractedSignalUS(1:UpSamplingRatio:end,1);
    
    
    

    