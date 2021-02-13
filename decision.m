function Output=decision(input)
% hard decision of the input sequences
if input<-2
    Output = -3;
elseif -2<=input&&input<0
    Output = -1;
elseif 0<=input&&input<2
    Output = 1;
elseif input>=2
    Output = 3;
else
    Output=100;
end
end