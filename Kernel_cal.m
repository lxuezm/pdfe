function [KernelSize, Maxkernel] = Kernel_cal(ChanLen1st,ChanLen2nd,ChanLen3rd)
KernelSize = ChanLen1st;% The length of first order
Kernel2ndSize = 0;% The length of second order
Kernel3rdSize = 0;% The length of third order
if ChanLen2nd ~= 0
    for k = 1 : ChanLen2nd
        for m = k : ChanLen2nd
            Kernel2ndSize = Kernel2ndSize + 1;
        end
    end
    KernelSize = KernelSize + Kernel2ndSize;
end
if ChanLen3rd ~= 0
    for k = 1 : ChanLen3rd
        for m = k : ChanLen3rd
            for n = m : ChanLen3rd
                Kernel3rdSize = Kernel3rdSize + 1;
            end
        end
    end
    KernelSize = KernelSize + Kernel3rdSize;
end
    Maxkernel=max([ChanLen1st,ChanLen2nd,ChanLen3rd]);
end
