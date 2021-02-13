function Volterra_in = Input_gen(data,ChanLen1st,ChanLen2nd)

Volterra_in_first = data;
t = 1;
if ChanLen2nd ~= 0
    for k = 1 : ChanLen2nd
        for m = k : ChanLen2nd
            Volterra_in_second(t) = data(floor(ChanLen1st/2)-floor(ChanLen2nd/2)+k)*data(floor(ChanLen1st/2)-floor(ChanLen2nd/2)+m);
            t=t+1;
        end
    end
    Volterra_in=[Volterra_in_first,Volterra_in_second];
else
    Volterra_in=Volterra_in_first;
end

end