function [output] = MLSD(input,BitPerSym,alpha)

h=modem.pammod('M', 2^BitPerSym, 'SymbolOrder', 'gray', 'InputType', 'Bit');
symbol = h.Constellation;
Len = length(input);
M = length(symbol);
Distance = zeros(M,Len);
Index = zeros(M,Len);

Distance(:,1) = abs(input(1) - symbol.').^2;
for ii = 2:Len
    for jj = 1:M
        Dis = zeros(1,M);
        for kk = 1:M
            Dis(kk) = Distance(kk,ii-1) + abs(input(ii) - symbol(jj) - alpha*symbol(kk)).^2;
        end
        [val,pos] = min(Dis);
        Distance(jj,ii) = val;
        Index(jj,ii) = pos;
    end
end

output = zeros(1,Len);
[~,pos] = min(Distance(:,Len));
output(Len) = symbol(pos);
for ii = Len:-1:2
    pos = Index(pos,ii);
    output(ii-1) = symbol(pos);
end

end

