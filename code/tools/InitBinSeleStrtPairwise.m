function [Sele, InitEC, InitED] = InitBinSeleStrtPairwise(HABC, HABD, HAB, Payload)

% Indic = false(2,numel(HAB(1, :))); % Indicator vector to select bins

% Indic(1,:) = true;
[Val, Idx] = max(HABD);

EDmin = 0; % 一个都不选
EDmax = sum(Val(:));

Val = 0;    % 在失真为mid时候，最大的嵌入容量。
            % 如果最大嵌入容量大于Payload，则在小于当前失真mid侧再次寻找
            % 如果最大嵌入容量小于Payload，则在大于大嵌失真mid侧再次寻找
            % 这个值应该是分组背包的返回值

A = HABC; % A表示价值
B = HABD; % B表示体积
            
l = EDmin;
r = EDmax;
InitSele = false(2,numel(HAB(1, :)));
while(l+1 < r)
    mid = ceil((l+r) / 2);
    [Val, Indic] = DP(A,B,mid);
            
    if (Val >= Payload)
        r = mid;
        InitSele = Indic;
    else
        l = mid;
    end
%     ec = HABC(Indic);
%     ed = HABD(Indic);
%     PSNR = 10*log10(512*512*255^2 / sum(ed));
%     [sum(ec), sum(ed), mid, PSNR]
end

% InitSele = Indic;
InitEC = sum(HABC(InitSele));
InitED = sum(HABD(InitSele));
% [InitEC InitED 10*log10(512*512*255^2 / InitED)]

Sele = zeros(1,numel(HAB)/2);
for i = 1:1:numel(HAB)/2
    if InitSele(1,i) == 1
        Sele(i) = 1;
    end
    if InitSele(2,i) == 1
        Sele(i) = 2;
    end
end

end

function [Val, Indic] = DP(A, B, M)

[~,N] = size(A);
A = A';
B = floor(B');
f = zeros(1, M);

% Indic = false(1,N)';
mark = zeros(N, M);

for i = 1 : 1 : N
    for j = M : -1 : 1
        for k = 1 : 1 : 2
            if B(i,k) < j
                if f(j) < f(j-B(i,k))+A(i,k)
                    f(j) = f(j-B(i,k))+A(i,k);
                    mark(i,j) = k;
                end
            elseif B(i,k) == j
                 if f(j) < A(i,k)
                    f(j) = A(i,k);
                    mark(i,j) = k;
                end
            end
        end
    end
end


maxs = 0;  
for j = M : -1 : 1
    if maxs <= f(j)
        maxs = f(j);
        vj = j;
    end
end
j = vj;
v = zeros(1,N);
for i = N : -1 : 1
    if 0 ~= j
        v(i) = mark(i,j);
        if 0 ~= v(i)
            j = j - B(i,v(i));
        end
    end
end

Indic = false(2,N);
for i = 1 : 1 : N
    if v(i) ~= 0
        Indic(v(i), i) = true;
    end
end
% A = A';
% B = B';
Val = sum(A(Indic'));
% [sum(A(Indic')) sum(B(Indic')) M]

end

% sum([795;1483;2273;3183;3998;3920;3614;2174;1796;1467;1268;1025])
