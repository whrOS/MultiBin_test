function [Tout] = GetTmax(I, a, b, TRatio, NLmax)
% Function : get the Tout,  such that  leads to number of the block
%           satisfying NL<Tout has a ratio of TRatio in all blocks

[A, B] = size(I);
N = floor((A-2)/a) * floor((B-2)/b); % Total number of blocks
Tmax = NLmax;

Hs = cell(1,Tmax); % T = 0 - Tmax-1
for i = 1:Tmax
    Hs{i} = zeros(1,512);
end

NL = zeros(floor((A-2)/a),floor((B-2)/b));
for i = 1:floor((A-2)/a)
    for j = 1:floor((B-2)/b)
        for ii = 1:a+1
            for jj = 1:b+2
                if ii == a+1 || jj == b+1 || jj == b+2
                    NL(i,j) = NL(i,j) + abs(I(a*(i-1)+ii,b*(j-1)+jj) - I(a*(i-1)+ii+1,b*(j-1)+jj));
                end
            end
        end
        for ii = 1:a+2
            for jj = 1:b+1
                if ii == a+1 || ii == a+2 || jj == b+1
                    NL(i,j) = NL(i,j) + abs(I(a*(i-1)+ii,b*(j-1)+jj) - I(a*(i-1)+ii,b*(j-1)+jj+1));
                end
            end
        end
        if NL(i,j) < Tmax
            X = I(a*(i-1)+1:a*i,b*(j-1)+1:b*j);
            X = X(:);
            [Y, In] = sort(X);
            % max
            if In(a*b) < In(a*b-1)
                dmax = Y(a*b-1) - Y(a*b);
            else
                dmax = Y(a*b) - Y(a*b-1);
            end
            % min
            if In(2) < In(1)
                dmin = Y(1) - Y(2);
            else
                dmin = Y(2) - Y(1);
            end
            Hs{NL(i,j)+1}(dmax+256) = Hs{NL(i,j)+1}(dmax+256) + 1;
            Hs{NL(i,j)+1}(dmin+256) = Hs{NL(i,j)+1}(dmin+256) + 1;
        end
    end
end
Tout = 0;
for i = 2:Tmax
    Hs{i} = Hs{i} + Hs{i-1};
    Ratio = sum(Hs{i}(:)) / (2*N);
    if Ratio >= TRatio
        Tout = i;
        break;
    end
end
if Tout == 0
    Tout = Tmax;
end
% Ratio

end

