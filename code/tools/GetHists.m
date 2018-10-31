function [MHs, k, m, HStep, NL, Tlog] = GetHists(I, a, b, Tmax)
% Inputs :
% I     : the cover image
% a,b   : one block size of a * b
% Tmax  : number of blocks with NL<Tmax larger than Ratio from Function:GetTmax()

% Outputs : 
% MHs   : cell(1,<=128), to get 128 histograms Hs by merging some small histograms
% k     : double, to merge first k histograms as one histogram -> parameter to embed
% m     : double, to merge each HStep histograms after first k ones as one histogram,
%        in total of m merged histograms -> parameter to embed
% HStep : double ,merge each HStep histograms to one histogram
% NL    : [floor((A-2)/a), floor((B-2)/b)], NoiseLevel of each block

[A, B] = size(I);
N = floor((A-2)/a) * floor((B-2)/b); % Total number of blocks

Hs = cell(1,Tmax);
for i = 1:Tmax
    Hs{i} = zeros(1,256);
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
                dmax = Y(a*b) - Y(a*b-1) - 1;
            else
                dmax = Y(a*b) - Y(a*b-1);
            end
            % min
            if In(2) < In(1)
                dmin = Y(2) - Y(1) - 1;
            else
                dmin = Y(2) - Y(1);
            end
            Hs{NL(i,j)+1}(dmax+1) = Hs{NL(i,j)+1}(dmax+1) + 1;
            Hs{NL(i,j)+1}(dmin+1) = Hs{NL(i,j)+1}(dmin+1) + 1;
        end
    end
end

MHs = cell(1,128);
for i = 1:128
    MHs{i} = zeros(1,256);
end

Tlog = [];

% hist 1
cnt = 1;
NLs = unique(NL(:));
for i = 1 : 1 : numel(NLs)
    if NLs(i) < Tmax
        T = NLs(i)+1;
        MHs{cnt} = MHs{cnt} + Hs{T};
        if sum(MHs{cnt}(:)) >= 2*N/128
            k = i;
            break;
        end
    end
end
Tlog = [Tlog, T];

NONum =  sum(NLs<Tmax); % number of histograms which are non-empty

HStep = ceil((NONum - k) / 127);
m = ceil((NONum - k - 127)/(HStep-1));

% HStep-merged hists
for i = k+1 : HStep : (k + m*HStep)
%     T = NLs(i)+1;
    cnt = cnt + 1;
    for j = 1 : 1 : HStep
        T = NLs(i+j-1)+1;
        MHs{cnt} = MHs{cnt} + Hs{T};
    end
    Tlog = [Tlog, T];
end

% no merged hists
for i = (k + m*HStep)+1 : 1 : NONum
    T = NLs(i)+1;
    cnt = cnt + 1;
    MHs{cnt} = Hs{T};
    Tlog = [Tlog, T];
end

% delete the unused histogram space
if cnt < 128
    for i = cnt+1 : 128
        MHs(i) = [];
    end
end

% cnt = 0;
% for i = 1 : numel(MHs)
%     cnt = cnt + sum(MHs{i}(:));
% end
% Ratio = cnt / (2*N)


end