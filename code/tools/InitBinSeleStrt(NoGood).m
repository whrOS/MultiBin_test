function [InitSele, InitEC] = InitBinSeleStrt(HABC, HABD, HAB, Payload)

[Val, Idx] = min(HAB);
SortIdx = 1 : 1 : numel(HAB(1, :));

Indic = false(2,numel(HAB(1, :))); % Indicator vector to select bins
for i = 1 : 1 : numel(HAB(1, :))
    Indic(Idx(i), i) = true;
end

CombArr = [SortIdx; Val; Idx; Indic; HABC; HABD; HAB];
CombArr = sortrows(CombArr',2)';

Idx   = CombArr(3, :);
Indic = ~~CombArr(4:5, :);
EC    = CombArr(6:7, :);

EC = EC(Indic);
if sum(EC(:)) < Payload
    % EC of hist-bin-selection strategy with min{A,B} not larger than payload
    SortedHAB = CombArr(10:11, :);
    Val2 = SortedHAB(~Indic)';
    CombArr(2, :) = Val2;
    CombArr = sortrows(CombArr',2)';
    
    Idx2   = CombArr(3, :);
    Indic2 = ~~CombArr(4:5, :);
    EC2    = CombArr(6:7, :);
    
    curEC = sum(EC2(Indic2));
    for i = 1 : 1 : numel(HAB(1, :))
        Indic2(:, i) = ~Indic2(:, i);
        newEC = sum(EC2(Indic2));
        if newEC <= curEC
            Indic2(:, i) = ~Indic2(:, i);
        else
            curEC = newEC;
            Idx2(i) = find(Indic2(:, i));
            if newEC >= Payload
                break;
            end
        end
    end
    
    EC = EC2(Indic2);
    Idx = Idx2;
end

ec       = 0;
InitSele = zeros(1, numel(HAB(1, :))) - 1;
for i = 1 : 1 : numel(HAB(1, :))
    InitSele(i) = Idx(i) - 1;
    ec = ec + EC(i);
    if ec >= Payload
        break
    end
end

InitEC = ec;

SortIdx = CombArr(1, :);
CombArr = [SortIdx; InitSele; EC'];
CombArr = sortrows(CombArr',1)';
InitSele = CombArr(2, :);

% EC = CombArr(3, :);
% Indic = false(2,numel(HABC(1, :))); % Indicator vector to select bins
% for i = 1 : 1 : numel(HABC(1, :))
%     Indic(Idx(i), i) = true;
% end
% sum(HABC(Indic))

% ec = 0;
% for i = 1 : 1 : numel(HAB(1, :))
%     if InitSele(i) == 0
%         ec = ec + HABC(1,i);
%     elseif InitSele(i) == 1
%         ec = ec + HABC(2,i);
%     elseif InitSele(i) == -1
% %         ec = ec;
%     end
% end

end