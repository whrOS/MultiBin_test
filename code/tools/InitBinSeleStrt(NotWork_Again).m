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
if sum(EC(:)) >= Payload
    % EC of hist-bin-selection strategy with min{A,B} larger than payload
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
    
else
    % EC of hist-bin-selection strategy with min{A,B} not larger than payload
    CombArr = sortrows(CombArr',1)';
    EC = CombArr(6:7, :);
    ED = CombArr(8:9, :);
    Indic = ~~CombArr(4:5, :);
    
    %
    EcExceed = Payload-sum(EC(Indic));
    EC1 = EC(Indic)';
    ED1 = ED(Indic)';
    EC2 = EC(~Indic)';
    ED2 = ED(~Indic)';
    
    deltaED = ED2-ED1;
    deltaEC = EC2-EC1;
    deltaEC(deltaEC==0) = -0.0001;
    deltaDC = (deltaED) ./ (deltaEC);
    
    Idx2 = deltaEC>-0.0001;
    
    %     CombArr = [SortIdx; Val; Idx; Indic; HABC; HABD; HAB];
    CombArr = [SortIdx; deltaDC; HABC; HABD; Indic];
    CombArr = CombArr(:,Idx2);
    CombArr = sortrows(CombArr',2)';
    
    ec = 0;
    for i = 1 : 1 : numel(CombArr(1,:))
        ec = ec + deltaEC(CombArr(1,i));
        %         aa = sum(EC(Indic));
        Indic(:,CombArr(1,i)) = ~Indic(:,CombArr(1,i));
        %         [sum(EC(Indic))-aa, deltaEC(CombArr(1,i))]
        if ec >= EcExceed
            break
        end
    end
%     sum(EC(Indic))
    % 直接调整容量增加部分，容量超出很多；
    % 接下来，反向调整容量减少部分，将容量将将大于Payload即可
    % 假设128个直方图都要参与？
    %     Indic = ~~CombArr(4:5, :);
    EcExceed = sum(EC(Indic)) - Payload;
    EC1 = EC(Indic)';
    ED1 = ED(Indic)';
    EC2 = EC(~Indic)';
    ED2 = ED(~Indic)';
    
    deltaED = ED1-ED2;
    deltaEC = EC1-EC2;
    deltaEC(deltaEC==0) = -0.0001;
    deltaDC = (deltaED) ./ (deltaEC);
    
    Idx2 = deltaEC>-0.0001;
    
    CombArr = [SortIdx; deltaDC; HABC; HABD; Indic];
    CombArr = CombArr(:,Idx2);
    CombArr = sortrows(CombArr',-2)';
    
    ec = 0;
    for i = 1 : 1 : numel(CombArr(1,:))
        ec = ec + deltaEC(CombArr(1,i));
%         aa = sum(EC(Indic));
        Indic(:,CombArr(1,i)) = ~Indic(:,CombArr(1,i));
%         [ec EcExceed aa - sum(EC(Indic))  deltaEC(CombArr(1,i))]
        if ec >= EcExceed
            ec = ec - deltaEC(CombArr(1,i));
            Indic(:,CombArr(1,i)) = ~Indic(:,CombArr(1,i));
            %             break
        end
    end
    InitSele = ~Indic(1,:);
    InitEC = sum(EC(Indic));
end

end