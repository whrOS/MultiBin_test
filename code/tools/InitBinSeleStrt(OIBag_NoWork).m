function [InitSelec, InitEC] = InitBinSeleStrt(HABC, HABD, HAB, Payload)

[Val, Idx] = min(HAB);
% SortIdx = 1 : 1 : numel(HAB(1, :));

Indic = false(2,numel(HAB(1, :))); % Indicator vector to select bins
for i = 1 : 1 : numel(HAB(1, :))
    Indic(Idx(i), i) = true;
end

EC = HABC(Indic);
ED = HABD(Indic);

InitSelec = false(1,numel(HAB(1, :)));

[ed, InitSelec] = OIBag(EC, ED, InitSelec, numel(HAB(1, :)), Payload);

InitEC = EC(InitSelec);

end


function [ed, Selec] = OIBag(EC, ED, Selec, k, payload)
    if k == 0
        ed = 0;
        return;
    end
    
    if payload < EC(k)
        ed = OIBag(EC, ED, Selec, k-1, payload);
        return;
    else
        ed1 = OIBag(EC, ED, Selec, k-1, payload - EC(k)) + ED(k);
        ed2 = OIBag(EC, ED, Selec, k-1, payload);
        if ed1 < ed2
            Selec(k) = true;
            ed = ed1;
        else
            Selec(k) = false;
            ed = ed2;
        end
        return;
    end
end
