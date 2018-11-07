function [HABC, HABD, HAB] = GetHRatioDisto(Hs)

HNum = numel(Hs);
AEC = zeros(1, HNum);AED = zeros(1, HNum);
BEC = zeros(1, HNum);BED = zeros(1, HNum);
% CEC = zeros(1, HNum);CED = zeros(1, HNum);
A = zeros(1, HNum);
B = zeros(1, HNum);
% C = zeros(1, HNum);

SetECO = 0.0001; % solution to the situation of ec = 0 => ed/ec = Inf

for i = 1 : 1 : HNum
   H = Hs{i};
   AEC(i)  = H(1);   AED(i)  = sum(H(2:end)) + 0.5*H(1);
   BEC(i)  = H(2);   BED(i)  = sum(H(3:end)) + 0.5*H(2);
%    CEC(i)  = H(4);   CED(i)  = sum(H(5:end)) + 0.5*H(4);
   if AEC(i) == 0
       AEC(i) = SetECO;
   end
   if BEC(i) == 0
       BEC(i) = SetECO;
   end
%    if CEC(i) == 0
%        CEC(i) = SetECO;
%    end
   A(i) = AED(i)/AEC(i);
   B(i) = BED(i)/BEC(i);
%    C(i) = CED(i)/CEC(i);
end

% HABC = [AEC; BEC; CEC];
% HABD = [AED; BED; CED];
% HAB = [A; B; C];
HABC = [AEC; BEC];
HABD = [AED; BED];
HAB = [A; B];

end
