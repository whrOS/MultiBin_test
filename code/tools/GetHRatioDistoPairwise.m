function [HABC, HABD, HAB] = GetHRatioDistoPairwise(Hs)

optMapEC0 = zeros(256,256); optMapEC0(1,:) = 1;  optMapEC0(:,1) = 1;  optMapEC0(1,1) = log2(3); optMapEC0(2,2) = 1;
optMapED0 = ones(256,256)*2;optMapED0(1,:) = 3/2;optMapED0(:,1) = 3/2;optMapED0(1,1) = 2/3;     optMapED0(2,2) = 1;

optMapEC1 = zeros(256,256); optMapEC1(2:end, 2:end) = optMapEC0(1:255,1:255); optMapEC1(2,1) = 1;  optMapEC1(1,2) = 1;
optMapED1 = zeros(256,256); optMapED1(2:end, 2:end) = optMapED0(1:255,1:255);
optMapED1(3:end,1) = 1; optMapED1(1,3:end) = 1; optMapED1(2,1) = 0.5; optMapED1(1,2) = 0.5;

HNum = numel(Hs);
AEC = zeros(1, HNum);
AED = zeros(1, HNum);
BEC = zeros(1, HNum);
BED = zeros(1, HNum);
A = zeros(1, HNum);
B = zeros(1, HNum);

SetECO = 0.0001; % solution to the situation of ec = 0 => ed/ec = Inf

for i = 1 : 1 : HNum
   H = Hs{i};
   AEC(i)  = sum(sum(H.*optMapEC0));
   AED(i)  = sum(sum(H.*optMapED0));
   BEC(i)  = sum(sum(H.*optMapEC1));
   BED(i)  = sum(sum(H.*optMapED1));
   if AEC(i) == 0
       AEC(i) = SetECO;
   end
   if BEC(i) == 0
       BEC(i) = SetECO;
   end
   A(i) = AED(i)/AEC(i);
   B(i) = BED(i)/BEC(i);
end

HABC = [AEC; BEC];
HABD = [AED; BED];
HAB = [A; B];

end
