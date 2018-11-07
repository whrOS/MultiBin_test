function [ec, PSNR, T, OptMapID] = EmbedProc2DCombined(I, a, b, Sele, NL, Tlog, Payload)
%%
[A, B] = size(I);
maps = {};
map = cell(3);
mapECs = {};
mapEDs = {};

[mapEDs, mapECs, maps] = getOneMap(mapEDs, mapECs, maps,map,[0,0]);
mapNum = numel(maps);

%%
Tmax = Tlog(find(Sele~=0, 1, 'last' ));

TSele = zeros(1, Tmax);
Cnt = 1;
for i = 1 : 1 : Tmax
    TSele(i) = Sele(Cnt);
    if i >= Tlog(Cnt)
        Cnt = Cnt + 1;
    end
end

Hs0 = cell(1,Tmax); % T = 0 - Tmax
Hs1 = cell(1,Tmax); % T = 0 - Tmax
for i = 1:Tmax
    Hs0{i} = zeros(256);
    Hs1{i} = zeros(256);
end
% EC = zeros(floor((A-2)/a),floor((B-2)/b));
% ED = zeros(floor((A-2)/a),floor((B-2)/b));
for i = 1:floor((A-2)/a)
    for j = 1:floor((B-2)/b)
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
            
            if TSele(NL(i,j)+1) == 1
                Hs0{NL(i,j)+1}(dmax+1, dmin+1) = Hs0{NL(i,j)+1}(dmax+1, dmin+1) + 1;
            elseif TSele(NL(i,j)+1) == 2
                Hs1{NL(i,j)+1}(dmax+1, dmin+1) = Hs1{NL(i,j)+1}(dmax+1, dmin+1) + 1;
            end
            
%             if TSele(NL(i,j)+1) == 1
%                 if dmax == 0
%                     EC(i,j) = EC(i,j) + 1;
%                     ED(i,j) = ED(i,j) + 0.5;
%                 else
%                     ED(i,j) = ED(i,j) + 1;
%                 end
%                 if dmin == 0
%                     EC(i,j) = EC(i,j) + 1;
%                     ED(i,j) = ED(i,j) + 0.5;
%                 else
%                     ED(i,j) = ED(i,j) + 1;
%                 end
%             elseif TSele(NL(i,j)+1) == 2
%                 if dmax == 1
%                     EC(i,j) = EC(i,j) + 1;
%                     ED(i,j) = ED(i,j) + 0.5;
%                 else
%                     if dmax > 1
%                         ED(i,j) = ED(i,j) + 1;
%                     end
%                 end
%                 if dmin == 1
%                     EC(i,j) = EC(i,j) + 1;
%                     ED(i,j) = ED(i,j) + 0.5;
%                 else
%                     if dmin > 1
%                         ED(i,j) = ED(i,j) + 1;
%                     end
%                 end
%             end
            
        end
    end
end

for i = 2:Tmax
    Hs0{i} = Hs0{i} + Hs0{i-1};
    Hs1{i} = Hs1{i} + Hs1{i-1};
end

%% 计算总的失真和容量
optMapEC0 = zeros(256,256); optMapEC0(1,:) = 1;  optMapEC0(:,1) = 1;  optMapEC0(1,1) = 2;
optMapED0 = ones(256,256)*2;optMapED0(1,:) = 3/2;optMapED0(:,1) = 3/2;optMapED0(1,1) = 1;

optMapEC1 = zeros(256,256); optMapEC1(2,:) = 1;  optMapEC1(:,2) = 1;  optMapEC1(2,2) = 2;
optMapED1 = ones(256,256)*2;optMapED1(2,:) = 3/2;optMapED1(:,2) = 3/2;optMapED1(2,2) = 1;
optMapED1(:,1) = 1; optMapED1(1,:) = 1; optMapED1(1,1) = 0; optMapED1(2,1) = 0.5; optMapED1(1,2) = 0.5;

EC0 = sum(sum(Hs0{Tmax} .* optMapEC0)); % 7045
ED0 = sum(sum(Hs0{Tmax} .* optMapED0)); % 9507.5
EC1 = sum(sum(Hs1{Tmax} .* optMapEC1)); % 2967
ED1 = sum(sum(Hs1{Tmax} .* optMapED1)); % 4697.5


%% 先记录非修按规定部分的失真
optMapEC0 = zeros(256,256); optMapEC0(1,:) = 1;  optMapEC0(:,1) = 1;  optMapEC0(1:3,1:3) = 0;
optMapED0 = ones(256,256)*2;optMapED0(1,:) = 3/2;optMapED0(:,1) = 3/2;optMapED0(1:3,1:3) = 0;

optMapEC1 = zeros(256,256); optMapEC1(2,:) = 1;  optMapEC1(:,2) = 1;  optMapEC1(2:4,2:4) = 0;
optMapED1 = ones(256,256)*2;optMapED1(2,:) = 3/2;optMapED1(:,2) = 3/2;optMapED1(2:4,2:4) = 0;
optMapED1(:,1) = 1; optMapED1(1,:) = 1; optMapED1(1,1) = 0; optMapED1(2,1) = 0.5; optMapED1(1,2) = 0.5;

ec0 = sum(sum(Hs0{Tmax} .* optMapEC0));
ed0 = sum(sum(Hs0{Tmax} .* optMapED0));
ec1 = sum(sum(Hs1{Tmax} .* optMapEC1));
ed1 = sum(sum(Hs1{Tmax} .* optMapED1));

ec0 = ec0 + ec1;
ed0 = ed0 + ed1;

%% 调优loacl
H2D = Hs0{Tmax}(1:3,1:3) + Hs1{Tmax}(2:4,2:4);
%%
OptEC = Payload;
OptED = 10000000;
OptMapID = 0;
log = [];
for m = 1 : mapNum
    ec = mapECs(m);
    ed = mapEDs(m);
    ec = ec{1}; optMapEC = flipud(rot90(ec'));
    ed = ed{1}; optMapED = flipud(rot90(ed'));
    ECMT = sum(sum(H2D .* optMapEC)) + ec0;
    EDMT = sum(sum(H2D .* optMapED)) + ed0;
    log = [log; [m ECMT EDMT]];
    if ECMT >= Payload && EDMT < OptED
%         OptEC = ECMT;
        OptED = EDMT;
        OptMapID = m;
    end
end

T = Tmax;
PSNR = 10*log10(A*B*255^2 / OptED);

end