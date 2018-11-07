% Embedding once with a given embedding payload for one image
%% Prepare
clc;
addpath(genpath('Origin Images')); addpath(genpath('tools')); addpath(genpath('result'));
Imgs = {'Lena', 'Baboon', 'Airplane', 'Barbara', 'Lake', 'Peppers', 'Boat', 'Elaine'};

%% Parameters
HNum = 128;
EdgInfoSize = HNum*log2(3) + 7 + 7 + 2 + 11;      % 128 + k + m + HStep + MapId

%
% Payload     = 38000;% + EdgInfoSize;  % payload + length of edge info
NLmax       = 1000;                 % only consider the blocks with NL < NLmax
TRatio      = 0.999;                % 90% blocks utilized to embed data


for IIdx = 1 : 1
    % Read cover image
    IName   = Imgs{IIdx};
    ISavStr1D = ['Proposed1D_2019_',IName,'.mat']
%     ISavStr2D = ['Proposed2D_2019_',IName,'.mat']
%     ISavStrPairwise = ['ProposedPairwise_2019_',IName,'.mat']
    I       = double(imread([IName,'.bmp']));
    [A, B]  = size(I);
    
    cnt = 0;
    Seles = {};
    R = zeros(8, 100000); % [EC PSNR a b T k m Hstep]
%     R2 = zeros(9, 100000); % [EC PSNR a b T k m Hstep]
    for Payload = 10000 : 1000 : 10000
        [Payload]
        maxPSNR = 0;
        for a = 2:5
            for b = 2:5
                %% Step 1 : Get <=128 merged histograms
                [Tmax] = GetTmax(I, a, b, TRatio, NLmax);
                [Hs, k, m, HStep, NL, Tlog] = GetHists(I, a, b, Tmax, HNum);
%                 [Hs, k, m, HStep, NL, Tlog] = GetHistsPairwise(I, a, b, Tmax, HNum);
                
                % test show
                IsShow = 0;
                for tmp = 1:IsShow
                    Bins = 20;
                    HNum = numel(Hs); H = zeros(HNum, 256);
                    for i = 1:1:HNum
                        H(i,:) = Hs{i};
                    end
                    x = 1:1:HNum;
                    y = 1:1:Bins;
                    xbins = 0:1:HNum;
                    ybins = 0:1:Bins;
                    bar3(H(x,y),1);
                    set(gca,'XTickLabel',xbins);
                    set(gca,'YTickLabel',ybins);
                    % axis([1 HistNum 1 Bins 0 max(H(:))+100])
                    % view([90,0])
                end
                
                %% Step 2 : Bins-selection strategy
                [HABC, HABD, HAB] = GetHRatioDisto(Hs); % HABC = [AEC BEC]; HABD = [AED BED]; HAB = [A B]
%                 [HABC, HABD, HAB] = GetHRatioDistoPairwise(Hs); % HABC = [AEC BEC]; HABD = [AED BED]; HAB = [A B]
                
                [v, ~] = max(HABC);
                MAXEC = sum(v);
                %         fprintf('\n-----Max Embedding Capacity : %d. \n', MAXEC);
                if MAXEC < Payload
                    fprintf('\n-----No Enough Embedding Capacity. \n');
                    continue
                end
                
                % Step 2.1 : DP to get minimized ED/EC
                [InitSele, InitEC, InitED] = InitBinSeleStrt(HABC, HABD, HAB, Payload);
%                 [InitSele, InitEC, InitED] = InitBinSeleStrtPairwise(HABC, HABD, HAB, Payload);
                
                %% Step 3 : Embedding Process
                [EC, PSNR, T] = EmbedProc(I, a, b, InitSele, NL, Tlog, Payload);
% % %                 [EC, PSNR, T] = EmbedProc_test(I, a, b, InitSele, NL, Tlog, Payload);
% % %                 
% % %                 InitSele(:) = 0; [~, ~, ~] = EmbedProc(I, a, b, InitSele, NL, Tlog, Payload);
%                 [EC, PSNR, T] = EmbedProcPairwise(I, a, b, InitSele, NL, Tlog, Payload);
                
                if PSNR > maxPSNR
                    maxPSNR = PSNR;
                    OptSele = InitSele;
                    OptA = a;
                    OptB = b;
                    OptT = T;
                    OptK = k;
                    OptM = m;
                    OptHStep = HStep;
                    OptNL = NL;
                    OptTlog = Tlog;
                end
            end
        end
        if maxPSNR == 0
            break
        end
        %% 2D
%         [EC, PSNR, T, MapId] = EmbedProc2DCombined(I, OptA, OptB, OptSele, OptNL, OptTlog, Payload);
        %%
        cnt = cnt + 1;
        Seles{cnt} = OptSele;
        R(:,cnt) = [Payload, maxPSNR, OptA, OptB, OptT, OptK, OptM, OptHStep]'; % [EC PSNR a b T k m Hstep MapId]
%         [EC, PSNR, T, MapId] = EmbedProc2D(I, OptA, OptB, OptSele, OptNL, OptTlog, Payload);
%         R2(:,cnt) = [Payload, maxPSNR, OptA, OptB, OptT, OptK, OptM, OptHStep, MapId]'; % [EC PSNR a b T k m Hstep MapId]
    end
%     cnt = sum(R(2,:)~=0);
%     res = R(:, 1:cnt);
%     save(ISavStrPairwise, 'res');
    cnt = sum(R(2,:)~=0);
    res = R(:, 1:cnt);
    save(ISavStr1D, 'res');
%     cnt = sum(R2(2,:)~=0);
%     res = R2(:, 1:cnt);
%     save(ISavStr2D, 'res');
end