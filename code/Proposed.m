% Embedding once with a given embedding payload for one image
%% Prepare
clc;
addpath(genpath('Origin Images')); addpath(genpath('tools')); addpath(genpath('result'));
Imgs = {'Lena', 'Baboon', 'Airplane', 'Barbara', 'Lake', 'Peppers', 'Boat', 'Elaine'};

%% Parameters
EdgInfoSize = 128 + 7 + 7 + 2;      % 128 + k + m + HStep
IIdx        = 3;                    % image index pof Imgs{}
%
Payload     = 10000 + EdgInfoSize;  % payload + length of edge info
TRatio      = 0.9;                  % 90% blocks utilized to embed data
a = 2;    b = 2;                    % block size : a * b

% Read cover image
IName   = Imgs{IIdx};
ISavStr = ['Proposed_2019_',IName,'.mat'];
I       = double(imread([IName,'.bmp']));
[A, B]  = size(I);

%% Step 1 : Get <=128 merged histograms
[Tmax] = GetTmax(I, a, b, TRatio);
[Hs, k, m, HStep, NL, Tlog] = GetHists(I, a, b, Tmax);

% test show
for tmp = 1
    HistNum = numel(Hs); H = zeros(HistNum, 256);    Bins = 20;
    for i = 1:1:HistNum  
        H(i,:) = Hs{i}; 
    end
    x = 1:1:HistNum; 
    y = 1:1:Bins; 
    xbins = 0:1:HistNum; 
    ybins = 0:1:Bins;
    bar3(H(x,y),1);
    set(gca,'XTickLabel',xbins); 
    set(gca,'YTickLabel',ybins); 
    % axis([1 HistNum 1 Bins 0 max(H(:))+100]) 
    % view([90,0])
end

%% Step 2 : Bins-selection strategy
% Step 2.1 : Initial bins-selection strategy
% Step 2.2 : Iteration to optimize the strategy
%% Step 3 : Embedding Process