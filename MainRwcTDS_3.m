clear variables
clc; close all
%%
c = 299792458;
NbrPnts = 1e3;

Freq = linspace(0.1,1,NbrPnts)*1e12; % Hz
Lambda0 = c./Freq;
FreqTHz = Freq *1e-12;

%% Experimental Relative Permittivity of Pure Water     (w)
FileName = 'ComplexIndexWater.txt';
[FreqTHzW, nW,kW,AlphaW, FreqTHzInterpW, nWInterp, kWInterp, AlphaWInterp] =...
    LoadExpFreqEps_(FileName,FreqTHz);

%% Experimental Relative Permittivity of Dry Leaf       (s)
FileName = 'ComplexIndexDryLeaf.txt';
[FreqTHzS, nS,kS,AlphaS, FreqTHzInterpS,nSInterp, kSInterp, AlphaSInterp] =...
    LoadExpFreqEps_(FileName,FreqTHz);

%% Constant

% Sweep Thickness of Sample
NbrePtsL = 50;
L =LineSpaceM(250,350,NbrePtsL)*1e-6; % Sample Thickness [m]
L_um = L*1e6;                         % Sample Thickness [µm]

FPCoeff = 1 ; % Fabry–Pérot Coefficient ==> if (== 0) no FP Effect

RowIndex = 1;                           % Start read row
c = 299792458;                          % Light Speed [m/s]
N_air = 1.00027;                        % Refractive Index of Air

n1 = 1; k1 = 0;             % n & k of Medium 1 (Placed Before the Sample)
n3 = 1; k3 = 0;             % n & k of Medium 3 (Placed After the Sample)

% Complex Index
N2Func = @(n2,k2) n2-1i*k2; % Handle Fnc of Complex Refractive Index
N1 = n1-1i*k1;              % Complex Refractive Index of Medium 1
N3 = n3-1i*k3;              % Complex Refractive Index of Medium 3
% N2 = N2Func(n2,k2);         % Complex Refractive Index of Sample

%% Load Experimental Data TDS/Weight
[Data, FileNameFull] = LoadTimeDomainSignal_1(RowIndex,'mat');
if isempty(Data)  || iscell(Data)==0;    return; end
% Data(2800:end)=[]; % Sample 3
Data(7688:end)=[]; % Sample 11
NbreFile = numel(Data);

% Load Weight Data
[DryTime0, DateTimeOut0, Weight, RWC01, RWC02, PathName, FileName,WDry] =...
    LoadQuintixData_1(40);
DateTimeOut0 = DateTimeOut0*3600*24;

%% Read Scan Info File (Start
ScanInfoFileName = 'ScanInfo.txt';
ScanInfoFilePath = fullfile(PathName,ScanInfoFileName);

[~,DstFileName] = fileparts(FileNameFull) ;

ScanInfoFileNameNew = strcat(DstFileName(15:24),'-Dst-',ScanInfoFileName);
ScanInfoFileNamePathNew = fullfile(PathName,ScanInfoFileNameNew);
if exist(ScanInfoFilePath,'file') % If 
    ScanInfo = dlmread(ScanInfoFilePath);
elseif exist(ScanInfoFileNamePathNew,'file')
    ScanInfo = dlmread(ScanInfoFileNamePathNew);
else
    error('ScanInfo File Missing')
end
TdsStartTime = ScanInfo(1); NbrRepeats = ScanInfo(2);
ScanAvg = ScanInfo(3); WfmRate = ScanInfo(4);
[TdsStartTimeS, ~] = LabVIEW2MatlabDate (TdsStartTime);
%%
DelayStartTDS = 0.8;        % Delay to start TDS measurement
TimerTDS = ScanAvg/WfmRate; % Time of an average acquisition

DryTimeTds = (0:TimerTDS:0+(NbreFile-1)*TimerTDS)'; % Dry Time Axis for TDS measurement (s)
TdsTimeTrigV0 = DryTimeTds + TdsStartTimeS+DelayStartTDS; % Dry Time Axis for TDS measurement (days)
DryTimeTds = DryTimeTds/60; % Dry Time Axis for TDS measurement (mn)

iTdsWeightTrigg = zeros(NbreFile,1);
for i1=1:NbreFile
    [~,iTdsWeightTrigg(i1)] = min(abs(DateTimeOut0-TdsTimeTrigV0(i1)));
end

DryTime = DryTime0(iTdsWeightTrigg);
RWC1 = RWC01(iTdsWeightTrigg);
RWC2 = RWC02(iTdsWeightTrigg);

FormatOut = 'HH:MM:SS.FFF dd-mmm-yyyy';
TdsDateInDaysFormated = datestr(TdsTimeTrigV0/(3600*24), FormatOut);
WeightDateInDaysFormated = datestr(DateTimeOut0(iTdsWeightTrigg)/(3600*24), FormatOut);

%% Separate Time, Ref and Dst vector for each File

FuncGetColumn = @(x,iRow,iColumn) x(iRow,iColumn);

TimeRefArray = cell2mat(cellfun(@(x)FuncGetColumn(x,':',1),Data,'Un',0));% Time [ps]
RefArray = cell2mat(cellfun(@(x)FuncGetColumn(x,':',2),Data,'Un',0));
TimeDstArray = cell2mat(cellfun(@(x)FuncGetColumn(x,':',3),Data,'Un',0));% Time [ps]
DstArray = cell2mat(cellfun(@(x)FuncGetColumn(x,':',4),Data,'Un',0));
NbrePts = cell2mat(cellfun(@(x) size(x,1),Data,'Un',0));
NbrePts = NbrePts(1);

%% Compute the magnitude and phase of transmitted field
[RefDstArrayOut0, FreqTHzOut, AbsTFOut0, PhaseTFOut0,...
    AbsDstRefOut0, PhaseDstRefOut0] = deal(cell(NbreFile,1));

% Fit Frequency Range for which the level of noise is very low
FreqMinUnwrapPhaseHz = 0.1 * 1e12;
FreqMaxUnwrapPhaseHZ = .6 * 1e12;

% Use or not the Fit for Phase Unwrapping:
% 0 :   No Fit /Unwrapping the phase starting from the noisy part at low
%       frequencies causes errror to propagate towards the phase at higher
%       frequencies
% 1 :   Phase Extrapolation (0-0.1 THz) (from the unwrapped phase at higer frquencies)
% 2 :   Use Linear Fit (0-0.1 THz)
PhaseLinearFit = 1;

% Returns the (NbrePoint*NFFTCoeff)-point DFT. Zero padding to the data
% Time. It imporve the apparent frequency resolution by forcing the
% transform to generate additional frequency points.
NFFTCoeff = 0;
clc

for i1=1:NbreFile
    [RefDstArrayOut0{i1}, FreqTHzOut{i1}, AbsTFOut0{i1}, PhaseTFOut0{i1},...
        AbsDstRefOut0{i1}, PhaseDstRefOut0{i1}] =...
        GetSpectrumPhase_1 (DstArray(:,i1), RefArray(:,i1),...
        TimeRefArray(:,i1),TimeDstArray(:,i1),FreqMinUnwrapPhaseHz,...
        FreqMaxUnwrapPhaseHZ, NbrePts,PhaseLinearFit,NFFTCoeff);
end

FreqRes = FreqTHzOut{1}(2)*1e3;

%% Plot
close all
clear hFig
Xlabel = {'Frequency (THz)','Dehydration  Time (min)'};
%          1              2          3                4       5
Ylabel = {'Weight (mg)', 'RWC (%)', 'Thickness (µm)','T (%)','Residual Error (%)'};
Legend = {'Gravimetry','THz-TDS'};

AllFreqTHz= cell2mat(cellfun(@(x)FuncGetColumn(x,':',1),FreqTHzOut','Un',0));

AllAbsTFOut = cell2mat(cellfun(@(x)FuncGetColumn(x,':',1),AbsTFOut0','Un',0));
AllPhaseTFOut = cell2mat(cellfun(@(x)FuncGetColumn(x,':',1),PhaseTFOut0','Un',0));

[FreqTHz, FreqIndex] = SubArray(FreqTHzOut{i1}, 0.1, 2);
Index = [8 9 10 11];
xMax = 120;

% Create a DisplayName 
i2 = 1;
for i1 = Index
    Name{i2} = sprintf('%0.1f GHz',FreqTHz(i1)*1e3);
    i2 = i2+1;
end

% Plot 3D T_DryTime_Freq
hFig(1) = figure;
hAxes = CreatAxisXY('Linear','Linear',12,'off','off','on');
surf(FreqTHz, DryTime, AllAbsTFOut(FreqIndex,:)'*100,'LineStyle','none')
AxisXYProperties(hAxes,Xlabel{1},Xlabel{2},Ylabel{4});xlim([0.2 2]);
view(17,22)
colormap('jet')

% Plot T(f1,f2) vs DryTime
hFig(2) = figure;
hAxes = CreatAxisXY('Linear','Linear',12,'off','off','on');
hPlotNotAvg = plot(DryTime, AllAbsTFOut(FreqIndex(Index),:)'*100);
AxisXYProperties(hAxes,Xlabel{2},Ylabel{4},'r');xlim([0 xMax]);legend(Name)

hFig(3) = figure;
hAxes = CreatAxisXY('Linear','Linear',12,'off','off','on');
plot(DryTime0, RWC01,'DisplayName','','LineStyle','-',...
    'LineWidth',0.5,'Marker','.','MarkerSize',5,...
    'MarkerFaceColor','none','MarkerEdgeColor','b');
AxisXYProperties(hAxes,Xlabel{2},Ylabel{2},'r');xlim([0 xMax])

error('stop')
SampleFolder = 'Sample 11';
FigName ={'3D_T_DryTime_Freq', '2D_T(f1,f2)_DryTime', '2D_RWC_DryTime'};
mkdir(SampleFolder);
for i1=1:numel(hFig)
    savefig(hFig(i1),fullfile(SampleFolder,FigName{i1}));
end

%% Increase the SNR by averaging the set of replicate measurements (not so, because it varies with time)
close all
clc; clear hFig
LMovMean = 400;

[AllAbsTFOutAvg,AllAbsTFOutAvgError, index,NbrePoint, NbreAvrg,...
    indexStart, indexStop] = DataAveraging(AllAbsTFOut,LMovMean);
[AllPhaseTFOutAvg,AllPhaseTFOutAvgError] = DataAveraging(AllPhaseTFOut,LMovMean);

[RWC1Avg, RWC1AvgError] = DataAveraging(RWC1',LMovMean);
[RWC2Avg, RWC2AvgError] = DataAveraging(RWC2',LMovMean);
DryTimeAvg = DataAveraging(DryTime',LMovMean);

AbsTFOut = mat2cell(AllAbsTFOutAvg,NbrePoint,ones(1,NbreAvrg))';
PhaseTFOut = mat2cell(AllPhaseTFOutAvg,NbrePoint,ones(1,NbreAvrg))';
NbreFileMean = NbreAvrg;
index = round(index);

[DryTimeAvg DryTime(index)];


hFig(1) = figure;
hAxes = CreatAxisXY('Linear','Linear',12,'off','off','on');
hPlotRWcAvg = plot(DryTimeAvg, RWC1Avg,'DisplayName','','LineStyle','-',...
    'LineWidth',0.5,'Marker','.','MarkerSize',5,...
    'MarkerFaceColor','none','MarkerEdgeColor','b');
Y = RWC1Avg;
ErrorRWC1Coeff = RWC1AvgError;
ErrorRWC1 = Y.*ErrorRWC1Coeff./100;
errorbar(DryTimeAvg,Y,ErrorRWC1,'Color',hPlotRWcAvg.Color) ;
AxisXYProperties(hAxes,Xlabel{2},Ylabel{2},'r');xlim([0 xMax])


hFig(2) = figure;
hAxes = CreatAxisXY('Linear','Linear',12,'off','off','on');
Y = AllAbsTFOutAvg(FreqIndex(Index),:)'*100;
hPlotAvg = plot(DryTimeAvg, Y,...
    'LineStyle','none','LineWidth',3,'Marker','.');
AxisXYProperties(hAxes,Xlabel{2},Ylabel{4},'r');xlim([0 xMax]);
ErrorTDSCoeff = AllAbsTFOutAvgError(FreqIndex(Index),:)';
ErrorTDS = Y.*ErrorTDSCoeff./100;
for i1=1:numel(Index)
    errorbar(DryTimeAvg,Y(:,i1),ErrorTDS(:,i1),'Color',hPlotAvg(i1).Color) ;
end

legend(strcat(Name,' -Avg-',num2str(LMovMean,'%d')))

iFreqSelect = 3
hFig(3) = figure;
hAxes = CreatAxisXY('Linear','Linear',12,'off','off','on');
plot(AllAbsTFOutAvg(FreqIndex(Index),:)'*100, RWC1Avg,'DisplayName','','LineStyle','none',...
    'LineWidth',0.5,'Marker','s','MarkerSize',5,...
    'MarkerFaceColor','none');
AxisXYProperties(hAxes,Ylabel{4},Ylabel{2},'r');

legend(hPlotAvg.DisplayName)

%%
FreqSelect = 0.3;
for i1=1:NbreFileMean
    [FreqTHz(i1), iFreqExtract(i1)] = SubArray(FreqTHzOut{index(i1)}, FreqSelect);
end

for i1=1:NbreFile
    [FreqTHz(i1), iFreqExtract(i1)] = SubArray(FreqTHzOut{i1}, FreqSelect);
end

close all
hold on
plot(AllAbsTFOut(iFreqExtract,:)')
plot(index,AllAbsTFOutAvg(iFreqExtract,:)','LineWidth',3,'Color','k',...
    'LineStyle','none','Marker','s')
%%
figure
for i1=1:10
    semilogy(FreqTHzOut{i1},AbsTFOut{i1},'LineWidth',i1*0.5)
    hold on
end

%% Plot Results Abs & Arg
% close all
% Xlabel = 'Frequency (THz)';
% Ylabel = 'Phase (rad)';
% YlabelT = 'T';
%
% iFile = 1;
% MaxX = 5;
%
% figure
% hAxes = CreatAxisXY('Linear','Linear',12,'off','off','on');
% plot(FreqTHzOut{iFile},PhaseTFOut{iFile}(:,1),'r','DisplayName','Maher')
% AxisXYProperties(hAxes,Xlabel,Ylabel); xlim([0 MaxX]);
% figure
% hAxes = CreatAxisXY('Linear','log',12,'off','off','off');
% semilogy(FreqTHzOut{iFile},AbsDstRefOut{iFile}(:,2)/max(AbsDstRefOut{iFile}(:,2)),'r','DisplayName','Maher')
% AxisXYProperties(hAxes,Xlabel,YlabelT);xlim([0 MaxX])
%
% figure
% hAxes = CreatAxisXY('Linear','Linear',12,'off','off','on');
% plot(RefDstArrayOut{iFile}(:,1),RefDstArrayOut{iFile}(:,2),'r','DisplayName','Ref')
% plot(RefDstArrayOut{iFile}(:,1),RefDstArrayOut{iFile}(:,3),'b','DisplayName','Dst')
% AxisXYProperties(hAxes,'Time (ps)','Amplitude (a.u.)');
% % return

%% Select the useful part of the Signal [0.1-2]THz

FreqMin = 0.2; FreqMax = 0.6;

[FreqTHz, iFreqExtract] = SubArray(FreqTHzOut{1}, FreqMin, FreqMax);
Lambda = Lambda0(iFreqExtract);
NbrPts = numel(iFreqExtract);
WOut = 2* pi * FreqTHzOut{1}(iFreqExtract,1)'*1e12;  % Angular Frequency [rad/s]
DeltaStart = 1e2;    DeltaMin = 1e-16;
h = 1e-3; hy = 1e-3;
iSample = (1:NbrPts);

[FreqTHzSelected, FreqIndex] = SubArray(FreqTHz, FreqMin,FreqMax);%28211
NbreFreq = numel(FreqIndex);
ResMethod2 = zeros(2,NbrPts,NbrePtsL,2);


options = optimoptions('fmincon','Algorithm','sqp',...
    'Display','off','TolFun',eps.^2,'TolX',eps.^2,'MaxFunEvals',500);

nk2EpsFnc = @(n,k) (n^2-k^2)-1i*(2*n*k);
Eps2nFnc = @(Eps) sqrt(Eps.^3);
EpsAa = 1.00027;

RhoW = 1000; % kg/m3 is the density of water
RhoS = 1740;
for i1=1:NbreFileMean
    AbsT_Meas = AbsTFOut{i1}(iFreqExtract,1)';
    ArgT_Meas = PhaseTFOut{i1}(iFreqExtract,1)';
    
    for i2=1:NbrePtsL
        % Analytical Estimation of n2 (n2_InitEst) & k2 (k2_InitEst)
        
        [n2_InitEst,k2_InitEst] = AnalyticalEstimation_N_K(AbsT_Meas,...
            ArgT_Meas,N1,N3,N_air,WOut,L(i2),c,k1,NbrPts);
        x = [n2_InitEst; k2_InitEst];
        % Extraction Method 2
                
        Fnc = @(x,y) ErrorFunctionExactTnk (x,y,N1,N3,N_air,WOut(iSample),L(i2),c, FPCoeff,...
            AbsT_Meas(iSample), ArgT_Meas(iSample));
        
        [ResMethod2(:,:,i2,i1), m] = ExtractionNK_M2(x,h, hy, Fnc, DeltaStart, DeltaMin,NbrPts);
        
        % Determine The Volumetric Fraction Of Water In a Leaf.

        for i3=1:NbreFreq
            
            [~, iFreqS] = SubArray(FreqTHzInterpS, FreqTHzSelected(i3));
            nSi = nSInterp(iFreqS); kSi = kSInterp(iFreqS);
            
            % [~, iFreqW] = SubArray(FreqTHzInterpW, FreqTHzSelected(FreqIndex));
            nWi= nWInterp(iFreqS); kWi = kWInterp(iFreqS);
            
            ni = ResMethod2(1,FreqIndex(i3),i2,i1); ki = ResMethod2(2,FreqIndex(i3),i2,i1);
            
            EpsAs = nk2EpsFnc(nSi,kSi);
            EpsAw = nk2EpsFnc(nWi,kWi);
            % Effective Medium Model
            EpsEffMMFnc = @(x)  x(1)*EpsAa^(1/3)+x(2)*EpsAs^(1/3)+x(3)*EpsAw^(1/3);
            
            EpsExp = nk2EpsFnc(ni,ki)^(1/3);
            f = @(x) ((abs(EpsExp)) - (abs(EpsEffMMFnc(x))))^2 + ((angle(EpsExp)) - (angle(EpsEffMMFnc(x))))^2;
            %     Air       Dry(Solid)      Water
            x0 = [0         0.2             0.2];
%             tic
            [AExp ,Fval ,exitflag] = fmincon(f,x0,[],[],[1 1 1],1,...
                [0.0 0.0 0],[.6 .6 1],[],options);
%             toc
            AExpAll(i3,i2,i1,:) = AExp;
            AwExp(i3,i2,i1) = AExp(3);
            AsExp(i3,i2,i1) = AExp(2);
            AaExp(i3,i2,i1) = AExp(1);

            WeightW1(i3,i2,i1) = (AwExp(i3,i2,i1) * RhoW)/(AwExp(i3,i2,i1)*RhoW + AsExp(i3,i2,i1)*RhoS)*100;
            WeightW(i3,i2,i1) = AwExp(i3,i2,i1)/AwExp(i3,i2,1)*100;
            EpsEffMM(i3,i2,i1) = EpsEffMMFnc(AExp);
            nLeaf(i3,i2,i1) = Eps2nFnc(EpsEffMM(i3,i2,i1));
        end
        
        nEffMM(:,i2,i1) = real(nLeaf(:,i2,i1));
%         DeltaEpsEffMM = ((sqrt(EpsEffMM(:,i2,i1))-1)).';
%         Roughness = 4e-6;
%         AlphaScat = (DeltaEpsEffMM*4*pi*cos(0).*Roughness./Lambda).^2*1/L(i2)^2;
%         kScat = abs(imag(AlphaScat.*Lambda/(4*pi)*0));
        kEffMM(:,i2,i1) = abs(imag(nLeaf(:,i2,i1)));%+kScat.';
        DeltaAll(i2,i1) = sum(Fnc(nEffMM(:,i2,i1)',kEffMM(:,i2,i1)'));
                
        [T, Phase] = ExactSolutionT(nLeaf(:,i2,i1).',N1,N3,N_air,WOut(iSample),L(i2),c, FPCoeff);
        AbsTExact(:,i2,i1) = T';
        ArgTExact(:,i2,i1) = Phase';
    end
    NbreFileMean-i1
end
tyf
%%
clc
for i1=1:NbreFileMean
    for i2=1:NbrePtsL
        for i3=1:NbreFreq
            WeightW1(i3,i2,i1) = (AwExp(i3,i2,i1) * RhoW)/(AwExp(i3,i2,i1)*RhoW + AsExp(i3,i2,i1)*RhoS)*100;
            WeightW(i3,i2,i1) = AwExp(i3,i2,i1)/max(AwExp(i3,i2,:))*100;

            EpsEffMM(i3,i2,i1) = EpsEffMMFnc(AExpAll(i3,i2,i1,:));
            nLeaf(i3,i2,i1) = Eps2nFnc(EpsEffMM(i3,i2,i1));
        end
        nEffMM(:,i2,i1) = real(nLeaf(:,i2,i1));
        DeltaEpsEffMM = ((sqrt(EpsEffMM(:,i2,i1))-1)).';
        kEffMM(:,i2,i1) = abs(imag(nLeaf(:,i2,i1)));
        DeltaAll(i2,i1) = sum(Fnc(nEffMM(:,i2,i1)',kEffMM(:,i2,i1)'));
        
        [T, Phase] = ExactSolutionT(nLeaf(:,i2,i1).',N1,N3,N_air,WOut(iSample),L(i2),c, FPCoeff);
        AbsTExact(:,i2,i1) = T';
        ArgTExact(:,i2,i1) = Phase';
    end
    NbreFileMean-i1
end
%%
close all
clc
iFreqSelect = 5;
% Diameter of Leaf Disc
dL = 7e-3 ; % 10 mm
% Leaf Disc Area
aL = pi*(dL/2)^2;


RWC_TdsLegend  =sprintf('THz TDS (%0.2f GHz)',FreqTHz(FreqIndex(iFreqSelect))*1e3)
WeightW2 = squeeze(mean(WeightW1(iFreqSelect,:,:)));

clear BestW res RWC_TdsFinal ThecknessLeaf AwExpFinal AaExpFinal AsExpFinal
[ThecknessLeaf,RWC_TdsFinal] = deal(zeros(NbreFreq,NbreFileMean));
for i1=1:NbreFreq
    WeightWBest = squeeze(WeightW1(i1,:,:));
    [RWC_TdsFinal(i1,:),iminBestW] = min(WeightWBest);
    AwExpT = squeeze(AwExp(i1,:,:));
    AwExpFinal(i1,:) = AwExpT(iminBestW+(0:40:(NbreFileMean-1)*40));
    
    AsExpT = squeeze(AsExp(i1,:,:));
    AsExpFinal(i1,:) = AsExpT(iminBestW+(0:40:(NbreFileMean-1)*40));
    
    AaExpT = squeeze(AaExp(i1,:,:));
    AaExpFinal(i1,:) = AaExpT(iminBestW+(0:40:(NbreFileMean-1)*40));
    ThecknessLeaf(i1,:) = L_um(iminBestW);
%     Vol0 = max(ThecknessLeaf(i1,2))*aL;
%     Voli = ThecknessLeaf(i1,:)*aL;
%     RWC_TdsFinal(i1,:) = RWC_TdsFinal(i1,:).*Voli/Vol0;
end

% [ThecknessLeaf,RWC_TdsFinal] = deal(zeros(NbreFreq,NbreFileMean));
% for i1=1:NbreFileMean
%     WeightWBest = DeltaAll(:,i1);
% %     [~,iminBestW(i1)] = min(abs(diff(WeightWBest)));
%     iminBestW(i1) = find(WeightWBest<0.005,1);
%     iminBestW(i1) = iminBestW(i1)-1
%     RWC_TdsFinal(:,i1) = squeeze(WeightW1(:,iminBestW(i1),i1));
%     ThecknessLeaf(i1) = L_um(iminBestW(i1));
%     DeltaAll3D(i1) = WeightWBest(iminBestW(i1))
% %     Vol0 = max(ThecknessLeaf(i1,2))*aL;
% %     Voli = ThecknessLeaf(i1,:)*aL;
% %     RWC_TdsFinal(i1,:) = RWC_TdsFinal(i1,:).*Voli/Vol0;
% end




% Plot 3D Graph RWC_Tds (%) vs Leaf Thickness (µm) vs Dehydration Time (min)
hAxes = CreatAxisXY('Linear','Linear',12,'off','off','on');
surf(L_um,DryTimeAvg, squeeze(WeightW1(iFreqSelect,:,:))','LineStyle','none')
view (-136,23);colormap(flipud(winter)); hold on

plot3(ones(1,NbreFileMean).*ThecknessLeaf(1,:),DryTimeAvg,RWC_TdsFinal(iFreqSelect,:),...
    'Marker','o','MarkerSize',6,'MarkerFaceColor','r','MarkerEdgeColor','r')
AxisXYProperties(hAxes,Ylabel{3},Xlabel{2},Ylabel{2});

% Plot 3D Graph Residual Error Texp-Tcal (%) vs Leaf Thickness (µm) vs Dehydration Time (min)
figure
hAxes = CreatAxisXY('Linear','Linear',12,'off','off','on');
surf(L_um,DryTimeAvg, DeltaAll','LineStyle','none')
view (-136,23);colormap(flipud(winter)); hold on
% plot3(ones(1,NbreFileMean).*ThecknessLeaf(iFreqSelect,:),DryTimeAvg,DeltaAll3D,...
%     'Marker','o','MarkerSize',6,'MarkerFaceColor','r','MarkerEdgeColor','r')
AxisXYProperties(hAxes,Ylabel{3},Xlabel{2},Ylabel{5});


% Plot Leaf Thickness (µm) vs Dehydration Time (min)
figure
hAxes = CreatAxisXY('Linear','Linear',12,'off','off','on');
plot(DryTimeAvg,ThecknessLeaf(1,:),'b','DisplayName',RWC_TdsLegend,...
    'LineStyle','none',...
    'LineWidth',0.5,'Marker','s','MarkerSize',5,...
    'MarkerFaceColor','b','MarkerEdgeColor','b')
AxisXYProperties(hAxes,Xlabel{2},Ylabel{3});
legend('show')
% Plot RWC obtained from TDS and Gravimetry
RWC_GravFinal = RWC2Avg;%/max(RWC2Avg)*100;
RWC_GravFinalError = RWC_GravFinal.*RWC2AvgError/100;

RWC_TdsFinal = RWC_TdsFinal(FreqIndex(iFreqSelect),:);
RWC_TdsFinal = RWC_TdsFinal+(RWC_GravFinal(1)-RWC_TdsFinal(1));
% RWC_TdsFinal = RWC_TdsFinal/(max(RWC_TdsFinal))*100;

ErrorTDSCoeff = AllAbsTFOutAvgError(iFreqExtract(iFreqSelect),:);
RWC_TdsFinalError = RWC_TdsFinal.*ErrorTDSCoeff/100;

figure;
hAxes = CreatAxisXY('Linear','Linear',12,'off','off','on');
hRWC_Grav = plot(DryTimeAvg, RWC_GravFinal,'DisplayName','Gravimetry',...
    'LineStyle','none','LineWidth',2,'Marker','v','MarkerSize',8,...
    'MarkerFaceColor','none','MarkerEdgeColor','r');
errorbar(DryTimeAvg,RWC_GravFinal,RWC_GravFinalError,'Color','r') ;

hold on

hRWC_Tds = plot(DryTimeAvg, RWC_TdsFinal,'DisplayName',RWC_TdsLegend,...
    'LineStyle','none','LineWidth',2,'Marker','s','MarkerSize',8,...
    'MarkerFaceColor','none','MarkerEdgeColor','b');
errorbar(DryTimeAvg,RWC_TdsFinal,RWC_TdsFinalError,'Color','b') ;
    
legend([hRWC_Grav hRWC_Tds])

AxisXYProperties(hAxes,Xlabel{2},Ylabel{2});


% Plot Residual Error (%) between Gravimetry and TDS vs Dehydration Time (min)
TotalError = abs(RWC_TdsFinal'-RWC_GravFinal)./RWC_GravFinal*100;
figure
hAxes = CreatAxisXY('Linear','Linear',12,'off','off','on');
plot(DryTimeAvg,TotalError,'b','DisplayName','',...
    'LineStyle','none',...
    'LineWidth',0.5,'Marker','s','MarkerSize',5,...
    'MarkerFaceColor','b','MarkerEdgeColor','b')
AxisXYProperties(hAxes,Xlabel{2},Ylabel{5});

%%
close all
    RWC_GravFinal = RWC1Avg;%/max(RWC1Avg)*100;
    RWC_GravFinalError = RWC_GravFinal.*RWC1AvgError/100;
    figure;
    hAxes = CreatAxisXY('Linear','Linear',12,'off','off','on');    
    
for i1=1:NbreFreq
    WeightWBest = squeeze(WeightW(i1,:,:));
    [RWC_TdsFinalMat(i1,:),iminBestW] = min(WeightWBest);
end
    
for i1=1:NbreFreq
    iFreqSelect = i1;
    % Plot RWC obtained from TDS and Gravimetry

    
    RWC_TdsFinal = RWC_TdsFinalMat(FreqIndex(iFreqSelect),:);
    % RWC_TdsFinal = RWC_TdsFinal-RWC_TdsFinal*36/100;
    RWC_TdsFinal = RWC_TdsFinal/(max(RWC_TdsFinal))*100;
    
    ErrorTDSCoeff = AllAbsTFOutAvgError(iFreqExtract(iFreqSelect),:);
    RWC_TdsFinalError = RWC_TdsFinal.*ErrorTDSCoeff/100;
    

    hRWC_Grav = plot(DryTimeAvg, RWC_GravFinal,'DisplayName','Gravimetry',...
        'LineStyle','none','LineWidth',2,'Marker','v','MarkerSize',8,...
        'MarkerFaceColor','none','MarkerEdgeColor','r');
    
    hold on
    
    hRWC_Tds = plot(DryTimeAvg, RWC_TdsFinal,'DisplayName',RWC_TdsLegend,...
        'LineStyle','none','LineWidth',2,'Marker','s','MarkerSize',8,...
        'MarkerFaceColor','none','MarkerEdgeColor','b');
    
    legend([hRWC_Grav hRWC_Tds])
    

end
    AxisXYProperties(hAxes,Xlabel{2},Ylabel{2});
return
%%
%
% plot(L_um,WeightWBest)
Gravimetric = RWC1Avg(iMeasSelected)/RWC1Avg(1)*100;

[~,iFreqBest] = SubArray(FreqTHz,0.6);
% iFreqBest = 6;
FreqTHz(iFreqBest)
RWCEMM = RWC_TdsFinal(iFreqBest,:)';
RWCEMMNorm = RWCEMM/max(RWCEMM)*100;
RWCExpNorm = RWC1Avg/RWC1Avg(1)*100;

hAxes = CreatAxisXY('Linear','Linear',12,'off','off','on');
plot(DryTimeAvg(indexSelected),RWCExpNorm,'rv','LineWidth',2,'MarkerSize',8)
hold on
plot(DryTimeAvg(indexSelected),RWCEMMNorm,'bs','LineWidth',2,'MarkerSize',8)
AxisXYProperties(hAxes,Xlabel,Ylabel{2}); 
legend({'Gravimetry', 'THz TDS'},...
    'FontName','Times New Roman','FontSize',10,'FontWeight','Bold')


%% Compare nExp_L & kExp_L with nEffMM_L & kEffMM_L for L(i)  
AlphaExpFnc = @(k) k.*4.*pi.*(FreqTHz*ones(1,size(k,2))*1e12)/c*1e-2;
close all
clc
iLSelected = 60;
iMeasSelected = [1 20 40];
L_um(iLSelected)
nEffMM_L = squeeze(nEffMM(:,iLSelected,iMeasSelected));
kEffMM_L = squeeze(kEffMM(:,iLSelected,iMeasSelected));
AlphaEffMM_L = AlphaExpFnc(kEffMM_L);

nExp_L = squeeze(ResMethod2(1,:,iLSelected,iMeasSelected));
kExp_L = squeeze(ResMethod2(2,:,iLSelected,iMeasSelected));
AlphaExp_L = AlphaExpFnc(kExp_L);

figure('Position',[117   567   560   420])
plot(FreqTHz,nEffMM_L); hold on
plot(FreqTHz,nExp_L,'o','LineWidth',2)
LegendStr = [strcat('EMM_ ',num2str(iMeasSelected','%03d'));...
    strcat('Exp_',num2str(iMeasSelected','%03d'))];
legend(cellstr(LegendStr),'Interpreter','none')


figure('Position',[681   567   560   420])
plot(FreqTHz,kEffMM_L); hold on
plot(FreqTHz,kExp_L,'o','LineWidth',2)
legend(cellstr(LegendStr),'Interpreter','none')

figure('Position',[1200   567   560   420])
plot(FreqTHz,AlphaEffMM_L); hold on
plot(FreqTHz,AlphaExp_L,'o','LineWidth',2)
legend(cellstr(LegendStr),'Interpreter','none')

% Compare TExp
AbsT_Meas = {AbsTFOut{iMeasSelected}};
AbsT_Meas = cell2mat(cellfun(@(x)FuncGetColumn(x,iFreqExtract,1),AbsT_Meas,'Un',0));

ArgT_Meas = {PhaseTFOut{iMeasSelected}};
ArgT_Meas = cell2mat(cellfun(@(x)FuncGetColumn(x,iFreqExtract,1),ArgT_Meas,'Un',0));
    
figure('Position',[117   46   560   420])
TExp_L = squeeze(AbsTExact(:,iLSelected,iMeasSelected));
ArgExp_L = squeeze(ArgTExact(:,iLSelected,iMeasSelected));

plot(FreqTHz,TExp_L,'LineWidth',2); hold on
plot(FreqTHz,AbsT_Meas,'o','LineWidth',2)
legend(cellstr(LegendStr),'Interpreter','none')

figure('Position',[681   46   560   420])
plot(FreqTHz,ArgExp_L,'LineWidth',2); hold on
plot(FreqTHz,ArgT_Meas,'o','LineWidth',2)
legend(cellstr(LegendStr),'Interpreter','none')


%%
iFreqTHzSelected = 3;
close all
clc
FreqTHz(iFreqTHzSelected)
plot(L_um,(squeeze(WeightW1(1:10,:,iFreqTHzSelected)))')
%%
iFreqTHzSelected = 1;
close all
clc
plot(L_um, squeeze(AwExp(:,:,50)))

%%
close all
clc
plot(FreqTHz, squeeze(ResMethod2(iFreqTHzSelected,:,:,1)))
%% WeightW
clc
iMeasSelected = 1:NbreFileMean;
indexSelected = iMeasSelected;
close all
clear BestW res
WeightW12 = WeightW;
for i1=1:NbreFreq
    WeightWBest = squeeze(WeightW12(i1,:,iMeasSelected));
    
    [RWC_TdsFinal(i1,:),iminBestW] = min(WeightWBest);
%     plot(L_um,WeightWBest')
    ThecknessLeaf(i1,:) = L_um(iminBestW);
end

%
% plot(L_um,WeightWBest)
Gravimetric = RWC1Avg(iMeasSelected)/RWC1Avg(1)*100;

[~,iFreqBest] = SubArray(FreqTHz,0.6);
% iFreqBest = 6;
FreqTHz(iFreqBest)
RWCEMM = RWC_TdsFinal(iFreqBest,:)';
RWCEMMNorm = RWCEMM/max(RWCEMM)*100;
RWCExpNorm = RWC1Avg/RWC1Avg(1)*100;

hAxes = CreatAxisXY('Linear','Linear',12,'off','off','on');
plot(DryTimeAvg(indexSelected),RWCExpNorm,'rv','LineWidth',2,'MarkerSize',8)
hold on
plot(DryTimeAvg(indexSelected),RWCEMMNorm,'bs','LineWidth',2,'MarkerSize',8)
AxisXYProperties(hAxes,Xlabel,Ylabel{2}); 
legend({'Gravimetry', 'THz TDS'},...
    'FontName','Times New Roman','FontSize',10,'FontWeight','Bold')
%%
ErrorTDS = abs(diff([RWCExpNorm RWCEMMNorm],1,2)./RWCExpNorm*100);
neg = RWCExp-RWC1(indexStart);
pos = RWCExp-RWC1(indexStop);
err = abs(neg)+abs(pos);
errorbar(DryTime0(indexSelected),RWCExpNorm,neg,pos) ;
% errorbar(DryingTime{2}(indexSelected),RWCEMMNorm,Error) ;
xlim([0 160])
ylim([40 100])
mean(ErrorTDS(1:17))
figure
hAxes = CreatAxisXY('Linear','Linear',12,'off','off','on');
plot(DryTime0(indexSelected),ErrorTDS,'bs','LineWidth',2,'MarkerSize',8)
AxisXYProperties(hAxes,Xlabel,'Error (%)'); 
xlim([0 160])
%%
RWC0 = RWC1;
clc
close all
i = 10;
plot(L_um,squeeze(WeightW12(:,:,i))); hold on
Gravimetric(i)
for i1=1:56
    a(i1) = L_um(find(DeltaAll(:,i1)<=0.05,1));
end
% plot(L_um,BestW(:,i),'LineWidth',2)
[RWCExp,RWCExpError, iRWCExp,~,~, indexStart, indexStop] =...
    DataAveraging(RWC0',LMovMean);
RWCExp = RWCExp(iMeasSelected);
indexStart = indexStart(iMeasSelected);
indexStop = indexStop(iMeasSelected);

%%
% surf(DryingTime{2}(indexSelected),L_um,DeltaAll)
figure
hAxes = CreatAxisXY('Linear','Linear',12,'off','off','on');
plot(DryingTime{2}(indexSelected),ThecknessLeaf(1,:),'rv','LineWidth',2,'MarkerSize',8)
AxisXYProperties(hAxes,Xlabel,'L (µm)'); 

   %%
   
   clc
   WeightWBest = squeeze(WeightW12(:,:,1));
   for i1=1:NbreFreq
  m{i1}= find(abs(WeightWBest(i1,:)-64)<.1);
   end
   m'
%%
[i,j] = SubArray(L_um,345)
WeightW(:,j,1)
%% Thickness
close all
LBest = L_um(iminBestW);
plot(LBest)


%% Index
close all
nLeafBest = squeeze(nLeaf(iFreqTHzSelected,:,:));
nLeafBest = nLeafBest(iminBestW*1:NbrePtsL:NbrePtsL*NbreFile)
plot(imag(nLeafBest),'.')
%%
ResX(iminBestW,:)
RWC_TdsFinal
BestL*1e6

return
%% Refractive index one file
clc;close all

iFreqTHzSelected = 1;
nEffMM = real(nLeaf);
kEffMM = -imag(nLeaf);



n = squeeze(ResMethod2(1,:,:,1));
k = squeeze(ResMethod2(2,:,:,1));

plot(FreqTHz,nEffMM(:,:))
hold on
plot(FreqTHz, n(:,100),'or')

figure
plot(FreqTHz,kEffMM(:,:))
hold on
plot(FreqTHz, k(:,100),'ok')
%%
clear rCoef
for i1=1:100
    for i2=1:100
        rCoef(i1,i2) = mean(std([nEffMM(:,i1) n(:,i2)],1,2));
    end
end
%%
[v,i] = min(rCoef,[],2)
plot(L_um(i),v)
