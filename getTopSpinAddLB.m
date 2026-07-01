%Function for MRSHIFT
%Function to read out experimental and processing parameters from Topspin

function [expParams]=getTopSpinAddLB(path01,l,k,scale,AddLB,LBmode)
% Load and process a TopSpin spectrum with optional LB handling.
%


% AddLB:
%   In default 'extra' mode:
%       0      = use original TopSpin LB only
%       nonzero = add extra exponential broadening via 1/AddLB
%
%   In 'total' mode:
%       interpreted as the desired total exponential LB
%       0 = true zero LB
%
% LBmode:
%   'extra'  default, backwards-compatible behavior
%   'total'  use AddLB as total LB instead of extra LB

format long
if ~exist('LBmode', 'var') || isempty(LBmode)
    LBmode = 'extra';
end

%For reading difference spectra created with difference_noesy_CPMG
if l == 1000
    matFile = fullfile(path01, num2str(k), 'pdata', num2str(l), 'diffNOESYCPMG.mat');

    S = load(matFile);
    expParams = S.expParams;

    return
end

%Get Acquisition parameter
acqus = fopen(fullfile(path01, num2str(k), 'acqus'));
B=textscan(acqus, '%s' ) ;
IndexB = strfind(B{1}, '$TD=' ); 
Index1 = find(not(cellfun( 'isempty' , IndexB))); 
np=B{1}{Index1+1} ;
np=str2num(np);
IndexC = strfind(B{1}, 'SW_h' ); 
Index2 = find(not(cellfun( 'isempty' , IndexC))); 
swh=B{1}{Index2+1} ;
swh=str2num(swh);
dwell=1/swh;
IndexC = strfind(B{1}, 'BF1' ); 
Index2 = find(not(cellfun( 'isempty' , IndexC))); 
BF1=str2num(B{1}{Index2+1});
IndexC = strfind(B{1}, 'SFO1' ); 
Index2 = find(not(cellfun( 'isempty' , IndexC))); 
SFO1=str2num(B{1}{Index2+1});
IndexC = strfind(B{1}, 'NUC1' ); 
Index2 = find(not(cellfun( 'isempty' , IndexC))); 
nucLine=B{1}{Index2+1};
    startIdx = strfind(nucLine, '<');
    endIdx = strfind(nucLine, '>');
    NUC1 = nucLine(startIdx+1:endIdx-1);
    IndexE = strfind(B{1}, '$TD=' ); 
Index4 = find(not(cellfun( 'isempty' , IndexE))); 
TD=str2num(B{1}{Index4+1});
   IndexE = strfind(B{1}, '$NS=' ); 
Index4 = find(not(cellfun( 'isempty' , IndexE))); 
NS=str2num(B{1}{Index4+1});
IndexE = strfind(B{1}, '$NBL=' ); 
Index4 = find(not(cellfun( 'isempty' , IndexE))); 
NBL=str2num(B{1}{Index4+1});
fclose(acqus);

%Get proc paramater:
procs = fopen(fullfile(path01, num2str(k), 'pdata', num2str(l), 'procs'));
B=textscan(procs, '%s' ) ;
IndexB = strfind(B{1}, 'OFFSET' ); 
Index1 = find(not(cellfun( 'isempty' , IndexB))); 
Offsetppm=str2num(B{1}{Index1+1});
% Offsetppm
IndexB = strfind(B{1}, '$SF' ); 
Index1 = find(not(cellfun( 'isempty' , IndexB))); 
SF=str2num(B{1}{Index1+1});
IndexC = strfind(B{1}, 'LB' ); 
Index2 = find(not(cellfun( 'isempty' , IndexC))); 
LB=str2num(B{1}{Index2+1});
% LB
% LB=0.8;
IndexD = strfind(B{1}, 'TDeff' ); 
Index3 = find(not(cellfun( 'isempty' , IndexD))); 
TDeff=str2num(B{1}{Index3+1});
IndexE = strfind(B{1}, '$SI=' ); 
Index4 = find(not(cellfun( 'isempty' , IndexE))); 
ZeroFill=str2num(B{1}{Index4+1});
IndexF = strfind(B{1}, 'PHC0' ); 
Index5 = find(not(cellfun( 'isempty' , IndexF))); 
PHC0=str2num(B{1}{Index5+1});
IndexG = strfind(B{1}, 'PHC1' ); 
Index6 = find(not(cellfun( 'isempty' , IndexG))); 
PHC1=str2num(B{1}{Index6+1});
IndexG = strfind(B{1}, 'WDW' ); 
Index6 = find(not(cellfun( 'isempty' , IndexG))); 
WDW=str2num(B{1}{Index6+1});
IndexG = strfind(B{1}, 'SSB' ); 
Index6 = find(not(cellfun( 'isempty' , IndexG))); 
SSB=str2num(B{1}{Index6+1});
fclose(procs);

%Get proc2s size
proc2s = fopen(fullfile(path01, num2str(k), 'pdata', num2str(l), 'proc2s'));
if proc2s~=-1
    B=textscan(proc2s, '%s' ); 
    IndexG = strfind(B{1}, '$SI=' ); 
    Index6 = find(not(cellfun( 'isempty' , IndexG))); 
    NBL=str2num(B{1}{Index6+1});
fclose(proc2s);
end

%Open the actual data, which is normally either the fid file (1D) or ser
%file (2D). Depending on the topspin version, it is saved as int32 or float
%64. the following should take care of the various options by checking that
%the number of points of the data matches what it should from TD and NBL!
serBruker = fopen(fullfile(path01, num2str(k), 'ser'));%Mac and windows have opposite slashes
if serBruker==-1
    serBruker = fopen(fullfile(path01, num2str(k), 'fid'));%Mac and windows have opposite slashes
end
data=fread(serBruker,'int32');

if length(data) ~= TD*NBL
    serBruker = fopen(fullfile(path01, num2str(k), 'ser'));%Mac and windows have opposite slashes
    if serBruker==-1
    serBruker = fopen(fullfile(path01, num2str(k), 'fid'));%Mac and windows have opposite slashes
    end
    data=fread(serBruker,'float64');
end
fclose(serBruker);

%reshape vector into matrix according to number of points
% % % np=size(data,1);
data=reshape(data,np,size(data,1)/np);

%Number of points to calculate
% In old "extra" mode, keep the previous behavior.
% In "total" mode, do not reinterpret AddLB as a decay time.
if strcmp(LBmode, 'extra') && AddLB~=0
    TDeff = round(20*AddLB/dwell/2)*2; % and ensure even number
end


%separate real and imaginary and make complex dataset
A = data(1:2:end,:);  % odd matrix
B = data(2:2:end,:); % even matrix
data=A+i*B;
clear A B

%remove the first 69 points (the filter) and cut to TDeff
datatemp=data(1:round(TDeff/2),:);
datatemp(end:ZeroFill,:)=0; %Do also zero filling according to SI F2
data=datatemp;
pts=size(data);
clear datatemp
%%


%first order phase correction: Getting the smoothest baseline by instead of
%removing the filter phase correct

PHC1b=27360;
tzero=round(PHC1b/360)+1; %Define t=0 point
%define first point as t=0 and use dwell time to get rest
t=dwell:dwell:dwell*(size(data,1));
t=t-t(tzero);

switch LBmode
    case 'extra'
        if AddLB~=0
            effectiveExpDecay = LB*2 + 1/AddLB;
            extraExpDecay = 1/AddLB;
        else
            effectiveExpDecay = LB*2;
            extraExpDecay = 0;
        end

    case 'total'
        effectiveExpDecay = AddLB*2;
        extraExpDecay = [];

    otherwise
        error('Unknown LBmode: %s', LBmode);
end
%Line Broadening Do this before the actual phase correction, otherwise the
%data axis doesnt match the time axis - but with the t=0 point in the right
%place
% figure(101)
% plot(t,data(:,1))

if WDW==0
elseif WDW==1
    data = data .* exp(-abs(t)' * effectiveExpDecay);
elseif WDW==3
    x = linspace(0, pi/SSB, round(TDeff/2));
    if SSB == 1
        window = sin(x);  % sin(0 to π) = 0 to 1
    elseif SSB == 2
        window = cos(x);  % cos(0 to π/2) = 1 to 0
    else
        uialert(app.UIFigure, 'This appodization window is not defined, modify in topspin or go to getTopSpinAddLB.m script and add it', 'Error');
    end
    window(end:ZeroFill)=0;
    data=data.*window';
elseif WDW==4
        x = linspace(0, pi/SSB, round(TDeff/2));
        if SSB == 1
            window = sin(x).^2;  % sin(0 to π) = 0 to 1 
        elseif SSB == 2
            window = cos(x).^2;  % cos(0 to π/2) = 1 to 0
        else
            uialert(app.UIFigure, 'This appodization window is not defined, modify in topspin or go to getTopSpinAddLB.m script and add it', 'Error');
        end
        window(end:ZeroFill)=0;
        data=data.*window';
else
    uialert(app.UIFigure, 'This appodization window is not defined, modify in topspin or go to getTopSpinAddLB.m script and add it', 'Error');
end

if strcmp(LBmode, 'extra') && AddLB~=0 && WDW~=1
    data = data .* exp(-abs(t)' * extraExpDecay);
elseif strcmp(LBmode, 'total') && WDW~=1
    warning('LBmode total currently only controls exponential LB cleanly for WDW==1.');
end


if PHC1b>0
dataL=data(1:tzero-1,:);
dataR=data(tzero:pts(1),:);
else
dataL=data(1:pts(1)+tzero-1,:);
dataR=data(pts(1)+tzero:pts(1),:);
end
data=cat(1,dataR,dataL);
clear dataR dataL




%Zeroth order phase correction from PHC0:
% data=data*exp(-1i*(PHC0)/180*pi);
data(1,:)=data(1,:)*1.0;%For nicer baseline
% figure(100)
% plot(t,data(:,1))
% size(data)



%FT last or first experimental point
% data(1,size(data,2))=data(1,size(data,2))/2;
SpecExp= circshift(fftshift(fft(data),1),-1,1);

corrPHC=1+PHC1/PHC0;%To match phase in topspin seems to require this
Phi=-corrPHC*PHC0+PHC1*[1:1:length(SpecExp)]'/length(SpecExp);

F1DappReal=real(SpecExp).*cos(Phi/180*pi)-imag(SpecExp).*sin(Phi/180*pi);
F1DappImag=imag(SpecExp).*cos(Phi/180*pi)+real(SpecExp).*sin(Phi/180*pi);

SpecExp=F1DappReal+1i*F1DappImag;
SpecExp=SpecExp/scale;


baxis=-1/(2*dwell)+1/(dwell*pts(1)):(1/(dwell*pts(1))):1/(2*dwell); % spectral width is given by 1/dwell. and maximal frequency by 1/2dwell (nyquist). from this and the used np calculate the step-size.

baxisppm=baxis/BF1;
XAxis=baxisppm-baxisppm(end)+Offsetppm;%In ppm with Offset

% figure(101)
% plot(XAxis, SpecExp)
% set(gca, 'XDir','reverse');

expParams.SpecExpTot2D=SpecExp;
expParams.pts=pts;
expParams.dwell=dwell;
expParams.BF1=BF1;
expParams.SFO1=SFO1;
expParams.SF=SF;
expParams.TDeff=TDeff;
expParams.NUC1=NUC1;
expParams.XAxis=XAxis;
expParams.SSB=SSB;
expParams.WDW=WDW;
expParams.LB=LB;
expParams.NS=NS;

end