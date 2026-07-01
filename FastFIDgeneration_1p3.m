%Function for MRSHIFT
%Generate spectrum from skeleton. Requires that the spin system was already treated
%with NSpin_Simulation to create skeleton from diagonalization. 
%Used for reploting, or for amp/T2 fits.

function [F]=FastFIDgeneration_1p3(params,expParams,AddLB,LBmode)
if ~exist('LBmode', 'var') || isempty(LBmode)
    LBmode = 'extra';
end
if ~exist('AddLB', 'var')
    AddLB = [];
end
w_CS_J  = params.w_CS_J;
T2      = params.T2;
GL      = params.GL;
HBlocks = params.block_sizes;
Skeleton = params.Skeleton;

if isfield(params, 'compressor')
    compressor = params.compressor;
else
    compressor = [];
end
% compressor

dwell=expParams.dwell;  np=expParams.pts;  SFO1=expParams.SFO1;    BF1=expParams.BF1;   SF=expParams.SF;   TDeff=expParams.TDeff;

FID=zeros(1,TDeff/2);

t  =0:expParams.dwell:expParams.dwell*(TDeff/2-1);


%For equivalent protons, only fit for one parameter, this is achieved in
%main script throguh compression of the w_CS_J matrix, here this step is
%reverted:

[decompressed_w_CS_J,decompressed_T2,decompressed_GL]=decompMethyl(w_CS_J,T2,GL,compressor);

% T2
w_CS_J=decompressed_w_CS_J;
T2=decompressed_T2;%Reconstruct T2
GL=decompressed_GL;%Reconstruct T2
% T2
wCS=diag(w_CS_J);

counter=1;%Counter for block structure
for m=1:numel(HBlocks)
wCS_temp=wCS(counter:HBlocks(m)+counter-1);
w_CS_J_temp=w_CS_J(counter:HBlocks(m)+counter-1,counter:HBlocks(m)+counter-1);
T2_temp=T2(counter:HBlocks(m)+counter-1);
GL_temp=GL(counter:HBlocks(m)+counter-1);
counter=counter+HBlocks(m);

isSpin1 = (abs(wCS_temp) > 0.5e14 & abs(wCS_temp) < 1.5e14);%This finds 14N (since spin 1, needs special care)
idx13C  = find(abs(wCS_temp) > 0.5e13 & abs(wCS_temp) < 1.5e13);%This finds 13C (since low natural abundance, needs special care)
numC = length(idx13C);


%Number of Spins
% nSpins=length(wCS_temp);
pass=1;
% Run multiple passes if 13C exists, otherwise 1 pass
for pass = 1:(numC + 1)
    

skel = Skeleton(pass, m);

T2_for_bins = T2_temp(skel.spinMap)'; 
GL_for_bins = GL_temp(skel.spinMap)';
        
        % Vectorized FID generation
        % freqAll is already in Hz, so multiply by 2*pi for radians
        freqRad = skel.freqAll * 2 * pi;
        
        
        
        % Apply the scaleFactor (for 13C abundance) and sum the contributions
        
FID = FID + skel.scaleFactor * (skel.Agrp' * exp(1i * freqRad * t - t ./ T2_for_bins - 0.5 * (GL_for_bins .* t).^2));
        
% T2_temp
% T2_for_bins
end

end %End blocks loop


FID(TDeff/2 + 1:np(1))=0; %Do also zero filling according to SI F2
FID(1)=FID(1)/2;%For FT artifact


ttot  =0:expParams.dwell:expParams.dwell*(np(1)-1);

x = linspace(0, pi/2, round(TDeff/2));

if expParams.WDW==0
elseif expParams.WDW==1
     switch LBmode
        case 'extra'
            if AddLB ~= 0
                effectiveExpDecay = expParams.LB*2 + 1/AddLB;
            else
                effectiveExpDecay = expParams.LB*2;
            end
    
        case 'total'
            effectiveExpDecay = AddLB*2;
    
        otherwise
            error('Unknown LBmode: %s', LBmode);
     end
         FID = FID .* exp(-ttot * effectiveExpDecay);

elseif expParams.WDW==3 && expParams.SSB == 2
    window = cos(x);  % cos(0 to π/2) = 1 to 0
    window(TDeff/2:np(1))=0;%Zero fill 
    FID=FID.*window;
elseif expParams.WDW==4 && expParams.SSB == 2
    window = cos(x).^2;  % cos(0 to π/2) = 1 to 0
    window(TDeff/2:np(1))=0;%Zero fill 
    FID=FID.*window;
else
    uialert(app.UIFigure, 'This appodization window is not defined, modify in topspin or go to getTopSpinAddLB.m script and add it', 'Error');
end

if strcmp(LBmode, 'extra') && ~isempty(AddLB) && AddLB~=0 && expParams.WDW~=1
    t(TDeff/2 + 1:np(1)) = 0; % Do also zero filling according to SI F2

    FID = FID .* exp(-t * (1/AddLB));

elseif strcmp(LBmode, 'total') && expParams.WDW~=1
    warning('LBmode total currently only controls exponential LB cleanly for WDW==1.');
end

  % figure(3)
  % plot(ttot,FID)
  
  
F=fftshift(fft(FID,np(1)));


F=real(F)';
F=circshift(F,-1,1);


function [decompressed_w_CS_J,decompressed_T2,decompressed_GL]=decompMethyl(w_CS_J,T2,GL,compressor)
decompressed_w_CS_J=w_CS_J;
decompressed_T2=T2;
decompressed_GL=GL;
    for k=1:length(compressor)
        p=compressor(k);
        decompressed_w_CS_J=[decompressed_w_CS_J(:,1:p-1) zeros(size(decompressed_w_CS_J,1),1) decompressed_w_CS_J(:,p:end)];
        decompressed_w_CS_J=[decompressed_w_CS_J(1:p-1,:) ; zeros(1,size(decompressed_w_CS_J,2)); decompressed_w_CS_J(p:end,:)];
        decompressed_w_CS_J(p,p)=decompressed_w_CS_J(p-1,p-1);
        decompressed_T2=[decompressed_T2(1:p-1) decompressed_T2(p-1) decompressed_T2(p:end)];
        decompressed_GL=[decompressed_GL(1:p-1) decompressed_GL(p-1) decompressed_GL(p:end)];
    end
    
    for k=1:length(compressor)
        p=compressor(k);
        decompressed_w_CS_J(:,p)=decompressed_w_CS_J(:,p)+decompressed_w_CS_J(:,p-1);
        decompressed_w_CS_J(p-1,p)=0;
        decompressed_w_CS_J(p,:)=decompressed_w_CS_J(p,:)+decompressed_w_CS_J(p-1,:);
        decompressed_w_CS_J(p,p-1)=0;
    end
end


end
