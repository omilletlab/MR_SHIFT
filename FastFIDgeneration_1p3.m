function [F]=FastFIDgeneration_1p3(params,expParams,AddLB,LBmode)
% FASTFIDGENERATION_1P3 Generate a processed spectrum from a previously
% calculated MRSHIFT transition skeleton.
%
% This is the fast simulation path used for spectrum replotting and repeated
% amplitude/T2/Gaussian-line-shape fitting. The expensive spin-Hamiltonian
% diagonalization is performed once by NSpin_Simulation_1p3; this function
% then reconstructs the FID directly from the stored transition groups.
%
% INPUTS
%   params
%       Spin-system and fitting structure containing:
%         w_CS_J      - Chemical-shift/J-coupling matrix. It is used here
%                       mainly to recover block layout and equivalent spins.
%         T2          - Current transverse-relaxation times.
%         GL          - Current Gaussian broadening parameters.
%         block_sizes - Number of spins in each independent simulation block.
%         Skeleton    - Precomputed grouped transitions from
%                       NSpin_Simulation_1p3.
%         compressor  - Optional equivalent-spin compression indices.
%
%   expParams
%       Acquisition and processing parameters:
%         dwell  - Sampling interval in seconds.
%         pts    - Number of points in the processed spectrum.
%         TDeff  - Effective acquired time-domain size.
%         WDW    - TopSpin window-function code.
%         LB     - Exponential line-broadening value.
%         SSB    - Sine-bell shift parameter.
%       SFO1, BF1, and SF are retained for interface compatibility, although
%       the skeleton already contains the required transition frequencies.
%
%   AddLB
%       Optional additional or total broadening parameter. Its interpretation
%       depends on LBmode.
%
%   LBmode
%       'extra' (default): combine AddLB with the original processing.
%       'total'          : request a total exponential broadening where the
%                          selected window supports it.
%
% OUTPUT
%   F
%       Real, processed, frequency-domain spectrum as a column vector.
%
% ALGORITHM
%   1. Expand compressed equivalent-spin parameters.
%   2. Split T2 and GL values according to the database block structure.
%   3. For every block and isotope pass, map grouped transitions to their
%      dominant physical spin and generate all damped sinusoids at once.
%   4. Sum the vectorized contributions into the complex FID.
%   5. Apply zero filling, apodization, and optional extra broadening.
%   6. Fourier transform and return the aligned real spectrum.

%% Optional broadening inputs
if ~exist('LBmode', 'var') || isempty(LBmode)
    LBmode = 'extra';
end
if ~exist('AddLB', 'var')
    AddLB = [];
end

%% Unpack spin-system and acquisition parameters
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
%% Restore compressed equivalent-spin parameters
[decompressed_w_CS_J,decompressed_T2,decompressed_GL]=decompMethyl(w_CS_J,T2,GL,compressor);

% T2
w_CS_J=decompressed_w_CS_J;
T2=decompressed_T2;%Reconstruct T2
GL=decompressed_GL;%Reconstruct T2
% T2
wCS=diag(w_CS_J);

%% Simulate each independent block of the spin system
% Blocks reduce the Hamiltonian dimension when the database specifies
% groups that can be simulated independently and summed in the FID.
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
    
% Retrieve grouped frequencies, amplitudes, abundance scale, and the local
% physical-spin map calculated by NSpin_Simulation.
skel = Skeleton(pass, m);

% Assign each grouped transition the T2 and GL of its dominant physical spin.
T2_for_bins = T2_temp(skel.spinMap)'; 
GL_for_bins = GL_temp(skel.spinMap)';
        

        % Vectorized FID generation: matrix dimensions are
        %   transitions x time points for the exponential term, and
        %   one amplitude per transition in skel.Agrp.
        % freqAll is already in Hz, so multiply by 2*pi for radians
        % Skeleton frequencies are stored in hertz; complex evolution requires
        % angular frequency in radians per second.
        freqRad = skel.freqAll * 2 * pi;
        
        
        
        % Apply the scaleFactor (for 13C abundance) and sum the contributions

        % Sum all grouped transitions in one matrix operation. The exponent combines
% oscillation, Lorentzian T2 decay, and Gaussian damping. scaleFactor applies
% the isotopic-abundance weight stored for this pass.        
FID = FID + skel.scaleFactor * (skel.Agrp' * exp(1i * freqRad * t - t ./ T2_for_bins - 0.5 * (GL_for_bins .* t).^2));

end

end %End blocks loop

%% Zero-fill and apply the selected TopSpin processing window
FID(TDeff/2 + 1:np(1))=0; %Do also zero filling according to SI F2
FID(1)=FID(1)/2;%For FT artifact


ttot  =0:expParams.dwell:expParams.dwell*(np(1)-1);

x = linspace(0, pi/2, round(TDeff/2));
% TopSpin WDW interpretation used here:
%   0 = no apodization
%   1 = exponential
%   3 = cosine when SSB equals 2
%   4 = cosine-squared when SSB equals 2
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
  
%% Transform the FID and format the returned spectrum  
F=fftshift(fft(FID,np(1)));


F=real(F)';
F=circshift(F,-1,1);

%% Local helper: expand compressed equivalent-spin entries
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
