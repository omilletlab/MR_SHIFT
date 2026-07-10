function [F,Skeleton]=NSpin_Simulation_1p3(params,expParams,AddLB,LBmode)
% NSPIN_SIMULATION_1P3 Simulate an NMR spectrum by explicit spin-Hamiltonian
% construction and matrix diagonalization.
%
% This function is the main spectral simulation engine used by MRSHIFT. It
% converts a database spin-system description into a processed frequency-
% domain spectrum and also returns a compact transition "Skeleton" that can
% be reused by faster fitting routines.
%
% INPUTS
%   params
%       Structure describing the spin system. The fields used here are:
%         w_CS_J      - Square matrix containing chemical shifts on the
%                       diagonal and scalar J couplings off-diagonal.
%         T2          - Transverse-relaxation time assigned to each spin.
%         GL          - Gaussian broadening parameter assigned to each spin.
%         block_sizes - Number of spins in each independently simulated
%                       block.
%         compressor  - Optional indices describing equivalent-spin entries
%                       that were compressed by the calling code.
%
%   expParams
%       Structure containing acquisition and processing parameters:
%         dwell  - Sampling interval in seconds.
%         pts    - Final number of frequency-domain points.
%         SFO1   - Spectrometer carrier frequency.
%         BF1    - Base frequency used to convert ppm offsets to hertz.
%         SF     - Processed spectrum reference frequency.
%         TDeff  - Effective acquired time-domain size.
%         WDW    - TopSpin window-function identifier.
%         LB     - Exponential line-broadening parameter.
%         SSB    - Sine-bell shift parameter.
%
%   AddLB
%       Optional additional or total broadening control. Its interpretation
%       depends on LBmode and on the selected window function.
%
%   LBmode
%       'extra' (default): combine AddLB with the experimental processing.
%       'total'          : use AddLB as the requested total exponential
%                          broadening where supported.
%
% OUTPUTS
%   F
%       Real, processed, frequency-domain spectrum as a column vector.
%
%   Skeleton
%       Structure array indexed by isotope pass and spin-system block.
%       Each element stores grouped transition frequencies, amplitudes,
%       dominant-spin assignments, and isotope-pass scaling. This avoids
%       repeating the costly Hamiltonian diagonalization during selected
%       T2/Gaussian-broadening fits.
%
% SIMULATION OVERVIEW
%   1. Restore any compressed equivalent-spin entries.
%   2. Simulate each independent spin block.
%   3. Build spin-1/2 and, where required, spin-1 operators.
%   4. Construct the chemical-shift and J-coupling Hamiltonian.
%   5. Diagonalize the Hamiltonian and calculate allowed transitions.
%   6. Group transitions that cannot be resolved on the acquired grid.
%   7. Construct the complex FID with Lorentzian/Gaussian decay.
%   8. Apply zero filling and the experimental apodization.
%   9. Fourier transform and return the real spectrum.
%
% IMPORTANT CONVENTIONS
%   Database entries use large diagonal marker values to identify
%   heteronuclei. The code distinguishes 14N, 13C, and the observed nucleus
%   from the magnitude of these markers before setting heteronuclear offsets
%   to zero for Hamiltonian construction.
%
%   Heteronuclear couplings are treated with the secular Iz*Iz term, whereas
%   homonuclear couplings use the full isotropic scalar-product operator.

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

if isfield(params, 'compressor')
    compressor = params.compressor;
else
    compressor = [];
end


dwell=expParams.dwell;  np=expParams.pts;  SFO1=expParams.SFO1;    BF1=expParams.BF1;   SF=expParams.SF;   TDeff=expParams.TDeff;
XAxis=-1/(2*dwell)+1/(dwell*np(1)):(1/(dwell*np(1))):1/(2*dwell);
F=zeros(1,np(1));

% Bruker complex data contain two real numbers per complex point, hence the
% use of TDeff/2 acquired complex samples.
FID=zeros(1,TDeff/2);
t  =0:expParams.dwell:expParams.dwell*(TDeff/2-1);


%For equivalent protons, only fit for one parameter, this is achieved in
%main script throguh compression of the w_CS_J matrix, here this step is
%reverted:

%% Restore compressed equivalent-spin parameters
[decompressed_w_CS_J,decompressed_T2,decompressed_GL]=decompMethyl(w_CS_J,T2,GL,compressor);

w_CS_J=decompressed_w_CS_J;
T2=decompressed_T2;%Reconstruct T2
GL=decompressed_GL;%Reconstruct T2

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


% Diagonal marker ranges encode isotope identity in the database.
%This finds 14N (since spin 1, needs special care)
isSpin1 = (abs(wCS_temp) > 0.5e14 & abs(wCS_temp) < 50e14);

%This finds 13C (since low natural abundance, needs special care)
idx13C  = find(abs(wCS_temp) > 0.5e13 & abs(wCS_temp) < 1.5e13);
has13C = ~isempty(idx13C);
abundance13C = 0.01108;
isotopeShiftHz = -2.3;
numC = length(idx13C);

% Any large diagonal marker denotes a heteronucleus rather than a physical
% chemical shift. obsIdx therefore identifies the detected nucleus spins.
isHet = abs(wCS_temp) > 1e3;        % 1e4 ppm is added in AddUndefinedPars, so this is the threshold
obsIdx = find(~isHet);              % observed nucleus spins (e.g. 1H)

% Put heteronucleus to 0 to avoid any problems. 
wCS_temp(isHet) = 0;

% Convert observed-nucleus chemical shifts from ppm into angular frequency,
% including acquisition-carrier and processed-reference corrections.
wCS_temp=wCS_temp*BF1-(SFO1-BF1)*1e6+(SF-BF1)*1e6;%From ppm to Hz
wCS_temp=wCS_temp*2*pi;

%Number of Spins
nSpins=length(wCS_temp);

%% Generate 13C-satellite and main-isotopologue contributions
% Run multiple passes if 13C exists, otherwise 1 pass
for pass = 1:(numC + 1)
    
w_CS_J_current=w_CS_J_temp;
if pass <= numC
        % --- SATELLITE PASSES ---
        scaleFactor = abundance13C;
        % Keep coupling for the "active" 13C, zero out all other carbons
        active_C_idx = idx13C(pass);
        other_C_indices = idx13C(idx13C ~= active_C_idx);
        w_CS_J_current(other_C_indices, :) = 0;
        w_CS_J_current(:, other_C_indices) = 0;
        isoShiftHz = isotopeShiftHz;
    else
        % --- MAIN PEAK PASS (ALL 12C) ---
        scaleFactor = 1 - (numC * abundance13C);
        % Zero out all carbon couplings
        w_CS_J_current(idx13C, :) = 0;
        w_CS_J_current(:, idx13C) = 0;
        isoShiftHz = 0;
end

%% Construct single-spin basis operators
% Sparse matrices limit memory use before Kronecker expansion into the
% complete Hilbert space.
%Pauli Matrices Spin 1/2: 
I_z = sparse([1/2 0; 0 -1/2]);  
I_y = sparse([0 -1i/2; 1i/2 0]);    
I_x = sparse([0 1/2; 1/2 0]);        
I_p = sparse([0 1; 0 0]);       
one=sparse([1,0;0,1]);

%Pauli Matrices Spin 1 (14N):
   I_zQ = sparse([1 0 0; 0 0 0; 0 0 -1]); 
I_yQ = 1/(sqrt(2)*1i) * sparse([0 1 0; -1 0 1; 0 -1 0]); 
I_xQ = 1/sqrt(2) * sparse([0 1 0; 1 0 1; 0 1 0]); 
I_pQ = sparse([0 sqrt(2) 0; 0 0 sqrt(2); 0 0 0]);
oneQ = speye(3);

%% Embed each single-spin operator in the full spin space
%Cell arrays of operators
Lx=cell(1,nSpins);  Ly=cell(1,nSpins);  Lz=cell(1,nSpins); Lp=cell(1,nSpins);
for n = 1:nSpins
    % Initialize the first spin's contribution
    if 1 == n
        if isSpin1(1)
            Lx_temp = I_xQ; Ly_temp = I_yQ; Lz_temp = I_zQ; Lp_temp = I_pQ;
        else
            Lx_temp = I_x; Ly_temp = I_y; Lz_temp = I_z; Lp_temp = I_p;
        end
    else
        if isSpin1(1)
            Lx_temp = oneQ; Ly_temp = oneQ; Lz_temp = oneQ; Lp_temp = oneQ;
        else
            Lx_temp = one; Ly_temp = one; Lz_temp = one; Lp_temp = one;
        end
    end

    % Loop through the remaining spins
    for k = 2:nSpins
        if k == n
            if isSpin1(k)
                Lx_temp = kron(Lx_temp, I_xQ);
                Ly_temp = kron(Ly_temp, I_yQ);
                Lz_temp = kron(Lz_temp, I_zQ);
                Lp_temp = kron(Lp_temp, I_pQ);
            else
                Lx_temp = kron(Lx_temp, I_x);
                Ly_temp = kron(Ly_temp, I_y);
                Lz_temp = kron(Lz_temp, I_z);
                Lp_temp = kron(Lp_temp, I_p);
            end
        else
            if isSpin1(k)
                Lx_temp = kron(Lx_temp, oneQ);
                Ly_temp = kron(Ly_temp, oneQ);
                Lz_temp = kron(Lz_temp, oneQ);
                Lp_temp = kron(Lp_temp, oneQ);
            else
                Lx_temp = kron(Lx_temp, one);
                Ly_temp = kron(Ly_temp, one);
                Lz_temp = kron(Lz_temp, one);
                Lp_temp = kron(Lp_temp, one);
            end
        end
    end
    Lx{n} = Lx_temp;
    Ly{n} = Ly_temp;
    Lz{n} = Lz_temp;
    Lp{n} = Lp_temp;
end

%% Assemble the spin Hamiltonian
% Preallocate Hamiltonian with the correct dynamic dimension
totalDim = size(Lx{1}, 1);
H = spalloc(totalDim, totalDim, (nSpins^2) + (nSpins^2));

%Chemical shift
for n = 1:nSpins
    % We apply the shift ONLY to the observed nucleus (obsIdx)
    % ismember(n, obsIdx) returns true if the current spin is the observer
    currentFreq = wCS_temp(n); 
    if ismember(n, obsIdx)
        currentFreq = currentFreq + (isoShiftHz * 2 * pi);
    end
    H = H + currentFreq * Lz{n};
end

%J-coupling
% Get row/column indices of non-zero entries
[row, col] = find(triu(w_CS_J_current, 1)); % triu(w_J,1) ensures k > n
% Loop only over non-zero entries
for idx = 1:numel(row)
    n = row(idx);
    k = col(idx);
     if isHet(n) || isHet(k)
        H = H + 2*pi*w_CS_J_current(n,k) * (Lz{n} * Lz{k});% heteronuclear: secular only
     else
        H = H + 2*pi*w_CS_J_current(n,k) * (Lx{n} * Lx{k} + Ly{n} * Ly{k} + Lz{n} * Lz{k});
     end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate spectra in freq domain%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Diagonalize the Hamiltonian and calculate spectral transitions
% Fp is the total raising operator for the observed nucleus. Transition
% intensities are evaluated after transforming Fp into the eigenbasis.

% Initialize a 4x4 matrix of zeros
Fp = zeros(totalDim, totalDim);
for ii = obsIdx(:).'
    Fp = Fp + Lp{ii}; 
end


% tic

% Matrix diagonalization is the dominant computational cost of this engine.
% Dividing by 2*pi expresses the eigenvalues in hertz rather than rad/s.
[V, D] = eig(full(H)/(2*pi),'vector');%By defining 'vector' eig already takes the diagonal of D


freq = D - D.'; % Outer subtraction
% Compute Amps(k,m) = |V(:,k)' * Fp * V(:,m)|^2 
Vinv=V';
Amps = abs(Vinv * Fp * V).^2; % Entire matrix at once
% Amps=Amps/(2^(nSpins-1));%Normalizes according to number of spins
% Normalization in case of quads
Amps = Amps / (totalDim / 2);
% toc
% tic

%% Assign line-shape parameters to each transition
% The observed-spin raising operators are transformed separately. For each
% transition, the spin with the largest contribution supplies T2 and GL.
results = cellfun(@(x) Vinv * x * V, Lp(obsIdx), 'UniformOutput', false);%Calculate all Ip operators in new frame. 

% Discard numerically negligible transitions to reduce subsequent work.
% This is an explicit approximation controlled by the relative threshold.
Amps=Amps.*(abs(Amps)>1e-5*max(max(abs((Amps)))));%Ignore small entries from diagonalization (Approximation!!!)
% Amps
AmpsIdx=find(Amps);
% length(AmpsIdx)
%Find for which spin the value is largest and take corresponding T2
for n=1:length(AmpsIdx)
    element_values = cellfun(@(x) full(x(AmpsIdx(n))), results);%Get the values for Ip in new frame
    [~, mi] = max(abs(element_values));
    max_index(n) = obsIdx(mi);   % map back to original spin index
end
% Get all T2 values at once
T2_terms = T2_temp(max_index)';
RG = GL_temp(max_index)';

freqAll=freq(AmpsIdx);



%% Merge transitions below the effective digital resolution
% Closely spaced transitions are grouped before FID generation. Their
% amplitudes are summed and their representative frequency is the
% amplitude-weighted centroid.
tolF = 2/(dwell*TDeff);%Tolerance for frequency differences. If below this value group together
fbin = round(freqAll(:) / tolF);  % integer bins

[ubins, ~, ic] = unique(fbin, 'stable');


A = Amps(AmpsIdx);
Agrp = accumarray(ic, A, [], @sum);


% amplitude-weighted representative frequency (centroid)
f_num = accumarray(ic, freqAll(:).*A, [], @sum);
funiq = f_num ./ Agrp;

% The strongest member of each merged group determines its T2 and GL,
% retaining the dominant physical line-shape assignment.
% ia = accumarray(ic, (1:numel(freqAll)).', [], @(v) v(1));
ia = zeros(length(ubins), 1);
for k = 1:length(ubins)
    group_indices = find(ic == k);
    [~, max_idx] = max(A(group_indices));
    ia(k) = group_indices(max_idx);
end
T2_terms = T2_terms(ia);
RG = RG(ia);
freqAll=funiq;


%% Generate and accumulate the time-domain signal
% Frequencies are converted back to rad/s for complex exponential evolution.
% Exponential and Gaussian factors produce the mixed line broadening.
freqAll=freqAll*2*pi;
allFIDtemp=exp(1i * freqAll * t).*exp(-t./T2_terms-0.5*(RG*t).^2);
% Weight by amplitudes and sum
FID = FID+scaleFactor *(Agrp'* allFIDtemp); 

%% Store the reusable transition skeleton
% These quantities are sufficient for downstream line-shape fitting without
% reconstructing and diagonalizing the Hamiltonian at every iteration.
% Capture the physical spin ID for each bin for fast fiting of T2 and Gaussian broadening
% max_index(ia) tells us which spin (1..nSpins) is the primary contributor to bin 'k'
Skeleton(pass, m).freqAll = funiq;
Skeleton(pass, m).Agrp = Agrp;
Skeleton(pass, m).spinMap = max_index(ia); 
Skeleton(pass, m).scaleFactor = scaleFactor;

end
clear allFFtemp max_index T2_temp T2_terms

end %End blocks loop

%% Apply zero filling and acquisition-processing windows
FID(TDeff/2 + 1:np(1))=0; %Do also zero filling according to SI F2
FID(1)=FID(1)/2;%For FT artifact

ttot  =0:expParams.dwell:expParams.dwell*(np(1)-1);

% figure(198)
% plot(ttot,real(FID),ttot,imag(FID))

% %%%%SHIMMING ISSUES!%%%%%
% %Add z-Shim polynomial distribution:
% Z1 =  0.0;           % Hz
% Z2 =  -0.8;           % Hz
% Z3 =  1.0;           % Hz
% Z4 =   -1.0;           % Hz
% X  = -0.0;    % X shim  -> linear in x
% Y  = 0.0;   % Y shim  -> linear in y
% XZ = 0;     % x*z
% YZ = 0;     % y*z
% XY = 2.0;     % x*y
% X2Y2 = 00;   % (x^2 - y^2)
% 
% % Spatial grid / weighting
% % ---------- Sample geometry (normalized cylinder) ----------
% % x,y in disk radius 1; z in [-1,1]
% Ns = 800;  % number of random spin packets (increase for smoother results)
% 
% % Uniform points in a cylinder:
% theta = 2*pi*rand(Ns,1);
% r     = sqrt(rand(Ns,1));       % sqrt for uniform area in disk
% x     = r .* cos(theta);
% y     = r .* sin(theta);
% z     = 2*rand(Ns,1) - 1;
% 
% % Weighting: uniform sample
% w = ones(Ns,1);                
% w = w / sum(w);
% 
% % ---------- Residual frequency map df(x,y,z) ----------
% df = X*x + Y*y + Z1*z ...
%    + Z2*(z.^2) + Z3*(z.^3) + Z4*(z.^4) ...
%    + XZ*(x.*z) + YZ*(y.*z) ...
%    + XY*(x.*y) + X2Y2*(x.^2 - y.^2);   % Hz
% 
% % Spatially averaged phase factor:
% phaseMat  = exp(1i * 2*pi * (df * t));      % [Nz x N]
% phase_avg = (w.' * phaseMat);             % [1 x N]
% 
% % Bad-shim FID
% % size(FID)
% % size(phase_avg)
% % size(t)
% FID = FID .* phase_avg;

% x parameterizes cosine and squared-cosine windows over acquired points.
x = linspace(0, pi/2, round(TDeff/2));

% Interpret the TopSpin WDW code: 0=no window, 1=exponential,
% 3=cosine for SSB=2, and 4=cosine-squared for SSB=2.
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
  %% Fourier transform and format the output spectrum
F=fftshift(fft(FID,np(1)));

F=real(F)';
F=circshift(F,-1,1);

%% Local helper: restore compressed equivalent-spin entries
% A compressed parameter represents multiple equivalent spins during
% fitting. This helper inserts duplicate rows/columns and relaxation values,
% then reconstructs their couplings while removing artificial self-links.
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