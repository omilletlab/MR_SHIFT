%Function for MRSHIFT

%Function for calculating the spectrum from matrix diagonalization

%Creates broadening by appodization of time domain signal with Exp and Gaussian decays

%Takes into account block structure from database

%decompresses methyl entries (Compression happens in main document)

%Takes into account experimental parameters, including processing (LB)

%Generates the skeleton for fast FID generation used for T2/GB fits


function [F,Skeleton]=NSpin_Simulation_1p3(params,expParams,AddLB,LBmode)
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

if isfield(params, 'compressor')
    compressor = params.compressor;
else
    compressor = [];
end

% compressor

dwell=expParams.dwell;  np=expParams.pts;  SFO1=expParams.SFO1;    BF1=expParams.BF1;   SF=expParams.SF;   TDeff=expParams.TDeff;
XAxis=-1/(2*dwell)+1/(dwell*np(1)):(1/(dwell*np(1))):1/(2*dwell);
F=zeros(1,np(1));

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
% T2_temp

isSpin1 = (abs(wCS_temp) > 0.5e14 & abs(wCS_temp) < 50e14);%This finds 14N (since spin 1, needs special care)

idx13C  = find(abs(wCS_temp) > 0.5e13 & abs(wCS_temp) < 1.5e13);%This finds 13C (since low natural abundance, needs special care)
has13C = ~isempty(idx13C);
abundance13C = 0.01108;
isotopeShiftHz = -2.3;
numC = length(idx13C);

isHet = abs(wCS_temp) > 1e3;        % 1e4 ppm is added in AddUndefinedPars, so this is the threshold
obsIdx = find(~isHet);              % observed nucleus spins (e.g. 1H)

% Put heteronucleus to 0 to avoid any problems. 
wCS_temp(isHet) = 0;


wCS_temp=wCS_temp*BF1-(SFO1-BF1)*1e6+(SF-BF1)*1e6;%From ppm to Hz
wCS_temp=wCS_temp*2*pi;

%Number of Spins
nSpins=length(wCS_temp);

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


% Initialize a 4x4 matrix of zeros
Fp = zeros(totalDim, totalDim);
for ii = obsIdx(:).'
    Fp = Fp + Lp{ii}; 
end


% tic

[V, D] = eig(full(H)/(2*pi),'vector');%By defining 'vector' eig already takes the diagonal of D
%This is the current bottleneck in terms of time for the calculation.
%Following the idea of Kuprov it is possible to project into a space span
%by the relevant states only and then using eigs it should be possible to
%reduce the time significantly
% D = diag(D); % Extract diagonal as vector

freq = D - D.'; % Outer subtraction
% Compute Amps(k,m) = |V(:,k)' * Fp * V(:,m)|^2 
Vinv=V';
Amps = abs(Vinv * Fp * V).^2; % Entire matrix at once
% Amps=Amps/(2^(nSpins-1));%Normalizes according to number of spins
% Normalization in case of quads
Amps = Amps / (totalDim / 2);
% toc
% tic

%Use T2 from dominating single spin transition (only for homonuclei relevant, thus Lp(obsIdx). Therefore:
results = cellfun(@(x) Vinv * x * V, Lp(obsIdx), 'UniformOutput', false);%Calculate all Ip operators in new frame. 


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
% T2_terms = 0.5./(pi*T2_temp(max_index)');%Old for Lorentzian
T2_terms = T2_temp(max_index)';
RG = GL_temp(max_index)';
%And all Lorentzians at once
% allFFtemp=T2_terms./((T2_terms.^2)+(XAxis-freq(AmpsIdx)).^2);
%Gauss Lorentzian mixture:

freqAll=freq(AmpsIdx);
% size(freqAll)
% T2_terms
% Amps(AmpsIdx)


%Old script: Take all eigenvalues 
% % % allFFtemp=(RG./(T2_terms.*sqrt(2*pi)).*exp(-0.5*((XAxis-freqAll)./T2_terms).^2)...
% % %     +1/pi*(1-RG).*(T2_terms)./(T2_terms.^2+(XAxis-freqAll).^2));
% % % % size(F)
% % % % size(allFFtemp)
% % % % Weight by amplitudes and sum
% % % F = F+sum(Amps(AmpsIdx) .* allFFtemp, 1); 

%Current approach. Group frequencies that are closer than the separation
%between two points
tolF = 2/(dwell*TDeff);%Tolerance for frequency differences. If below this value group together
fbin = round(freqAll(:) / tolF);  % integer bins

[ubins, ~, ic] = unique(fbin, 'stable');


A = Amps(AmpsIdx);
Agrp = accumarray(ic, A, [], @sum);


% amplitude-weighted representative frequency (centroid)
f_num = accumarray(ic, freqAll(:).*A, [], @sum);
funiq = f_num ./ Agrp;


% pick T2, GL from member with largest amplitude of each group 
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

% allFFtemp=(RG./(T2_terms.*sqrt(2*pi)).*exp(-0.5*((XAxis-freqAll)./T2_terms).^2)...
%     +1/pi*(1-RG).*(T2_terms)./(T2_terms.^2+(XAxis-freqAll).^2));
% % Weight by amplitudes and sum
% F = F+sum(Agrp  .* allFFtemp, 1); 

freqAll=freqAll*2*pi;
allFIDtemp=exp(1i * freqAll * t).*exp(-t./T2_terms-0.5*(RG*t).^2);
% Weight by amplitudes and sum
% FID = FID+scaleFactor *sum(Agrp  .* allFIDtemp, 1); 
FID = FID+scaleFactor *(Agrp'* allFIDtemp); 
% T2_terms
% toc

% Capture the physical spin ID for each bin for fast fiting of T2 and Gaussian broadening
% max_index(ia) tells us which spin (1..nSpins) is the primary contributor to bin 'k'
Skeleton(pass, m).freqAll = funiq;
Skeleton(pass, m).Agrp = Agrp;
Skeleton(pass, m).spinMap = max_index(ia); 
Skeleton(pass, m).scaleFactor = scaleFactor;

end
clear allFFtemp max_index T2_temp T2_terms

end %End blocks loop

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