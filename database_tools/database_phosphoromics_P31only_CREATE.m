%% Build the 31P metabolite simulation database
% This script assembles the metabolite definitions used by the fitting
% workflow and saves each definition as a separate variable in a MAT-file.
%
% Database conventions:
%   - Chemical shifts are expressed in ppm.
%   - Scalar J-couplings are expressed in Hz.
%   - T2 and T2err values are expressed in seconds.
%   - w_CS_J stores the spin-system topology. Diagonal entries are chemical
%     shifts; off-diagonal entries are J-couplings between spin sites.
%   - CSerr and T2err define allowed fitting uncertainty/ranges.
%   - Pivot identifies the spin site used as the reference when applicable.
%   - notes provides a human-readable metabolite identifier.

% Use to create the file: 'simulation_database_31P.mat'
%Used in the fitting script

%Creates a database with entries as topological matrices following the GISSMO formalism.


% Initialize database structure
simDB = struct();

%Single P site
simDB.OPMe3 = struct('w_CS_J', 53.22,'CSerr', 0.02, 'T2', 0.45, 'T2err', 0.1,'notes', 'OPMe3');
simDB.TMP = struct('w_CS_J', 22.872,'CSerr', 0.02, 'T2', 0.8, 'T2err', 0.1,'notes', 'TMP');
simDB.MPHO = struct('w_CS_J', 4.69,'CSerr', 0.02, 'T2', 0.6,'T2err', 0.1, 'notes', 'MPHO');
simDB.PG6 = struct('w_CS_J', 4.585,'CSerr', 0.04, 'T2', 0.45,'T2err', 0.05, 'notes', 'PG6');
simDB.G6Pc = struct('w_CS_J', 4.392,'CSerr', 0.015, 'T2', 0.4,'T2err', 0.05, 'notes', 'G6Pc');
simDB.G6Pt = struct('w_CS_J', 4.355,'CSerr', 0.015, 'T2', 0.4,'T2err', 0.05, 'notes', 'G6Pt');%Dont actually know which is which (cis vs trans)
simDB.S7P = struct('w_CS_J', 4.263,'CSerr', 0.02, 'T2', 0.3,'T2err', 0.05, 'notes', 'S7P');
simDB.M6P = struct('w_CS_J', 4.232,'CSerr', 0.02, 'T2', 0.3,'T2err', 0.05, 'notes', 'M6P');
simDB.G3P = struct('w_CS_J', 4.219,'CSerr', 0.02, 'T2',  0.5,'T2err', 0.05, 'notes', 'G3P');
simDB.PG3 = struct('w_CS_J', 4.04,'CSerr', 0.03, 'T2', 0.55,'T2err', 0.05, 'notes', 'PG3');
simDB.PEA = struct('w_CS_J', 3.748,'CSerr', 0.02, 'T2', 0.22,'T2err', 0.05, 'notes', 'PEA');
simDB.GMP = struct('w_CS_J', 3.704,'CSerr', 0.005, 'T2', 0.33,'T2err', 0.05, 'notes', 'GMP');
simDB.AMP = struct('w_CS_J', 3.698,'CSerr', 0.02, 'T2', 0.33,'T2err', 0.05, 'notes', 'AMP');
simDB.CMP = struct('w_CS_J', 3.62,'CSerr', 0.02, 'T2', 0.4,'T2err', 0.05, 'notes', 'CMP');
simDB.PCHO = struct('w_CS_J', 3.251,'CSerr', 0.02, 'T2', 0.5,'T2err', 0.1, 'notes', 'PCHO');
simDB.PO43 = struct('w_CS_J', 2.45,'CSerr', 0.02, 'T2', 0.38,'T2err', 0.05, 'notes', 'PO43');
simDB.R1P = struct('w_CS_J', 2.22,'CSerr', 0.02, 'T2', 0.35,'T2err', 0.05, 'notes', 'R1P');
simDB.Gal1P = struct('w_CS_J', 2.358,'CSerr', 0.02, 'T2', 0.55,'T2err', 0.05, 'notes', 'Gal1P');
simDB.GPE = struct('w_CS_J', 0.49,'CSerr', 0.1, 'T2', 0.35,'T2err', 0.05, 'notes', 'GPE');
simDB.GPC = struct('w_CS_J', -0.22,'CSerr', 0.02, 'T2', 0.45,'T2err', 0.05, 'notes', 'GPC');
simDB.PCR = struct('w_CS_J', -3.293,'CSerr', 0.02, 'T2', 0.2,'T2err', 0.05, 'notes', 'PCR');
simDB.PEP = struct('w_CS_J', -0.81,'CSerr', 0.02, 'T2', 0.7,'T2err', 0.1, 'notes', 'PEP');
% simDB.Aux1 = struct('w_CS_J', 0, 'T2', 0.2, 'notes', 'Aux1');

%Two Sites (If only one T2 value given the autofit will give one T2 for
%both sites of metabolite)
simDB.ADP = struct('w_CS_J', [-6.392 22.57; 0 -10.874],'CSerr', 0.03, 'T2', [0.2],'Pivot',1,  'notes', 'ADP');
simDB.UDP = struct('w_CS_J', [-6.45 22.67; 0 -11.0],'CSerr', 0.01, 'T2', [0.2],'Pivot',1,  'notes', 'UDP');
simDB.GDP = struct('w_CS_J', [-6.392 22.47; 0 -10.864],'CSerr',  [0.005 0.005], 'T2', [0.2],'Pivot',2,  'notes', 'GDP');
simDB.FADH2 = struct('w_CS_J', [-10.687 20.5; 0 -11.35],'CSerr', 0.03, 'T2', [0.15],'Pivot',1, 'notes', 'FADH2');
simDB.UDPGA = struct('w_CS_J', [-11.2 20.18; 0 -12.88],'CSerr', 0.01, 'T2', [0.15],'Pivot',2,  'notes', 'UDPGA');
simDB.UDPGal = struct('w_CS_J', [-11.245 20.85; 0 -12.76],'CSerr', 0.01, 'T2', [0.15],'Pivot',2,  'notes', 'UDPGal');
simDB.UDPG = struct('w_CS_J', [-11.322 20.77; 0 -12.945],'CSerr', 0.005, 'T2', [0.15], 'T2err', [0.02],'Pivot',2,  'notes', 'UDPG');
simDB.UDPNAcGal = struct('w_CS_J', [-11.45 21.18; 0 -13.01],'CSerr', 0.01, 'T2', [0.15],'Pivot',2,  'notes', 'UDPNAcGal');
simDB.NADplus = struct('w_CS_J', [-11.407 20.5; 0 -11.73],'CSerr', 0.02, 'T2', [0.15],'Pivot',2,  'notes', 'NADplus');
simDB.CDPCh = struct('w_CS_J', [-11.489 21.2; 0 -12.264],'CSerr', 0.03, 'T2', [0.2 0.2],'Pivot',2, 'notes', 'CDPCh');
simDB.UDPNAcGlu = struct('w_CS_J', [-11.5 21; 0 -13.17],'CSerr', 0.01, 'T2', [0.15],'Pivot',2, 'notes', 'UDPNAcGlu');
simDB.Aux2a = struct('w_CS_J',[-10.69 21.5; 0 -11.34],'CSerr', 0.02, 'T2', [0.15 0.15],  'notes', 'Aux2a');
simDB.Aux2b = struct('w_CS_J',[-10.69 21.5; 0 -11.34],'CSerr', 0.02, 'T2', [0.15 0.15],  'notes', 'Aux2b');



% Two Sites Singlets
simDB.BPG = struct('w_CS_J', [3.72 0; 0 3.345],'CSerr', 0.02, 'T2', [0.3],'Pivot',2, 'notes', 'BPG');
simDB.FBPc = struct('w_CS_J', [4.76 0; 0 3.817],'CSerr', 0.02, 'T2', [0.3],'Pivot',1, 'notes', 'FBPc');%alpha (old c)
simDB.FBPt = struct('w_CS_J', [3.87 0; 0 3.804],'CSerr', 0.02, 'T2', [0.2 0.2],'T2err', [0.02],'Pivot',2, 'notes', 'FBPt');%beta (old t)
%Two sites magnetic equivalent
simDB.CDPEA = struct('w_CS_J', [-11.325],'CSerr', 0.01, 'T2', [0.2], 'Pivot',1,'notes', 'CDPEA');
simDB.NADH = struct('w_CS_J', [-11.36],'CSerr', 0.01, 'T2', [0.2], 'notes', 'NADH');
simDB.NADPH = struct('w_CS_J', [-11.34 20 0;0 -11.34 0;0 0 3.505],'CSerr', 0.01, 'T2', [0.3],'Pivot',3, 'notes', 'NADPH');


%Three Sites (s)(d)(d)
simDB.NADPplus = struct('w_CS_J', [3.43 0 0; 0 -11.38 20;0 0 -11.75],'CSerr', 0.02, 'T2', [0.2 0.2 0.2],'Pivot',3,  'notes', 'NADPplus');
simDB.CoA = struct('w_CS_J', [3.876 0 0; 0 -10.9 20.64;0 0 -11.38],'CSerr', 0.02, 'T2', [0.2 0.2 0.2],'T2err', [0.02],'Pivot',2, 'notes', 'CoA');
%Three Sites (d)(d)(t)
simDB.ATP = struct('w_CS_J', [-6.08 0 20.05; 0 -11.18 19.75;0 0 -21.96],'CSerr', [0.03 0.03 0.1], 'T2', [0.16 0.14 0.14],'Pivot',3,  'notes', 'ATP');
simDB.UTP = struct('w_CS_J', [-6.11 0 20.05; 0 -11.28 19.75;0 0 -22.01],'CSerr', [0.02 0.02 0.02], 'T2', [0.16 0.14 0.14],'Pivot',1,  'notes', 'UTP');


fnames = fieldnames(simDB);  % Get all field names in data

save('database_phosphoromics_P31only.mat', '-struct', 'simDB');
