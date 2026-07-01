function [FitOrder,fitRules,Expnr] = createAutoFitRules()

                
Expnr=[];
%Define which experiments to run with autofit
Expnr=[491,481,471,461,451,441,431,421,411,401,391,381,371,361,351,341,331,321,311,301,291,281,271,261,251,241,231,221,211,201,191,181,171,161,151,141,131,121,111,101,91,81,71,61,51,41,31,21,11];
%801,791,781,771,761,751,741,731,721,711,701,691,681,671,661,651,641,631,621,611,601,591,581,571,561,551,541,531,521,511,501,

% Define original fit order
FitOrder = {'PO43','OPMe3', 'TMP', 'MPHO', 'PG6',  'G6Pc', 'G6Pt', ...
             'G3P','M6P', 'S7P','PG3', 'AMP', 'GMP','CMP','PEA','BPG', 'PCHO', 'R1P',...
            'Gal1P','GPE','GPC', 'PCR','PEP','FBPt','FBPc', 'GMP','ATP',{'ADP'},'UDP','UTP','GDP',...
            'NADp', 'UDPNAcGlu', 'UDPG', 'UDPNAcGal', 'UDPGal','UDPGA',  'NADPp','CDPCh','FADH2','CoA',...
            'CDPEA','NADH'};
% FitOrder = {'PO43','OPMe3', 'TMP', 'MPHO', 'PG6',  'G6Pc', 'G6Pt', ...
%              'G3P','M6P', 'S7P','PG3', 'AMP', 'GMP','CMP','PEA','BPG', 'PCHO', 'R1P',...
%             'Gal1P','GPE','GPC', 'PCR','PEP','FBPt','FBPc', 'GMP','ATP',{'ADP'},...
%             'NADH', 'UDPNAcGlu', 'UDPG', 'UDPNAcGal', 'UDPGal','UDPGA',  'NADPH','CDPCh','FADH2','CoA'};
           
 % FitOrder = {'OPMe3', 'TMP','AMP', 'GMP','PEA','BPG','ATP','ADP','UDP','UTP','GDP'};
  % FitOrder = {'OPMe3', 'TMP','ATP','ADP','UDP','UTP','GDP'};
 % FitOrder = {'OPMe3', 'TMP','AMP', 'GMP','PEA','BPG'};


    fitRules = struct();

    % Pivot chemical-shift dependency rules
    fitRules.ATP.useLBPrefit = true;
    fitRules.ATP.nonPivotCorrectionScale = 0.5;

    fitRules.ADP.shiftParents = {'ATP'};
    fitRules.ADP.shiftScale   = 1.0;
    fitRules.ADP.useLBPrefit = true;
    fitRules.ADP.fitWithLB = 0.0;
    fitRules.ADP.fitWithLBMode = 'total';% total Means, the lb value define in fitwithLB is the absolute value, removing the experiment lb. to add to experiment use 'extra'


    fitRules.GDP.startCSFrom = 'ADP';
    fitRules.GDP.startCSOffset = [0 0; 0 0.01];
    fitRules.GDP.nonPivotShiftParent = 'TMP';%The non-pivot spins are very stable. 
    fitRules.GDP.fitWithLB = 0.0;
    fitRules.GDP.fitWithLBMode = 'total';% total Means, the lb value define in fitwithLB is the absolute value, removing the experiment lb. to add to experiment use 'extra'

    fitRules.CoA.shiftParents = {'ATP'};
    fitRules.CoA.shiftScale   = 0.25;
    fitRules.CoA.nonPivotShiftParentBySpin = {'TMP', '', ''};%Spin1: link to TMP, Spin2 is pivot anyway, Spin3 link to pivot (no exception).

    fitRules.UDP.shiftParents = {'ADP'};
    fitRules.UDP.shiftScale   = 1.0;

    fitRules.NADH.shiftParents = {'ATP', 'ADP'};
    fitRules.NADH.shiftScale   = 0.25;

    fitRules.UDPNAcGlu.shiftParents = {'ATP', 'ADP'};
    fitRules.UDPNAcGlu.shiftScale   = 0.25;

    fitRules.UDPG.shiftParents = {'ATP', 'ADP'};
    fitRules.UDPG.shiftScale   = 0.25;

    fitRules.UDPNAcGal.shiftParents = {'ATP', 'ADP'};
    fitRules.UDPNAcGal.shiftScale   = 0.25;

    fitRules.UDPGal.shiftParents = {'ATP', 'ADP'};
    fitRules.UDPGal.shiftScale   = 0.25;

    fitRules.UTP.shiftParents = {'ATP', 'ADP'};
    fitRules.UTP.shiftScale   = 1.0;

    fitRules.CoA.shiftParents = {'ATP', 'ADP'};
    fitRules.CoA.shiftScale   = 1.0;

    fitRules.FADH2.shiftParents = {'ATP', 'ADP'};
    fitRules.FADH2.shiftScale   = 1.0;

    fitRules.CDPEA.shiftParents = {'ATP', 'ADP'};
    fitRules.CDPEA.shiftScale   = 0.25;

    fitRules.NADH.shiftParents = {'NADp'};
    fitRules.NADH.shiftScale   = 1;

    % 
    fitRules.AMP.fitWithLB = 0.0;
    fitRules.AMP.fitWithLBMode = 'total';% total Means, the lb value define in fitwithLB is the absolute value, removing the experiment lb. to add to experiment use 'extra'

    fitRules.GMP.shiftParents = {'AMP'};
    fitRules.GMP.shiftScale   = 1.0;
    fitRules.GMP.copyT2From = 'AMP';
    fitRules.GMP.fitWithLB = 0.0;
    fitRules.GMP.fitWithLBMode = 'total';% total Means, the lb value define in fitwithLB is the absolute value, removing the experiment lb. to add to experiment use 'extra'

    fitRules.NADPH.shiftParents = {'NADH'};
    fitRules.NADPH.shiftScale   = 1.0;

    fitRules.BPG.nonPivotShiftParent = 'AMP';%The non-pivot spins are more related to the Parent than to its pivot

    fitRules.GPE.useLBPrefit = true;

    if ~exist('Expnr', 'var')
        Expnr=[];
    end
end