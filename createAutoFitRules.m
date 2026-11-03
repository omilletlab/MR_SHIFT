function [FitOrder,fitRules,Expnr, additionalFits,mainFitSteps,autoFitOptions] = createAutoFitRules()
% CREATEAUTOFITRULES Define the default automated-fitting strategy used by
% MRSHIFT.
%
% This function centralizes the configuration that controls:
%   - which experiments are processed automatically;
%   - the order in which metabolites are introduced into the fit;
%   - metabolite-specific chemical-shift and line-shape dependencies;
%   - the sequence of the main fitting stages;
%   - additional refinement passes after the main fit.
%
% OUTPUTS
%   FitOrder
%       Ordered cell array of metabolite names. Earlier entries are fitted
%       before later entries so that reliable metabolites can provide
%       starting values or constraints for dependent metabolites.
%
%   fitRules
%       Structure whose fields are metabolite names. Each metabolite field
%       may contain one or more optional rules, including:
%
%         shiftParents
%             Previously fitted metabolites whose shift changes are used to
%             estimate the starting shift of the current metabolite.
%
%         shiftScale
%             Scaling applied to the parent-derived shift correction.
%
%         startCSFrom / startCSOffset
%             Initialize a complete chemical-shift matrix from another
%             metabolite and then add the supplied matrix offset.
%
%         nonPivotCorrectionScale
%             Scale applied to non-pivot shift corrections.
%
%         nonPivotShiftParent / nonPivotShiftParentBySpin
%             Override the usual pivot-based correction for all non-pivot
%             spins or for selected spins individually.
%
%         skipNonPivotShiftCorrection
%             Disable the usual non-pivot shift correction.
%
%         copyT2From
%             Reuse T2 values from a previously fitted metabolite.
%
%         useLBPrefit
%             Include the metabolite in the preliminary line-broadening fit.
%
%         fitWithLB / fitWithLBMode
%             Override line broadening during fitting. In 'total' mode,
%             fitWithLB is the requested absolute broadening; in 'extra'
%             mode it is added to the experimental broadening.
%
%   Expnr
%       Optional list of experiment numbers to process automatically. An
%       empty array leaves experiment selection to the calling workflow.
%
%   additionalFits
%       Ordered structure array describing refinement passes. Common fields:
%         names    - Metabolites included in the pass.
%         config   - Parameter subset fitted in that pass.
%         strategy - Fitting strategy used by the caller.
%         label    - Human-readable pass label.
%         MaxIter  - Optional optimizer iteration limit.
%
%   mainFitSteps
%       Ordered names of the principal fitting stages.
%
%   autoFitOptions
%       Reserved structure for global autofit options. It is initialized
%       here so new options can be added without changing the function
%       interface.

autoFitOptions = struct();

                
Expnr=[];
%Define which experiments to run with autofit
% Expnr=[261, 271, 381, 391, 581, 591, 641, 651, 861, 871, 931, 941, 1121, 1131, 1171, 1181, 1221, 1231, 2221, 2231, 2241, 2251, 2261, 2271, 2281, 2901, 2911, 2921, 2931, 2941, 2951, 2961, 2971, 3691, 3701, 3711, 3721, 3731, 3751, 3761, 3771, 4141, 4151, 4161, 4171, 4181, 4191, 4201, 4211];
Expnr=[881, 951, 1141, 1191, 1241, 2171, 2181, 2191, 2201, 2211, 2781, 2791, 2801, 2811, 3581, 3591, 3601, 3611, 3621, 4261, 4271, 4281, 4291];

% Define fit order
FitOrder = {'PO43','OPMe3', 'TMP', 'MPHO', 'PG6',  'G6Pc', 'G6Pt', ...
             'G3P','M6P', 'S7P','PG3', 'AMP', 'GMP','CMP','PEA','BPG', 'PCHO', 'R1P',...
            'Gal1P','GPE','GPC', 'PCR','PEP','FBPt','FBPc','ATP',{'ADP'},'UDP','UTP','GDP',...
            'NADplus', 'UDPNAcGlu', 'UDPG', 'UDPNAcGal', 'UDPGal','UDPGA',  'NADPplus','CDPCh','FADH2','CoA',...
            'CDPEA','NADH'};


    fitRules = struct();

    fitRules.ATP.useLBPrefit = true;
    fitRules.ATP.nonPivotCorrectionScale = 0.5;

    fitRules.ADP.shiftParents = {'ATP'};
    fitRules.ADP.shiftScale   = 1.0;
    fitRules.ADP.useLBPrefit = true;
    fitRules.ADP.fitWithLB = 0.0;
    fitRules.ADP.fitWithLBMode = 'total';% total Means, the lb value define in fitwithLB is the absolute value, removing the experiment lb. to add to experiment use 'extra'


    fitRules.GDP.startCSFrom = 'ADP';
    fitRules.GDP.startCSOffset = [0 0; 0 0.01];%This means that i am using the ADP CS matrix as a whole and adding this as offset. instead of the regular addition of an offset.
    fitRules.GDP.skipNonPivotShiftCorrection = true;
    fitRules.GDP.fitWithLB = 0.0;
    fitRules.GDP.fitWithLBMode = 'total';% total Means, the lb value define in fitwithLB is the absolute value, removing the experiment lb. to add to experiment use 'extra'
    fitRules.GDP.copyT2From = 'ADP';
    fitRules.GDP.shiftScale   = 1.0;


    fitRules.CoA.shiftParents = {'ATP'};
    fitRules.CoA.shiftScale   = 0.25;
    fitRules.CoA.nonPivotShiftParentBySpin = {'TMP', '', ''};%Spin1: link to TMP, Spin2 is pivot anyway, Spin3 link to pivot (no exception).

    fitRules.UDP.shiftParents = {'ADP'};
    fitRules.UDP.shiftScale   = 1.0;
    fitRules.UDP.copyT2From = 'ADP';


    fitRules.NADplus.shiftParents = {'ADP', 'ATP'};
    fitRules.NADplus.shiftScale   = 0.25;

    fitRules.UDPNAcGlu.shiftParents = {'ADP', 'ATP'};
    fitRules.UDPNAcGlu.shiftScale   = 0.25;

    fitRules.UDPG.shiftParents = {'ADP', 'ATP'};
    fitRules.UDPG.shiftScale   = 0.25;

    fitRules.UDPNAcGal.shiftParents = {'ADP', 'ATP'};
    fitRules.UDPNAcGal.shiftScale   = 0.25;

    fitRules.UDPGal.shiftParents = {'ADP', 'ATP'};
    fitRules.UDPGal.shiftScale   = 0.25;

    fitRules.UTP.shiftParents = {'ATP', 'ADP'};
    fitRules.UTP.shiftScale   = 1.0;

    fitRules.CoA.shiftParents = {'ADP', 'ATP'};
    fitRules.CoA.shiftScale   = 1.0;

    fitRules.FADH2.shiftParents = {'ADP', 'ATP'};
    fitRules.FADH2.shiftScale   = 1.0;

    fitRules.CDPEA.shiftParents = {'ADP', 'ATP'};
    fitRules.CDPEA.shiftScale   = 0.25;
    fitRules.CDPEA.copyT2From = 'ADP';


    fitRules.NADH.shiftParents = {'NADplus'};
    fitRules.NADH.shiftScale   = 1;
    fitRules.NADH.copyT2From = 'NADplus';


    % 
    fitRules.AMP.fitWithLB = 0.0;
    fitRules.AMP.fitWithLBMode = 'total';% total Means, the lb value define in fitwithLB is the absolute value, removing the experiment lb. to add to experiment use 'extra'

    % fitRules.GMP.shiftParents = {'AMP'};
    % fitRules.GMP.shiftScale   = 1.0;
    fitRules.GMP.copyT2From = 'AMP';
    fitRules.GMP.fitWithLB = 0.0;
    fitRules.GMP.fitWithLBMode = 'total';% total Means, the lb value define in fitwithLB is the absolute value, removing the experiment lb. to add to experiment use 'extra'

    fitRules.NADPplus.shiftParents = {'NADplus'};
    fitRules.NADPplus.shiftScale   = 1.0;

    fitRules.BPG.nonPivotShiftParent = 'AMP';%The non-pivot spins are more related to the Parent than to its pivot

    fitRules.GPE.useLBPrefit = true;

% =====================================================
% Main fitting routine strategy
% =====================================================
% The script executes these stages from left to right:
%   LBPrefit   - preliminary line-broadening estimation;
%   Amp        - amplitudes only;
%   AmpCS      - amplitudes and chemical shifts;
%   AmpT2      - amplitudes and transverse relaxation;
%   NonPivotCS - final refinement of non-pivot chemical shifts.
mainFitSteps = {'LBPrefit', 'Amp', 'AmpCS', 'AmpT2', 'NonPivotCS'};

% =====================================================
% Additional fitting passes
% =====================================================
% Each structure element defines one post-main-fit refinement pass. Their
% order is significant because later passes use the results of earlier ones.
additionalFits(1).names = {'FBPt'};
additionalFits(1).config = 'AmpT2';
additionalFits(1).strategy = 'PivotStandard';
additionalFits(1).label = {'FBPt'};

additionalFits(2).names = {'UDPG'};
additionalFits(2).config = 'CS';
additionalFits(2).strategy = 'PivotStandard';
additionalFits(2).label = {'UDPG'};

additionalFits(3).names = {'ATP'};
additionalFits(3).config = 'CST2GL';
additionalFits(3).strategy = 'PivotScan';
additionalFits(3).label = {'ATP'};

additionalFits(4).names = {'ADP'};
additionalFits(4).config = 'AmpT2';
additionalFits(4).strategy = 'PivotScan';
additionalFits(4).label = {'ADP'};


fitNames = flattenNameList(FitOrder);

copyT2Names = fitNames(cellfun(@(n) ...
    isfield(fitRules, n) && isfield(fitRules.(n), 'copyT2From'), fitNames));

additionalFits(5).names = setdiff(fitNames, [{'ATP','ADP'}, copyT2Names], 'stable');
additionalFits(5).config = 'AmpT2GL';
additionalFits(5).MaxIter = 80;
additionalFits(5).strategy = 'GaussianBroadening';
additionalFits(5).label = {'Gaussian broadening'};

additionalFits(6).names = copyT2Names;
additionalFits(6).config = 'Amp';
additionalFits(6).strategy = 'PivotStandard';
additionalFits(6).label = {'Rest Amp'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    if ~exist('Expnr', 'var')
        Expnr=[];
    end

    function namesOut = flattenNameList(namesIn)

    namesOut = {};

    if isempty(namesIn)
        return
    end

    if ischar(namesIn)
        namesOut = {namesIn};
        return
    end

    for n = 1:numel(namesIn)
        if iscell(namesIn{n})
            namesOut = [namesOut, flattenNameList(namesIn{n})]; %#ok<AGROW>
        else
            namesOut{end+1} = char(namesIn{n}); %#ok<AGROW>
        end
    end
    end

end