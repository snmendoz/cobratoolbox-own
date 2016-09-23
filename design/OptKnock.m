function [optKnockSol,bilevelMILPproblem] = OptKnock(model,selectedRxnList,options,constrOpt,prevSolutions,verbFlag,solutionFileNameTmp)
%OptKnock Run OptKnock in the most general form
%
% OptKnock(model,selectedRxnList,options,constrOpt,prevSolutions,verbFlag,solutionFileNameTmp)
%
%INPUTS
% model            Structure containing all necessary variables to described a
%                  stoichiometric model
%   rxns            Rxns in the model
%   mets            Metabolites in the model
%   S               Stoichiometric matrix (sparse)
%   b               RHS of Sv = b (usually zeros)
%   c               Objective coefficients
%   lb              Lower bounds for fluxes
%   ub              Upper bounds for fluxes
%   rev             Reversibility of fluxes
%
% selectedRxnList  List of reactions that can be knocked-out in OptKnock
%
% options          OptKnock options
%   targetRxn       Target flux to be maximized. Desire product formation
%                   reaction for which gene deletions will be suggested to
%                   lead to overproduction of that product. Type: String.
%
%OPTIONAL INPUTS
% options             OptKnock options (Type: structure)
%   numDel             # of deletions allowed (Default: 5)
%   numDelSense        Direction of # of deletions constraint (G/E/L)
%                      (G: Greater than; E: Equal to; L: Lower than)
%                      (Default: L)
%   vMax               Max flux (Default: 1000). See NOTES below.
%   solveOptKnock      Solve problem within Matlab (Default: true)
%
% constrOpt           Explicitly constrained reaction options (Type: structure)
%   rxnList             Reaction list (Type: cell array)
%   values              Values for constrained reactions (Type: double array)
%   sense               Constraint senses for constrained reactions (G/E/L)
%                       (G: Greater than; E: Equal to; L: Lower than)
%                       (Type: char array)
%                     Example: constrOpt=struct('rxnList',{{'EX_for_e';'EX_etoh_e'}},'values',[1;5],'sense',['G';'G']);
%
% prevSolutions       A cell array of nx1 for which each row i contains a
%                     m_ix1 array of reactions. Each array of reactions i 
%                     is a previously calculated solution of optKnock.
%                     Example: prevSolutions={{'PFL';'RPI'};{'PFL';'TKT2'}}
%
% verbFlag            Verbose flag (Type: bolean) (Default: false)
%
% solutionFileNameTmp File name for storing temporary solutions
%
%OUTPUTS
% optKnockSol           OptKnock solution structure
%   rxnList              Reaction KO list
%   fluxes               Flux distribution
% bilevelMILPproblem    optKnock problem structure
%
%NOTES
% OptKnock uses bounds of -vMax to vMax or 0 to vMax for reversible and
% irreversible reactions. If you wish to constrain a reaction, use
% constrOpt.
%
% Markus Herrgard 3/28/05
% Richard Que (04/27/10) - Added some default parameters.
% Sebastian Mendoza (22/09/2016) - revision of code, check the correct 
% functioning, improvement of comments, unnecessary code lines were removed.

% Set these for MILP callbacks
global MILPproblemType;
global selectedRxnIndIrrev;
global rxnList;
global irrev2rev;
global solutionFileName;
global biomassRxnID;
global OptKnockKOrxnList;
global OptKnockObjective;
global OptKnockGrowth;
global solID;

% empty variable handling
if (isempty(model)); error('OptKnock: No model specified'); end

if (~exist('options','var') || isempty(options) ); error('OptKnock: No target reaction specified'); 
else
    if ~isfield(options,'vMax'), options.vMax = 1000; end
    if ~isfield(options,'numDel'), options.numDel = 5; end
    if ~isfield(options,'numDelSense'), options.numDelSense = 'L'; end
    if ~isfield(options,'solveOptKnock'), options.solveOptKnock = true; end
end

if ~exist('constrOpt','var')
    constrOpt.values = [];
    constrOpt.sense = [];
    constrOpt.rxnList = [];
end

if (nargin < 5); prevSolutions = []; end
if (nargin < 6); verbFlag = false; end
if (nargin < 7)
    solutionFileName = 'optKnockSolutions.mat';
else
    solutionFileName = solutionFileNameTmp;
end

% Convert to irreversible rxns
[modelIrrev,matchRev,rev2irrev,irrev2rev] = convertToIrreversible(model);

% Create the index of the previous KO's suggested by OptKnock to avoid obtaining the same
% solution again
selPrevSolIrrev = [];
for i = 1:length(prevSolutions)
    prevSolRxnList = prevSolutions{i};
    if ~isempty(prevSolRxnList)
        selPrevSol = ismember(model.rxns,prevSolRxnList);
        selPrevSolIrrev(:,i) = selPrevSol(irrev2rev);
    end
end

[~,nRxns] = size(modelIrrev.S);

% Create matchings for reversible reactions in the set selected for KOs
% This is to ensure that both directions of the reaction are knocked out
selSelectedRxn = ismember(model.rxns,selectedRxnList);
selSelectedRxnIrrev = selSelectedRxn(irrev2rev);
selectedRxnIndIrrev = find(selSelectedRxnIrrev);
cnt = 0;
nSelected = length(selectedRxnIndIrrev);
selRxnCnt = 1;
while selRxnCnt <= nSelected
    rxnID = selectedRxnIndIrrev(selRxnCnt);
    if (matchRev(rxnID)>0)
        cnt = cnt + 1;
        selectedRxnMatch(cnt,1) = selRxnCnt;
        selectedRxnMatch(cnt,2) = selRxnCnt+1;
        selRxnCnt = selRxnCnt + 1;
    end
    selRxnCnt = selRxnCnt + 1;
end

% Set inner constraints for the LP
constrOptIrrev = setConstraintsIrrevModel(constrOpt,model,modelIrrev,rev2irrev);

% Set objectives for linear and integer parts
cLinear = zeros(nRxns,1);
cInteger = zeros(sum(selSelectedRxnIrrev),1);

% Set the correct objective coefficient
targetRxnID = find(ismember(model.rxns,options.targetRxn));
targetRxnIDirrev = rev2irrev{targetRxnID}(1);
cLinear(targetRxnIDirrev) = 1;

% Create the constraint matrices for the bilevel MILP
bilevelMILPproblem = createBilevelMILPproblem(modelIrrev,cLinear,cInteger,selSelectedRxnIrrev,...
    selectedRxnMatch,constrOptIrrev,[],options,selPrevSolIrrev);

% Initialize initial solution x0
bilevelMILPproblem.x0 = [];

% Maximize
bilevelMILPproblem.osense = -1;

% If verbose, print number of constraints, integer and continous variables 
if (verbFlag)
    [nConstr,nVar] = size(bilevelMILPproblem.A);
    nInt = length(bilevelMILPproblem.intSolInd);
    fprintf('MILP problem with %d constraints %d integer variables and %d continuous variables\n',...
        nConstr,nInt,nVar);
end

% Set model for MILP problem
bilevelMILPproblem.model = modelIrrev;

% Set these for CPLEX callbacks
MILPproblemType = 'OptKnock';
rxnList = model.rxns;
biomassRxnID = find(modelIrrev.c==1);
solID = 0;
OptKnockObjective = [];
OptKnockGrowth = [];
OptKnockKOrxnList = {};

% Solve problem
if (options.solveOptKnock)
    optKnockSol = solveCobraMILP(bilevelMILPproblem,'printLevel',verbFlag);
    if (~isempty(optKnockSol.cont))
        % Find fluxes of network with gene deletions
        optKnockSol.fluxes = convertIrrevFluxDistribution(optKnockSol.cont(1:length(matchRev)),matchRev);
    end
    if (~isempty(optKnockSol.int))
        % Figure out the KO reactions
        optKnockRxnInd = selectedRxnIndIrrev(optKnockSol.int < 1e-4);
        optKnockSol.rxnList = model.rxns(unique(irrev2rev(optKnockRxnInd)));
    end
else
    optKnockSol.rxnList = {};
    optKnockSol.fluxes = [];
end