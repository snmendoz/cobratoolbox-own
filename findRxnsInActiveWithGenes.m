function Reaclist = findRxnsInActiveWithGenes(model,genes)
% findRxnsInActiveWithGenes finds all reactions for which the provided genes
% show sufficient evidence of their absence. (i.e. make the corresponding GPRs always be zero)
%
% [Reaclist] = findRxnsActiveWithGenes(model,genes)
%
%INPUT
% model             COBRA model structure
% genes             A list of gene identifiers
%
%OUTPUTS
% Reaclist          Cell array of reactions which are supported by the
%                   provided genes 
%
% Only reactions which do have a GPR are considered for this function.
% Reactions without GPRs are ignored, as we don't have evidence.
%
% 04/03/15 Thomas Pfau

%get the positions of the provided genes (anything not in the model will
%simply be ignored)

genepos = find(ismember(model.genes,genes));
Reaclist = {};
%Set up the x vector for evaluation of the GPR rules
x = ones(size(model.genes));
x(genepos) = 0;

%Evaluate all reactions (ignoring those which have no GPR i.e. where we
%have no idea
for i=1:numel(model.rules)
    if ~isempty(model.rules{i})        
        res = eval(model.rules{i});
        if ~res 
            Reaclist{end+1} = model.rxns{i};
        end
    end
end

