%addpath(genpath('C:\Program Files\MATLAB\R2015b\toolbox\cobra'));
%initCobraToolbox

load('geni_names.mat'); % Name of genes
load('reaction_expression.mat'); % Reaction expressions
load('yeastmm.mat');
load('pos_genes_in_react_expr.mat'); % Position of each of the genes in the reaction expression var ^^
load('ixs_geni_sorted_by_length.mat');
load('msbdataAltered.mat'); 
load('growthRatesMSB.mat');

changeCobraSolver('pdco', 'LP');
changeCobraSolver('pdco', 'QP');

growthLoc = 3487; %(BIOMASS rxn index))

growthReactName = model.rxns(growthLoc);
model = changeObjective(model, growthReactName); %Set the objetive to be the biomass

genes = model.genes;
genes_in_dataset = msbdata.geneName;

NumberOfObjectives = 1 % Number of objectives
NumberOfGenes = numel(genes); % Number of variables

%% Flux distribution control
%%
GeneExpressionArray = ones(numel(genes),1);  % We start from the all-one configuration
%[v1_control, f_out_control] = evaluate_objective(1, GeneExpressionArray,NumberOfObjectives,NumberOfGenes,model,genes,reaction_expression,pos_genes_in_react_expr,ixs_geni_sorted_by_length)


%% Flux distribution disease
%%
maxcorcof = 0;
bestGamma = 0;
corrcoefs = [];

gamma = [1, 2, 4, 8, 16]; % Exploring the different gamma values 

for s = 1 : numel(gamma) % Try each of the gamma values in metrade
    lowerPercentileBounded = [1]; % We initially explored different low bound values and settled on 1
    for lowbound = 1: numel(lowerPercentileBounded)
    bounds = [-50];
    lowerPercentBounded = lowerPercentileBounded(lowbound);
   for j = 1 : numel(bounds)
    model = setMediaConditions(model, bounds(j));
    fluxes = model.rxnNames;
    for t=3:(width(msbdata))
        expr_profile = table2array(msbdata(:,t));
        deletedGeneIndex = find(strcmp(msbdata.commonName,msbdata.Properties.VariableNames{t}),1);
        deletedGeneExpressionIndexInArray = inf;
        pos_genes_in_dataset = zeros(numel(genes),1);
        for i=1:numel(genes)
            position = find(strcmp(genes{i},genes_in_dataset));
            if ~isempty(position)
                pos_genes_in_dataset(i) = find(strcmp(genes{i},genes_in_dataset));
                if pos_genes_in_dataset(i) == deletedGeneIndex
                    deletedGeneExpressionIndexInArray = i;
                end
                GeneExpressionArray(i) = expr_profile(pos_genes_in_dataset(i));
            end
        end
        
        NumberOfGenes = numel(genes);
      
        [v1, f_out] = evaluate_objective(gamma(s),GeneExpressionArray,NumberOfObjectives,NumberOfGenes, model,genes,reaction_expression,pos_genes_in_react_expr,ixs_geni_sorted_by_length, lowerPercentBounded, deletedGeneExpressionIndexInArray);
        f_out

        try
        if length ( v1(:,1) ) > 1 
            fluxes = [fluxes,array2table(v1,'VariableNames',{msbdata.Properties.VariableNames{t}})];
        end
        catch 
        end
 
      end
    %-------------------------------------------------------------------------------------------
    % Extract just the growth rate from fluxes 
    justGrowth =  fluxes(growthLoc,2:end);
    justGrowthColumn = table2array(justGrowth);
    justGrowthColumn = array2table(justGrowthColumn.');
    justGrowthColumn = [justGrowth.Properties.VariableNames', justGrowthColumn];
    justGrowthColumn.Properties.VariableNames{1} = 'commonName';
    growthRatesMSB.commonName = cellstr(growthRatesMSB.commonName);
    joined = innerjoin(justGrowthColumn, growthRatesMSB);
    joinedForComp = table2array(joined(:,2:end));
    
    [x p] = corrcoef(joinedForComp(:,1) , joinedForComp(:,2));
    x = x(2)
    corrcoefs = [corrcoefs , x];
    if (x < maxcorcof) % We keep the gamma value with the best correlation and save the flux rates
        maxcorcof = x;
        bestFlux = fluxes; % The flux rates
        bestGamma = gamma(s);
        bestLowerBound = lowerPercentileBounded(lowbound)
    end
                end
          
    end
end

%Here we are organising the data and plotting the correlated results 
    justGrowthColumn = table2array(justGrowth);
    justGrowthColumn = array2table(justGrowthColumn.');
    justGrowthColumn = [justGrowth.Properties.VariableNames', justGrowthColumn];
    justGrowthColumn.Properties.VariableNames{1} = 'commonName';
    growthRatesMSB.commonName = cellstr(growthRatesMSB.commonName);
    joined = innerjoin(justGrowthColumn, growthRatesMSB);
    joinedForComp = table2array(joined(:,2:end));
    scatter(joinedForComp(:,1),joinedForComp(:,2));
    corr(joinedForComp(:,1),joinedForComp(:,2))
    xlabel('Biomass growth rate captured in FBA');
    ylabel ('doubling time change log2(Strand/W.T)');
    p = polyfit(joinedForComp(:,1),joinedForComp(:,2),1)
    f = polyval(p,joinedForComp(:,1));
    hold on
    plot(joinedForComp(:,1),f,'--r', 'LineWidth', 3)

%Setting the initial media conditions, adjusting the uptake rates
function model = setMediaConditions(model,bound)
exchangeReactions = {'ammonium exchange'
    'sulphate exchange'
    'biotin exchange'
    '(R)-pantothenate exchange'
    'folic acid exchange'
    'myo-inositol exchange'
    'nicotinate exchange'
    '4-aminobenzoate exchange'
    'pyridoxine exchange'
    'H+ exchange'
    'riboflavin exchange'
    'thiamine(1+) exchange'
    'sulphate exchange'
    'potassium exchange'
    'phosphate exchange'
    'sulphate exchange'
    'sodium exchange'
    'L-alanine exchange'
    'L-arginine exchange'
    'L-asparagine exchange'
    'L-aspartate exchange'
    'L-cysteine exchange'
    'L-glutamate exchange'
    'L-glutamine exchange'
    'glycine exchange'
    'L-histidine exchange'
    'L-isoleucine exchange'
    'L-leucine exchange'
    'L-lysine exchange'
    'L-methionine exchange'
    'L-phenylalanine exchange'
    'L-proline exchange'
    'L-serine exchange'
    'L-threonine exchange'
    'L-tryptophan exchange'
    'L-tyrosine exchange'
    'L-valine exchange'
    'oxygen exchange'
    'adenine exchange'
    'uracil exchange'};
    for i = 1 : numel(exchangeReactions)
         model.lb(model.rxns{find(ismember(model.rxnNames, exchangeReactions{i}),1)}) = bound;
    end
end
