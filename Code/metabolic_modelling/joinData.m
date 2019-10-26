load('msbdataAltered.mat');
load('growthRatesMSB.mat');
load('minNormFluxWithRxnRef.mat');


msbdata(:,2:3) = [];
expressionData = table2array(msbdata(:, 2:end));
expressionData = array2table(expressionData.'); %transposes data
expressionData.Properties.RowNames = editStructureToBecomeHeader(msbdata.Properties.VariableNames(2:end)); %adds the old column names as row names
expressionData.Properties.VariableNames = editStructureToBecomeHeader(msbdata{:,1}); % add old row names as column names

metexpressionData = table2array(metabolic_genes(:, 2:end));
metexpressionData = array2table(metexpressionData.'); %transposes data
metexpressionData.Properties.RowNames = editStructureToBecomeHeader(genesInMsbAndMeta.Properties.VariableNames(2:end)); %adds the old column names as row names
metexpressionData.Properties.VariableNames = editStructureToBecomeHeader(genesInMsbAndMeta{:,1}); % genes_in_dataset = msbdata.geneName;add old row names as column names

fluxData = table2array(bestFlux(:, 2:end));
fluxData = array2table(fluxData.');%transposes data
fluxData.Properties.RowNames = editStructureToBecomeHeader(bestFlux.Properties.VariableNames(2:end)); %adds the old column names as row names
fluxData.Properties.VariableNames = editStructureToBecomeHeader(bestFlux{:,1});% add old row names as column names

fluxAndExpressionData = join( fluxData, expressionData, 'keys', 'RowNames');

% Readies the growth rate data
growthRates = growthRatesMSB(3:end, 2);
growthRateRowNames = editStructureToBecomeHeader(growthRatesMSB{3:end, 1});
growthRates.Properties.RowNames = growthRateRowNames;

%Since there were some min-norms with no solutions these are removed 
knockOutsWithNoFluxData = setdiff(growthRateRowNames, fluxAndExpressionData.Properties.RowNames);

for x = 1: numel(knockOutsWithNoFluxData)
    index = find(strcmp(growthRates.Properties.RowNames , knockOutsWithNoFluxData{x}),1)
    growthRates(index,:) = [];
end

%Joining the rest of the data 
completeDataSet = join(   growthRates,fcolumnsluxAndExpressionData , 'keys', 'RowNames');
metabolic_data = join(   growthRates, metexpressionData , 'keys', 'RowNames');

justGeneExpressionDataset = completeDataSet(:,3496:end);
justGeneExpressionDataset = [completeDataSet(:,1), justGeneExpressionDataset];

justFluxDataset = completeDataSet(:,1:3495);


%Makes adjustments to the structure of the header removing disallowed
%characters 
function res = editStructureToBecomeHeader (colHeaders)
    colHeaders = strrep(colHeaders,'-','_');
    colHeaders = strrep(colHeaders,"'",'_');
    colHeaders = strrep(colHeaders,'-','_');
    colHeaders = strrep(colHeaders,'(','_');
    colHeaders = strrep(colHeaders,')','_');
    colHeaders = strrep(colHeaders,'[','_');
    colHeaders = strrep(colHeaders,',','_');
    colHeaders = strrep(colHeaders,' ','_');
    colHeaders = strrep(colHeaders,':','_');
    colHeaders = strrep(colHeaders,'___','_');
   colHeaders = strrep(colHeaders,'__','_');
    
    for x = 1 : numel(colHeaders)
        word = char(colHeaders(x));
        if ~isletter(word(1)) 
            colHeaders = insertBefore(colHeaders,1,"AMEND");
        end
    end
    
    
    colHeaders = cellstr(colHeaders);

    res = colHeaders;
end



