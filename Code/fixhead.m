%Fixes the issues with the headers so that they simply have the gene name
%which is knocked out in the mutation 
function [fixheaders] = fixhead(dataset)
    %First we removed the excess information from the columns and capitalise 
    for k=3:c
         dataset.Properties.VariableNames{k} = upper(extractBefore(char(dataset.Properties.VariableNames(k)),"_"));
    end
    fixheaders = dataset;
end