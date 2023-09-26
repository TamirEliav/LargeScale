function addFieldsToWorkspace(inputStruct)
    % Get the field names of the input struct
    fieldNames = fieldnames(inputStruct);

    % Add each field of the struct as a variable in the caller's workspace
    for i = 1:numel(fieldNames)
        varName = fieldNames{i};
        varValue = inputStruct.(varName);
        assignin('caller', varName, varValue);
    end
end
