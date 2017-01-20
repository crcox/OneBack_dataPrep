function [beta,xyz,ijk] = import_beta_txt(filename,nSession,sessionSize)
%IMPORTFILE Import numeric data from a text file as column vectors.
%   [beta,ijk,xyz] = IMPORTFILE(FILENAME) Reads MRI data from text file
%   FILENAME.
%
%   See also TEXTSCAN.

    %% Initialize variables.
    delimiter = ' ';

    %% Format string for each line of text:
    %   beta: double (%f)
    %    xyz: double (%f)
    %	 ijk: integer (%d)

    %% Open the text file.
    fileID = fopen(filename,'r');
    tLines = fgets(fileID); frewind(fileID);
    numCols = numel(strfind(tLines,delimiter)) + 1;
    formatSpec = ['%d%d%d%f%f%f',repmat('%f',1, numCols-6)];

    %% Read columns of data according to format string.
    % This call is based on the structure of the file used to generate this
    % code. If an error occurs for a different file, try regenerating the code
    % from the Import Tool.
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', false, 'EmptyValue' ,NaN, 'HeaderLines', 0, 'ReturnOnError', false);

    %% Close the text file.
    fclose(fileID);

    %% Allocate imported array to column variable names
    ijk = cell2mat(dataArray(:, 1:3));
    xyz = cell2mat(dataArray(:, 4:6));
    n = nSession * sessionSize;
    x = dataArray(:, 7:(n+6));
    beta = mat2cell(cell2mat(x),size(ijk,1),repmat(sessionSize,1,nSession));
end