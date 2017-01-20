function [ STIMULI, TARGETS ] = import_stimulus_txt(filename)
%IMPORT_STIMULUS_TXT Import trial data, and compose target structure.
%   [STIMULI] = IMPORTFILE(FILENAME) Reads data from stimlist.csv
%   and returns it as a structure with a named fields for each column.
%
%   [STIMULI,TARGETS] = IMPORTFILE(FILENAME) Reads data from stimlist.csv
%   and, additionally, composes a TARGET structure. The TARGET structure
%   will be a nSubject x nSession x nTarget structured array (22 x 4 x 8).
%
%    See also TEXTSCAN.

%% Initialize variables.
delimiter = ',';

%% Read columns of data as strings.
% For more information, see the TEXTSCAN documentation.
formatHeader = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s';
formatBody = '%s%d%d%s%s%f%d%d%d%d%d%d%d%d';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataHeader = textscan(fileID, formatHeader, 1, 'Delimiter', delimiter, 'ReturnOnError', false);
dataHeader = cellfun(@(x) x{:}, dataHeader, 'unif', 0);
dataArray = textscan(fileID, formatBody, 'Delimiter', delimiter, 'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Compose temporary structure for ease of reference
STIMULI = cell2struct(dataArray,dataHeader,2);

%% Initialize STIMULI struct, subject x session
subjects = unique(STIMULI.subject);
sessions = unique(STIMULI.session);
nSubject = numel(subjects);
nSession = numel(sessions);
targetLabels = dataHeader(end-7:end);
TARGETS = struct(...
    'subject',[],...
    'session',[],...
    'label',repmat(targetLabels',1,nSession,nSubject),...
    'type','category',...
    'sim_source',[],...
    'sim_metric',[],...
    'target',[]...
);
TARGETS = permute(TARGETS,[3,2,1]);

for iSubj = 1:nSubject
    subj = subjects(iSubj);
    for iSession = 1:nSession
        sess = sessions(iSession);
        z = (STIMULI.subject==subj) & (STIMULI.session==sess);
        for iTarg = 1:numel(targetLabels);
            tlab = targetLabels{iTarg};
            TARGETS(iSubj,iSession,iTarg).subject = subj;
            TARGETS(iSubj,iSession,iTarg).session = sess;
            TARGETS(iSubj,iSession,iTarg).label = tlab;
            switch tlab
                case 'multinomial'
                    TARGETS(iSubj,iSession,iTarg).target = STIMULI.(tlab)(z);
                otherwise
                    TARGETS(iSubj,iSession,iTarg).target = logical(STIMULI.(tlab)(z));
            end
        end
    end
end