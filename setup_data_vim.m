function setup_data(varargin)
    %% Setup
    OVERWRITE = false;
    DATA_DIR = '/media/chris/Data/MRI/PictureOneBack/data';
    SRC_DIR = '/home/chris/src/OneBack_dataPrep';
    SUBJ_DIR = arrayfun(@(i) sprintf('MRH026_4%02d', i), 1:22, 'unif', 0);
    TXT_DIR = cellfun(@(sdir) fullfile(DATA_DIR,sdir,'analysis_RSA_corr'), SUBJ_DIR, 'unif', 0)';
    STIM_DIR = fullfile(SRC_DIR,'stimuli');
    COORD_DIR = fullfile(SRC_DIR,'coords');
    SIM_DIR = fullfile(SRC_DIR,'similarity');
    SUBJECTS = 1:22; % subjects 11 and 12 do not have labeled coordinates.

    AverageOverSessions = true;

    if nargin > 1
      args = reshape(varargin, 2, [])';
      for iArg = 1:size(args,1)
        key = args{iArg,1};
        value = args{iArg,2};
        switch lower(key)
        case 'average'
          AverageOverSessions = value;
        end
      end
    end

    %% Define Output directory
    if AverageOverSessions == 1
        bdir = 'avg';
    else
        bdir = 'full';
    end

    DATA_DIR_OUT = fullfile(...
        '/home/chris/MRI/PictureOneBack/data',...
        bdir...
    );

    if ~exist(DATA_DIR_OUT,'dir')
        mkdir(DATA_DIR_OUT)
    end

    fprintf('Looking for source MRI data in:\n\t%s\n',DATA_DIR);
    fprintf('Looking for stimulus labels and trial orders in:\n\t%s\n',STIM_DIR);
    fprintf('Looking for looking for (basal) label-to-coordinate mapping in:\n\t%s\n',COORD_DIR);
    fprintf('Looking for similarity structures in:\n\t%s\n\n',SIM_DIR);

    fprintf('Average over sessions: %d\n', AverageOverSessions);
    fprintf('Processed data will be written to:\n\t%s\n', DATA_DIR_OUT);

    %% Read presentation order and stim labels from file
    % All subjects have the same order
    file_stim_key = fullfile(STIM_DIR, 'labels.txt');

    fid = fopen(file_stim_key);
    stim_key = textscan(fid,'%s');
    fclose(fid);

    nItems = numel(stim_key{1});
    nsessions = numel(unique(stim_order{1}));

    %% Generate sort indexes
    % We do not need them

    %% load similarity structure
    % NEXT (judge similarity in kind based on word)
    % The rows in the embedding are already in an order that matches the order
    % in the stim_key.
    embedding = csvread(fullfile(SIM_DIR,'NEXT_CK_KIND_5D.csv'));
    if usejava('jvm')
      plot_similarity_decompositions(embedding);
    end
    S = corr(embedding','type','Pearson');

    %% Check dimensionality of decomposition
    addpath('~/src/WholeBrain_RSA/src');
    tau = 0.2;
    [~,r] = sqrt_truncate_r(S, tau);
    fprintf('-----\n');
    fprintf('To approximate S with error tolerance %.2f, %d dimensions are required.\n', tau, r);
    fprintf('-----\n');

    %% Prep metadata structure
    metadata = struct(...
        'subject',num2cell(SUBJECTS),...
        'targets',struct(),...
        'stimuli',stim_key(1),...
        'filters',struct(),...
        'coords',struct(),...
        'cvind',[],...
        'nrow',[],...
        'ncol',[],...
        'sessions',[],...
        'AverageOverSessions',[]...
    );

    %% Define metadata
    NSUBJ = numel(SUBJECTS);
    NCOND = 2;
    animate = [zeros(50,1);ones(50,1)];

    % TARGETS
    TARGETS = struct(...
        'label', {'visual','semantic'},...
        'type', {'similarity','similarity'},...
        'sim_source',{[],'NEXT'},...
        'sim_metric',{[],'correlation'},...
        'target',{animate,S}...
    );

    % FILTERS
    FILTERS = generate_outlier_filters(X);

    % CV
    nSchemes = 10;
    nFolds = 10;
    SCHEMES = zeros(nItems, nSchemes);
    for iScheme = 1:nSchemes
        c = cvpartition(animate,'KFold', nFolds);
        for iFold = 1:nFolds
            SCHEMES(:,iScheme) = SCHEMES(:,iScheme) + (test(c, iFold) * iFold);
        end
    end

    for iSubj = 1:NSUBJ
        % DATA
        fpath = fullfile(TXT_DIR{iSubj}, 'all+tlrc.txt');
        [beta,xyz,ijk] = import_beta_txt(fpath);

        if AverageOverSessions
          X = mean(cat(3,beta{:}),3)';
        else
          X = cell2mat(beta(:));
        end
        clear beta;

        COORDS = struct('orientation','mni','labels',ELECTRODE{iSubj},'ijk',[],'ind',[],'xyz',XYZ{iSubj});

        % ---------
        metadata(iSubj).AverageOverSessions = AverageOverSessions;
        metadata(iSubj).BoxCarSize = BoxCarSize;
        metadata(iSubj).WindowStartInMilliseconds = WindowStartInMilliseconds;
        metadata(iSubj).WindowSizeInMilliseconds = WindowSizeInMilliseconds;
        metadata(iSubj).filters = FILTERS;
        metadata(iSubj).coords = COORDS;
        if AverageOverSessions == 1;
            for iTarget = 1:numel(TARGETS)
                switch lower(TARGETS(iTarget).type)
                case 'category'
                    TARGETS(iTarget).target = repmat(TARGETS(iTarget).target,nsessions,1);
                case 'similarity'
                    TARGETS(iTarget).target = repmat(TARGETS(iTarget).target,nsessions,nsessions);
                end
                SCHEMES = repmat(SCHEMES,nsessions,1);
            end
            metadata(iSubj).sessions = [];
            metadata(iSubj).nrow = 100;
        else
            metadata(iSubj).sessions = stim_order{1};
            metadata(iSubj).nrow = 400;
        end
        metadata(iSubj).targets = TARGETS;
        metadata(iSubj).cvind = SCHEMES;
        metadata(iSubj).ncol = 0; % will be set later
    end

    %% Load And Process Data and Coordinates
    fmt = 's%02d.mat';
    mode = 'raw';
    for iSubj=1:NSUBJ
        sdir = sprintf('Pt%02d',iSubj);
        sfile = filelist{iSubj,1};
        if isempty(sfile)||isempty(filelist{iSubj,4});
            fprintf('Skipping subject %d because of missing data.\n',iSubj);
            continue;
          else
            fprintf('Beginning subject %d.\n',iSubj);
        end
        dpath_out = fullfile(DATA_DIR_OUT, sprintf(fmt,iSubj));
        if exist(dpath_out,'file') && ~OVERWRITE;
            fprintf('Skipping subject %d because output already exists.\n',iSubj)
            continue
        end
        spath = fullfile(DATA_DIR,sdir,sfile);
        fprintf('Loading %s...\n', spath);
        Pt = load(spath);

        nChunks = numel(filelist{iSubj,2});
        if nChunks > 1
          nTicks = 0;
          for iChunk = 1:nChunks
            cvar = filelist{iSubj,2}{iChunk};
            nTicks = nTicks + size(Pt.(cvar).DATA,1);
          end
          interval = Pt.(cvar).DIM(1).interval;
          electrodeLabels = Pt.(cvar).DIM(2).label;
          Pt.LFP(1) = init_source_struct(nTicks,electrodeLabels,interval);
          for iChunk = 1:nChunks
            cvar = filelist{iSubj,2}{iChunk};
            sessions = filelist{iSubj,3}{iChunk};
            nSessions = numel(sessions);
            tagfmt = filelist{iSubj,4};
            if iChunk == 1
              a = 1;
              b = size(Pt.(cvar).DATA,1);
              Pt.LFP.DATA(a:b,:) = Pt.(cvar).DATA;
              Pt.LFP.DIM(1).scale(a:b) = Pt.(cvar).DIM(1).scale;
              psize = b;
              pscale = max(Pt.(cvar).DIM(1).scale);
              Pt = rmfield(Pt,cvar);
            else
              a = psize + 1;
              b = psize + size(Pt.(cvar).DATA,1);
              Pt.LFP.DATA(a:b,:) = Pt.(cvar).DATA;
              Pt.LFP.DIM(1).scale(a:b) = Pt.(cvar).DIM(1).scale + pscale;

              for iSession = 1:nSessions
                tag = sprintf(tagfmt,iSession);
                Pt.(tag) = Pt.(tag) + psize;
              end

              psize = b;
              pscale = max(Pt.(cvar).DIM(1).scale);
              Pt = rmfield(Pt,cvar);
            end
          end
        else
          cvar = filelist{iSubj,2}{1};
          Pt.LFP = Pt.(cvar);
          Pt = rmfield(Pt,cvar);
        end

        ecoord = ELECTRODE{iSubj};
        edata = cellstr(Pt.LFP.DIM(2).label);

        zd = ismember(edata, ecoord);
        zc = ismember(ecoord, edata);

        Pt.LFP.DATA = Pt.LFP.DATA(:,zd);

        onsetIndex = cell(1,4);
        tagfmt = filelist{iSubj,4};
        for iSession = 1:4
          tagname = sprintf(tagfmt,iSession);
          onsetIndex{iSession} = Pt.(tagname);
        end

        Hz = 1 / Pt.LFP.DIM(1).interval; % ticks per second
        window_start = (WindowStartInMilliseconds / 1000) * Hz;
        window_size = (WindowSizeInMilliseconds / 1000) * Hz; % in ticks (where a tick is a single time-step).

        % Will return a session -by- electrode cell array, each containing a
        % trial -by- time matrix.
        M = arrangeElectrodeData(Pt.LFP.DATA, onsetIndex, [window_start, window_size]);

        % Sort and average time-points
        nElectrodes = size(M,2);
        for iElectrode = 1:nElectrodes
            for iSession = 1:4
                M{iSession,iElectrode} = M{iSession,iElectrode}(stim_sort_ix{iSession},:);
                if BoxCarSize > 1
                    M{iSession,iElectrode} = boxcarmean(M{iSession,iElectrode},BoxCarSize,'KeepPartial',0);
                end
            end
            % Average Sessions
            if AverageOverSessions
                tmp = cat(3,M{:,iElectrode});
                M{1,iElectrode} = mean(tmp,3);
            end
        end
        if AverageOverSessions
            M(2:end,:) = [];
        end
        X = cell2mat(M);
        [~,reduxFilter] = removeOutliers(X);
        metadata(iSubj).filters(end+1) = struct('label','rowfilter','dimension',1,'filter',reduxFilter.words);
        metadata(iSubj).filters(end+1) = struct('label','colfilter','dimension',2,'filter',reduxFilter.voxels);
        metadata(iSubj).ncol = size(X,2);
        metadata(iSubj).samplingrate = Hz;

        save(dpath_out, 'X');
    end

    %% Save metadata
    save(fullfile(DATA_DIR_OUT,'metadata.mat'));
end
