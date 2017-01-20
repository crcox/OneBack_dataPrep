function setup_data(varargin)
    %% Setup
%    OVERWRITE = false;
    DATA_DIR = '/media/chris/Data/MRI/PictureOneBack/data';
    SRC_DIR = '/home/chris/src/OneBack_dataPrep';
    SUBJ_DIR = arrayfun(@(i) sprintf('MRH026_4%02d', i), 1:22, 'unif', 0);
    BASELINE = 'active';
    STIM_DIR = fullfile(SRC_DIR,'stimuli');
    COORD_DIR = fullfile(SRC_DIR,'coords');
    SIM_DIR = fullfile(SRC_DIR,'similarity');
    SUBJECTS = 1:22; % subjects 11 and 12 do not have labeled coordinates.
    nSubj = numel(SUBJECTS);
    AverageOverSessions = true;

    if nargin > 1
      args = reshape(varargin, 2, [])';
      for iArg = 1:size(args,1)
        key = args{iArg,1};
        value = args{iArg,2};
        switch lower(key)
        case 'average'
          AverageOverSessions = value;
        case 'baseline'
          BASELINE = value;
        end
      end
    end
    
    switch BASELINE
      case 'active'
        TXT_DIR = cellfun(@(sdir) fullfile(DATA_DIR,sdir,'analysis_RSA_corr'), SUBJ_DIR, 'unif', 0)';
      case 'passive'
        TXT_DIR = cellfun(@(sdir) fullfile(DATA_DIR,sdir,'analysis_RSA'), SUBJ_DIR, 'unif', 0)';
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
    file_stim_key = fullfile(STIM_DIR, 'stimlist.csv');

    [ STIMULI, TARGETS ] = import_stimulus_txt(file_stim_key);
    TARGETS = rmfield(permute(TARGETS(:,1,:),[1,3,2]), {'subject','session'});

    nItems = numel(unique(STIMULI.word));
    nSessions = numel(unique(STIMULI.session));

    %% Generate sort indexes
    % We do not need them

    %% load similarity structure
    % NEXT (judge similarity in kind based on word)
    % The rows in the embedding are already in an order that matches the order
    % in the stim_key.
    embedding = csvread(fullfile(SIM_DIR,'NEXT_US_KIND_IMAGE_3D.csv'));
    if usejava('jvm')
      plot_similarity_decompositions(embedding);
    end
    S = corr(embedding','type','Pearson');
    t = struct(...
        'label','semantic',...
        'type','similarity',...
        'sim_source','NEXT',...
        'sim_metric','correlation',...
        'target',S...
    );
    TARGETS = [TARGETS,repmat(t,nSubj,1)];
    
    t = struct(...
        'label','semantic',...
        'type','embedding',...
        'sim_source','NEXT',...
        'sim_metric','euclidean',...
        'target',embedding...
    );
    TARGETS = [TARGETS,repmat(t,nSubj,1)];

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
        'stimuli',{unique(STIMULI.word)},...
        'filters',struct(),...
        'coords',struct(),...
        'cvind',[],...
        'nrow',[],...
        'ncol',[],...
        'sessions',[],...
        'AverageOverSessions',[]...
    );

    %% Setup CV indexes
    nschemes = 10;
    nfolds = 9;
    SCHEMES = zeros(nItems, nschemes);
    z1 = strcmpi({TARGETS(1,:).label},'multinomial');
    multinomial = TARGETS(1,z1).target;
    for iScheme = 1:nschemes
        c = cvpartition(multinomial,'KFold', nfolds);
        for iFold = 1:nfolds
            SCHEMES(:,iScheme) = SCHEMES(:,iScheme) + (test(c, iFold) * iFold);
        end
    end

    %% Define Metadata and Process Data
    for iSubj = 1:nSubj
        disp(iSubj)
        % DATA
        dpath_in = fullfile(TXT_DIR{iSubj}, 'all+tlrc.txt');
        dpath_out = fullfile(DATA_DIR_OUT,sprintf('MRH026_4%02d_%s_%s.mat',iSubj,BASELINE,bdir));
        mpath_out = fullfile(DATA_DIR_OUT,sprintf('metadata_%s_%s.mat',BASELINE,bdir));

        [beta,xyz,ijk] = import_beta_txt(dpath_in, nSessions, nItems);

        % DROP ALL-ZERO VOXELS
        z = cell(size(beta));
        for ib = 1:numel(beta);
          z{ib} = any(beta{ib}, 2);
        end
        z1 = all(cell2mat(z),2);
        for ib = 1:numel(beta)
          beta{ib} = beta{ib}(z1,:);
        end
        
        % COMPOSE MATRIX (average or concatenate)
        if AverageOverSessions == 1
          X = mean(cat(3,beta{:}),3)';
        else
          X = cell2mat(beta(:)')';
        end 

        % FILTERS
        FILTERS = generate_outlier_filters(X);
        FILTERS(3) = struct('label','animate','dimension',1,'filter',TARGETS(iSubj,1).target);
        FILTERS(4) = struct('label','notanimate','dimension',1,'filter',~TARGETS(iSubj,1).target);

        COORDS = struct('orientation','mni','ijk',ijk,'ind',[],'xyz',xyz);

        % ---------
        metadata(iSubj).AverageOverSessions = AverageOverSessions;
        metadata(iSubj).coords = COORDS;
        cvind = SCHEMES;
        if AverageOverSessions == 0;
            for iTarget = 1:numel(TARGETS(iSubj,:))
                switch lower(TARGETS(iTarget).type)
                case 'category'
                    TARGETS(iSubj,iTarget).target = repmat(TARGETS(iSubj,iTarget).target,nSessions,1);
                case 'similarity'
                    TARGETS(iSubj,iTarget).target = repmat(TARGETS(iSubj,iTarget).target,nSessions,nSessions);
                end
            end
            for iFilter = 1:numel(FILTERS)
              if FILTERS(iFilter).dimension == 1 && ~strcmp(FILTERS(iFilter).label, 'rowfilter')
                FILTERS(iFilter).filter = repmat(FILTERS(iFilter).filter(:),nSessions,1);
              end
            end
            cvind = repmat(cvind,nSessions,1);
            z = STIMULI.subject == iSubj+400;
            metadata(iSubj).sessions = STIMULI.session(z);
        end
        metadata(iSubj).filters = FILTERS;
        metadata(iSubj).targets = TARGETS(iSubj,:);
        metadata(iSubj).cvind = cvind;
        metadata(iSubj).nrow = size(X,1);
        metadata(iSubj).ncol = size(X,2);
        %% Save data
        save(dpath_out, 'X');
    end

    %% Save metadata
    save(mpath_out,'metadata');
end
