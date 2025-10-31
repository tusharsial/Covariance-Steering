function createDirectoryStructure(plot_dir, main_dir, di_dir, cw_dir, di_subdirs, cw_subdirs)
    
% Create Plot directory
    if ~exist(plot_dir, 'dir')
        mkdir(plot_dir);
        fprintf('Created directory: %s\n', plot_dir);
    end
    
% Create main directory
    if ~exist(main_dir, 'dir')
        mkdir(main_dir);
        fprintf('Created directory: %s\n', main_dir);
    end
    
% Create DI directory and subdirectories
    if ~exist(di_dir, 'dir')
        mkdir(di_dir);
        fprintf('Created directory: %s\n', di_dir);
    end
    for i = 1:length(di_subdirs)
        subdir = fullfile(di_dir, di_subdirs{i});
        if ~exist(subdir, 'dir')
            mkdir(subdir);
            fprintf('Created directory: %s\n', subdir);
        end
    end
    
% Create Noisy CW directory and subdirectories
    if ~exist(cw_dir, 'dir')
        mkdir(cw_dir);
        fprintf('Created directory: %s\n', cw_dir);
    end
    for i = 1:length(cw_subdirs)
        subdir = fullfile(cw_dir, cw_subdirs{i});
        if ~exist(subdir, 'dir')
            mkdir(subdir);
            fprintf('Created directory: %s\n', subdir);
        end
    end
    
    fprintf('Directory structure created successfully!\n');
end