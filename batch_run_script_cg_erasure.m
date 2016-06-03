% this is the script to batch run cg_erasure on the matrices.

%mat_folder = '/u/antor/u50/zhu36/ResearchProjects/ErasureCodedComputation/data';
mat_folder = '.';
addpath(mat_folder);

% This set of matrices are nice cases selected by me.
matrix_file_list = {'Ltridiag500', 'mhdb416', 'nos3'};

% The number of faults.
num_fail_list = [0, 1, 0.2];  % 0.2 means 20% percentage faults.
%%
rand('seed', 0);
randn('seed', 0);

tol = 1e-10;
for i = 1:numel(matrix_file_list)
    matrix_file = matrix_file_list{i};
    for j = 1:length(num_fail_list)
        num_fail = num_fail_list(j);
        log_file = [matrix_file '_num_fail=' num2str(num_fail) '.out'];
        result_file = [matrix_file '_num_fail=' num2str(num_fail) '_results.mat'];
        diary(log_file); run_script_cg_erasure; diary off;
        % save the results.
        save(result_file);
    end
end

rmpath(mat_folder);

%% make the figures
tol = 1e-10;
for i = 1:numel(matrix_file_list)
    matrix_file = matrix_file_list{i};
    for j = 1:length(num_fail_list)
        num_fail = num_fail_list(j);
        log_file = [matrix_file '_num_fail=' num2str(num_fail) '.out'];
        result_file = [matrix_file '_num_fail=' num2str(num_fail) '_results.mat'];
        load(result_file)
        
        residual_figure(matrix_file, resvec, fail_point, num_fail);
    end
end