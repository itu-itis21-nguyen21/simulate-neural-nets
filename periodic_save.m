function [dN, n, file_index] = periodic_save(dN, n, file_index, M)
    % Periodic save operation
    % Determine the range to save
    save_range = 1:n-M-1;
    

    % Save dN to file excluding the last M columns
    save_file = sprintf('./chunks/dN_chunk_%d.mat', file_index);
    dN_to_save = dN(:, save_range);
    save(save_file, 'dN_to_save', '-v7.3');
    disp(['Saved ', save_file]);

    % Retain the last M columns for history
    dN(:, 1:M+1) = dN(:, n-M-1:n-1);

    % Reset dN for the next chunk
    dN(:, M+2:end) = 0;

    % Reset n
    n = M + 1;

    % Increment file index
    file_index = file_index + 1;
end