% Initialize an empty matrix for the concatenated data
dN = [];

% Specify the directory where the chunk files are stored
chunk_dir = './chunks/'; % Adjust this to the correct directory if needed

% Loop through all chunk files
file_index = 1;
while true
    % Construct the filename
    chunk_file = sprintf('%sdN_chunk_%d.mat', chunk_dir, file_index);
    
    % Check if the file exists
    if isfile(chunk_file)
        % Load the chunk
        loaded_data = load(chunk_file, 'dN_to_save');
        
        % Append the chunk to the final matrix
        dN = [dN, loaded_data.dN_to_save]; %#ok<AGROW>
        
        % Increment the file index
        file_index = file_index + 1;
    else
        % Exit the loop if no more files are found
        break;
    end
end

% Save the concatenated matrix to a single file
save('simdN.mat', 'dN', '-v7.3');
disp('Concatenation complete. Final matrix saved as simdN.mat.');
