function Y = generate_stable_diagonal_matrix(a, b, n)
    % Generates an n x n diagonal matrix with eigenvalues
    % having real parts less than 0.

    % Preallocate the diagonal vector
    eigenvalues = zeros(n, 1);

    % Generate random eigenvalues
    for i = 1:n
        % Generate a random real part less than 0
        real_part = -(a + (b-a)*rand()); % Random negative real part

        % Generate a random imaginary part
        %imag_part = rand() * 2 - 1; % Random value in [-1, 1]
        imag_part = 0;

        % Combine real and imaginary parts
        eigenvalues(i) = real_part + 1i * imag_part;
    end

    % Create the diagonal matrix
    Y = diag(eigenvalues);
end