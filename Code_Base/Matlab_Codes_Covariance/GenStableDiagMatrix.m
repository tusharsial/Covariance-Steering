% Function for generating an n x n random stable diagonal matrix with 
% eigenvalues having real parts less than 0.

function X = GenStableDiagMatrix(w1, w2, nx)
    
% Preallocate the diagonal vector
eigenvalues = zeros(nx, 1);

% Generate random eigenvalues
for i = 1:nx
     % Generate a random real part less than 0
     real_part = -(w1 + (w2-w1)*rand()); % Random negative real part

     % Generate a random imaginary part
     %imag_part = rand() * 2 - 1; % Random value in [-1, 1]
     imag_part = 0;

     % Combine real and imaginary parts
     eigenvalues(i) = real_part + 1i * imag_part;
 end

 % Create the diagonal matrix
 X = diag(eigenvalues);
end