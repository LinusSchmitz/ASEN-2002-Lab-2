% Least Squares Function
function [VelAvg,VelErr,A,B] = LSFunc(x,y)
    
    %Find number of data points in the vectors
    N = length(x);
    
    % Find linear best fit coefficients A and B
    % Create H matrix
    H = [ones(N,1) x];

%     del = N*sum(x.^2)-(sum(x))^2;
%     A = (sum(x.^2)+sum(y.^2)-sum(x)*sum(x.*y))/del;
%     B = (N*sum(x.*y)-sum(x)*sum(y))/del;
%     
%     SigmaY = sqrt(sum((y-A-B.*x).^2))/sqrt(N-2);

    SigmaY = 0.1;
    
    % Create W matrix (hint: type <help diag> in command line)
    W = ones(N,1); % set W as identity matrix
    W = (1/SigmaY^2)*W;% recompute the W matrix using SigmaY
    W = diag(W);
    
    % Solve for P matrix
    p1 = H' * W;
    p2 = p1 * H;
    P = inv(p2);

    % Solve for x_hat matrix and extract A and B parameters
    x_hat = P * H' * W * y;
    A = x_hat(1);
    B = x_hat(2);

    % extract uncertainty in A and uncertainty in B from P matrix
   % A_error = sqrt(P(1));
    %B_error = sqrt(P(4));

    VelAvg = A + B * x; % best fit line

    VelErr = H * P * H'; % predicted error
    VelErr = diag(VelErr);
    
end