% DATA 303 ASSIGNMENT - Danielle Hall - 15/10/25
%--------------------------------------------------------------------------

%========================= QUESTION 1 =====================================
clc, clearvars % clearing previous stuff

% Initialising my matricies
A = [1 2; 0 1; 1 4; 1 2; 1 0];
b = [4 1 5 6 1]'; % transposed
agmt = [A b]; % add b as third column of A

% Householder Reflections
% reflection 1
[newA,H] = house303(agmt,1,1); % row 1 column 1
% reflection 2
[newA2,H2] = house303(newA,2,2); % row 2 column 2

% extract R and Q^Tb
R = newA2(1:2,1:2)
Qtb = newA2(1:2,3)

% solves for x1 x2
x = R \ Qtb

% the residual vetor is the 'leftover' at the bottom of the third column
residual = newA2(3:5,3)
length_residual = norm(residual)
%========================= QUESTION 1 END =================================

%========================= QUESTION 2 =====================================

%--------------------------------------------------------------------------
function [newA,H] = house303(A,ii,kk,ForceExact)
% This subroutine does one Householder reflection on the matrix A
% so that all elements in column k which lie below row i are set to
% zero.    The routine outputs the reflected matrix newA and the
% Householder matrix H used to perform that reflection.
%
% USAGE:  [newA,H] = house303(A,i,k)
%
format compact;
if ~exist('ForceExact'),  ForceExact = 1;  end
[m,n] = size(A);             % get the size (m by n) of the matrix A
Acolumn = zeros(m,1);        % initialize the lower part of column k
Acolumn(ii:m) = A(ii:m,kk);  % copy column k from row i to last row
lencol = norm(Acolumn,2);    % get the length of Acolumn
w = Acolumn;                 % set the Householder vector to Acolumn

% now add +/- length of Acolumn to w in the ii^th element.
if Acolumn(ii) > 0
    w(ii) = w(ii) + lencol;    
else
    w(ii) = w(ii) - lencol;
end  

len_w = norm(w,2);          
if len_w == 0                % avoid divide by zero (do nothing and exit)
    newA = A;  H = eye(m);  return;  
end
w = w/len_w;                  % normalize w
H = eye(m) - (2*w)*w';
newA = H * A;
if ForceExact,   newA(ii+1:m,kk) = zeros(m-ii,1);  end
return;

end
%--------------------------------------------------------------------------
