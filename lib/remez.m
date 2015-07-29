function [A, b] = remez(fun, fun_der,interval,order,varargin)
% REMEZ: An algorithm to calculate the minimax polynomial to a given function over an interval.
%
% Inputs:
%   -fun: A string which computes the function that we want to model.
%   -fun_der: A string which computes the derivative of fun.
%   -interval: The start and end points of the interval over which we would like to model fun.
%   -order: Order of the polynomial that we want to find.
% Optional Inputs (given in 'varible', value pairs):
%   -maxIter:
%   -thresh:
%   -findZeroThresh:
%   -findZeroMaxIter:
% Outputs:
%   -A: Coefficients for the polynomial, p(x-interval(1)), and the maximum error A(end).
%       Coefficients are in order of increasing powers of (x-interval(1)).
%   -b: Coefficients for the polynomial, p(x).
%
% Originally written by Sherif A. Tawfik, Faculty of Engineering, Cairo University.
% Modified by Nicholas J. Fraser, University of Sydney.

% Specify some default settings.
maxIter = 10;
thresh = 2^-30;
findZeroThresh = 2^-52;
findZeroMaxIter = 1000;

% Create the input parser.
p = inputParser;
% Add values to parse.
addOptional(p, 'maxIter', maxIter, @isnumeric);
addOptional(p, 'thresh', thresh, @isnumeric);
addOptional(p, 'findZeroThresh', findZeroThresh, @isnumeric);
addOptional(p, 'findZeroMaxIter', findZeroMaxIter, @isnumeric);
parse(p, varargin{:});

% Fetch parsed values.
maxIter = p.Results.maxIter;
thresh = p.Results.thresh;
findZeroThresh = p.Results.findZeroThresh;
findZeroMaxIter = p.Results.findZeroMaxIter;

powers=ones(order+2,1)*([0:order]);% the powers of the polynomial repeated in rows (order +2) times
coeff_E =(-1).^[1:order+2];  
coeff_E=coeff_E(:);     % the coefficients of the E as a column array
t=1:order;  
t=t(:); % the powers of the polynomial starting from 1 in a column. This is used when differntiation
%the polynomial

y=linspace(interval(1),interval(2),order+2); % the first choice of the (order+2) points

for i=1:maxIter,
    y=y(:); % make the points array a column array
    h=(y-interval(1))*ones(1,order+1); % repeat the points column minus the start of the interval 
                                       %(order +1) times 
    coeff_h=h.^powers;   % raise the h matrix by the power matrix elementwise
    M=[coeff_h coeff_E]; % the matrix of the LHS of the linear system of equations   
    N= feval(fun,y); % the column vector of the RHS of the linear system of equations   
    A=M\N;  % solution of the linear system of equations, first (order +1) element are the    
            % coefficients of the polynomial. Last element is the value of
            % the error at these points 
    A1=A(1:end-1); % the coefficients only         
    A_der=A(2:end-1).*t;   % the coeffcients of the derivative of the polynomial                   
    
    z(1)=interval(1);  % z(1) is the start point of the interval
    z(order+3)=interval(2);   % z(order+3) is the end point of the interval
    % in between we fill in with the roots of the error function
    for k=1: order+1        
            z(k+1)=findzero(@err, y(k), y(k+1), findZeroMaxIter, findZeroThresh, fun, A1, interval(1));
    end

    % between every two points in the array z, we seek the point that
    % maximizes the magnitude of the error function. If there is an extreme
    % point (local maximum or local minimum) between such two points of the
    % z array then the derivative of the error function is zero at this
    % extreme point. We thus look for the extreme point by looking for the
    % root of the derivative of the error function between these two
    % points. If the extreme point doesn't exist then we check the value of the  error function
    % at the two current points of z and pick the one that gives maximum
    % magnitude
    
    for k=1:order+2,
        if sign(err(z(k),fun_der,A_der,interval(1) ))~=sign(err(z(k+1),fun_der,A_der,interval(1))) % check for a change in sign
            y1(k)=findzero(@err, z(k), z(k+1), findZeroMaxIter, findZeroThresh, fun_der, A_der, interval(1)); % the extreme point that we seek
            v(k)=abs(err(y1(k),fun,A1,interval(1))); % the value of the error function at the extreme point
        else  % if there is no change in sign therefore there is no extreme point and we compare the endpoints of the sub-interval
            v1=abs(err(z(k),fun,A1,interval(1))); % magnitude of the error function at the start of the sub-interval
            v2=abs(err(z(k+1),fun,A1,interval(1))); % magnitude of the error function at the end of the sub-interval
            % pick the larger of the two
            if v1>v2
                y1(k)=z(k);
                v(k)=v1;
            else
                y1(k)=z(k+1);
                v(k)=v2;
            end
        end            
    end
    [mx ind]=max(v); % search for the point in the extreme points array that gives maximum magnitude for the error function
    % if the difference between this point and the corressponding point in
    % the old array is less than a certain threshold then quit the loop
    if abs(y(ind)-y1(ind)) < thresh, break; end
    % compare it also with the following point if it is not the last point
    if ind<length(y) & abs(y(ind+1)-y1(ind))  < thresh, break; end
    % replace the old points with the new points
    y=y1;
end

if i == maxIter, warning('Remez did not converge after %d iterations.', i); end

% Convert coefficients to ones not biased by x0.
% TODO: Make this its own function?
% TODO: Calculate this using recursion?
c = A(1:end-1);
b = zeros(size(c));
previousRow = zeros(size(c));
x0 = interval(1);
for i=1:length(c),
    nextRow = zeros(size(c));
    for j=i:-1:1,
        d = i-j;
        if j == 1,
            nextRow(j) = 1;
        else
            nextRow(j) = previousRow(j) + previousRow(j-1);
        end
        b(j) = b(j) + c(i)*nextRow(j)*(-x0)^d;
    end
    previousRow = nextRow;
end

