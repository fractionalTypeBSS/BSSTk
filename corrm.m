function M = corrm(x,tau)

if nargin==0, error('Not enough input arguments.'); end
if nargin>2, error('Too many input arguments.'); end
if ndims(x)>2, error('Inputs must be 2-D.'); end

if min(size(x)) == 1 % Handle special case
    x = x(:);  
end

[m,n] = size(x);

if nargin < 2 
  tau=0;
end

tau=fix(abs(tau));
if tau>m 
  error('Choose tau smaller than vector size');
end

xc = x - repmat(sum(x)/m,m,1);  % Remove mean

R=x(1:m-tau,:);   % make time delayed vector
L=x(1+tau:m,:);

M=L'*R / (m-tau); % calculate correlation
