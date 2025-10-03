function vq = sincInterp(x,v,xq,w)
%%SINCINTERP    Ideal band limited interpolation
%
%   Vq = SINCINTERP(X, V, Xq, W) interpolates to find Vq, the values of the
%   underlying function V=F(X) at the query points Xq. 
%  
%   X must be a vector. The length of X is equal to N.
%   V must have length N, and Vq is the same size as Xq.
%   W is the bandlimit in Hz

%   2022 NIST/CTL Steve Blandino

%   This file is available under the terms of the NIST License.

%% Input processing
if isvector(v)
    v = v(:); % Ensure v is a column vector
end

if size(v, 1) ~= numel(x)
    error('The number of rows in v must match the length of x.');
end
%% Interpolation

[x, id] = sort(x, 'ascend');
v = v(id,:);
[Ts,T] = ndgrid(xq,x);

vq= sinc((Ts - T)*w)*v;

end