%NDIM_PROD Display of polydimensional matrices.
% str = ndim_disp( A, A_ind, is_disp)
%
% Example 
%   A = rand(2,3,4); % make 3 dim matrix
%   ndim_disp( A, 'a b c');% display matrix
%   str = ndim_disp( A, 'a b c',0);% get display string without view
%
%   See also ndim_prod, ndim_prod_123.
%
%   (c) Leonid Lazarev 2025

% history, mmat_disp 2016-01-09

function str = ndim_disp( A, A_ind, is_disp_str, is_disp_A )

max_Anumel = 40;
in1 = inputname(1);
if nargin<2 
    A_ind = [];in2 = 'ind';
else 
    in2 = inputname(2); 
    if isempty(in2), in2 = 'ind';end
end
if nargin<3, is_disp_str = 1;end
if nargin<4, if  numel(A)<max_Anumel, is_disp_A =1; else, is_disp_A = 0; end;end

str = ['size(',in1,') = [',num2str(size(A)),']'];
if ~isempty(A_ind),  str = [str,', ',in2,' = ''',A_ind,''''];
end
if is_disp_str, disp(str); end
if is_disp_A, disp(A);end




