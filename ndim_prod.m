%NDIM_PROD Product of polydimensional matrices with named indexes. 
%
%   [C,C_indexes] = NDIM_PROD(A,A_ind,B,B_ind,C_ind) - is the product of matrices A and B with named indexes. 
%
%     The index names denotes like this 'a b c' or 'ab1,bcd,cdef35,678' with space or comma delimiters.
%       They must be unique ('n m n m' is wrong!). Number of indexes may be
%       any, from zero ('','a','a b','a b c ',...). It may be equal or more then number
%       dimensions of marixes, without leading and trailing singular dimensions.  
%       All coinciding names of indexes in A_indexes, B_indexes treated as common
%       indexes (for 'a b c', 'b c d e' common 'b', 'c').
%       Size of matrices A,B must coincide on all common indexes! Not coinciding
%       indexes treated as individual indexes (for 'a b c', 'b c d e' individual 'a', 'd', 'e').
%       Common indexes may be summarizing or static.
%       Summarizing indexes must be eclosed in bracket individually like this 'a (b) (c)'
%       or 'a,(b),(c)'. All other common indexes are static
%       (for 'a (b) c', 'b c d e' remaining 'c').
%       Default order of result indexes in matrix C is
%       [remaining(A,B),individual(A),individual(B)]  (for 'a (b) c', 'b c d e' C_indexes = 'с a d e').
%       In order to permute result matrix or simply to insist their correctness default
%       parameter C_indexes may be included.
%       Braced index in B_ind are ignored.
%
%   [C,C_indexes] = NDIM_PROD(A,A_ind,B,B_ind,D,D_ind,C_ind) - is the product Ax(BxD). 
%   [C,C_indexes] = NDIM_PROD(A,A_ind,B,B_ind,D,D_ind,E,E_ind,C_ind) - is the product Ax(Bx(DxE)). 
%   [C,C_indexes] = NDIM_PROD(A,A_ind,B,B_ind,D,D_ind,E,E_ind,...,...,C_ind) - product of any number of matrixes. 
%
%       Production of many matrixes may be calculated at once. Produnction
%       is fulfilled form right lo left. Any time summarizing indexes takes from
%       index string of left matrix. 
%       Any matrix may have singularities, i.e. dimensions with 1 size in any positions, f.e. rand(1,2,3,1,1,1,3).
%       Scalars, columns and rows also may be used as a matrix with any
%       number of indexses. 
%
% RUS
%       Произведение матриц с именованными индексами. Индексы перечисляются в одной строке через пробел
%       или запятую 'a b c' или 'a,b,c'. Не допускается повторения имен индексов у одной матрицы
%       ('n m n m' - так нельзя!). Число индесов в строке может быть
%       любое, начиная с нуля ('','a','a b','a b c d e f'). Оно должно быть
%       равно или больше мерности матрицы, не считая ведущих сингулярных
%       размерностей. Имена могут состоять  из любых символов кроме скобок. 
%       Все совпадающие по именам индексы у A,B считаются общими (для 'a b c','b c d e' общие 'b','c').
%       Размер матриц A,B вдоль общих индексов должен совпадать! Невоспадающие по именам индексы
%       считаются индивидуальными (for 'a b c','b c d e' индивидуальные 'a','d','e').
%       Общие индексы могут быть суммируюшими или стационарными.
%       Суммирующие индексы нужно заключить в овальные скобки 'a (b) (c)'
%       или 'a,(b),(c)'. Все остальные общие индексы считаются стационарными
%       (для 'a (b) c','b c d e' сьтационарный 'c').
%       Рекомендуется всегда последним аргументом ставить порядок индексов в результирующей матрице. 
%       По умолчанию порядок следования индексов в результирующей матрице
%       такой [стационарные(A,B),индивидыальные(A),индивидуальные(B)]  (для 'a (b) c','b c d e' C_indexes = 'с a d e').
%       Скобки в индексах последней матрицы (B_ind) не принимаются во
%       внимание.
%       В матрицах в любой позиции могут быть сингулярные индексы, т.е.
%       единичной размерности. Поэтому, например,строку rand(1,3) можно
%       описать размыми сопосбами, например, 'a' - как столбец, 'a1 b' - как строку,
%       'a1 b c1' - как 3D матрицу. %
%
%       Можно за раз перемножать три и более матриц. Перемножаться они
%       будут справа на лево, исползуя обозначения суммирующих индексов в
%       левом из сомножителей.
%
%  Examples:
%       A = rand(2,3,4); % make 3 dim matrix
%       B = rand(3,4,5); % make 3 dim matrix
%       [C1,C1in] = ndim_prod(A,'a b c',B,'b c d','a b c d');     size(C1),C1in% direct product without summation
%       [C2,C2in] = ndim_prod(A,'a (b) c',B,'b c d','c a d');    size(C2),C2in% product with summarizing on b
%       [C3,C3in] = ndim_prod(A,'a (b) (c)',B,'b c d','a d');       size(C3),C3in% product with summarizing on b,c
%       [C4,C4in] = ndim_prod(A,'a b c',B,'d e f','a b c d e f');   size(C4),C4in% no com, no sum 
%
%   Examples 2:
%       A = rand(2,3); % make 2 dim matrix
%       B = rand(3,4,5); % make 3 dim matrix
%       D = rand(4,5,6,2); % make 4 dim matrix
%       E = rand(5,6,2,2,2); % make 5 dim matrix
%       [C,C_ind] = ndim_prod(A,'a (b)',B,'b (c) (d)',D,'c d (e) f',E,'d e f j h','a f j h') ;
%       size(C)
%
%   See also ndim_disp ndim_prod_123.
%
% (c) Leonid Lazarev 2025

%   [C,C_indexes] = NDIM_PROD(A,A_ind,C_ind)
%       summiarize matrix A and/or permute
% RUS   сумма матрицы и/или изменение порядка индексов
%
%   Example 3:
%       A = rand(4,5,6,2); % make 4 dim matrix
%       [C,C_ind] = NDIM_PROD(A,'a b c d','d c b a'); % equal to permute(A,[4 3 2 1])
%       [C,C_ind] = NDIM_PROD(A,'a (b) c (d)','c a');% equal to permute( squeeze( sum(A,[2,4] ))
%       [C,C_ind] = NDIM_PROD(A,'a (b) c (d)');% equal to squeeze(sum(A,[2,4] ))
%   set global param equal 0
%       global is_check_ndim_prod; is_check_ndim_prod=0;
%   for speed-up after checking code
%
% можно вызывать все внутренние функции из командной строки или другой функции
% [X,Y] = ndim_prod( FUN_NAME, PAR1, PAR2, ... )

%   See also (do).
%
% (c) Leonid Lazarev 2025

% remake of 2012-12-28 mmat_xx

function [C,indC] = ndim_prod(varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEMO +
if nargin==0 % CHECK
   
    strAA = {[],'1','2 1','1 3 2','4 2 3 1','3 2 5 4 1'};
     strAAx = {[],'(1)','(2) (1)','(1) (3) (2)','(4) (2) (3) (1)','3 2 5 4 1'};
   j1 = 3;
    for j  = j1:6
        for jj = (j+1):6
%             [j jj ] = deal(3    , 4);
            disp([ j jj ])
            strA = strAA{j};
            strAx = strAAx{j};
            strB = strAA{jj};
            if ~isempty(strA), Ai = str2num(strA);else Ai = [];end
              if ~isempty(strB), Bi = str2num(strB);else Bi = [];end
            if ~isempty(Ai) A = rand(Ai);else A=5;end
            if ~isempty(Bi) B = rand(Bi);else B=4;end
            
            [C1,C1in] = ndim_prod(A,strAx,B,strB);
            [indA,indAX] = str2ind(strAx);
             [indB] = str2ind(strB);
             
            C2 = mmat_xx(A,[strA,' X ',ind2str(indAX)],B,strB,C1in);
            
            max(C2(:)-C1(:))
%             ndim_ind2str_for_mmat_xx
%             size(C1)
%             C1in
        end
    end
    return; %------->
end
% SUB FUNCTION +
if ischar(varargin{1})
    varargout = feval(varargin{1},varargin{2:end});
    if length( varargout )>0, C = varargout{1};end
    if length( varargout )>1, Cind = varargout{2};end
    return; %------->
end
% NDIM_PROD +
%%%%%%%%%%%%%%%%%%%%%%%%%%
lv = length(varargin);
if lv<4, error('ndim_prod: No enought param!');end
% если нет выходного списка индексов
if mod(lv,2)==0% нет indC % вариант без indC нужен в рекурсии, нельзя делать обязательным
    indC = [];
    mat_num = lv/2;
    par = varargin; % no indC
    % если он есть
else % есть indC
    indC = varargin{end};% output indexes
    mat_num = (lv-1)/2; % num of matrix
    par = varargin(1:end-1);% indata without indC
end
if mat_num==2 % если 2 матрицы м=то умножаем
    [C, indC] = As_x_Bs(par{:},indC); %  ndim_prod(A,indA,B,indB) - product
else % >2
    % рекурсия без первой матрицы
    [T, indT] = ndim_prod(par{3:end}); % BxCxD...
    % первая матрица на произведения остальных
    [C, indC] = As_x_Bs(varargin{1:2},T,indT,indC);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROD AxB
% произведение двух матриц , строковая индексация
function [C, strC] = As_x_Bs(A,strA,B,strB,strC)

if nargin<5, strC = [];end

if isempty(strC) ,Ci = []; is_per = 0;
else ,[Ci] = str2ind(strC); is_per = 1;
end

[Ai,sumi] = str2ind(strA); % cell of indexes
[Bi] = str2ind(strB);

% comi,sumi,remAi,remBi,
% comAj,comBj,sumAj,sumBj,remAj,remBj
[~,sumAj] = ndim_ismember(sumi,Ai);% X in A
[is_in,sumBj] = ndim_ismember(sumi,Bi);% X in B
if ~all(is_in), error(['WRONG index string, no (sum) index in B']);
end
[Comi] = intersect(Ai,Bi);% com + sum
[comi] = setdiff(Comi,sumi);
[remAi,remAj] = setdiff(Ai,Comi); % work well fo []
[remBi,remBj] = setdiff(Bi,Comi); % individual 
[~,comAj] = ndim_ismember(comi,Ai);% common j
[~,comBj] = ndim_ismember(comi,Bi);

Ci_default = [comi;remAi;remBi]; % default order of indexes   , colunm of cells
if ~isempty(Ci) % check with default
    if ~all( ndim_ismember(Ci,Ci_default) ), error(['Result INDEX is WRONG, must be ',ind2str(Ci_default),' in some order. Now ',ind2str(Ci),'.']); end
else,    Ci = Ci_default;
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[C] = ndim_prod_123(A,comAj,sumAj,remAj,B,comBj,sumBj,remBj);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% permute
if is_per
    [~,per] = ndim_ismember(Ci,Ci_default);
    C = ndim_permute(C,per);
%     strC the same
else, strC = ind2str(Ci_default); %C the same 
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NDIM_PROD_IND2STR cell of indexes to string
%
%   [str] = NDIM_PROD_IND2STR(ind,is_X,is_1)
%       Function reverse to str2ind , create string from cell of indexes
%       and flags of summarizing and singularity. Space as delimeters.

% revers function to str2ind
function [str] = ind2str(ind,is_X)
if nargin<2, is_is_X = 0;else, is_is_X = 1;end
str = '';
for j =1:length(ind)
    s1 = ind{j};
    if is_is_X, if is_X(j), s1 =['(',s1,')']; end;end
    if j<length(ind), s1 = [s1,' '];end
    str = [str s1];
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STR2IND transform string to cell of indexes
%
%   [хind,indX] = NDIM_PROD_STR2IND(str)
%       Transform index string for ndim_prod into cell of all indexes,
%       cell of summarizing indexes (in bracked), bracked omitted, and cell of possibly
%       singular indexes (like 1a), sign 1 remain in names.
%       Delimeters are spaces or commas.

% remake of mmat_str2ind %2011-04-12
function [ind,indX] = str2ind(str)

if isempty(str), ind = [];indX=[]; return; % % for '' for scalar
end
str = replace(str,',',' ');% меняем запятые на пробелы
ind = split(strtrim(str));% режем по пробелам  со скобками
indX = [];
% выделяем те что в скобках, cкобки убираем
for j = 1:length(ind)
    s = ind{j};
    if s(1)=='(' && s(end)==')' % если скобки с двух сторон, убираем и ставим флаг
        ind{j} = s(2:end-1);
        indX{end+1} = ind{j};% ставим флаг
    elseif s(1)=='(' || s(end)==')', error(['WRONG index string, unclosed brace: ' , s]); % если закрыть забыли
    end
end
% проверяем уникальность
if length(ind)> length(unique(ind)),  error(['Wrong notation: repetition in indexes: ',str]);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% permute with cases pe=[],pe=[s1]
% same as in ndim_prod_123
function B = ndim_permute(A,pe)
if ~isempty(pe)
    if length(pe)==1, B = A;  
    else,  B = permute(A,pe);
    end
else
    if numel(A)==1, B = A;  else, error('WRONG permute size');
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ismember + ismember([],A)
function [is, j] = ndim_ismember(a,A)
if isempty(a), is = 1; j =[];
else,[is, j] = ismember(a,A);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % intersect + intersect([],[])
% function [C,iA,iB] = ndim_intersect(A,B)
% if isempty(A) || isempty(B), C=[];iA=[];iB=[];
% else, [C,iA,iB] = intersect(A,B);
% end
% end



