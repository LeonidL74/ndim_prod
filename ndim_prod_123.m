% NDIM_PROD_123 
% [C] = ndim_prod_123(A,comAj,sumAj,remAj,B,comBj,sumBj,remBj )
%       Production of matrixes A and B with positioned static, summirising and
%       individual idexes. 
% Here
% comAj, comBj - positions of static indexes of A and B  
% sumAj, sumBj - positions of summarizing indexes of A and B
% remAj, remBj - positions of individual indexes
% Resulting matrix C have indexes in following order: com,remAj,remBj
% Let D(scalar)=0,D(row or column)=1, D(A)=ndims(A) - real ndims
% Let csrAj = [comAj,sumAj,remAj] - nion of positions 
% Let D_need = max(csrAj) - needed number of dimensions
% Then must be
% length(csrAj)==length(unique(csrAj)) - all positions must be unique
% D_need >= D, - needed D must be equal or larger then real D. 
% max(csrAj)=length(csrAj) - they must be full, i.e. sirt(csrAj) = 1,2,3,4,...,max
%
% RUS
% Произведение матриц с указанием номеров стационарных,
% суммирующих и индивидуальных индексов
% допустимы любые мерности: число,строка, столбец, матрица любой мерности
% здесь
% comAj, comBj - номера стационарных индексов в A и B
% sumAj, sumBj - номера суммирующих индексов в A и B
% remAj, remBj - номера индивидуальных индексов
% Резултирующая матрици имеет порядок индексов: com,remAj,remBj
% полагается что мерность D(число)=0,D(строка, вектор)=1, D(A)=ndims(A)
% пусть D_need = max([comAj,sumAj,remAj]) - требуемая мерность
% должы выполняться равенство
% Пусть D(scalar)=0,D(row or column)=1, D(A)=ndims(A) - реальная мерность
% Let csrAj = [comAj,sumAj,remAj] - объединение позиций
% Let D_need = max(csrAj) - требуемая мерность
% Тогда должно выполняться
% length(csrAj)==length(unique(csrAj)) - все позиции уникальны
% D_need >= D, - требуемая мерность не меньше реальной
% max(csrAj)=length(csrAj) - позиции полны, т.е. sirt(csrAj) = 1,2,3,4,...,max
%
%
% Examples 
%     ndim_prod_123( rand(2,3,4),[],[2 3],[1],rand(3,4,5),[],[1 2],3)
%     ndim_prod_123( rand(1,2,3,4),[1],1+[2 3],1+[1],rand(1,3,4,5),[1],1+[1 2],1+3)
%
%   See also ndim_prod, ndim_disp
%
%   (c) Lazarev Leonid

function [C] = ndim_prod_123(A,comAj,sumAj,remAj,B,comBj,sumBj,remBj )

if nargin==0% check %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ndim_prod_123( rand(2,3,4),[],[2 3],[1],rand(3,4,5),[],[1 2],3)
    ndim_prod_123( rand(1,2,3,4),[1],1+[2 3],1+[1],rand(1,3,4,5),[1],1+[1 2],1+3)
    ndim_prod_123( rand(1,2,3,4),[1],1+[2 3],1+[1],rand(1,3,4,5),[1],1+[1 2],1+3)
    
    % return
    % все размерности одинаковы
    r = 1;
    AA = {4,rand(r,1),rand(r,r),rand(r,r,r),rand(r,r,r,r)};
    for nc = 0:4;
        for ns = 0:4;
            for j = (nc+ns+1):length(AA)
                for jj = (nc+ns+1):length(AA)
                    %    [nc,ns,j,jj] = deal(   0   ,  2 ,    3   ,  4);
                    disp( [nc,ns,j, jj] )
                    [C] = ndim_prod_123( AA{j},1:nc,nc+(1:ns),(nc+ns)+(1:ndim_ndims(AA{j})-nc-ns),AA{jj},1:nc,nc+(1:ns),(nc+ns)+(1:ndim_ndims(AA{jj})-nc-ns) )  ;
                end
            end
        end
    end
    return
end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

max_numel_com_sum_rem_rem = 10000000; % если требуется большая по размеру матрица, то переходим к циклу по com
% to row
comAj = comAj(:).';sumAj = sumAj(:).';remAj = remAj(:).';comBj = comBj(:).';sumBj = sumBj(:).';remBj = remBj(:).';

% CHECK %%%%%%%%%%%%%
[dA,~,siA] = check_indj(A,[comAj,sumAj,remAj]);% проверка своих индексов с учетом возможности сингулярных размерностей
[dB,~,siB] = check_indj(B,[comBj,sumBj,remBj]);

is_com = ~isempty([comAj,comBj]);% наличие стационарных индексов
if is_com % проверяем совпадение размерностей по ним
    if length(comAj)~=length(comBj), error('WRONG com indexes,length(comAj)~=length(comBj)');end
    if any(siA(comAj)~=siB(comBj)), error('WRONG size, any(siA(comAj)~=siB(comBj))');end
end
is_sum = ~isempty([sumAj,sumBj]);% наличие суммирующих индексов
if is_sum % проверяем совпадение размерностей по ним
    if length(sumAj)~=length(sumBj), error('WRONG com indexes,length(comAj)~=length(comBj)');end
    if any(siA(sumAj)~=siB(sumBj)), error('WRONG size, any(siA(comAj)~=siB(comBj))');end
end
%%%%%%%%%%%%%%%%%%%%%

is_remA = ~isempty(remAj);% наличие индивидуальных индексов у A
is_remB = ~isempty(remBj);% наличие индивидуальных индексов у B

sicom = siA(comAj);% размер по стационарным размерностям
sisum = siA(sumAj);% размер по суммирующим размерностям
siremA = siA(remAj);% размер по индивидуальным размерностям у А
siremB = siB(remBj);% размер по индивидуальным размерностям у В

% сжимаем A,B до 1D,2D,3D
% p_com==[] for empty sicom, not 1 as prod([])!
if ~isempty(comAj), p_com = prod(sicom); else, p_com = [];end % произведение размерностей общих
if ~isempty(sumAj), p_sum = prod(sisum); else, p_sum = [];end% произведение размерностей суммы
if ~isempty(remAj), p_remA = prod(siremA);else, p_remA = [];end% произведение индивидуальных
if ~isempty(remBj), p_remB = prod(siremB);else, p_remB = [];end% произведение индивидуальных

pA = ndim_permute(A,[comAj,sumAj,remAj]);% в порядке com-sum-rem
pB = ndim_permute(B,[comBj,sumBj,remBj]);%

% сжимаем размерности если необходимо
if max([length(comAj),length(sumAj),length(remAj),length(remBj)])>1 % если есть не одиночные размеры
    rA = ndim_reshape( pA,[p_com, p_sum,p_remA]);% сжатые матрицы
    rB = ndim_reshape( pB,[p_com, p_sum,p_remB]);
    is_reshape = 1;% потом их распакуем
else
    rA = pA;rB = pB;  
    is_reshape = 0;% как есть, не нужно сжимать
end

% row -> column
if dA == 1, rA = rA(:);end% column % если строка переводим в столбец
if dB == 1, rB = rB(:);end% column

% is_com is_sum is_remA is_remB
switch 100*(10*is_com + is_sum) + 10*is_remA + is_remB
    case 0000, C = rA*rB; %a a, a x a
    case 0001, C = rA*rB; %a r, a x r
    case 0010, C = rA*rB; %r a, r x a
    case 0011, C = rA*rB.'; %r r, r x r столбец на строку
    case 0100, C = rA.'*rB; %s s, s x s, строка на столбец
    case 0101, C = rB.'*rA; % s sr , rs x s, матрица на столбец
    case 0110, C = rA.'*rB; % sr s, rs x s, матрица на столбец
    case 0111, C =  rA.'*rB;% sr sr, rs x sr, матрица на матрицу
    case 1000, C = rA.*rB; % c c, c x c% общие вектор x вектор
    case 1001, C = (repmat(rA(:),[1,p_remB])).*rB; % c crB, crB x crB% дополняем до size(B)
    case 1010, C = rA.*(repmat(rB(:),[1,p_remA]));%% cr c, crA x crA дополняем до size(A)
    case 1011, C = (rA(:,:,ones(1,p_remB))).*(permute( rB(:,:,ones(1,p_remA)),[1 3 2]));% crA crB, crArB x crArB добавляем c,r,+r
    case 1100, C = sum(rA.*rB,2); % cs cs cs x cs общие вектор x вектор
    case 1101, C = permute( sum( (rA(:,:,ones(1,p_remB))).*rB,2), [1 3 2]);%% cs csr, csrB x csrB добавляем c,r,+r% csrArB -> c1rArB -> crArB
    case 1110, C = permute( sum((rB(:,:,ones(1,p_remA))).*rA,2), [1 3 2]);% % csr cs, csrA csrA добавляем c,r,+r
    case 1111,% csr csr, csrArB x csrArB
        if p_com*p_sum*p_remA*p_remB <=max_numel_com_sum_rem_rem % если не много элементов
            C = sum( (permute( rA(:,:,:,ones(1,p_remB)),[1,3,4,2] )).*(permute( rB(:,:,:,ones(1,p_remA)),[1,4,3,2] )),4);%csrArB -> crArBs %csrBrA -> crArBs
        else % в цикле если много элементов
            C = zeros(p_com,p_remA,p_remB);% crr
            for jcom = 1:p_com, 
                C(jcom,:,:) = sum((permute( rA(jcom,:,:,ones(1,p_remB)),[1,3,4,2])).*(permute( rB(jcom,:,:,ones(1,p_remA)),[1,4,3,2])),4);% jsrArB -> jrArBs, % jsrBrA -> jrArBs
            end
        end
end

if is_reshape, C = ndim_reshape(C,[sicom,siremA,siremB]);% распаковываем
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NDIM_CHECK_INDJ 
% проверка совместимости матрицы и номеров индексов, и их уникальность
% + размеры матрицы с учетом требуемых индексов
% indj = [comj,remj,sumj] - общие-несуммирующие, индивидуальные, суммирующие
function [dA,dA_real,siA,is_ok] = check_indj(A,indj,is_terminate)

is_ok = 1;
dA = length(indj); % требуемая мерность
dA_real = ndim_ndims(A);% реальная мерность
si = size(A); % размерности
% проверка не мало ли индексов
if max(indj)~=length(indj), if is_terminate, error( ['WRONG INDEX, not full, max(indj)~=length(indj), ', num2str(indj)]);else, is_ok = 0;end;end
if isempty(A), if is_terminate, error(['matrix is empty []']);else, is_ok = 0;end;end
if dA_real>dA, if is_terminate, error(['dA_real>dA - real dim lage then on ind ',num2str(remAj)]);else, is_ok = 0;end;end
% проверка на уникальность
if length(indj)>length(unique(indj)),if is_terminate, error( ['WRONG INDEX, not unique ', num2str(indj)]);else, is_ok = 0;end;end
% все числа не больше чем требуется мерность
if any (indj>dA),if is_terminate, error( ['WRONG INDEX, index exseed D=',num2str(dA),' ind= ', num2str(indj)]);else, is_ok = 0;end;end
if dA>length(si), siA = [si , ones(1,dA-length(si))]; % [s s 1 1 1 ]% если требуется больше чем дает матлаб и dA>2
elseif dA==length(si), siA = si; % если совпадает, то ok
else% dA<length(si) % если меньше но dA_real<=dA, то dA  =0,1
    if dA == 1,  siA = numel(A);  else, siA =[];%dA = 0;
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% permute with cases pe=[],pe=[s1]
function B = ndim_permute(A,pe)
if ~isempty(pe)
    if length(pe)==1, B = A;  else,  B = permute(A,pe);
    end
else
    if numel(A)==1, B = A;  else, error('WRONG permute size');
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reshape with cases re=[],re=[s1]
function B = ndim_reshape(A,re)
if ~isempty(re)
    if length(re)==1, re = [re 1]; end
    B = reshape(A,re);
else,% []
    if numel(A)==1, B = A;
    else, error('WRONG reshape size');
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NDIM_NDIMS - Number of dimensions with 0 for scalar, 1 for column or ow
%   N = NDIMS(X) returns the number of dimensions in the matrix X.
%   The number of dimensions in an array may be 0, 1, 2 or greater.
%   Trailing singleton dimensions are ignored.
%   Put simply, it is NDIMS(X) - (NDIMS(X)==2)*SUM(SIZE(X)==1).
%
% Examples 
% ndim_ndims(rand(2,2,2)) ,%  returns 3
% ndim_ndims(rand(2,2))   ,%returns 2
% ndim_ndims(rand(2,1))   ,%returns 1
% ndim_ndims(rand(1,2))   ,%returns 1
% ndim_ndims(rand)   ,%returns 0

% NDIM_NDIMS + + +
function D = ndim_ndims(A)
if isempty(A), D=[];
else, D = ndims(A);
    if D==2, D = D - sum(size(A)==1);
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

