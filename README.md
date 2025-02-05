# ndim_prod
Product of multi-dimensional matrixes on Matlab

Any number of dimensions, of common-static indexes, of common-summarizing indexes, of individual indexes may be used.

Example

A = rand(2,3,4);

B = rand(3,4,5,6);

C = rand(5,6,7,8,9);

D = ndim_prod( A,'a b (c)',B,'b c (d) (e)',C,'d e f g h','a b f g h');
