# ndim_prod
Product of multi-dimensional matrixes on Matlab

Any number of dimensions

Any number of common-static indexes

Any number of common-summarizing indexes

Any number of individael indexes

Example

A = rand(2,3,4);

B = rand(3,4,5,6);

C = rand(5,6,7,8,9);

D = ndin_prod( A,'a b (c)',B,'b c (d) (e)',C,'d e f g h','a b f g h');
