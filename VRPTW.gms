sets
v      /s1,cu1*cu9/
s(v)   /s1/
cu(v)  /cu1*cu9/
k      /k1*k6/
alias(v,i,j,n)
;

table c(i,j)
         s1       cu1     cu2     cu3     cu4     cu5    cu6     cu7     cu8     cu9
s1       0        5       5       5       5       5      5       5       5       5
cu1      5        0       5       5       5       5      5       5       5       5
cu2      5        5       0       5       5       5      5       5       5       5
cu3      5        5       5       0       5       5      5       5       5       5
cu4      5        5       5       5       0       5      5       5       5       5
cu5      5        5       5       5       5       0      5       5       5       5
cu6      5        5       5       5       5       5      0       5       5       5
cu7      5        5       5       5       5       5      5       0       5       5
cu8      5        5       5       5       5       5      5       5       0       5
cu9      5        5       5       5       5       5      5       5       5       0
;

table t(i,j)
         s1     cu1     cu2     cu3     cu4     cu5       cu6     cu7     cu8     cu9
s1       0      2.3     2.2     2.1     2.8     2.4       3.3     1.8     2.6     3.6
cu1      2.3    0       0.8     2.3     3.8     3.9       4.5     2.1     4.9     1.4
cu2      2.2    0.8     0       2.8     4.2     4.2       4.9     2.6     4.6     1.5
cu3      2.1    2.3     2.8     0       1.5     1.8       2.2     0.3     4.4     3.5
cu4      2.8    3.8     4.2     1.5     0       0.8       0.8     1.7     4.3     5
cu5      2.4    3.9     4.2     1.8     0.8     0         0.9     1.9     3.5     5.2
cu6      3.3    4.5     4.9     2.2     0.8     0.9       0       2.4     4.3     5.8
cu7      1.8    2.1     2.6     0.3     1.7     1.9       2.4     0       4.2     3.4
cu8      2.6    4.9     4.6     4.4     4.3     3.5       4.3     4.2     0       6.1
cu9      3.6    1.4     1.5     3.5     5       5.2       5.8     3.4     6.1     0
;


*t(i,j)=t(i,j)*0.5;

parameter
demand(j) /cu1 9,cu2 5,cu3 6,cu4 6,cu5 2,cu6 9,cu7 5,cu8 9,cu9 2/
ss(i)     /cu1 0.25,cu2 0.25,cu3 0.25,cu4 0.25,cu5 0.25,cu6 0.25,cu7 0.25,cu8 0.25,cu9 0.25/
a(i)      /cu1 2,cu2 1,cu3  3,cu4 3,cu5 3,cu6 3,cu7 2,cu8 2,cu9 2/
b(i)      /cu1 9,cu2 7.25,cu3 7,cu4 8,cu5 10,cu6 6,cu7 8,cu8 9,cu9 8/
tmax
q
;
tmax =15;
q=40;

*b(i) =b(i) +5;



free variable
z ;
binary variables
sigma(k) , y(i,j,k) ;
positive variable
x(i) ;

equations
of
c1
c2
c3
c4
c5
c6
c7
c8
c9
c10
c11
*c12
;

of..              z =e= sum( (k,i,j)$(cu(j)), c(i,j) * y(i,j,k) ) + 1000*sum(k, sigma(k) ) + sum( j, x(j)     ) ;

c1(k)..           sum( (i,j)$(cu(j) ), demand(j)*y(i,j,k) ) =l= q * sigma(k) ;

c2(k)..           sum(j$(cu(j)) ,y('s1',j,k) ) =e= sigma(k) ;

c3(k)..           sum(i$(cu(i)),y(i,'s1',k) ) =e= sigma(k) ;

c4(k,j)$(cu(j)).. sum(i,y(i,j,k)) =e= sum(i,y(j,i,k)) ;

c5(j)$(cu(j))..   sum((k,i), y(i,j,k)) =e= 1 ;

c6..    x("s1") =e= 0 ;

c7(i)$(cu(i))..   x(i) =g= a(i) ;

c8(i)$(cu(i))..   x(i) =l= b(i) ;

c9(K,I,J)$( cu(j) and ( ord(i) ne ord(j) ) )..    x(j) =g= x(i)+(ss(i) + t(i,j))*y(i,j,k) - tmax*(1-y(i,j,k)) ;

c10(k,i)$( cu(i)) ..                              x(i) + (ss(i)+t(i,'s1'))*y(i,'s1',k) =l= tmax + tmax*(1-y(i,'s1',k));
c11(k,i,j)$(ord(i) eq ord(j))..                   y(i,j,k)=e=0;

*c12(K,I)..  x('cu7') =e= x(i)+(ss(i) + t(i,'cu7'))*y(i,'cu7',k) - tmax*(1-y(i,'cu7',k)) ;
y.fx(i,i,k)=0;


model ali /all/;
options mip=cplex , optcr=0;
solve ali us mip min z ;
display z.l , y.l,x.l ;

parameter
aa
bb;

aa=  sum( (k,i,j)$(cu(j)), c(i,j) * y.l(i,j,k) ) + 1000*sum(k, sigma.l(k) ) ;
bb =  sum( j$( cu(j) ), x.l(j)     )      ;

display aa,bb;
