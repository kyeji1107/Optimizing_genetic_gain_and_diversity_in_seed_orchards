##Box 1
#Parameters and sets
param nc > 0;
set S := 1..nc by 1;
param ns > 0;
param b {i in S};
param l {i in S} >= 0, <= 1;
param u {i in S} >= 0, <= 1;
param c {i in S, j in S} >= 0, <= 1;

#Variables
var p {i in S} >= 0, <= 1;
var r {i in S} binary; 

#Objective function
maximize gain: sum{i in S} b[i] * p[i];

#Contraints
s.t. c1: sum {i in S,j in S} 
     c[i,j] * p[i] * p[j] <= 1/(2*ns);
s.t. c2: sum {i in S} p[i] = 1;
s.t. c3a {i in S}: p[i] >= r[i] * l[i];
s.t. c3b {i in S}: p[i] <= r[i] * u[i];


##Box 2
#Constraints
s.t. c4 {i in S}: 
     p[i] == 0 or l[i] <= p[i] <= u[i];


##Box 3
#Parameters and sets
param cu >= 0, <= 1;

#Constraints
s.t. c5 {i in S, j in S: i < j}:
     if (p[i] * p[j]) > 0 then c[i,j] <= cu;


##Box 4
#Variables
var y {i in S, j in S} binary;

#Constraints
s.t. c6 {i in S, j in S : i < j and c[i,j] > cu}: 
p[j] <= 1 - y[i,j];
s.t. c7 {i in S, j in S : i < j and c[i,j] > cu}: 
p[i] <= y[i,j];


##Box 5
#Variables
var f {i in S} >= 0, <= 1;
var m {i in S} >= 0, <= 1;

#Contraints
s.t. c8 {i in S}: p[i] = (f[i]+m[i])/2;


##Box 6
#Parameters and sets
param nd > 0;

#Constraints
s.t. c9: {i in (nc-nd+1)..nc}: f[i] = 0;
s.t. c11 {i in 1..(nc-nd)}: p[i] = (f[i]+m[i])/2;
s.t. c12 {i in (nc-nd+1)..nc}: p[i] = m[i];
