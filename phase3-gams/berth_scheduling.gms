* Sensitivity-analysis results will be written (appended) to this file.
$set excelFile Sensitivity_Results.xlsx


Sets
    i   Ships           /1*8/
    j   Docks           /1*4/
    k   Objectives      /1*3/
    bin BinarySet       /1*2/;

* Define set iq such that i <> q
Alias (i,q),(i,ii);

Scalar
    big_M  "Big M constant" /80/

Parameters
    C(j)      "Dock capacity"   /1 57, 2 40, 3 52, 4 30/
    CA(i)     "Containers on each ship"  /1 6, 2 18, 3 20, 4 5, 5 11, 6 19, 7 20, 8 13/
    L(j)      "Max length at each dock"    /1 118, 2 236, 3 120, 4 138/
    LS(i)     "Length of each ship"      /1 114, 2 126, 3 115, 4 102, 5 102, 6 131, 7 119, 8 138/
    P(i)      "Priority of each ship"   /1 0.963797, 2 0.773138, 3 0.441551, 4 0.981939, 5 0.190926, 6 0.348681, 7 0.970487, 8 0.445607/
    O(j)      "Lag time at each dock"   /1 6.922094, 2 9.727825, 3 9.564813, 4 5.748308/
    Un(i,j)   "Unload feasibility (or cost/indicator) between ship i and dock j"
    a(i)      "Arrival time of each ship"  /1 4.76553, 2 1.00459, 3 3.653, 4 6.05418, 5 1.86919, 6 2.11799, 7 4.62295, 8 6.60118/
    CIC(j)    "Cost for extra capacity at dock j"  /1 5.87154, 2 5.880447, 3 8.864247, 4 6.64032/
    ext(j)    "Max extra allowable capacity at dock j" /1 29, 2 15, 3 49, 4 17/;



Table Un(i,j)
       1             2             3             4
    1  9.306945068   8.320459113   9.041338719   8.655545623
    2  8.499993353   9.905633818   9.993113985   8.089112765
    3  9.720322075   9.206381222   8.763211972   8.567236436
    4  9.349929694   8.913662302   9.371722971   9.32369264
    5  8.265956289   9.535675628   9.964826498   9.938776321
    6  9.226653641   8.088521266   8.008110288   8.267945054
    7  9.882004543   8.605721124   8.732291203   9.796392489
    8  8.62872761    9.097964368   8.872061915   8.129988352;
* Similarly for other parameters if 2D or more convenient

Binary Variables
    X(i,j)     "Assignment of ship i to dock j"
    y(i,q,j)   "Ordering of ships i and q at dock j (i <> q)"
    f(i,bin)   "Binary indicator for ship i (e.g., for disjunctive constraints)";

Positive Variables
    W(i)       "Waiting time or related measure for ship i",
    d(i,j)     "Delay or slack variable for ship i at dock j",
    u(i,k)     "Auxiliary continuous variable for ship i and objective k",
    s(i)       "Start time or similar continuous variable for ship i";

Integer Variables
    IC(j)      "Extra capacity used at dock j";

Free Variable z;

Equations
obj
c1
c2eq
c2le
c3
c4
c5
c6
c7
c8
c9
c10
c11
c12a
c12b
c13
c14a
c14b
c15a
c15b
c15c
c16;
* c1: Each ship assigned to exactly one dock
c1(i).. sum(j, X(i,j)) =e= 1;

* c2: Number of ships assigned to each dock
* For j = 3
c2eq(j)$ (ord(j) = 3)..
    sum(i, X(i,j)) =e= 3;

* For j <> 3
c2le(j)$ (ord(j) <> 3)..
    sum(i, X(i,j)) =l= 2;

* c3: Total containers at each dock within capacity + extra
c3(j)..
    sum(i, X(i,j)*CA(i)) =l= C(j) + IC(j);

* c4: Extra capacity must be within allowed extension
c4(j)..
    IC(j) =l= ext(j);

* c5: Ship length must not exceed dock capacity
c5(i,j)..
    LS(i)*X(i,j) =l= L(j);

* c6: Maximum wait time based on priority
c6(i)..
    W(i) =l= 80*(1 - P(i));

* c7: Wait time equals start time minus arrival
c7(i)..
    W(i) =e= s(i) - a(i);

* c8: Non-overlapping ship scheduling on same dock
c8(i,q,j)$(ord(i) <> ord(q))..
    s(i) + Un(i,j) + O(j) =l= s(q) + big_M*(1 - y(i,q,j));

* c9: At most one ordering per pair
c9(i,q,j)$(ord(i) <> ord(q))..
    y(i,q,j) + y(q,i,j) =l= 1;

* c10: If both ships assigned to same dock, one must precede the other
c10(i,q,j)$(ord(i) <> ord(q))..
    y(i,q,j) + y(q,i,j) =g= X(i,j) + X(q,j) - 1;

* c11: Conditional linking between y and X
c11(i,q,j)$(ord(i) <> ord(q))..
    y(i,q,j) + y(q,i,j) =l= X(i,j);

* c12: Waiting time zero if alone at dock (dock 3 special case)
c12a(i,q,j)$(ord(i) <> ord(q) and ord(j) <> 3)..
    W(i) =l= 80 * (1 - y(i,q,j));

c12b(i,j)$(ord(j) = 3)..
    W(i) =l= 80 * (2 - sum(ii$(ord(ii) <> ord(i)), y(i,ii,j)));

* c13: W(i) = sum of u(i,k)
c13(i)..
    W(i) =e= sum(k, u(i,k));

* c14: u(i,k) = {10, 30} * f(i,bin)
c14a(i,bin,k)$(ord(bin) = 1 and ord(k) = 1)..
    10 * f(i,bin) =l= u(i,k);

c14b(i,bin,k)$(ord(bin) = 2 and ord(k) = 2)..
    30 * f(i,bin) =l= u(i,k);


* c15: upper bounds on u(i,k) based on f(i,bin)
c15a(i,k)$ (ord(k) = 1)..
    u(i,k) =l= 10;
c15b(i,k)$(ord(k) = 2)..
    u(i,k) =l= 30 * f(i,'1');
c15c(i,k)$ (ord(k) = 3)..
    u(i,k) =l= 30*f(i,'2');

* c16: Binary order constraint
c16(i)..
    f(i,'1') =g= f(i,'2');

obj..
    z =e=
        25 * sum(i, u(i,'1'))
      + 10 * sum(i, u(i,'2'))
      +  5 * sum(i, u(i,'3'))
      +  5 * sum((i,j), Un(i,j) * X(i,j))
      +     sum(j, CIC(j) * IC(j));

Model shipdock /all/;
Solve shipdock using mip minimizing z;
display w.l;

sets
 sen_anal_1 /1*11/,
 sen_anal_2 /1*21/,
 sen_anal_3 "can be used for the fourth sen analysis as well" /1*36/,
 sen_anal_5 /1*51/;
parameter
 O_sen_anal(sen_anal_1,j) "Modified Lag Time at Each Dock" /1.2 5/,
 UN_sen_anal (sen_anal_2,i,j) "Modified Unloading Time for each ship at each dock" /1.7.3 5/,
 C_sen_anal (sen_anal_3, j) "Modified Dock Capacity (it can be used for the third and forth sen analysis)" /1.2 5
                                                                                                            1.3 5/,
 a_sen_anal (sen_anal_5, i) "Modified Arrival Time of each ship"/1.6 1/,
 z_sen_anal_1 (sen_anal_1),
 z_sen_anal_2 (sen_anal_2),
 z_sen_anal_3 (sen_anal_3),
 z_sen_anal_4 (sen_anal_3),
 z_sen_anal_5 (sen_anal_5);

* --- Sensitivity scenario result parameters ---
Parameter
    X_sa1(sen_anal_1,i,j)    "X for each scenario (analysis 1)",
    y_sa1(sen_anal_1,i,q,j)    "y for each scenario (analysis 1)",
    f_sa1(sen_anal_1,i,bin)    "f for each scenario (analysis 1)",
    W_sa1(sen_anal_1,i)    "W for each scenario (analysis 1)",
    d_sa1(sen_anal_1,i,j)    "d for each scenario (analysis 1)",
    u_sa1(sen_anal_1,i,k)    "u for each scenario (analysis 1)",
    s_sa1(sen_anal_1,i)    "s for each scenario (analysis 1)",
    IC_sa1(sen_anal_1,j)    "IC for each scenario (analysis 1)",
    X_sa2(sen_anal_2,i,j)    "X for each scenario (analysis 2)",
    y_sa2(sen_anal_2,i,q,j)    "y for each scenario (analysis 2)",
    f_sa2(sen_anal_2,i,bin)    "f for each scenario (analysis 2)",
    W_sa2(sen_anal_2,i)    "W for each scenario (analysis 2)",
    d_sa2(sen_anal_2,i,j)    "d for each scenario (analysis 2)",
    u_sa2(sen_anal_2,i,k)    "u for each scenario (analysis 2)",
    s_sa2(sen_anal_2,i)    "s for each scenario (analysis 2)",
    IC_sa2(sen_anal_2,j)    "IC for each scenario (analysis 2)",
    X_sa3(sen_anal_3,i,j)    "X for each scenario (analysis 3)",
    y_sa3(sen_anal_3,i,q,j)    "y for each scenario (analysis 3)",
    f_sa3(sen_anal_3,i,bin)    "f for each scenario (analysis 3)",
    W_sa3(sen_anal_3,i)    "W for each scenario (analysis 3)",
    d_sa3(sen_anal_3,i,j)    "d for each scenario (analysis 3)",
    u_sa3(sen_anal_3,i,k)    "u for each scenario (analysis 3)",
    s_sa3(sen_anal_3,i)    "s for each scenario (analysis 3)",
    IC_sa3(sen_anal_3,j)    "IC for each scenario (analysis 3)",
    X_sa4(sen_anal_3,i,j)    "X for each scenario (analysis 4)",
    y_sa4(sen_anal_3,i,q,j)    "y for each scenario (analysis 4)",
    f_sa4(sen_anal_3,i,bin)    "f for each scenario (analysis 4)",
    W_sa4(sen_anal_3,i)    "W for each scenario (analysis 4)",
    d_sa4(sen_anal_3,i,j)    ,
    u_sa4(sen_anal_3,i,k)    "u for each scenario (analysis 4)",
    s_sa4(sen_anal_3,i)    "s for each scenario (analysis 4)",
    IC_sa4(sen_anal_3,j)    "IC for each scenario (analysis 4)",
    X_sa5(sen_anal_5,i,j)    "X for each scenario (analysis 5)",
    y_sa5(sen_anal_5,i,q,j)    "y for each scenario (analysis 5)",
    f_sa5(sen_anal_5,i,bin)    "f for each scenario (analysis 5)",
    W_sa5(sen_anal_5,i)    "W for each scenario (analysis 5)",
    d_sa5(sen_anal_5,i,j)    "d for each scenario (analysis 5)",
    u_sa5(sen_anal_5,i,k)    "u for each scenario (analysis 5)",
    s_sa5(sen_anal_5,i)    "s for each scenario (analysis 5)",
    IC_sa5(sen_anal_5,j)    "IC for each scenario (analysis 5)";
* --- end parameter block ---
scalars
 base_O /5/,
 base_UN /5/,
 base_C /5/,
 base_a /1/;

loop(sen_anal_1,
  O_sen_anal(sen_anal_1, '2') = base_O + 0.5 * (ord(sen_anal_1) - 1);
);



loop(sen_anal_2,
  UN_sen_anal(sen_anal_2,'7','3') = base_UN + 0.5 * (ord(sen_anal_2) - 1);
);



*loop (sen_anal_2, UN_sen_anal(ord(sen_anal_2) + 1, 7, 3) = UN_sen_anal(sen_anal_2, 7, 3) +0.5);
loop(sen_anal_3,
  C_sen_anal(sen_anal_3, '2') = base_C + (ord(sen_anal_3) - 1);
);


*loop (sen_anal_3, C_sen_anal(ord(sen_anal_3) + 1, 2) = C_sen_anal(sen_anal_3, 2) +1);
loop(sen_anal_3,
  C_sen_anal(sen_anal_3, '3') = base_C + (ord(sen_anal_3) - 1);
);


*loop (sen_anal_3, C_sen_anal(ord(sen_anal_3) + 1, 3) = C_sen_anal(sen_anal_3, 3) +1);

loop(sen_anal_5,
  a_sen_anal(sen_anal_5,'6') = base_a + 0.1*(ord(sen_anal_5) -1);
);


*loop (sen_anal_5, a_sen_anal(ord(sen_anal_5) + 1, 6) = a_sen_anal(sen_anal_5, 6) +0.1);
loop (sen_anal_1,
      O('2') = O_sen_anal(sen_anal_1, '2');
      solve shipdock using MIP minimizing z;

      X_sa1(sen_anal_1,i,j)   = X.l(i,j);
      y_sa1(sen_anal_1,i,q,j) = y.l(i,q,j);
      f_sa1(sen_anal_1,i,bin) = f.l(i,bin);
      W_sa1(sen_anal_1,i)     = W.l(i);
      u_sa1(sen_anal_1,i,k)   = u.l(i,k);
      s_sa1(sen_anal_1,i)     = s.l(i);
      IC_sa1(sen_anal_1,j)    = IC.l(j);
*d_sa1(sen_anal_1,i,j)   = d.l(i,j);

      z_sen_anal_1 (sen_anal_1) = z.l;
);
* --- Export results for sensitivity analysis 1 ---
execute_unload 'Sensitivity_sa1.gdx'
    X_sa1, y_sa1, f_sa1, W_sa1, d_sa1, u_sa1, s_sa1, IC_sa1, z_sen_anal_1
    ;


loop (sen_anal_2,
      UN('7','3') = UN_sen_anal(sen_anal_2, '7','3');
      solve shipdock using MIP minimizing z;

      X_sa2(sen_anal_2,i,j)   = X.l(i,j);
      y_sa2(sen_anal_2,i,q,j) = y.l(i,q,j);
      f_sa2(sen_anal_2,i,bin) = f.l(i,bin);
      W_sa2(sen_anal_2,i)     = W.l(i);
*d_sa2(sen_anal_2,i,j)   = d.l(i,j);
      u_sa2(sen_anal_2,i,k)   = u.l(i,k);
      s_sa2(sen_anal_2,i)     = s.l(i);
      IC_sa2(sen_anal_2,j)    = IC.l(j);

      z_sen_anal_2(sen_anal_2) = z.l;
);
* --- Export results for sensitivity analysis 2 ---
execute_unload 'Sensitivity_sa2.gdx'
    X_sa2, y_sa2, f_sa2, W_sa2, d_sa2, u_sa2, s_sa2, IC_sa2, z_sen_anal_2
    ;


loop (sen_anal_3,
      C('2') = C_sen_anal(sen_anal_3, '2');
      solve shipdock using MIP minimizing z;

      X_sa3(sen_anal_3,i,j)   = X.l(i,j);
      y_sa3(sen_anal_3,i,q,j) = y.l(i,q,j);
      f_sa3(sen_anal_3,i,bin) = f.l(i,bin);
      W_sa3(sen_anal_3,i)     = W.l(i);
*d_sa3(sen_anal_3,i,j)   = d.l(i,j);
      u_sa3(sen_anal_3,i,k)   = u.l(i,k);
      s_sa3(sen_anal_3,i)     = s.l(i);
      IC_sa3(sen_anal_3,j)    = IC.l(j);

      z_sen_anal_3(sen_anal_3) = z.l;
);
* --- Export results for sensitivity analysis 3 ---
execute_unload 'Sensitivity_sa3.gdx'
    X_sa3, y_sa3, f_sa3, W_sa3, d_sa3, u_sa3, s_sa3, IC_sa3, z_sen_anal_3
    ;



loop (sen_anal_3,
      C('3') = C_sen_anal(sen_anal_3, '3');
      solve shipdock using MIP minimizing z;

      X_sa4(sen_anal_3,i,j)   = X.l(i,j);
      y_sa4(sen_anal_3,i,q,j) = y.l(i,q,j);
      f_sa4(sen_anal_3,i,bin) = f.l(i,bin);
      W_sa4(sen_anal_3,i)     = W.l(i);
*d_sa4(sen_anal_3,i,j)   = d.l(i,j);
      u_sa4(sen_anal_3,i,k)   = u.l(i,k);
      s_sa4(sen_anal_3,i)     = s.l(i);
      IC_sa4(sen_anal_3,j)    = IC.l(j);

      z_sen_anal_4(sen_anal_3) = z.l;
);
* --- Export results for sensitivity analysis4 ---
execute_unload 'Sensitivity_sa4.gdx'
    X_sa4, y_sa4, f_sa4, W_sa4, d_sa4, u_sa4, s_sa4, IC_sa4, z_sen_anal_4
    ;


loop (sen_anal_5,
      a('6') = a_sen_anal(sen_anal_5, '6');
      solve shipdock using MIP minimizing z;

      X_sa5(sen_anal_5,i,j)   = X.l(i,j);
      y_sa5(sen_anal_5,i,q,j) = y.l(i,q,j);
      f_sa5(sen_anal_5,i,bin) = f.l(i,bin);
      W_sa5(sen_anal_5,i)     = W.l(i);
*      d_sa5(sen_anal_5,i,j)   = d.l(i,j);
      u_sa5(sen_anal_5,i,k)   = u.l(i,k);
      s_sa5(sen_anal_5,i)     = s.l(i);
      IC_sa5(sen_anal_5,j)    = IC.l(j);

      z_sen_anal_5(sen_anal_5) = z.l;
);
* --- Export results for sensitivity analysis 5 ---
execute_unload 'Sensitivity_sa5.gdx'
    X_sa5, y_sa5, f_sa5, W_sa5, d_sa5, u_sa5, s_sa5, IC_sa5, z_sen_anal_5
    ;


display O_sen_anal, UN_sen_anal, C_sen_anal, a_sen_anal,
 z_sen_anal_1, z_sen_anal_2, z_sen_anal_3, z_sen_anal_4,z_sen_anal_5;
