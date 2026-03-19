1. keeping y = 0.d0 while changing x can give us a clear cut-off distance. also A=3.35 gave the cut-off for us at 2.15 equivalent to GWU's 4.33
2. Due to some reason y was kept at 1.5 hence the behavior of ewald sum changed.
3. Changing ka=0.5 made the swap_vel graph much steeper as opposed to the original ka=0.8.
4. Trying to match the cut-off distance to 2.15 in both kappas-by varing the other parameters judiciously. Then will use them both to see what simulations we get.

5. Set of parameters for ka=0.5 (to make cut-off distance 2.15):
   A : 1.745
   b : 0.65

6. Set of parameters for ka=0.8 (to make cut-off distance 2.15):
   A : 3.325
   b : 0.65
   
