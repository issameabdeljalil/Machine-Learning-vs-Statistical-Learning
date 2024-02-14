/* define a store a SAS/IML function that computes the Iman-Conover transformation */
proc iml;
/* Input: X is a data matrix with k columns
          C is a (k x k) rank correlation matrix
   Output: A matrix W. The i_th column of W is a permutation of the i_th
          columns of X. The rank correlation of W is close to the 
          specified correlation matrix, C.
*/          
start ImanConoverTransform(X, C);
   N = nrow(X);
   S = J(N, ncol(X));
   /* T1: Create normal scores of each column */
   do i = 1 to ncol(X);
      ranks = ranktie(X[,i], "mean");          /* tied ranks */
      S[,i] = quantile("Normal", ranks/(N+1)); /* van der Waerden scores */
   end;
   /* T2: apply two linear transformations to the scores */
   CS = corr(S);        /* correlation of scores */
   Q = root(CS);        /* Cholesky root of correlation of scores */
   P = root(C);         /* Cholesky root of target correlation */
   T = solve(Q,P);      /* same as  T = inv(Q) * P; */
   Y = S*T;             /* transform scores: Y has rank corr close to target C */
 
   /* T3: Permute or reorder data in the columns of X to have the same ranks as Y */
   W = X;
   do i = 1 to ncol(Y);
      rank = rank(Y[,i]);          /* use ranks as subscripts, so no tied ranks */
      tmp = W[,i]; call sort(tmp); /* sort column by ranks */
      W[,i] = tmp[rank];           /* reorder the column of X by the ranks of M */
   end;
   return( W );
finish;
store module=(ImanConoverTransform);  /* store definition for later use */
quit;
