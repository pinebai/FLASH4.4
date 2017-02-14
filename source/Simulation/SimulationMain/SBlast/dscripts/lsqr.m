function x = lsqr(A,b,opt)

[Am,An] = size(A);
if (Am < An)
   printf('lsqr.m: Not enough rows in A.\n');
   break;
end

[bm,bn] = size(b);
if ((bm != Am)||(bn != 1))
   printf('lsqr.m: b has incorrect size.\n');
   break;
end

if (opt==0)
   [Q,R] = qr(A);
   x = R \ (Q'*b);
else
   [U,S,V] = svd(A);

   [m,n] = size(S);
   D = zeros(n,m);
   for i=1:n
      if (S(i,i)<eps)
         D(i,i) = 0;
      else
         D(i,i) = 1 / S(i,i);
      end    
   end

   x = (V*(D*U')) * b;
end
