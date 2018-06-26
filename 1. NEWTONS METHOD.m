%Initlialize Variables
clear
clc
format long
syms('X1','X2','X3','P');

%Formulated Function
Q = [0.02778,0.00387,0.00021; 0.00387,0.01112,-0.00020; 0.00021,-0.00020,0.00115];
C = [0.1073; 0.0737; 0.0627];
X = [X1; X2; X3];
A = [1,1,1];
b = 1;
DELTA = 3.5;
diff = 1;
count = 0;

F = (DELTA/2)*transpose(X)*Q*X - transpose(C)*X + P*(A*X-b);
G = gradient(F);
H = hessian(F);
%F = P*(X1 + X2 + X3 - 1) - C(1)*X1 - C(2)*X2 - C(3)*X3 + (DELTA/2)*( Q(1,1)*X1^2 + Q(2,2)*X2^2 + Q(3,3)*X3^2 ) + DELTA*Q(1,2)*X1*X2 + DELTA*Q(1,3)*X1*X3 + DELTA*Q(2,3)*X2*X3;

%Multi-Variate Newtons Method   
step = 1; 
X0 = [0.8; 0.5; 0.3; 0.2];

fprintf("NEWTONS METHOD \n")
fprintf("          P                X1                X2                 X3\n")
tic;
while (diff >= 0.00000001)
   
   %Multi-Variate Newton Iterate
   G1 = eval(subs(G, {P, X1, X2, X3}, {X0(1,1), X0(2,1), X0(3,1), X0(4,1)}));
   H1 = eval(subs(H, {P, X1, X2, X3}, {X0(1,1), X0(2,1), X0(3,1), X0(4,1)}));
   XI = X0 - step * H1\G1;
   
   %Display Iterate
   disp(double(transpose(XI)))
   
   %Stopping Criteria
   diff = abs(norm(XI)-norm(X0));
   
   %Prepare Iterate 
   X0 = XI;
   count = count + 1; 

   if count > 100
       break
   end
end
toc;