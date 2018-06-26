%Initlialize Variables
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

F = (DELTA/2)*transpose(X)*Q*X - transpose(C)*X + P*(A*X - b);
G = gradient(F);
H = hessian(F);

%Backtracking
damp = 0.0005;
adj = 0.5;
beta = 0.75;

%Multi-Variate Newtons Method   
step = 1; 
X0 = [0.1; 0.4; 0.4; 0.4];

fprintf("NEWTONS METHOD \n")
fprintf("          P                X1                X2                 X3\n")
tic;
while (diff >= 0.00000001)
    
   %Multi-Variate Newton Iterate
   G1 = eval(subs(G, {P, X1, X2, X3}, {X0(1,1), X0(2,1), X0(3,1), X0(4,1)}));
   H1 = eval(subs(H, {P, X1, X2, X3}, {X0(1,1), X0(2,1), X0(3,1), X0(4,1)}));
   
   %Direction Vector
   d = -H1\G1;
   
   %Calc Iterate
   XI = X0 + step * d; 
   
   %Function Evaluation
   F1 = eval(subs(F, {P, X1, X2, X3}, {X0(1,1), X0(2,1), X0(3,1), X0(4,1)}));
   F2 = eval(subs(F, {P, X1, X2, X3}, {XI(1,1), XI(2,1), XI(3,1), XI(4,1)}));
   
   %Backtracking
    while (F2 > F1 + damp*step*dot(d,G1))
        if(dot(d,G2) < beta*dot(d,G1))
            step = step * adj;
            XB = X0 + step*d;
            F2 = eval(subs(F, {P, X1,  X2, X3}, {XB(1,1), XB(2,1), XB(3,1), XB(4,1)})); 
        else
            break
        end
    end
  
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
