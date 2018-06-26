%Initilize Variables
format long
clear; clc;
syms('X1','X2','X3','P');
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

%Initial Points
H0 = eye(4);
X0 = [0.8; 0.5; 0.3; 0.2];
step = 1;


fprintf("\n\nBFGS QUASI NEWTON METHOD \n")
fprintf("       P                X1                X2                 X3\n")
tic;
while (diff >= 0.00000001)
    
    d = - H0\eval(subs(G, {P, X1, X2, X3}, {X0(1,1), X0(2,1), X0(3,1), X0(4,1)}));
    XI = X0 + step*d;
    
    %Gradient at X0
    G1 = eval(subs(G, {P, X1, X2, X3}, {X0(1,1), X0(2,1), X0(3,1), X0(4,1)}));
    G2 = eval(subs(G, {P, X1, X2, X3}, {XI(1,1), XI(2,1), XI(3,1), XI(4,1)}));
    
    %BFGS Inputs
    S = XI - X0;
    Y = G2 - G1;
    H1 = H0 + ((Y*transpose(Y))/(transpose(Y)*S)) - (H0*S*transpose(S)*H0)/(transpose(S)*H0*S);
    
    %Stopping Criteria
    diff = abs(norm(XI)-norm(X0));

    %Re-Assign new Values for Next Iterations
    H0 = H1;
    X0 = XI;
    count = count + 1;
       
   %Display Iterate
   disp(double(transpose(XI))) 
 
   if count > 50
       break
   end
end
toc;
