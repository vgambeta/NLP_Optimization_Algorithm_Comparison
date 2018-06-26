%Initlialize Variables
clear
% clc
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

% Steepest Descent Method
step = 1;
X0 = [0.8; 0.5; 0.3; 0.2];

fprintf("\n\nSTEEPEST DESCENT METHOD \n")
fprintf("          P                X1                X2                 X3\n")
while (diff >= 0.00000001)
    
   %Steepest Descent Iterate
   d = -eval(subs(G, {P, X1, X2, X3}, {X0(1,1), X0(2,1), X0(3,1), X0(4,1)}));
   XI = X0 + step * d; 
   
   %Display Iterate
   disp(double(transpose(XI)))
   
   %Stopping Criteria
   diff = abs(norm(XI)-norm(X0));
   
   %Prepare Iterate
   X0 = XI;
   count = count + 1;
   
   if count > 25
       break
   end
end

