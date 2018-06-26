%Initlialize Variables
clear
clc
format SHORT G
syms('X1','X2','X3','X4','X5','X6','P','F','F1','F2');

Q =   [  1.28462907e-03,  -4.85979426e-05,  -2.68182984e-05, -3.74681934e-05,  -2.05155057e-05,   4.42920732e-05;
         -4.85979426e-05,   6.61936638e-04,   1.11035038e-05, -1.57122599e-05,  4.99880926e-06,  -1.05281564e-05;
         -2.68182984e-05,   1.11035038e-05,   1.03117358e-04, 2.94318841e-05,   1.93242416e-05,   1.79852927e-05;
         -3.74681934e-05,  -1.57122599e-05,   2.94318841e-05, 9.15032741e-05,   2.32157567e-05,   1.16020161e-05;
         -2.05155057e-05,   4.99880926e-06,   1.93242416e-05, 2.32157567e-05,   5.36922808e-05,   9.82900409e-06;
          4.42920732e-05,  -1.05281564e-05,   1.79852927e-05, 1.16020161e-05,   9.82900409e-06,   2.69267858e-04];
C = [-0.0008727288011368683; 0.0032457521609619735; 0.0012820815841300232; 0.000687296745627091; 0.0009273986806007257; 0.0006356549214420348];
X = [X1; X2; X3; X4; X5; X6];
A = [1,1,1,1,1,1];
b = 1;
DELTA = 3.5;
diff = 1;
count = 0;

F = (DELTA/2)*transpose(X)*Q*X - transpose(C)*X + P*(A*X - b);
G = gradient(F);

%Backtracking
damp = 0.0005;
adj = 0.5;
beta = 0.75;

%Initial Points
H0 = eye(7);
X0 = [0.8; 0.2; 0.1; 0.2; 0.1; 0.2; 0.2];
step = 1;

tic;
fprintf("\n\nBFGS QUASI NEWTON METHOD \n")
fprintf("          P          BTE         NVDA          TMO         LLY         JNJ         GOLD\n")
while (diff >= 0.00000001)   
   
    d = -H0\eval(subs(G, {P, X1, X2, X3, X4, X5, X6}, {X0(1,1), X0(2,1), X0(3,1), X0(4,1), X0(5,1), X0(6,1), X0(7,1)}));
    XI = X0 + step*d;
    
    %Gradient at X0
    G1 = eval(subs(G, {P, X1, X2, X3, X4, X5, X6}, {X0(1,1), X0(2,1), X0(3,1), X0(4,1), X0(5,1), X0(6,1), X0(7,1)}));
    G2 = eval(subs(G, {P, X1, X2, X3, X4, X5, X6}, {XI(1,1), XI(2,1), XI(3,1), XI(4,1), XI(5,1), XI(6,1), XI(7,1)}));
  
    %BFGS Inputs
    S = XI - X0;
    Y = G2 - G1;
    H1 = H0 + ((Y*transpose(Y))/(transpose(Y)*S)) - (H0*S*transpose(S)*H0)/(transpose(S)*H0*S);   
      
    F1 = eval(subs(F, {P, X1, X2, X3, X4, X5, X6}, {X0(1,1), X0(2,1), X0(3,1), X0(4,1), X0(5,1), X0(6,1), X0(7,1)}));
    F2 = eval(subs(F, {P, X1, X2, X3, X4, X5, X6}, {XI(1,1), XI(2,1), XI(3,1), XI(4,1), XI(5,1), XI(6,1), XI(7,1)}));

    %Backtracking
    if (F2 > F1 + damp*step*dot(d,G1) )
        if(dot(d,G2) < beta*dot(d,G1))
            step = step * adj;
            XB = X0 + step*d;
            F2 = eval(subs(F, {P, X1, X2, X3, X4, X5, X6}, {XB(1,1), XB(2,1), XB(3,1), XB(4,1), XB(5,1), XB(6,1), XB(7,1)}));
        else
            break
        end
    end

    %Stopping Criteria
    diff = abs(norm(XI)-norm(X0));

    %Re-Assign new Values for Next Iterations
    H0 = H1;
    X0 = XI;
    count = count + 1;
       
   %Display Iterate
   disp(count)
   disp( round(double(transpose(XI)),5) ) 
  
   
   if count > 150
       break
   end
end
toc;