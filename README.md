# NLP_Optimization_Algorithm_Comparison
Optimization reformulation using a Langrangian and implementation of three different optimization techniques. 

Non-linear programming is the study of solving non-linear problems. Solutions are estimated based on various different search algorithms. Three methods are explord here; Newton's Method, Steepest Descent Method (SDM) and BFGS - Quasi Newton Method.

Backtracking is implemented for each method which seeks to find an improvement to the solvers solution time by determining an optimal step length for the search algorithm at eash step.  

This computational project set out to compare these three methods of optimization. It was determined that the newtons method, superior due to the quadratic objective function, functioned well with and without the additional backtracking line search for determining step length. The steepest descent method, although reliable under certain circumstances does not converge with given objective function even with the backtracking to determine the optimal step length.

These files contain alogrithms which are the basic foundational implementation of a solution method used to solve a non-linear problem subject to only a single linear equality constraint. The reformualation and basic optimization cannot be conudcted for a problem that includes ineuqality constraints.
