README

PGF_equation_generator will only generate the equations. Use the function episolve to obtain ODE solutions for the time evolution of the susceptible, infected and recovered population. 

Included files: combinator.m, episolve.m, inf_neighbors.m, PGF_equation_generator, subgraphs.mat 

subgraphs.mat contains various subgraphs that can be included in the code generation algorithm. Cn and Ln denotes complete subgraphs and cycles of n nodes respectively. The subgraphs are stored as adjacency matrices and must be undirected and unweighted. 
PGF_equation_generator.m contains more detailed instructions and example calls. 

http://arxiv.org/abs/1405.6234 is an arxiv version of the paper on which this code has been based. 

A poisson random network example 

PGF_equation_generator(4, 'C2');
[S, I, R, T] = episolve();
plot(T,I,'r','linewidth',2);
xlabel('Time')
ylabel('Infectious incidence')
grid on
box on