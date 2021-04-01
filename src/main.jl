using GL2017Replication

p = ModelParameters()
grid = ModelGrid()
policy = ModelSolutions()
# Run all
solve_steady_state!(p,grid,policy)