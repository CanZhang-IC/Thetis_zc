from firedrake_adjoint import *
from firedrake import *
from pyadjoint import Block
from pyadjoint.overloaded_function import overload_function
from prepare_cable.centretoland import two_length, dlength_dx


backend_two_length = two_length

class two_lengthBlock(Block):
    def __init__(self,turbine_locations,**kwargs):
        super(two_lengthBlock,self).__init__()
        
        self.kwargs = kwargs
        self.turbine_locations = turbine_locations

        for location in self.turbine_locations:
            self.add_dependency(location)

    def __str__(self):
        return "two_lengthBlock"

    def prepare_evaluate_adj(self,inputs, adj_inputs, relevant_dependencies):
        dldx = dlength_dx(inputs[:])
        return dldx


    def evaluate_adj_component(self,inputs,adj_inputs,block_variable,idx,prepared):
        adj_input = adj_inputs[0]
        result_adj = adj_input*prepared[idx]
        return result_adj.item()

    def recompute_component(self,inputs,block_variable,idx,prepared):
        return backend_two_length(inputs[:])
    

two_length = overload_function(two_length,two_lengthBlock)



