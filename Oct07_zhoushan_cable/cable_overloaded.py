from firedrake_adjoint import *
from firedrake import *
from pyadjoint import Block
from pyadjoint.overloaded_function import overload_function
from order_in_length_out import cablelength, dlength_dx
import Hybrid_Code


backend_cablelength = cablelength

class CablelengthBlock(Block):
    def __init__(self,turbine_locations,substation_location,order_w,**kwargs):
        super(CablelengthBlock,self).__init__()
        
        self.kwargs = kwargs
        self.turbine_locations = turbine_locations
        self.substation_location = substation_location
        self.order = order_w
        for location in self.turbine_locations:
            self.add_dependency(location)
        # for sublocation in self.substation_location:
        #     self.add_dependency(sublocation)
        # for order in order_w:
        #     self.add_dependency(order)

    def __str__(self):
        return "CablelengthBlock"

    def prepare_evaluate_adj(self,inputs, adj_inputs, relevant_dependencies):
        ####Document###
        # the first two derivates are for landlocation, 
        # which is not considered as the dependencies.
        # However, in the Jacobian, its derivates are considered,
        # so here the first two values are ignored .
        turbine_index = len(self.turbine_locations)
        #substation_index = len(self.turbine_locations)+len(self.substation_location)
        dldx = dlength_dx(inputs[:],self.substation_location,self.order)[2:]
        return dldx


    def evaluate_adj_component(self,inputs,adj_inputs,block_variable,idx,prepared):
        adj_input = adj_inputs[0]
        result_adj = adj_input*prepared[idx]
        return result_adj.item()

    def prepare_recompute_component(self, inputs, relevant_outputs):
        t_location_float = [] 
        for i in inputs[:]:
            t_location_float.append(float(i.values()))
        l_location_float = [float(self.substation_location[0].values()),float(self.substation_location[1].values())]
        # print(t_location_float)
        # print(l_location_float)
        cableclass = Hybrid_Code.CableCostGA(t_location_float,l_location_float)
        order_w = cableclass.compute_cable_cost_order()
        order_con = [i for j in order_w for i in j]
        # print(order_con)
        for i,order in enumerate(self.order):
            order.assign(order_con[i])
        relevant_outputs = self.order
        return relevant_outputs
         

    def recompute_component(self,inputs,block_variable,idx,prepared):
        turbine_index = len(self.turbine_locations)
        # substation_index = len(self.turbine_locations)+len(self.substation_location)
        return backend_cablelength(inputs[:],self.substation_location,prepared)
    

cablelength = overload_function(cablelength,CablelengthBlock)



