classdef ConstantBC < NodalBC
  
  properties
    value;
  end
  
  methods
    
    function this = ConstantBC(problem, boundary, variable, value)
      this = this@NodalBC(problem, boundary, variable);
      this.value = value;
    end
    
    function c = computeNodalConstraint(this, n)
      c = this.value;
    end
    
  end
  
end