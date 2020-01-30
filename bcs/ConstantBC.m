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
      c = this.u-this.value;
    end
    
    function gc = computeNodalConstraintGradient(this, n)
      gc = 1;
    end
    
  end
  
end