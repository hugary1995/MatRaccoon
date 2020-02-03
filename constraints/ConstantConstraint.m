classdef ConstantConstraint < EqualityConstraint
  
  properties
    var_id;
    value;
  end
  
  methods
    
    function this = ConstantConstraint(problem, variable, value)
      this = this@EqualityConstraint(problem);
      this.var_id = this.coupledValue(variable);
      this.value = value;
    end
    
    function g = computeQpConstraint(this, ivar)
      u = this.values(this.var_id);
      g = this.elem.test(this.i, this.qp)*(u(this.qp)-this.value);
    end
    
    function h = computeQpConstraintGradient(this, ivar, jvar)
      h = this.elem.test(this.i, this.qp)*this.elem.test(this.j, this.qp);
    end
    
  end
  
end