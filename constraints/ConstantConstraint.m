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
    
    function value = computeQpConstraint(this)
      u = this.values(this.var_id);
      value = 0.5*(u(this.qp)-this.value)^2;
    end
    
    function value = computeQpConstraintGradient(this, ivar)
      u = this.values(this.var_id);
      value = this.elem.test(this.i, this.qp)*(u(this.qp)-this.value);
    end
    
    function value = computeQpConstraintHessian(this, ivar, jvar)
      value = this.elem.test(this.i, this.qp)*this.elem.test(this.j, this.qp);
    end
    
  end
  
end