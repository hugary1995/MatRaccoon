classdef CoupledValueConstraint < EqualityConstraint
  
  properties
    var_id;
    v_var_id;
  end
  
  methods
    
    function this = CoupledValueConstraint(problem, variable, coupled_variable)
      this = this@EqualityConstraint(problem);
      this.var_id = this.coupledValue(variable);
      this.v_var_id = this.coupledValue(coupled_variable);
    end
    
    function value = computeQpConstraint(this)
      u = this.values(this.var_id);
      v = this.values(this.v_var_id);
      value = 0.5*(u(this.qp)-v(this.qp))^2;
    end
    
    function value = computeQpConstraintGradient(this, ivar)
      u = this.values(this.var_id);
      v = this.values(this.v_var_id);
      if ivar == this.var_id
        value = this.elem.test(this.i, this.qp)*(u(this.qp)-v(this.qp));
      elseif ivar == this.v_var_id
        value = -this.elem.test(this.i, this.qp)*(u(this.qp)-v(this.qp));
      end
    end
    
    function value = computeQpConstraintHessian(this, ivar, jvar)
      if ivar == this.var_id && jvar == this.var_id
        value = this.elem.test(this.i, this.qp)*this.elem.test(this.j, this.qp);
      elseif ivar == this.var_id && jvar == this.v_var_id
        value = -this.elem.test(this.i, this.qp)*this.elem.test(this.j, this.qp);
      elseif ivar == this.v_var_id && jvar == this.var_id
        value = -this.elem.test(this.i, this.qp)*this.elem.test(this.j, this.qp);
      elseif ivar == this.v_var_id && jvar == this.v_var_id
        value = this.elem.test(this.i, this.qp)*this.elem.test(this.j, this.qp);
      end
    end
    
  end
  
end