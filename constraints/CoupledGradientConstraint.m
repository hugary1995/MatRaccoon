classdef CoupledGradientConstraint < EqualityConstraint
  
  properties
    var_id;
    v_var_id;
  end
  
  methods
    
    function this = CoupledGradientConstraint(problem, variable, coupled_variable)
      this = this@EqualityConstraint(problem);
      this.var_id = this.coupledValue(variable);
      this.v_var_id = this.coupledGradient(coupled_variable);
    end
    
    function value = computeQpConstraint(this)
      u = this.values(this.var_id);
      grad_v = this.gradients(this.v_var_id);
      value = 0.5*(u(this.qp)-grad_v(this.qp).x)^2;
    end
    
    function value = computeQpConstraintGradient(this, ivar)
      u = this.values(this.var_id);
      grad_v = this.gradients(this.v_var_id);
      value = 0;
      if ivar == this.var_id
        value = this.elem.test(this.i, this.qp)*(u(this.qp)-grad_v(this.qp).x);
      elseif ivar == this.v_var_id
        value = -this.elem.grad_test(this.i, this.qp).x*(u(this.qp)-grad_v(this.qp).x);
      end
    end
    
    function value = computeQpConstraintHessian(this, ivar, jvar)
      value = 0;
      if ivar == this.var_id && jvar == this.var_id
        value = this.elem.test(this.i, this.qp)*this.elem.test(this.j, this.qp);
      elseif ivar == this.var_id && jvar == this.v_var_id
        value = -this.elem.test(this.i, this.qp)*this.elem.grad_test(this.j, this.qp).x;
      elseif ivar == this.v_var_id && jvar == this.var_id
        value = -this.elem.grad_test(this.i, this.qp).x*this.elem.test(this.j, this.qp);
      elseif ivar == this.v_var_id && jvar == this.v_var_id
        value = this.elem.grad_test(this.i, this.qp).x*this.elem.grad_test(this.j, this.qp).x;
      end
    end
    
  end
  
end