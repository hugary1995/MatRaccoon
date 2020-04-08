classdef EqualGradientConstraint < EqualityConstraint
  
  properties
    var_id;
    v_var_id;
  end
  
  methods
    
    function this = EqualGradientConstraint(problem, variable, coupled_variable)
      this = this@EqualityConstraint(problem);
      this.var_id = this.coupledGradient(variable);
      this.v_var_id = this.coupledGradient(coupled_variable);
    end
    
    function value = computeQpConstraint(this)
      grad_u = this.gradients(this.var_id);
      grad_v = this.gradients(this.v_var_id);
      value = 0.5*(grad_u(this.qp).x-grad_v(this.qp).x)^2;
    end
    
    function value = computeQpConstraintGradient(this, ivar)
      grad_u = this.gradients(this.var_id);
      grad_v = this.gradients(this.v_var_id);
      value = 0;
      if ivar == this.var_id
        value = this.elem.grad_test(this.i, this.qp).x*(grad_u(this.qp).x-grad_v(this.qp).x);
      elseif ivar == this.v_var_id
        value = -this.elem.grad_test(this.i, this.qp).x*(grad_u(this.qp).x-grad_v(this.qp).x);
      end
    end
    
    function value = computeQpConstraintHessian(this, ivar, jvar)
      value = 0;
      if ivar == this.var_id && jvar == this.var_id
        value = this.elem.grad_test(this.i, this.qp).x*this.elem.grad_test(this.j, this.qp).x;
      elseif ivar == this.var_id && jvar == this.v_var_id
        value = -this.elem.grad_test(this.i, this.qp).x*this.elem.grad_test(this.j, this.qp).x;
      elseif ivar == this.v_var_id && jvar == this.var_id
        value = -this.elem.grad_test(this.i, this.qp).x*this.elem.grad_test(this.j, this.qp).x;
      elseif ivar == this.v_var_id && jvar == this.v_var_id
        value = this.elem.grad_test(this.i, this.qp).x*this.elem.grad_test(this.j, this.qp).x;
      end
    end
    
  end
  
end