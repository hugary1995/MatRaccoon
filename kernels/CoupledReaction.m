classdef CoupledReaction < Kernel
  
  properties
    u_var_id;
    v_var_id;
  end
  
  methods
    
    function this = CoupledReaction(problem, u, v)
      this = this@Kernel(problem);
      this.u_var_id = this.coupledValue(u);
      this.v_var_id = this.coupledValue(v);
    end
    
    function f = computeQpObjective(this)
      u = this.values(this.u_var_id);
      v = this.values(this.v_var_id);
      f = u(this.qp)*v(this.qp);
    end
    
    function g = computeQpGradient(this, ivar)
      if ivar == this.u_var_id
        value = this.values(this.v_var_id);
      elseif ivar == this.v_var_id
        value = this.values(this.u_var_id);
      end
      g = this.elem.test(this.i, this.qp)*value(this.qp);
    end
    
    function h = computeQpHessian(this, ivar, jvar)
      if ivar == jvar
        h = 0;
      else
        h = this.elem.test(this.i, this.qp)*this.elem.test(this.j, this.qp);
      end
    end
  end
  
end