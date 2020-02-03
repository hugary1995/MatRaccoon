classdef FunctionMaterial < Material
  
  properties
    fun;
  end
  
  methods
    
    function this = FunctionMaterial(problem, fun)
      this = this@Material(problem);
      this.fun = fun;
    end
    
    function c = computeQpMaterial(this)
      p = this.elem.mapped_q_points(this.qp);
      c = this.fun(p.x, p.y);
    end
    
  end
  
end