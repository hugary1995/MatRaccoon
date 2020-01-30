classdef POINT1 < Elem

properties
end

methods
  
  function this = initQRule(this)
    this.q_points = Vector(0, 0);
    this.weights = 1;
  end
  
  function this = initJxW(this)
    this.JxW = 1;
  end
  
  function this = initTest(this)
    this.test = 1;
    this.mapped_q_points = Vector(this.nodes.x, this.nodes.y);
  end
  
  function this = initGradTest(this)
    this.grad_test = Vector(0, 0);
  end
  
end

end