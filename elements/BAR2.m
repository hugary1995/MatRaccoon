classdef BAR2 < Elem

properties
end

methods
  
  function this = initQRule(this)
    this.q_points = [Vector(-sqrt(3)/3,0), Vector(sqrt(3)/3,0)];
    this.weights = [1, 1];
  end
  
  function this = initJxW(this)
    this.JxW = 0.5*[this.nodes(2).x - this.nodes(1).x, this.nodes(2).x - this.nodes(1).x];
  end
  
  function this = initTest(this)
    a = 0.5*(1+sqrt(3)/3);
    b = 0.5*(1-sqrt(3)/3);
    this.test = [a, b; b, a];
    
    qp1 = Vector(this.test(1,1)*this.nodes(1).x+this.test(2,1)*this.nodes(2).x,...
      this.test(1,1)*this.nodes(1).y+this.test(2,1)*this.nodes(2).y);
    qp2 = Vector(this.test(1,2)*this.nodes(1).x+this.test(2,2)*this.nodes(2).x,...
      this.test(1,2)*this.nodes(1).y+this.test(2,2)*this.nodes(2).y);
    this.mapped_q_points = [qp1, qp2];
  end
  
  function this = initGradTest(this)
    dN1_qp1 = Vector(-0.5/this.JxW(1),0);
    dN1_qp2 = Vector(-0.5/this.JxW(2),0);
    dN2_qp1 = Vector(0.5/this.JxW(1),0);
    dN2_qp2 = Vector(0.5/this.JxW(2),0);
    
    this.grad_test = [dN1_qp1, dN1_qp2; dN2_qp1, dN2_qp2];
  end
  
end

end