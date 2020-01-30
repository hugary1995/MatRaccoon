classdef (Abstract) Elem < handle
  
  properties
    id;
    nodes;
    q_points;
    mapped_q_points;
    weights;
    JxW;
    test;
    grad_test;
  end
  
  methods
    
    function this = Elem(nodes)
      persistent counter_elem;
      if isempty(counter_elem)
        counter_elem = 1;
      else
        counter_elem = counter_elem+1;
      end
      this.id = counter_elem;
      this.nodes = nodes;
      this.initQRule();
      this.initJxW();
      this.initTest();
      this.initGradTest();
    end
    
  end
  
  methods (Abstract)
    
    initQRule(this)
    
    initJxW(this)
    
    initTest(this)
    
    initGradTest(this)
    
  end
  
end