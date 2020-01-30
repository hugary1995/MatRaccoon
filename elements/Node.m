classdef Node < Vector
  
  properties
    id;
  end
  
  methods
    
    function this = Node(x,y)
      persistent counter_node;
      if isempty(counter_node)
        counter_node = 1;
      else
        counter_node = counter_node+1;
      end
      this = this@Vector(x,y);
      this.id = counter_node;
    end
    
  end
  
end