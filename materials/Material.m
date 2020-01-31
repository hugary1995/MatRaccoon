classdef (Abstract) Material < handle
  
  properties
    problem;
    data;
    
    var_ids;
    values;
    gradients;
    
    qp;
    
    elem;
  end
  
  methods
    
    function this = Material(problem)
      this.problem = problem;
      this.data = cell(1, length(problem.mesh.elems));
      
      this.values = containers.Map('KeyType', 'uint32', 'ValueType', 'any');
      this.gradients = containers.Map('KeyType', 'uint32', 'ValueType', 'any');
    end
    
    function id = coupledValue(this, variable)
      id = this.problem.getVariableId(variable);
      this.values(id) = [];
      this.var_ids = unique([this.var_ids, id]);
    end
    
    function id = coupledGradient(this, variable)
      id = this.problem.getVariableId(variable);
      this.gradients(id) = [];
      this.var_ids = unique([this.var_ids, id]);
    end
    
    function reinitElem(this, e)
      this.elem = e;
      
      for jvar = this.values.keys()
        this.values(jvar{1}) = this.reinitCoupledValue(jvar{1});
      end
      
      for jvar = this.gradients.keys()
        this.gradients(jvar{1}) = this.reinitCoupledGradient(jvar{1});
      end
    end
    
    function v = reinitCoupledValue(this, var_id)
      v = zeros(1,length(this.elem.q_points));
      for i_ = 1:length(this.elem.nodes)
        dof = this.problem.globalDoF(this.elem.nodes(i_), var_id);
        dof_value = this.problem.solution(dof);
        for qp_ = 1:length(this.elem.q_points)
          v(qp_) = v(qp_)+this.elem.test(i_, qp_)*dof_value;
        end
      end
    end
    
    function v = reinitCoupledGradient(this, var_id)
      v(1:length(this.elem.q_points)) = Vector(0, 0);
      for i_ = 1:length(this.elem.nodes)
        dof = this.problem.globalDoF(this.elem.nodes(i_), var_id);
        dof_value = this.problem.solution(dof);
        for qp_ = 1:length(this.elem.q_points)
          v(qp_) = v(qp_)+this.elem.grad_test(i_, qp_)*dof_value;
        end
      end
    end
    
    function computeMaterial(this)
      c = cell(1, length(this.elem.q_points));
      for qp_ = 1:length(this.elem.q_points)
        this.qp = qp_;
        c{qp_} = this.computeQpMaterial();
      end
      this.data{this.elem.id} = c;
    end
    
  end
  
  methods (Abstract)
    
    computeQpMaterial(this)
    
  end
  
end