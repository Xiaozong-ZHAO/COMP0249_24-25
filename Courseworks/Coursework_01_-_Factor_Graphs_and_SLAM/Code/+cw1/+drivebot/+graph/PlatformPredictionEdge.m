classdef PlatformPredictionEdge < g2o.core.BaseBinaryEdge
    %PLATFORMPREDICTIONEDGE Factor representing the process model from x_k to x_{k+1}.
    %
    % The state update is:
    %   x_{k+1} = x_k + M * [v_x; v_y; v_theta]
    % where
    %   M = dT * [cos(psi)  -sin(psi)  0
    %             sin(psi)   cos(psi)  0
    %             0          0         1]
    %
    % The 'measurement' in this context is actually the mean control input
    % (v_x, v_y, v_theta). The error is:
    %   e = M^{-1} * (x_{k+1} - x_k) - u,
    % so that if the prediction is perfect, e ~ 0.
    %
    % Vertex slot 1: x_k (the previous platform state).
    % Vertex slot 2: x_{k+1} (the next platform state).
    
    properties(Access = protected)
        dT  % The time step between consecutive states.
    end

    methods(Access = public)
        function obj = PlatformPredictionEdge(dT)
            % Constructor. dT is the nominal time step, which may be
            % overridden if the vertices have time stamps.
            assert(dT >= 0, 'Time step must be non-negative');
            obj = obj@g2o.core.BaseBinaryEdge(3);
            obj.dT = dT;
        end
        
        function initialEstimate(obj)
            %INITIAL ESTIMATE of x_{k+1} given x_k and control input u_k.
            x_k = obj.edgeVertices{1}.x;       % [x; y; theta]
            u_k = obj.measurement;            % control input
            dt  = obj.getDeltaT();            % might override with vertex time fields
            
            psi = x_k(3);
            M = dt * [cos(psi), -sin(psi), 0;
                      sin(psi),  cos(psi), 0;
                      0,         0,        1];
            
            x_next = x_k + M * u_k;           % predict
            x_next(3) = g2o.stuff.normalize_theta(x_next(3));
            obj.edgeVertices{2}.setEstimate(x_next);
        end
        
        function computeError(obj)
            %COMPUTEERROR: error e = M^{-1}( x_{k+1} - x_k ) - u_k.
            x_k   = obj.edgeVertices{1}.x;
            x_k1  = obj.edgeVertices{2}.x;
            u_k   = obj.z;        % 'measurement' is the control input
            dt    = obj.getDeltaT();
            
            psi   = x_k(3);
            
            % M^-1 = (1/dt)*[ cos(psi)   sin(psi)  0
            %                 -sin(psi)  cos(psi)  0
            %                  0         0         1 ]
            % i.e. the transpose of R(psi) scaled by 1/dt.
            M_inv = (1/dt) * [ cos(psi),  sin(psi),  0;
                              -sin(psi),  cos(psi),  0;
                               0,         0,         1];
            
            obj.errorZ = M_inv * (x_k1 - x_k) - u_k;
            obj.errorZ(3) = g2o.stuff.normalize_theta(obj.errorZ(3));
        end
        
        function linearizeOplus(obj)
            %LINEARIZEOPLUS: compute Jacobians wrt x_k (vertex1) and x_{k+1} (vertex2).
            x_k   = obj.edgeVertices{1}.x;
            x_k1  = obj.edgeVertices{2}.x;
            dt    = obj.getDeltaT();
            
            psi = x_k(3);
            
            % A = M_inv = (1/dt)*R(psi)^T
            R_trans = [ cos(psi),  sin(psi), 0;
                       -sin(psi),  cos(psi), 0;
                        0,         0,        1];
            A = (1/dt) * R_trans;
            
            % Residual e = A*(x_k1 - x_k) - u => wrt x_{k+1} => partial e / partial x_{k+1} = A
            J2 = A;
            
            % Wrt x_k => partial e / partial x_k = -A + dA/d(psi)* (x_k1 - x_k)
            dx = (x_k1 - x_k);
            
            % derivative of R_trans wrt psi
            dR_dpsi = [ -sin(psi), cos(psi), 0;
                        -cos(psi), -sin(psi),0;
                         0,        0,        0];
            dA_dpsi = (1/dt) * dR_dpsi;
            
            J1 = -A; 
            J1(:,3) = J1(:,3) + dA_dpsi*dx;    % chain rule for d( A(psi)*dx )/d(psi)
            
            obj.J{1} = J1;
            obj.J{2} = J2;
        end
    end
    
    methods (Access = private)
        function dt = getDeltaT(obj)
            %GETDELTAT tries to read time fields in vertices, otherwise uses obj.dT.
            v1 = obj.edgeVertices{1};
            v2 = obj.edgeVertices{2};
            
            dt = obj.dT;
            if isfield(v1, 'time') && isfield(v2, 'time')
                dtOverride = abs(v2.time - v1.time);
                if dtOverride > 0
                    dt = dtOverride;
                end
            end
            % If dt=0, fallback to a tiny positive to avoid division by zero
            dt = max(dt, 1e-9);
        end
    end
end
