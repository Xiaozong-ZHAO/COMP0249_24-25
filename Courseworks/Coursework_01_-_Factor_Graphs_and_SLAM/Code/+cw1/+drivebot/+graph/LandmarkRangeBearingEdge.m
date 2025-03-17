classdef LandmarkRangeBearingEdge < g2o.core.BaseBinaryEdge
    %LANDMARKRANGEBEARINGEDGE Factor for range-bearing observations of a landmark.
    %
    % Observing the range r and bearing beta of a landmark (lx, ly) from the
    % robot pose x_{k+1} = [x, y, theta]. The measurement model is:
    %   r_pred   = sqrt((lx - x)^2 + (ly - y)^2)
    %   beta_pred= wrap( atan2(ly-y, lx-x) - theta )
    %
    % The error is: e = [r_pred - r_meas;  wrap(beta_pred - beta_meas)]
    % with angle wrapping handled by g2o.stuff.normalize_theta or wrapToPi.
    %
    % Vertex slot 1: x_{k+1} (robot pose).
    % Vertex slot 2: l_{j}   (landmark position).
    
    methods(Access = public)
        
        function obj = LandmarkRangeBearingEdge()
            % Constructor. The measurement dimension is 2: [range; bearing].
            obj = obj@g2o.core.BaseBinaryEdge(2);
        end
        
        function initialEstimate(obj)
            %INITIAL ESTIMATE of the landmark position given the robot pose and measured [r, beta].
            xPose = obj.edgeVertices{1}.estimate();    % [x, y, theta]
            rPose = xPose(1:2);                        % robot (x,y)
            thetaR= xPose(3);
            
            zMeas = obj.measurement;                   % [r, beta]
            rMeas = zMeas(1);
            bMeas = zMeas(2);
            
            % Landmark in world coords: l = robotXY + r*[cos(psi); sin(psi)],
            % where psi=thetaR + bMeas. Normalize angle to handle wrap.
            psi   = g2o.stuff.normalize_theta(thetaR + bMeas);
            lX    = rPose(1) + rMeas*cos(psi);
            lY    = rPose(2) + rMeas*sin(psi);
            
            obj.edgeVertices{2}.setEstimate([lX; lY]);
        end
        
        function computeError(obj)
            %COMPUTEERROR: e = [r_pred - r_meas; wrap(beta_pred - beta_meas)].
            xPose = obj.edgeVertices{1}.estimate();  % [x, y, theta]
            lPose = obj.edgeVertices{2}.estimate();  % [lx, ly]
            
            x     = xPose(1);
            y     = xPose(2);
            th    = xPose(3);
            lx    = lPose(1);
            ly    = lPose(2);
            
            zMeas = obj.measurement;                % [r, beta]
            rMeas = zMeas(1);
            bMeas = zMeas(2);
            
            dx    = (lx - x);
            dy    = (ly - y);
            
            rPred = sqrt(dx^2 + dy^2);
            bPred = wrapToPi(atan2(dy, dx) - th);    % or g2o.stuff.normalize_theta
            
            e_r   = rPred - rMeas;
            e_b   = g2o.stuff.normalize_theta(bPred - bMeas);
            
            obj.errorZ = [e_r; e_b];
        end
        
        function linearizeOplus(obj)
            %LINEARIZEOPLUS: compute Jacobians w.r.t. robot pose and landmark coords.
            xPose = obj.edgeVertices{1}.x; % [x, y, theta]
            lPose = obj.edgeVertices{2}.x; % [lx, ly]
            
            x   = xPose(1);
            y   = xPose(2);
            th  = xPose(3);
            lx  = lPose(1);
            ly  = lPose(2);
            
            dx  = (lx - x);
            dy  = (ly - y);
            r   = sqrt(dx^2 + dy^2);
            rr  = dx^2 + dy^2;  % r^2
            
            % Partial derivatives w.r.t. xPose = (x, y, th)
            % e_r = (rPred - rMeas),   rPred= sqrt(dx^2 + dy^2) => ∂rPred/∂x= dx/r * ∂dx/∂x => -dx/r
            d_er_dx  = -dx / r;
            d_er_dy  = -dy / r;
            d_er_dth =  0;
            
            % e_beta = (atan2(dy,dx) - th) - bMeas => partial wrt x => partial(atan2(dy,dx))/∂x + partial(-th)/∂x
            % dx=lx-x => ∂dx/∂x=-1 => so partial(atan2(dy,dx))/∂x= dy/rr
            d_eb_dx  =  dy / rr;
            d_eb_dy  = -dx / rr;
            d_eb_dth = -1;
            
            J1 = [
                d_er_dx,  d_er_dy,  d_er_dth;
                d_eb_dx,  d_eb_dy,  d_eb_dth
            ];
            
            % Partial derivatives w.r.t. lPose = (lx, ly)
            % e_r wrt lx => +dx / r,   e_r wrt ly => +dy / r
            d_er_dlx = dx / r;
            d_er_dly = dy / r;
            
            % e_beta wrt lx => we consider partial(atan2(dy,dx)) wrt dx => -dy/rr, but dx= (lx-x) => ∂dx/∂lx=1 => so total is -dy/rr
            d_eb_dlx = -dy / rr;
            d_eb_dly =  dx / rr;
            
            J2 = [
                d_er_dlx,  d_er_dly;
                d_eb_dlx,  d_eb_dly
            ];
            
            obj.J{1} = J1;
            obj.J{2} = J2;
        end
    end
end