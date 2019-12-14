%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Position, Velocity and Acceleration analysis of a single Pendulum with Revolute
% Joint with ground. The motion is simulated along with its position,
% velocity and acceleration. There's an option to save the video of the
% motion.
% Input: input all the attributes in the 'input.txt' file
% Output: Plot for position, velocity and acceleration, option to save
% video
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 disp('The two bodies are connected by a Revolute Joint.')

A1 = @(th)[0,  0,      1;                        % Rotation matrix as a function of angle
    sin(th),  cos(th), 0;
    -cos(th), sin(th), 0];

theta = @(t)pi/4*cos(2*t);                       % theta, followed by its time derivatives
dtheta = @(t)-pi/2*sin(2*t);
ddtheta = @(t)-pi*cos(2*t);
 if strcmp(dynamicAnalysis, 'Yes')
        disp('Running an Inverse Dynamic Analysis')
 end
% variables derived from above input
for n = 1:length(ts)
    t = ts(n);
    omega = [dtheta(t),0,0].';                      % angular velocity obtained from f(t) = Acos(omega*t)
    theta0 = theta(t);
    r1 = A1(theta0)*[L,0,0].';
    r2 = zeros(3,1);
    ft = cos(theta(t) + pi/2);                       % function for driving constraint, and its derivatives
    dft =-sin(theta(t) + pi/2)*dtheta(t);
    ddft = -sin(theta(t) + pi/2)*ddtheta(t) - cos(theta(t) + pi/2)*dtheta(t)^2;
    
    % calculating p and pdot from input
    p1 = getp(A1(theta0));                           % calculate p from A
    unitycheck = p1(1)^2 + p1(2:4).'*p1(2:4);        % verifying that the A.'*A=eye(3) condition is met
    e1 = p1(2:4);
    E1 = [-e1, tilda(e1)+p1(1)*eye(3)];
    
    p1dot = 0.5*E1.'*omega;                          % using omega and E to get pdot
    orthocheck = p1dot.'*p1;
    
    p2 = getp(eye(3));                               % euler paramters for G-RF
    e2 = p2(2:4);
    p2dot = zeros(4,1);
    
    % calling 3 CD GCon to form the spherical joint (SJ)
    gcon_inp.p1 = p1;
    gcon_inp.p1dot = p1dot;
    gcon_inp.p2 = p2;
    gcon_inp.p2dot = p2dot;
    gcon_inp.ft = 0;
    gcon_inp.dft = 0;
    gcon_inp.ddft = 0;
    gcon_inp.r1 = r1;
    gcon_inp.r2 = r2;
    gcon_inp.s1_p = s1_p;
    gcon_inp.s2_q = s2_q;
    gcon_inp.i = i;
    gcon_inp.j = j;
    
    for w = 1:3
        c = zeros(3,1);
        c(w) = 1;
        gcon_inp.c = c;
        [phi(w,1),nu(w,1),gamma(w,1),dphi_dr(w,:),dphi_dp(w,:)] = cons_CD(gcon_inp);
    end
    
    
    % calling 2 DP1 GCon to add to SJ to form Revolute joint
    % defining such that x'-y' plane is perpendicular to X-axis
    a_2 = [1,0,0].';                                % x-axis in G-RF
    a_1 = [1,0,0].';                                % x-axis in L-RF
    b_1 = [0,1,0].';                                % y-axis in L-RF
    
    gcon_inp.a_2 = a_2;
    gcon_inp.a_1 = a_1;
    w = w+1;
    [phi(w,1),nu(w,1),gamma(w,1),dphi_dr(w,:),dphi_dp(w,:)] = cons_DP1(gcon_inp);
    gcon_inp.a_1 = b_1;
    w = w+1;
    [phi(w,1),nu(w,1),gamma(w,1),dphi_dr(w,:),dphi_dp(w,:)] = cons_DP1(gcon_inp);
    
    
    % calling another DP1 for the driving constraint
    a_1 = [0,1,0].';
    a_2 = [0,0,-1].';
    gcon_inp.a_2 = a_2;
    gcon_inp.a_1 = a_1;
    gcon_inp.ft = ft;
    gcon_inp.dft = dft;
    gcon_inp.ddft = ddft;
    w = w+1;
    [phi(w,1),nu(w,1),gamma(w,1),dphi_dr(w,:),dphi_dp(w,:)] = cons_DP1(gcon_inp);
    
    % euler parameter normalization constraint
    w = w+1;
    [phi(w,1),nu(w,1),gamma(w,1),dphi_dr(w,:),dphi_dp(w,:)] = getpnorm(p1,p1dot);
    
    % assembling the Jacobian matrix by combining dphi_dr and dphi_dp
    Phi_q = [dphi_dr,dphi_dp];
    q(:,n) = [r1;p1];
    qbar(:,n) = A1(theta0)\q(1:3,n);      % in L-RF
    qdot(:,n) = Phi_q\nu;
    q_dot(:,n) = A1(theta0)\qdot(1:3,n); % in L-RF
    qddot(:,n) = Phi_q\gamma;
    q_ddot(:,n) = A1(theta0)\qddot(1:3,n);% in L-RF
    
    % Inverse Dynamic Analysis
    if strcmp(dynamicAnalysis, 'Yes')
        
    rddot = qddot(1:3,n);
    pddot = qddot(4:end,n);
    nb = 1;
    Mmatrix  = M*eye(3);
    Jbar = diag(J)*eye(3);
    Gmatrix = [p1(2:4),tilda(p1)+p1(1)*eye(3)];
    Gdot = [p1dot(2:4),tilda(p1dot)+p1dot(1)*eye(3)];
    nBar = zeros(3,1);
    tauHat = 2*Gmatrix'*nBar + 8*Gdot'*Jbar*Gdot*p1;
    Jp = 4*Gmatrix'*Jbar*Gmatrix;
    Fvector = [0 0 M*9.81]';
    
    % Left hand side (LHS) matrix :
    %   [phi_r' zeros(3*nb,nb);
    %    phi_p'        P1]
    
    LHS = [dphi_dr(w,:)' zeros(3,1)
        dphi_dp(w,:)', p1];
    
    % RHS matrix :
    %       -[  M*rddot - F;
    %         J^p*pddot - tauHat]
    
    RHS = -[Mmatrix*rddot - Fvector;
        Jp*pddot - tauHat];
    
    % lagrange multipliers Lambdas
    lambda = LHS\RHS;
    lambdaVector = lambda(1)*ones(7,1);
    % for reaction forces
    phi_r_i = dphi_dr(:,(3*1-2):3*1);
    rForce = -phi_r_i'*lambdaVector;
    
    % for reaction torque
    phi_p_i = dphi_dp(:,(4*1-3):4*1);
    rTorque = phi_p_i'*lambdaVector;
    
    % this reaction torque is in r-p formulation, we need to convert this to
    % r-w formulation
    magnitude = sqrt(rTorque'*rTorque);
    pf = rTorque/magnitude;
    torque = pf(2)/sqrt(1-pf(1)^2)*magnitude;
    Torque(n) = torque;
    end
end
output = {Phi_q, q,qbar, qdot q_dot qddot q_ddot ts};

if strcmp(dynamicAnalysis,'Yes')
    figure(2)
    plot(ts,Torque)
    grid on
    title('Reaction torque on a pendulum')
    xlabel('Time(s)')
    ylabel('Torque')
    box on
    print('-dpng','-painters','ReactionTorque.png');
 end
