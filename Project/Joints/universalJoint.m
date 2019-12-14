%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Position, Velocity and Acceleration of Universal Joint
% Input: input all the attributes in the 'input.txt' file
% Output: Values for Position, velocity and Acceleration
% Universal Joint has 3 Co-ordinate Difference Constraints and one DP1
% Constraints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
warning off;
 disp('The two bodies are connected by a Universal Joint.')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reads the input file from input.txt and saves them as variables in
% workspace

A1 = @(th)[0,  0,      1;                % Rotation matrix as a function of angle
    sin(th),  cos(th), 0;
    -cos(th), sin(th), 0];

theta = @(t)pi/4*cos(2*t);                       % theta, followed by its time derivatives
dtheta = @(t)-pi/2*sin(2*t);
ddtheta = @(t)-pi*cos(2*t);

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
    
    % calling 1 DP1 GCon to add to SJ to form Universal joint
    % defining such that x'-y' plane is perpendicular to X-axis
    a_2 = [1,0,0].';                                % x-axis in G-RF
    a_1 = [1,0,0].';                                % x-axis in L-RF
    b_1 = [0,1,0].';                                % y-axis in L-RF
    
    gcon_inp.a_2 = a_2;
    gcon_inp.a_1 = a_1;
    w = w+1;
    [phi(w,1),nu(w,1),gamma(w,1),dphi_dr(w,:),dphi_dp(w,:)] = cons_DP1(gcon_inp);
    w = w+1;
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
    
end
output = {Phi_q, q,qbar, qdot q_dot qddot q_ddot ts};
disp('The results have been saved as Output')
