function varargout = cons_D(gcon_inp,varargin)
% function that calculates phi, nu,gamma, dphi_dr, dphi_dp for GCon D


%%%%%%%%%%%% enter vararg conditions here %%%%%%%%%%%%%%
% extracting input variables from structure
% s1_p = gcon_inp.s1_p;
% p1 = gcon_inp.p1;
% p1dot = gcon_inp.p1dot;
% s2_q = gcon_inp.s2_q;
% p2 = gcon_inp.p2;
% p2dot = gcon_inp.p2dot;
% r1 = gcon_inp.r1;
% r1dot = gcon_inp.r1dot;
% r2 = gcon_inp.r2;
% r2dot = gcon_inp.r2dot;
% ft = gcon_inp.ft;
% dft = gcon_inp.dft;
% ddft = gcon_inp.ddft;
% i = gcon_inp.i;
% j = gcon_inp.j;

A1 = Amatrix(p1);
A2 = Amatrix(p2);

% calculating all the B matrices required
B_p2s2q = Bmatrix(gcon_inp.p2,gcon_inp.s2_q);
B_p2dots2q = Bmatrix(gcon_inp.p2dot,gcon_inp.s2_q);   
B_p1s1p = Bmatrix(gcon_inp.p1,gcon_inp.s1_p);
B_p1dots1p = Bmatrix(gcon_inp.p1dot,gcon_inp.s1_p);

d12 = (gcon_inp.r2 + A2*gcon_inp.s2_q - gcon_inp.r1 - A1*gcon_inp.s1_p);
d12dot = gcon_inp.r2dot + B_p2s2q*gcon_inp.p2dot - gcon_inp.r1dot - B_p1s1p*gcon_inp.p1dot;
phi = d12.'*d12 - gcon_inp.ft;

nu = gcon_inp.dft;

gamma = -2*d12.'*B_p2dots2q*gcon_inp.p2dot + 2*d12.'*B_p1dots1p*gcon_inp.p1dot - 2*(d12dot.'*d12dot) + gcon_inp.ddft;

ddr1 = -2*d12.'*eye(3);               % DP1 has no r terms
ddr2 = 2*d12.'*eye(3);
ddp1 = -2*d12.'*B_p1s1p;
ddp2 = 2*d12.'*B_p2s2q;

% checking if either bodies are ground
if gcon_inp.i == 0
    dphi_dr = ddr2;
    dphi_dp = ddp2;
elseif gcon_inp.j == 0
    dphi_dr = ddr1;
    dphi_dp = ddp1;
else
    dphi_dr = [ddr1,ddr2];
    dphi_dp = [ddp1,ddp2];
end

arg = 1;
if nargin > 1
    if any(strcmp('phi',varargin{1}))
        varargout{arg} = phi;
        arg = arg+1;
    end
    if any(strcmp('nu',varargin{1}))
        varargout{arg} = nu;
        arg = arg+1;
    end
    if any(strcmp('gamma',varargin{1}))
        varargout{arg} = gamma;
        arg = arg+1;
    end
    if any(strcmp('dphi_dr',varargin{1}))
        varargout{arg} = dphi_dr;
        arg = arg+1;
    end
    if any(strcmp('dphi_dp',varargin{1}))
        varargout{arg} = dphi_dp;
    end
else
    varargout = {phi,nu,gamma,dphi_dr,dphi_dp};
end