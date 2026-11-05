function fig = DC_Motor_Speed_Control_GUI()
%DC_MOTOR_SPEED_CONTROL_GUI  Interactive P / PI / PD / PID speed control of a DC motor.
%
%   Run with:   DC_Motor_Speed_Control_GUI
%
%   Model (linearized, constant field):
%       Va = La*di/dt + Ra*i + Ke*w        J*dw/dt + B*w = Kt*i        dtheta/dt = w
%       D(s)     = (La*s + Ra)*(J*s + B) + Kt*Ke
%       P_i(s)   = (J*s + B)/D(s)          P_w(s) = Kt/D(s)          P_theta(s) = Kt/(s*D(s))
%
%   Unity negative feedback on SPEED, error e = w_ref - w, controller output = Va:
%       X(s)/Wref(s) = P_X(s)*C(s) / (1 + P_w(s)*C(s))
%
%   Writing C = Nc/Dc and clearing fractions, the common factor D(s) cancels
%   EXACTLY (done analytically, which is numerically safer than minreal):
%       Dcl(s)        = D(s)*Dc(s) + Kt*Nc(s)        (closed-loop characteristic polynomial)
%       I(s)/Wref(s)  = (J*s + B)*Nc(s) / Dcl(s)
%       W(s)/Wref(s)  =  Kt*Nc(s)       / Dcl(s)
%       Th(s)/Wref(s) =  Kt*Nc(s)       / (s*Dcl(s))
%   These equal P_X*C/(1 + P_w*C) identically.
%
%   MODEL PREDICTIVE CONTROL (MPC) - added controller option "MPC"
%   Same linearized motor, state x = [ia; w; theta], input u = Va:
%       dx/dt = A x + Bc u,  A = [-Ra/La -Ke/La 0; Kt/J -B/J 0; 0 1 0],  Bc = [1/La; 0; 0]
%   Exact zero-order-hold discretization with sampling interval Ts (computed from the
%   CURRENT motor parameters via expm):   x(k+1) = Ad x(k) + Bd u(k).
%   At every sampling instant k (receding horizon):
%       1. take the current simulated state x(k)
%       2. predict  w(k+j|k) = C Ad^j x(k) + sum_m C Ad^(j-1-m) Bd u(k+m),  j = 1..Np
%          (inputs after the control horizon Nc are held at u(k+Nc-1))
%       3. minimise  J = Q*sum_{j=1..Np}(w(k+j|k) - wref)^2 + R*sum_{j=0..Nc-1} u(k+j|k)^2
%          -> QP in U = [u(k)..u(k+Nc-1)]:  H = Q*G'*G + R*I,  f = Q*G'*(Phi*x(k) - wref)
%       4. apply ONLY u(k) = first element of the optimum, advance the plant one Ts
%       5. repeat at k+1
%   Default parameters: Q_MPC = 1, R_MPC = 1, Np = 20, Nc = 10 (all four editable in section E),
%   Ts = 5 ms (reduced automatically only if needed to keep 1001 samples for short runs).
%   By default MPC applies NO constraints (the classical controllers are not limited
%   either, so the comparison is like-for-like). The two tick boxes in the "Voltage /
%   current limits" panel (right column) switch on ENFORCEMENT inside the MPC
%   optimisation: "MPC enforces V limits" adds  Vmin <= u(k+j) <= Vmax  and "MPC enforces
%   I limit" adds  -Imax <= ia(k+j|k) <= Imax, j = 1..Np (predicted at the sampling
%   instants). Either one makes the QP constrained and then requires quadprog
%   (Optimization Toolbox). Without enforcement the QP is solved with a Cholesky
%   factorisation of H (base MATLAB). MPC never silently falls back to another
%   controller: any failure (including infeasible limits) is reported in the Status panel.
%
%   CONTROL MODE (radio buttons "Speed Control" / "Position Control")
%   Speed mode (original): reference wref [rad/s], feedback = w, e = wref - w, output Va.
%   Position mode (added): reference thref [rad], feedback = THETA, e = thref - theta, output Va.
%   Direct position feedback (no inner speed loop). With C = Nc/Dc the exact closed loop is
%       Dcl(s)        = s*D(s)*Dc(s) + Kt*Nc(s)          (extra integrator of P_theta)
%       Th(s)/Thref(s) = Kt*Nc(s) / Dcl(s)
%       W(s)/Thref(s)  = s*Kt*Nc(s) / Dcl(s)
%       I(s)/Thref(s)  = s*(J*s + B)*Nc(s) / Dcl(s)
%   which equal P_X*C/(1 + P_theta*C). Stability = all roots of Dcl strictly in the left
%   half plane (for P control this needs Kp < about 2.48 with the default motor).
%   As in speed mode the derivative term acts on the full error, so a reference step gives
%   PD/PID an initial current jump of about Kd*thref/La (ideal derivative, no filtering).
%   MPC in position mode uses the SAME code: output matrix C = [0 0 1] and the cost
%       J = Q*sum_j (theta(k+j|k) - thref)^2 + R*sum_j u(k+j|k)^2.
%   The GUI fields Q, R, Np keep their values when the mode changes; Q weights the speed
%   error in speed mode and the position error in position mode (the label says which).
%   Initial conditions in both modes: ia(0) = 0, w(0) = 0, theta(0) = 0.
%   Position-mode defaults (illustrative, chosen to be stable with the default motor):
%   thref = 1 rad; Kp/Ki/Kd = P 1/0/0, PI 1/2/0, PD 1/0/0.005, PID 1/2/0.005.
%
%   KALMAN FILTER (KF) + MPC  (added; everything above is unchanged)
%   Two selectable MPC feedback configurations (section F; they apply to the MPC controller):
%     Mode A  "MPC using Direct States"           : MPC receives the ACTUAL simulated state x(k).
%     Mode B  "MPC using Kalman Filter Estimates" : MPC receives the KF estimate xhat(k|k).
%   Naming used in code, plots and GUI:
%     x    ACTUAL state of the plant  [ia; w; theta]   (the plant ALWAYS evolves with this state)
%     z    NOISY MEASUREMENT          z = C x + v,  C = [1 0 0; 0 1 0]  (ia and w are measured,
%                                     theta is NOT measured; the reference is never a measurement).
%                                     OPTIONAL tick box "Position sensor": C = [1 0 0; 0 1 0; 0 0 1].
%     NOTE: with the default C the pair (Ad,C) is NOT observable - theta does not influence ia or w -
%     so xhat_theta is dead-reckoned from the w estimate: it keeps the initial error e0_theta and drifts
%     (random walk). That is harmless in Speed Control (MPC does not use theta) but in Position Control
%     the MPC then regulates the ESTIMATED position; tick "Position sensor" to make theta observable.
%     xh   ESTIMATE of the KF         (never written back into the plant)
%   Plant used in both modes (identical, same ZOH matrices Ad, Bd as MPC and KF):
%       x(k+1) = Ad x(k) + Bd u(k) + w(k),    w(k) ~ N(0, Qtrue),  Qtrue = diag(sigq.^2)
%       z(k)   = C x(k) + v(k),               v(k) ~ N(0, Rtrue),  Rtrue = diag(sigr.^2)
%   Closed-loop order at every sampling instant k (KF mode):
%       1. z(k) = C x(k) + v(k)                               (noisy sensor reading of the ACTUAL state)
%       2. predict   xh^-(k) = Ad xh(k-1) + Bd u(k-1),   P^-(k) = Ad P(k-1) Ad' + Q_KF
%          (k = 0: xh^-(0) = x(0) + e0 and P^-(0) = P0 are the prior, e0 = initial estimate error)
%       3. correct   K = P^- C' / (C P^- C' + R_KF),  xh = xh^- + K (z - C xh^-),
%                    P = (I-KC) P^- (I-KC)' + K R_KF K'   (Joseph form = (I-KC)P^- when K is optimal,
%                    but keeps P symmetric / positive semidefinite; P is also symmetrised and any
%                    tiny negative eigenvalue from round-off is clamped to 0)
%       4. MPC computes u(k) from xh(k) (Mode B) or from x(k) (Mode A) - same QP, same Q_MPC, R_MPC,
%          Np, Nc, Ts, constraints
%       5. x(k+1) = Ad x(k) + Bd u(k) + w(k)  (ACTUAL plant);  u(k) is stored and used in step 2 at k+1.
%   Unit conventions (all GUI noise values are STANDARD DEVIATIONS, covariances are their squares):
%     Section G (noise really injected into the simulated MPC loop; also used for BOTH modes A and B
%     for the process noise, so they are compared under identical disturbances):
%         sigr = [sigma_ia [A], sigma_w [rad/s], sigma_theta [rad]] measurement-noise std. dev.
%                (sigma_theta only used if the position sensor is ticked) -> Rtrue = diag(sigr.^2)
%         sigq = [sigma_ia [A], sigma_w [rad/s], sigma_theta [rad]] PER-SAMPLE process-noise std. dev.
%                                                                    -> Qtrue = diag(sigq.^2)
%         Level menus None / Low / Medium / High / Custom fill these fields with presets (editable).
%         Random numbers: one private RandStream(mt19937ar, seed); the unit-normal sequences are drawn
%         once (independent of the controller configuration) and then scaled by the std. devs., so
%         Mode A and Mode B, and different noise levels, see the SAME underlying noise realisation.
%     Section H (what the Kalman Filter ASSUMES):  Q_KF = qScale * Qtrue,  R_KF = rScale * Rtrue
%         (variances; scale 1 = filter knows the true noise, scale <> 1 = deliberately mismatched
%         filter). Q_KF / R_KF are floored at 1e-12 so they stay positive definite when a noise
%         level is 'None'. P0 = diag([P0_ia [A^2], P0_w [(rad/s)^2], P0_theta [rad^2]]) is the initial
%         estimation-error covariance (VARIANCES); e0 = [ia, w, theta] initial estimate error so that
%         xh^-(0) = x(0) + e0 (a non-zero default shows the filter converging).
%   Q_KF / R_KF are completely separate from Q_MPC / R_MPC (different variables, different panels).
%   Noise is injected only in the MPC simulation; the classical P/PI/PD/PID transfer-function
%   simulations are deterministic and unchanged. With process noise 'None' and Mode A the MPC loop is
%   identical to the original deterministic one.
%   Failure handling: a non-finite KF value stops the run with an error; an ill-conditioned innovation
%   covariance skips that single correction (prediction only) and is COUNTED and reported in the Status
%   panel; if the QP is infeasible or the solver fails (e.g. an ENFORCED limit that the noisy / estimated
%   state already violates), the run stops with an error message - constraints are never relaxed or
%   dropped silently.
%   Experiments (section I): 1) Direct vs KF MPC, 2) measurement-noise Low/Medium/High,
%   3) process-noise Low/Medium/High, 4) combined 3 x 3 noise grid; results in a table window and
%   overlaid in the plots.
%
%   PERFORMANCE METRICS (panel under the plots; recomputed after every simulation)
%   All metrics are computed by ONE function, computeMotorMetrics, from the signals the
%   simulation actually produced (time, w, theta, ia and the applied armature voltage Va).
%   Tracking error: e = wref - w (speed mode) or e = thref - theta (position mode).
%     1. Steady-state error  e_ss = ref - y_ss, y_ss = time-weighted mean of y over the
%        final DEF.ssWindow (10 %) of the run. The output counts as settled only if its
%        spread inside that window is <= DEF.settleTol (0.1 %) of max(|ref|,|y|); otherwise
%        the value is labelled "NOT settled" (end-of-run error, not a steady-state error).
%     2. ISE = trapz(t, e.^2) on the simulation time vector     [(rad/s)^2.s or rad^2.s]
%     3. Peak armature current  max|ia| over all samples                         [A]
%     4. Control effort  Ju = integral of Va^2 dt  and  V_RMS = sqrt(Ju/T)  [V^2.s, V]
%        Va is the voltage applied to the plant: MPC -> u(k) held over each Ts (exact
%        integral sum(u.^2.*dt)); classical -> derived from the SAME closed loop,
%        Va/Ref = Nc*D/Dcl (speed) or s*D*Nc/Dcl (position), simulated with lsim.
%        PD/PID: the ideal derivative acts on a reference STEP, so Va contains a Dirac
%        impulse Kd*ref*delta(t) (it is what makes ia jump by Kd*ref/La at t = 0).
%        Its energy is infinite, so Ju, V_RMS and max|Va| are reported as Inf and the
%        values for t > 0 are given separately. Nothing is invented or clipped.
%     5. Limit violations: number of samples with Va > Vmax, Va < Vmin or |ia| > Imax
%        (a sample exactly ON a limit is not a violation; a 1e-6 relative numerical
%        tolerance, DEF.limTolRel, absorbs solver round-off), plus max|Va| and max|ia|.
%   MONITORING vs ENFORCEMENT: the limits entered in the GUI are always MONITORED for every
%   controller. They are ENFORCED only by MPC, only for the tick boxes that are ticked, and
%   only at the sampling instants. P/PI/PD/PID are linear transfer-function simulations that
%   are NOT saturated: their voltage and current are reported, never limited.
%   Defaults of the limits (+-12 V, 1.5 A) are ILLUSTRATIVE placeholders, like Kt and Ke.
%
%   NOTE: Kt = Ke = 0.0395 are ILLUSTRATIVE PLACEHOLDER values. Replace them with
%   datasheet values for your motor.
%
%   ESTIMATION METRICS (KF mode only; otherwise shown as n/a, never 0):
%       RMSE_x = sqrt(mean((x_actual - x_estimated).^2)) over all samples, for ia, w and theta.
%   Speed-tracking ISE and steady-state speed error are the existing ISE / e_ss of the panel in Speed
%   Control mode (in Position Control mode they refer to the position, as before).
%
%   Requires: MATLAB R2019b or newer (uifigure/uigridlayout/uitable with table data)
%             and the Control System Toolbox (tf, lsim) for P/PI/PD/PID.
%             MPC itself needs only base MATLAB (expm, chol); quadprog is needed
%             only if a voltage and/or current limit is ENFORCED in MPC (tick boxes).
%             The performance metrics need only base MATLAB (trapz).

% ------------------------------------------------------------------------
% Toolbox check
% ------------------------------------------------------------------------
if ~hasControlToolbox()
    msg = ['The Control System Toolbox is required (functions tf and lsim) ' ...
           'but is not available in this MATLAB installation or license.'];
    tmpFig = uifigure('Name','Control System Toolbox required');
    uialert(tmpFig, msg, 'Toolbox Missing', 'Icon', 'error');
    error('DCMotorGUI:NoToolbox', '%s', msg);
end

% ------------------------------------------------------------------------
% Defaults and constants
% ------------------------------------------------------------------------
DEF.La    = 0.12;       % H
DEF.Ra    = 7;          % Ohm
DEF.J     = 2.12e-6;    % kg.m^2
DEF.B     = 6.04e-6;    % N.m.s/rad
DEF.Kt    = 0.0395;     % N.m/A   (illustrative placeholder)
DEF.Ke    = 0.0395;     % V.s/rad (illustrative placeholder)
DEF.wref  = 10;         % rad/s
DEF.T     = 5;          % s
DEF.names = {'P','PI','PD','PID'};
DEF.gains = [0.05 0    0;       % P   : Kp Ki Kd
             0.05 0.20 0;       % PI
             0.05 0    0.001;   % PD
             0.05 0.20 0.001];  % PID

% Control-mode defaults
DEF.mode     = 'Speed Control';
DEF.thref    = 1;       % rad  (position-mode reference)
DEF.gainsPos = [1 0    0;       % P   : Kp Ki Kd   (position-mode gains, independent of speed mode)
                1 2    0;       % PI
                1 0    0.005;   % PD
                1 2    0.005];  % PID

% MPC defaults (documented in the header comment)
DEF.all      = [DEF.names {'MPC'}];  % dropdown entries
DEF.Q        = 1;       % tracking-error coefficient on (w - wref)^2   [1/(rad/s)^2]
DEF.R        = 1;       % control-effort coefficient on Va^2           [1/V^2]
DEF.Np       = 20;      % prediction horizon [samples] (20 x 5 ms = 100 ms)
DEF.Nc       = 10;      % control horizon [samples] (editable in the GUI, 1 <= Nc <= Np)
DEF.TsMax    = 5e-3;    % sampling interval [s]
DEF.MaxNp    = 200;     % practical upper limit for Np
DEF.MaxSteps = 20000;   % practical limit on MPC samples per run (=> T <= 100 s at 5 ms)

% Noise and Kalman-filter defaults (documented in the header comment). All noise values are
% STANDARD DEVIATIONS; covariances (variances) are built from their squares.
DEF.fbNames   = {'MPC using Direct States','MPC using Kalman Filter Estimates'};
DEF.levelNames = {'None','Low','Medium','High'};
%                        sigma_ia[A] sigma_w[rad/s] sigma_theta[rad] (position sensor, optional)
DEF.measPresets = [0      0      0;        % None
                   0.002  0.02   0.002;    % Low
                   0.01   0.1    0.01;     % Medium
                   0.05   0.5    0.05];    % High
DEF.measTheta = false;             % position sensor (theta measured) off by default
%                        sigma_ia[A] sigma_w[rad/s] sigma_theta[rad]   (per sample)
DEF.procPresets = [0      0      0;       % None
                   0.001  0.005  2.5e-4;  % Low
                   0.005  0.02   1e-3;    % Medium
                   0.02   0.1    5e-3];   % High
DEF.measLevel = 'Medium';
DEF.procLevel = 'Low';
DEF.seed      = 1;                 % random-number seed (reproducible, identical for both modes)
DEF.qScale    = 1;                 % Q_KF = qScale*Qtrue  (1 = filter knows the true process noise)
DEF.rScale    = 1;                 % R_KF = rScale*Rtrue  (1 = filter knows the true measurement noise)
DEF.P0        = [2.5e-3 1 1e-2];   % diag of P0 [A^2, (rad/s)^2, rad^2]  (= e0.^2 by default)
DEF.e0        = [0.05 1 0.1];      % initial estimate error [A, rad/s, rad]: xh^-(0) = x(0) + e0
DEF.covFloor  = 1e-12;             % floor on the diagonal of Q_KF and R_KF
DEF.expNames  = {'1. Direct-state MPC vs Kalman Filter MPC', ...
                 '2. Measurement-noise variation (Low / Medium / High)', ...
                 '3. Process-noise variation (Low / Medium / High)', ...
                 '4. Combined noise variation (3 x 3 grid)'};

% Voltage / current limits (GUI fields). ILLUSTRATIVE placeholders - replace with your motor's
% ratings. They are always MONITORED; MPC enforces them only if the tick boxes are ticked.
DEF.Vmin     = -12;     % lower armature-voltage limit [V]
DEF.Vmax     =  12;     % upper armature-voltage limit [V]
DEF.Imax     = 1.5;     % maximum permitted |armature current| [A]

% Performance-metric settings (documented in the header comment)
DEF.ssWindow  = 0.10;   % steady-state window = final 10 % of the run
DEF.settleTol = 1e-3;   % window spread must be <= 0.1 % of max(|ref|,|y|) to call it settled
DEF.limTolRel = 1e-6;   % relative tolerance for "exceeds a limit" (solver / round-off)

NPTS = 1001;                    % number of time samples per simulation

COL.bg     = [0.94 0.94 0.94];  % window background
COL.panel  = [0.98 0.98 0.98];  % section background
COL.blue   = [0.00 0.32 0.60];  % blue section headers
COL.ok     = [0.00 0.40 0.10];
COL.warn   = [0.78 0.42 0.00];
COL.err    = [0.75 0.00 0.00];
COL.info   = [0.20 0.20 0.20];

% ------------------------------------------------------------------------
% Build the interface
% ------------------------------------------------------------------------
fig = uifigure('Name','DC Motor Speed / Position Control - P, PI, PD, PID, MPC + Kalman Filter', ...
               'Position',[40 10 1400 1040], 'Color',COL.bg);

mainGrid = uigridlayout(fig,[1 2]);
mainGrid.ColumnWidth   = {450,'1x'};
mainGrid.RowHeight     = {'1x'};
mainGrid.Padding       = [8 8 8 8];
mainGrid.ColumnSpacing = 8;
mainGrid.BackgroundColor = COL.bg;

% ---- Left column: simulation settings -----------------------------------
leftGrid = uigridlayout(mainGrid,[12 1]);
leftGrid.Layout.Column = 1;
%                        A   Mode  B   C   D   E   F  G   H   I   Act Status
leftGrid.RowHeight     = {180,74,104,76,246,132,76,280,260,110,120,110};
leftGrid.Padding       = [0 0 0 0];
leftGrid.RowSpacing    = 4;
leftGrid.BackgroundColor = COL.bg;
leftGrid.Scrollable    = 'on';

% A. Motor parameters
motorGrid = createSection(leftGrid,1,'A. Motor parameters', ...
                          {26,26,26,40}, {110,'1x',110,'1x'}, COL);
addLabel(motorGrid,'La [H]',2,1);
edLa = addNumField(motorGrid,DEF.La,2,2,0,Inf,'off');
addLabel(motorGrid,'Ra [Ohm]',2,3);
edRa = addNumField(motorGrid,DEF.Ra,2,4,0,Inf,'off');
addLabel(motorGrid,'J [kg.m^2]',3,1);
edJ  = addNumField(motorGrid,DEF.J,3,2,0,Inf,'off');
addLabel(motorGrid,'B [N.m.s/rad]',3,3);
edB  = addNumField(motorGrid,DEF.B,3,4,0,Inf,'on');
addLabel(motorGrid,'Kt [N.m/A] *',4,1);
edKt = addNumField(motorGrid,DEF.Kt,4,2,0,Inf,'off');
addLabel(motorGrid,'Ke [V.s/rad] *',4,3);
edKe = addNumField(motorGrid,DEF.Ke,4,4,0,Inf,'off');
noteMotor = addLabel(motorGrid,['* Kt and Ke are illustrative placeholder values - ' ...
                                'replace them with your motor''s datasheet values.'],5,1);
noteMotor.Layout.Column = [1 4];
noteMotor.WordWrap = 'on';
noteMotor.FontAngle = 'italic';
noteMotor.FontSize = 11;
noteMotor.FontColor = COL.err;

% Control Mode (two mutually exclusive radio buttons)
modeGrid = createSection(leftGrid,2,'Control Mode',{28},{'1x'},COL);
bgMode = uibuttongroup(modeGrid,'BorderType','none','BackgroundColor',COL.panel, ...
                       'SelectionChangedFcn',@onModeChanged);
bgMode.Layout.Row = 2;  bgMode.Layout.Column = 1;
rbSpeed = uiradiobutton(bgMode,'Text','Speed Control','FontSize',12, ...
                        'Position',[10 3 190 22]);
rbPos   = uiradiobutton(bgMode,'Text','Position Control','FontSize',12, ...
                        'Position',[220 3 190 22]);
bgMode.SelectedObject = rbSpeed;

% B. Reference input (only the field for the selected mode is visible)
[refGrid,hdrRef] = createSection(leftGrid,3,'B. Reference input (speed)', ...
                        {26,26}, {170,'1x'}, COL);
lblRefSpd = addLabel(refGrid,'Reference speed [rad/s]',2,1);
edWref = addNumField(refGrid,DEF.wref,2,2,-Inf,Inf,'on');
lblRefPos = addLabel(refGrid,'Reference position [rad]',2,1);
edThref = addNumField(refGrid,DEF.thref,2,2,-Inf,Inf,'on');
lblRefPos.Visible = 'off';  edThref.Visible = 'off';
addLabel(refGrid,'Simulation duration [s]',3,1);
edTsim = addNumField(refGrid,DEF.T,3,2,0,1000,'off');

% C. Controller selection
ctrlGrid = createSection(leftGrid,4,'C. Controller selection', ...
                         {28}, {170,'1x'}, COL);
addLabel(ctrlGrid,'Selected controller',2,1);
ddCtrl = uidropdown(ctrlGrid,'Items',DEF.all,'Value','P','FontSize',12, ...
                    'ValueChangedFcn',@onControllerChanged);
ddCtrl.Layout.Row = 2;  ddCtrl.Layout.Column = 2;

% D. Independent controller gain table
[tblGrid,hdrGain] = createSection(leftGrid,5,'D. Controller gains - speed mode (editable)', ...
                        {150,48}, {'1x'}, COL);
tblGains = uitable(tblGrid, ...
                   'Data',makeGainTable(DEF), ...
                   'ColumnEditable',[false true true true], ...
                   'ColumnWidth',{105,100,100,100}, ...
                   'RowName',{}, ...
                   'ColumnSortable',false, ...
                   'FontSize',13, ...
                   'CellEditCallback',@onTableEdited);
tblGains.Layout.Row = 2;
% Separate gain table for position mode, in the same cell (only one is visible).
tblGainsPos = uitable(tblGrid, ...
                   'Data',makeGainTable(DEF,DEF.gainsPos), ...
                   'ColumnEditable',[false true true true], ...
                   'ColumnWidth',{105,100,100,100}, ...
                   'RowName',{}, ...
                   'ColumnSortable',false, ...
                   'FontSize',13, ...
                   'Visible','off', ...
                   'CellEditCallback',@onTableEdited);
tblGainsPos.Layout.Row = 2;
noteTbl = addLabel(tblGrid,['Double-click a cell to edit; press Enter or click another cell ' ...
                            'to commit. Each controller keeps its own gains, separately for speed and position mode. Unused gains are ' ...
                            'ignored: P uses Kp; PI uses Kp, Ki; PD uses Kp, Kd; PID uses all.'],3,1);
noteTbl.WordWrap = 'on';
noteTbl.FontSize = 11;
noteTbl.FontColor = [0.3 0.3 0.3];

% E. MPC parameters (read at every simulation run)
mpcGrid = createSection(leftGrid,6,'E. MPC parameters (used when MPC is selected)', ...
                        {26,26,22}, {110,'1x',110,'1x'}, COL);
lblQ = addLabel(mpcGrid,'Q_MPC (speed)',2,1);
edQ  = addNumField(mpcGrid,DEF.Q,2,2,0,Inf,'off');
lblR = addLabel(mpcGrid,'R_MPC (effort)',2,3);
edR  = addNumField(mpcGrid,DEF.R,2,4,0,Inf,'off');
addLabel(mpcGrid,'Np [samples]',3,1);
edNp = addNumField(mpcGrid,DEF.Np,3,2,1,DEF.MaxNp,'on');
edNp.RoundFractionalValues = 'on';
edNp.ValueDisplayFormat    = '%d';
addLabel(mpcGrid,'Nc [samples]',3,3);
edNc = addNumField(mpcGrid,DEF.Nc,3,4,1,DEF.MaxNp,'on');
edNc.RoundFractionalValues = 'on';
edNc.ValueDisplayFormat    = '%d';
noteMpc = addLabel(mpcGrid,sprintf('Q_MPC/R_MPC: cost weights (NOT the Kalman Q_KF/R_KF). Np, Nc integers, Nc <= Np. Ts <= %g ms (fixed)',DEF.TsMax*1e3),4,1);
noteMpc.Layout.Column = [1 4];
noteMpc.FontSize = 10;
noteMpc.FontColor = [0.3 0.3 0.3];

% F. MPC feedback configuration (Direct states / Kalman Filter estimates)
fbGrid = createSection(leftGrid,7,'F. MPC feedback configuration', ...
                       {28}, {110,'1x'}, COL);
addLabel(fbGrid,'MPC receives',2,1);
ddFeedback = uidropdown(fbGrid,'Items',DEF.fbNames,'Value',DEF.fbNames{1},'FontSize',12, ...
                        'ValueChangedFcn',@onFeedbackChanged);
ddFeedback.Layout.Row = 2;  ddFeedback.Layout.Column = 2;

% G. Noise injected into the simulated MPC loop (standard deviations)
noiseGrid = createSection(leftGrid,8,'G. Measurement / process noise (simulated, std. dev.)', ...
                          {26,26,26,26,26,26,44}, {110,'1x',110,'1x'}, COL);
addLabel(noiseGrid,'Meas. noise level',2,1);
ddMeasLevel = uidropdown(noiseGrid,'Items',[DEF.levelNames {'Custom'}],'Value',DEF.measLevel, ...
                         'FontSize',12,'ValueChangedFcn',@onMeasLevelChanged);
ddMeasLevel.Layout.Row = 2;  ddMeasLevel.Layout.Column = 2;
addLabel(noiseGrid,'Random seed',2,3);
edSeed = addNumField(noiseGrid,DEF.seed,2,4,0,4294967295,'on');
edSeed.RoundFractionalValues = 'on';  edSeed.ValueDisplayFormat = '%d';
addLabel(noiseGrid,'sigma_r ia [A]',3,1);
edSigRi = addNumField(noiseGrid,DEF.measPresets(3,1),3,2,0,Inf,'on');
addLabel(noiseGrid,'sigma_r w [rad/s]',3,3);
edSigRw = addNumField(noiseGrid,DEF.measPresets(3,2),3,4,0,Inf,'on');
addLabel(noiseGrid,'sigma_r theta [rad]',4,1);
edSigRt = addNumField(noiseGrid,DEF.measPresets(3,3),4,2,0,Inf,'on');
cbMeasTheta = uicheckbox(noiseGrid,'Text','Position sensor (measure theta)','Value',DEF.measTheta,'FontSize',11, ...
                         'Tooltip',['Unticked (default): only ia and w are measured, C = [1 0 0; 0 1 0]; theta is ' ...
                                    'then unobservable. Ticked: theta is also measured (sigma_r theta).']);
cbMeasTheta.Layout.Row = 4;  cbMeasTheta.Layout.Column = [3 4];
addLabel(noiseGrid,'Proc. noise level',5,1);
ddProcLevel = uidropdown(noiseGrid,'Items',[DEF.levelNames {'Custom'}],'Value',DEF.procLevel, ...
                         'FontSize',12,'ValueChangedFcn',@onProcLevelChanged);
ddProcLevel.Layout.Row = 5;  ddProcLevel.Layout.Column = 2;
addLabel(noiseGrid,'sigma_q ia [A]',6,1);
edSigQi = addNumField(noiseGrid,DEF.procPresets(2,1),6,2,0,Inf,'on');
addLabel(noiseGrid,'sigma_q w [rad/s]',6,3);
edSigQw = addNumField(noiseGrid,DEF.procPresets(2,2),6,4,0,Inf,'on');
addLabel(noiseGrid,'sigma_q theta [rad]',7,1);
edSigQt = addNumField(noiseGrid,DEF.procPresets(2,3),7,2,0,Inf,'on');
noteNoise = addLabel(noiseGrid,['Standard deviations (variance = square). sigma_r: sensor noise (ia, w; theta only ' ...
    'with the position sensor). sigma_q: per-sample process noise on [ia, w, theta], added to the plant in ' ...
    'BOTH MPC modes. Same seed = same noise realisation. Used only by MPC; 0 = no noise.'],8,1);
noteNoise.Layout.Column = [1 4];  noteNoise.WordWrap = 'on';  noteNoise.FontSize = 10;
noteNoise.VerticalAlignment = 'top';  noteNoise.FontColor = [0.3 0.3 0.3];

% H. Kalman Filter parameters (separate from the MPC weights)
kfGrid = createSection(leftGrid,9,'H. Kalman Filter Parameters', ...
                       {26,26,26,26,26,58}, {110,'1x',110,'1x'}, COL);
addLabel(kfGrid,'Q_KF scale',2,1);
edQkfScale = addNumField(kfGrid,DEF.qScale,2,2,0,Inf,'off');
addLabel(kfGrid,'R_KF scale',2,3);
edRkfScale = addNumField(kfGrid,DEF.rScale,2,4,0,Inf,'off');
addLabel(kfGrid,'P0 ia [A^2]',3,1);
edP0i = addNumField(kfGrid,DEF.P0(1),3,2,0,Inf,'off');
addLabel(kfGrid,'P0 w [(rad/s)^2]',3,3);
edP0w = addNumField(kfGrid,DEF.P0(2),3,4,0,Inf,'off');
addLabel(kfGrid,'P0 theta [rad^2]',4,1);
edP0t = addNumField(kfGrid,DEF.P0(3),4,2,0,Inf,'off');
addLabel(kfGrid,'e0 ia [A]',5,1);
edE0i = addNumField(kfGrid,DEF.e0(1),5,2,-Inf,Inf,'on');
addLabel(kfGrid,'e0 w [rad/s]',5,3);
edE0w = addNumField(kfGrid,DEF.e0(2),5,4,-Inf,Inf,'on');
addLabel(kfGrid,'e0 theta [rad]',6,1);
edE0t = addNumField(kfGrid,DEF.e0(3),6,2,-Inf,Inf,'on');
noteKf = addLabel(kfGrid,['Q_KF = (Q_KF scale) x diag(sigma_q^2),  R_KF = (R_KF scale) x diag(sigma_r^2) '  ...
    '(variances; scale 1 = filter matched to the noise in G, <> 1 = mismatched filter; large Q_KF/R_KF ' ...
    '= trust the measurements more). P0 = diag of the initial error covariance (variances). ' ...
    'e0 = initial estimate error: xhat(0) = x(0) + e0.'],7,1);
noteKf.Layout.Column = [1 4];  noteKf.WordWrap = 'on';  noteKf.FontSize = 10;
noteKf.VerticalAlignment = 'top';  noteKf.FontColor = [0.3 0.3 0.3];

% I. Experiments
expGrid = createSection(leftGrid,10,'I. MPC / Kalman Filter experiments', ...
                        {28,28}, {'1x'}, COL);
ddExp = uidropdown(expGrid,'Items',DEF.expNames,'Value',DEF.expNames{1},'FontSize',12);
ddExp.Layout.Row = 2;  ddExp.Layout.Column = 1;
btnExp = uibutton(expGrid,'push','Text','Run Experiment (uses current settings)', ...
                  'FontSize',12,'ButtonPushedFcn',@(~,~) runExperiment());
btnExp.Layout.Row = 3;  btnExp.Layout.Column = 1;

% Actions
btnGrid = createSection(leftGrid,11,'Actions',{34,34},{'1x','1x'},COL);
btnRun = uibutton(btnGrid,'push','Text','Run Selected Controller', ...
                  'FontSize',13,'FontWeight','bold', ...
                  'BackgroundColor',COL.blue,'FontColor',[1 1 1], ...
                  'ButtonPushedFcn',@(~,~) runSelected(''));
btnRun.Layout.Row = 2;  btnRun.Layout.Column = [1 2];
btnAll = uibutton(btnGrid,'push','Text','Show All Controllers', ...
                  'FontSize',12,'ButtonPushedFcn',@(~,~) runAll());
btnAll.Layout.Row = 3;  btnAll.Layout.Column = 1;
btnReset = uibutton(btnGrid,'push','Text','Reset Defaults', ...
                    'FontSize',12,'ButtonPushedFcn',@onReset);
btnReset.Layout.Row = 3;  btnReset.Layout.Column = 2;

% Status
statGrid = createSection(leftGrid,12,'Status',{'1x'},{'1x'},COL);
statusLabel = uilabel(statGrid,'Text','Ready.','WordWrap','on', ...
                      'VerticalAlignment','top','FontSize',12,'FontColor',COL.info);
statusLabel.Layout.Row = 2;

% ---- Right column: three vertically stacked plots + metrics / limits panels ----
rightGrid = uigridlayout(mainGrid,[5 1]);
rightGrid.Layout.Column = 2;
rightGrid.RowHeight     = {'1x','1x','1x','1x',330};
rightGrid.ColumnWidth   = {'1x'};
rightGrid.Padding       = [0 0 0 0];
rightGrid.RowSpacing    = 8;
rightGrid.BackgroundColor = COL.bg;

axI  = uiaxes(rightGrid);  axI.Layout.Row  = 1;
axTh = uiaxes(rightGrid);  axTh.Layout.Row = 2;
axW  = uiaxes(rightGrid);  axW.Layout.Row  = 3;
axV  = uiaxes(rightGrid);  axV.Layout.Row  = 4;
prepareAxes(axI, 'Armature Current i(t)',   'Current i(t) [A]',            DEF.T,[0 0 0]);
prepareAxes(axTh,'Angular Position \theta(t)','Position \theta(t) [rad]', DEF.T,[0 0 0]);
prepareAxes(axW, 'Angular Speed \omega(t)', 'Speed \omega(t) [rad/s]',    DEF.T,[0 0 0]);
prepareAxes(axV, 'Applied Armature Voltage V_a(t)', 'Voltage V_a(t) [V]',  DEF.T,[0 0 0]);

% ---- Performance Metrics panel (left) and voltage / current limits (right) ----
metricsShown = false;   % true only while the panel shows the result of the CURRENT inputs
metCont = uigridlayout(rightGrid,[1 2]);
metCont.Layout.Row      = 5;
metCont.ColumnWidth     = {'1x',245};
metCont.RowHeight       = {'1x'};
metCont.Padding         = [0 0 0 0];
metCont.ColumnSpacing   = 8;
metCont.BackgroundColor = COL.bg;

mg = createSection(metCont,1,'Performance Metrics', ...
                   {22,22,22,22,22,22,22,22,22,56},{135,'1x',135,'1x'},COL);
mg.Parent.Layout.Column = 1;
mg.RowSpacing = 3;
addLabel(mg,'Result for',2,1);
MH.valCase = addLabel(mg,'',2,2);   MH.valCase.Layout.Column = [2 4];
MH.valCase.FontWeight = 'bold';
addLabel(mg,'Steady-state error',3,1);
MH.valSS   = addLabel(mg,'',3,2);   MH.valSS.Layout.Column = [2 4];
addLabel(mg,'ISE',4,1);                    MH.valISE  = addLabel(mg,'',4,2);
addLabel(mg,'Peak armature current',4,3);  MH.valIpk  = addLabel(mg,'',4,4);
addLabel(mg,'Control effort Ju',5,1);      MH.valJu   = addLabel(mg,'',5,2);
addLabel(mg,'RMS applied voltage',5,3);    MH.valVrms = addLabel(mg,'',5,4);
addLabel(mg,'Max |Va| observed',6,1);      MH.valVmax = addLabel(mg,'',6,2);
addLabel(mg,'Max |ia| observed',6,3);      MH.valImax = addLabel(mg,'',6,4);
addLabel(mg,'Voltage violations',7,1);
MH.valVviol = addLabel(mg,'',7,2);  MH.valVviol.Layout.Column = [2 4];
addLabel(mg,'Current violations',8,1);
MH.valIviol = addLabel(mg,'',8,2);  MH.valIviol.Layout.Column = [2 4];
addLabel(mg,'Est. RMSE ia',9,1);           MH.valRmseI = addLabel(mg,'',9,2);
addLabel(mg,'Est. RMSE omega',9,3);        MH.valRmseW = addLabel(mg,'',9,4);
addLabel(mg,'Est. RMSE theta',10,1);       MH.valRmseT = addLabel(mg,'',10,2);
addLabel(mg,'MPC feedback',10,3);          MH.valFb    = addLabel(mg,'',10,4);
MH.valNote = addLabel(mg,'',11,1);  MH.valNote.Layout.Column = [1 4];
MH.valNote.WordWrap = 'on';  MH.valNote.FontSize = 11;  MH.valNote.VerticalAlignment = 'top';

lg = createSection(metCont,1,'Voltage / current limits', ...
                   {22,22,22,22,22,46},{100,'1x'},COL);
lg.Parent.Layout.Column = 2;
lg.RowSpacing = 3;
addLabel(lg,'Vmin [V]',2,1);       edVmin = addNumField(lg,DEF.Vmin,2,2,-Inf,Inf,'on');
addLabel(lg,'Vmax [V]',3,1);       edVmax = addNumField(lg,DEF.Vmax,3,2,-Inf,Inf,'on');
addLabel(lg,'Max |ia| [A]',4,1);   edImax = addNumField(lg,DEF.Imax,4,2,0,Inf,'off');
cbEnfV = uicheckbox(lg,'Text','MPC enforces V limits','Value',false,'FontSize',11, ...
                    'Tooltip',['Unticked: voltage limits are only MONITORED. Ticked: MPC adds ' ...
                               'Vmin <= u <= Vmax to its optimisation (needs quadprog).']);
cbEnfV.Layout.Row = 5;  cbEnfV.Layout.Column = [1 2];
cbEnfI = uicheckbox(lg,'Text','MPC enforces I limit','Value',false,'FontSize',11, ...
                    'Tooltip',['Unticked: the current limit is only MONITORED. Ticked: MPC adds ' ...
                               '|ia| <= Imax at the predicted sampling instants (needs quadprog).']);
cbEnfI.Layout.Row = 6;  cbEnfI.Layout.Column = [1 2];
noteLim = addLabel(lg,['Illustrative placeholders: enter your motor''s ratings. Limits are always ' ...
                       'monitored; only ticked boxes are enforced (MPC only).'],7,1);
noteLim.Layout.Column = [1 2];  noteLim.WordWrap = 'on';  noteLim.FontSize = 10;
noteLim.FontAngle = 'italic';   noteLim.VerticalAlignment = 'top';
noteLim.FontColor = [0.3 0.3 0.3];

% Any input change makes the numbers on screen out of date -> flag them (never silently stale).
staleList = {edLa,edRa,edJ,edB,edKt,edKe,edWref,edThref,edTsim,edQ,edR,edNp,edNc, ...
             edVmin,edVmax,edImax,cbEnfV,cbEnfI,edSeed,cbMeasTheta,edQkfScale,edRkfScale, ...
             edP0i,edP0w,edP0t,edE0i,edE0w,edE0t};
for kStale = 1:numel(staleList)
    staleList{kStale}.ValueChangedFcn = @(~,~) markMetricsStale();
end

% Editing a sigma field switches its level menu to 'Custom' (values are then user-defined).
edSigRi.ValueChangedFcn = @(~,~) onSigmaEdited('meas');
edSigRw.ValueChangedFcn = @(~,~) onSigmaEdited('meas');
edSigRt.ValueChangedFcn = @(~,~) onSigmaEdited('meas');
edSigQi.ValueChangedFcn = @(~,~) onSigmaEdited('proc');
edSigQw.ValueChangedFcn = @(~,~) onSigmaEdited('proc');
edSigQt.ValueChangedFcn = @(~,~) onSigmaEdited('proc');
expFig = [];            % handle of the experiment-results window (if open)

highlightSelectedRow();
applyModeUI();
resetMetricsPanel('Run a controller to see its performance metrics.');
drawnow;

% ---- Initial default simulation so all three axes are populated ---------
runSelected('');

% ========================================================================
%                          NESTED CALLBACKS / HELPERS
%   (nested so they share the UI handles defined above)
% ========================================================================

    function onModeChanged(~,~)
        % Switching mode only changes which reference field / gain table is shown and how Q
        % is labelled. Motor parameters, gains and MPC values are NOT modified.
        applyModeUI();
        resetMetricsPanel();      % definitions / units of the new mode, no stale numbers
        setStatus(sprintf(['%s selected. Motor parameters, controller gains and MPC settings are ' ...
                           'unchanged. Click "Run Selected Controller" or "Show All Controllers".'], ...
                          curMode()),'info');
    end

    function applyModeUI()
        isPos = strcmp(curMode(),'Position Control');
        lblRefSpd.Visible = onoff(~isPos);   edWref.Visible  = onoff(~isPos);
        lblRefPos.Visible = onoff(isPos);    edThref.Visible = onoff(isPos);
        tblGains.Visible    = onoff(~isPos);
        tblGainsPos.Visible = onoff(isPos);
        if isPos
            hdrRef.Text  = '  B. Reference input (position)';
            hdrGain.Text = '  D. Controller gains - position mode (editable)';
            lblQ.Text    = 'Q_MPC (position)';
        else
            hdrRef.Text  = '  B. Reference input (speed)';
            hdrGain.Text = '  D. Controller gains - speed mode (editable)';
            lblQ.Text    = 'Q_MPC (speed)';
        end
    end

    function md = curMode()
        md = bgMode.SelectedObject.Text;
    end

    function k = modeKey()
        if strcmp(curMode(),'Position Control'), k = 'position'; else, k = 'speed'; end
    end

    function t = modeTag()
        % '' in speed mode so the original speed-mode messages are unchanged.
        if strcmp(modeKey(),'position'), t = ' (position control)'; else, t = ''; end
    end

    function tb = activeTable()
        if strcmp(modeKey(),'position'), tb = tblGainsPos; else, tb = tblGains; end
    end

    function onControllerChanged(~,~)
        % Changing the selection NEVER touches the gain table.
        highlightSelectedRow();
        resetMetricsPanel();      % the metrics on screen belonged to the previous controller
        if strcmp(ddCtrl.Value,'MPC')
            setStatus(['MPC selected. Q_MPC, R_MPC, Np, Nc (section E), the feedback configuration (F), ' ...
                       'noise (G) and Kalman Filter settings (H) are read when you run; the PID gain ' ...
                       'table is unchanged. Click "Run Selected Controller" to simulate.'],'info');
        else
            setStatus(sprintf(['%s selected. Gains in the table are unchanged. ' ...
                               'Click "Run Selected Controller" to simulate.'], ddCtrl.Value),'info');
        end
    end

    function onFeedbackChanged(~,~)
        % Mode A / Mode B selector: nothing is modified, the shown metrics belong to the other mode.
        resetMetricsPanel();
        setStatus(sprintf(['"%s" selected (applies to the MPC controller). Click "Run Selected ' ...
                           'Controller" to simulate.'],ddFeedback.Value),'info');
    end

    function onMeasLevelChanged(~,~)
        % Presets fill the (editable) sigma_r fields; 'Custom' leaves them untouched.
        L = find(strcmp(DEF.levelNames,ddMeasLevel.Value),1);
        if ~isempty(L)
            edSigRi.Value = DEF.measPresets(L,1);
            edSigRw.Value = DEF.measPresets(L,2);
            edSigRt.Value = DEF.measPresets(L,3);
        end
        markMetricsStale();
    end

    function onProcLevelChanged(~,~)
        L = find(strcmp(DEF.levelNames,ddProcLevel.Value),1);
        if ~isempty(L)
            edSigQi.Value = DEF.procPresets(L,1);
            edSigQw.Value = DEF.procPresets(L,2);
            edSigQt.Value = DEF.procPresets(L,3);
        end
        markMetricsStale();
    end

    function onSigmaEdited(which)
        if strcmp(which,'meas'), ddMeasLevel.Value = 'Custom'; else, ddProcLevel.Value = 'Custom'; end
        markMetricsStale();
    end

    function onTableEdited(src,evt)
        % Validate the single edited cell; revert it if invalid.
        editRow = evt.Indices(1);
        editCol = evt.Indices(2);
        newVal  = evt.NewData;
        isGood  = isnumeric(newVal) && isscalar(newVal) && isreal(newVal) && ...
                  isfinite(newVal) && newVal >= 0;
        if ~isGood
            curTbl = src.Data;
            curTbl{editRow,editCol} = evt.PreviousData;   % revert only this cell
            src.Data = curTbl;
            colNames = src.Data.Properties.VariableNames;
            msgTxt = sprintf(['Invalid entry for %s of the %s controller. Gains must be ' ...
                              'finite, non-negative numbers. The previous value was restored.'], ...
                             colNames{editCol}, curTbl.Controller{editRow});
            setStatus(msgTxt,'error');
            uialert(fig,msgTxt,'Invalid Gain','Icon','warning');
        else
            markMetricsStale();   % a valid gain change makes the shown metrics out of date
        end
        drawnow;
    end

    function runSelected(prefixTxt)
        % Simulate the controller chosen in the dropdown and update the plots.
        if nargin < 1, prefixTxt = ''; end
        try
            setStatus('Simulating...','info');
            drawnow;
            mot  = readMotorAndReference();
            lim  = readLimits();
            nm   = ddCtrl.Value;
            if strcmp(nm,'MPC')
                mpc = readMPC();
                res = simulateMPC(mot,mpc,NPTS);
                plotSingle(nm,mot,res);
                showMetrics(nm,mot,res,lim,'');
                if res.stable
                    stText = [mpcStabilityText(res) mpcKfText(res)];
                    if res.kfSkipped > 0 || (res.useKF && strcmp(res.mode,'position') && ~res.measTheta)
                        setStatus(sprintf('%sMPC simulated (%s). %s',prefixTxt,mpcString(res),stText),'warn');
                    else
                        setStatus(sprintf('%sMPC simulated (%s). %s',prefixTxt,mpcString(res),stText),'ok');
                    end
                else
                    wTxt = sprintf(['MPC (%s) gives an UNSTABLE closed loop (largest discrete-time pole ' ...
                        'magnitude %.4f >= 1). The plotted response is NOT a successful control result. ' ...
                        'Change Q, R or Np.'], mpcString(res), res.rho);
                    setStatus(['WARNING: ' wTxt],'warn');
                    uialert(fig,wTxt,'Unstable Closed Loop','Icon','warning');
                end
            else
                [kp,ki,kd] = readGains(nm);
                validateGains(nm,kp,ki,kd);
                res  = simulateController(nm,mot,kp,ki,kd,NPTS);
                plotSingle(nm,mot,res);
                showMetrics(nm,mot,res,lim,'');
                gTxt = gainString(nm,kp,ki,kd);
                if res.stable
                    setStatus(sprintf(['%s%s controller%s simulated (%s). Closed loop is stable ' ...
                        '(largest pole real part %.4g 1/s).'], prefixTxt, nm, modeTag(), gTxt, res.maxRe),'ok');
                else
                    wTxt = sprintf(['%s controller%s (%s) gives an UNSTABLE closed loop ' ...
                        '(largest pole real part %.4g 1/s). The plotted response is diverging and is ' ...
                        'NOT a successful control result. Reduce the gains.'], nm, modeTag(), gTxt, res.maxRe);
                    setStatus(['WARNING: ' wTxt],'warn');
                    uialert(fig,wTxt,'Unstable Closed Loop','Icon','warning');
                end
            end
        catch ME
            resetMetricsPanel('Simulation failed - there is no valid result, so no metrics are shown.');
            setStatus(['Error: ' ME.message],'error');
            uialert(fig,ME.message,'Simulation Error','Icon','error');
        end
        drawnow;
    end

    function runAll()
        % Overlay P, PI, PD, PID and MPC with identical motor/reference/duration.
        try
            setStatus('Simulating all controllers...','info');
            drawnow;
            mot = readMotorAndReference();
            lim = readLimits();
            nCtl = numel(DEF.all);
            results = cell(1,nCtl);
            failList = {};
            for idx = 1:nCtl
                nm = DEF.all{idx};
                if strcmp(nm,'MPC')
                    % An MPC failure is reported but does not hide the classical results.
                    try
                        results{idx} = simulateMPC(mot,readMPC(),NPTS);
                    catch ME2
                        failList{end+1} = sprintf('MPC: %s',ME2.message); %#ok<AGROW>
                    end
                    continue;
                end
                try
                    [kp,ki,kd] = readGains(nm);
                    validateGains(nm,kp,ki,kd);
                    results{idx} = simulateController(nm,mot,kp,ki,kd,NPTS);
                catch ME2
                    error('DCMotorGUI:SimFailed','%s controller: %s',nm,ME2.message);
                end
            end
            sfx = modeTag();
            prepareAxes(axI, ['Armature Current i(t) - All Controllers' sfx], 'Current i(t) [A]',            mot.T,[0 0 0]);
            prepareAxes(axTh,['Angular Position \theta(t) - All Controllers' sfx],'Position \theta(t) [rad]', mot.T,[0 0 0]);
            prepareAxes(axW, ['Angular Speed \omega(t) - All Controllers' sfx], 'Speed \omega(t) [rad/s]',   mot.T,[0 0 0]);
            prepareAxes(axV, ['Applied Armature Voltage V_a(t) - All Controllers' sfx],'Voltage V_a(t) [V]', mot.T,[0 0 0]);
            unstableList = {};
            for idx = 1:nCtl
                nm = DEF.all{idx};
                r  = results{idx};
                if isempty(r), continue; end   % MPC failed: reported below
                if ~r.stable
                    unstableList{end+1} = nm; %#ok<AGROW>
                    continue;       % unstable responses are not drawn (they would hide the others)
                end
                cc = ctrlColor(nm);
                plot(axI, r.t,r.i, 'Color',cc,'LineWidth',1.8,'DisplayName',nm);
                plot(axTh,r.t,r.th,'Color',cc,'LineWidth',1.8,'DisplayName',nm);
                plot(axW, r.t,r.w, 'Color',cc,'LineWidth',1.8,'DisplayName',nm);
                plotVoltage(axV,r,cc,nm);
            end
            if strcmp(mot.mode,'position')
                plot(axTh,[0 mot.T],[mot.thref mot.thref],'--','Color',[0.15 0.15 0.15], ...
                     'LineWidth',1.4,'DisplayName','Reference position');
            else
                plot(axW,[0 mot.T],[mot.wref mot.wref],'--','Color',[0.15 0.15 0.15], ...
                     'LineWidth',1.4,'DisplayName','Reference speed');
            end
            legend(axI, 'Location','best');
            legend(axTh,'Location','best');
            legend(axW, 'Location','best');
            legend(axV, 'Location','best');
            % Metrics are per controller: show those of the controller selected in section C.
            selIdx = find(strcmp(DEF.all,ddCtrl.Value),1);
            if ~isempty(selIdx) && ~isempty(results{selIdx})
                showMetrics(ddCtrl.Value,mot,results{selIdx},lim, ...
                            'from the all-controllers run (controller selected in section C)');
            else
                resetMetricsPanel(sprintf(['%s did not simulate in the all-controllers run - ' ...
                                           'no metrics available.'],ddCtrl.Value));
            end
            if isempty(unstableList) && isempty(failList)
                setStatus(sprintf('All five controllers simulated%s; all closed loops are stable.%s', ...
                                  modeTag(),mpcEnforceNote()),'ok');
            else
                wTxt = '';
                if ~isempty(unstableList)
                    wTxt = sprintf(['Unstable closed loop for: %s. These controllers are NOT plotted ' ...
                        'because their responses diverge; this is not a successful control result. ' ...
                        'Reduce their gains (or change Q, R, Np for MPC) and run again. '], ...
                        strjoin(unstableList,', '));
                end
                if ~isempty(failList)
                    wTxt = [wTxt 'Simulation failed and is NOT plotted - ' strjoin(failList,'; ')];
                end
                setStatus(['WARNING: ' wTxt],'warn');
                uialert(fig,wTxt,'Controller Warning','Icon','warning');
            end
        catch ME
            resetMetricsPanel('Simulation failed - there is no valid result, so no metrics are shown.');
            setStatus(['Error: ' ME.message],'error');
            uialert(fig,ME.message,'Simulation Error','Icon','error');
        end
        drawnow;
    end

    function onReset(~,~)
        % Restore every default (motor, reference, duration, ALL gains, MPC Q/R/Np/Nc, feedback mode, noise, KF, selection = P).
        edLa.Value   = DEF.La;
        edRa.Value   = DEF.Ra;
        edJ.Value    = DEF.J;
        edB.Value    = DEF.B;
        edKt.Value   = DEF.Kt;
        edKe.Value   = DEF.Ke;
        edWref.Value = DEF.wref;
        edThref.Value = DEF.thref;
        edTsim.Value = DEF.T;
        tblGains.Data    = makeGainTable(DEF);
        tblGainsPos.Data = makeGainTable(DEF,DEF.gainsPos);
        bgMode.SelectedObject = rbSpeed;            % default mode = Speed Control
        applyModeUI();
        edQ.Value    = DEF.Q;
        edR.Value    = DEF.R;
        edNp.Value   = DEF.Np;
        edNc.Value   = DEF.Nc;
        ddFeedback.Value = DEF.fbNames{1};
        ddMeasLevel.Value = DEF.measLevel;   onMeasLevelChanged();
        ddProcLevel.Value = DEF.procLevel;   onProcLevelChanged();
        edSeed.Value = DEF.seed;
        cbMeasTheta.Value = DEF.measTheta;
        edQkfScale.Value = DEF.qScale;   edRkfScale.Value = DEF.rScale;
        edP0i.Value = DEF.P0(1);  edP0w.Value = DEF.P0(2);  edP0t.Value = DEF.P0(3);
        edE0i.Value = DEF.e0(1);  edE0w.Value = DEF.e0(2);  edE0t.Value = DEF.e0(3);
        ddExp.Value = DEF.expNames{1};
        edVmin.Value = DEF.Vmin;
        edVmax.Value = DEF.Vmax;
        edImax.Value = DEF.Imax;
        cbEnfV.Value = false;                       % limits back to "monitor only"
        cbEnfI.Value = false;
        ddCtrl.Value  = 'P';
        highlightSelectedRow();
        resetMetricsPanel();
        drawnow;
        runSelected('Defaults restored. ');
    end

    function mot = readMotorAndReference()
        % Always read the CURRENT values of the edit fields.
        mot.La   = edLa.Value;
        mot.Ra   = edRa.Value;
        mot.J    = edJ.Value;
        mot.B    = edB.Value;
        mot.Kt   = edKt.Value;
        mot.Ke   = edKe.Value;
        mot.wref  = edWref.Value;
        mot.thref = edThref.Value;
        mot.mode  = modeKey();                  % 'speed' or 'position'
        if strcmp(mot.mode,'position'), mot.ref = mot.thref; else, mot.ref = mot.wref; end
        mot.T    = edTsim.Value;
        validateInputs(mot);
    end

    function [kp,ki,kd] = readGains(nm)
        % Read the LATEST committed gains of controller 'nm' from the current mode's table.
        tb  = activeTable();            % gain table of the CURRENT mode
        tbl = tb.Data;
        if ~istable(tbl)
            error('DCMotorGUI:BadTable','The gain table has an unexpected data type.');
        end
        row = find(strcmp(tbl.Controller,nm),1);
        if isempty(row)
            error('DCMotorGUI:BadTable','Controller "%s" was not found in the gain table.',nm);
        end
        kp = tbl.Kp(row);
        ki = tbl.Ki(row);
        kd = tbl.Kd(row);
    end

    function mpc = readMPC()
        % Always read the CURRENT MPC edit fields (Q, R, Np), then validate them.
        mpc.Q     = edQ.Value;           % Q_MPC (tracking-error weight)
        mpc.R     = edR.Value;           % R_MPC (control-effort weight)
        mpc.Np    = edNp.Value;
        mpc.Nc    = edNc.Value;          % now editable (validated: integer, 1 <= Nc <= Np)
        mpc.TsMax = DEF.TsMax;
        % Feedback configuration, noise (standard deviations) and Kalman Filter settings.
        if strcmp(ddFeedback.Value,DEF.fbNames{2}), mpc.fb = 'kf'; else, mpc.fb = 'direct'; end
        mpc.noise.sigR = [edSigRi.Value; edSigRw.Value; edSigRt.Value];  % [A; rad/s; rad]
        mpc.measTheta  = cbMeasTheta.Value;      % true: position sensor, C = eye(3); false: C = [1 0 0; 0 1 0]
        mpc.noise.sigQ = [edSigQi.Value; edSigQw.Value; edSigQt.Value];  % [A; rad/s; rad] per sample
        mpc.noise.seed = edSeed.Value;
        mpc.kf.qScale  = edQkfScale.Value;        % Q_KF = qScale*diag(sigQ.^2)   (variances)
        mpc.kf.rScale  = edRkfScale.Value;        % R_KF = rScale*diag(sigR.^2)   (variances)
        mpc.kf.P0      = [edP0i.Value; edP0w.Value; edP0t.Value];   % diag of P0 (variances)
        mpc.kf.e0      = [edE0i.Value; edE0w.Value; edE0t.Value];   % initial estimate error
        mpc.kf.floor   = DEF.covFloor;
        mpc.MaxSteps = DEF.MaxSteps;
        % Limits are ENFORCED inside the MPC optimisation only for ticked boxes; otherwise
        % the corresponding bound is infinite here (the limit is then only monitored).
        lim = readLimits();
        if cbEnfV.Value, mpc.Vmin = lim.Vmin; mpc.Vmax = lim.Vmax;
        else,            mpc.Vmin = -Inf;     mpc.Vmax = Inf;      end
        if cbEnfI.Value, mpc.Imax = lim.Imax;
        else,            mpc.Imax = Inf;      end
        validateMPC(mpc,DEF.MaxNp);
    end

    function lim = readLimits()
        % Current values of the voltage / current limit fields (always finite, validated).
        lim.Vmin = edVmin.Value;
        lim.Vmax = edVmax.Value;
        lim.Imax = edImax.Value;
        validateLimits(lim);
    end

    function txt = mpcEnforceNote()
        % Extra sentence for the all-controllers status when MPC runs with enforced limits.
        parts = {};
        if cbEnfV.Value, parts{end+1} = 'voltage'; end
        if cbEnfI.Value, parts{end+1} = 'current'; end
        if isempty(parts)
            txt = '';
        else
            txt = sprintf([' MPC enforced its %s limit(s); P, PI, PD and PID are NOT limited ' ...
                           '(their limits are only monitored).'],strjoin(parts,' and '));
        end
    end

    function M = collectMetrics(res,mot,lim)
        % Gather the REAL simulated signals and call the single metrics function.
        opts.ssWindow  = DEF.ssWindow;
        opts.settleTol = DEF.settleTol;
        opts.tolRel    = DEF.limTolRel;
        Va = [];
        if isfield(res,'u'), Va = res.u; end
        opts.vKind = 'sampled';
        if isfield(res,'uHold'), opts.vKind = res.uHold; end
        opts.impulseV = 0;
        if isfield(res,'uImpulse'), opts.impulseV = res.uImpulse; end
        opts.minTau = NaN;
        if isfield(res,'poles') && ~isempty(res.poles) && all(isfinite(res.poles))
            pm = max(abs(res.poles));
            if pm > 0, opts.minTau = 1/pm; end
        end
        M = computeMotorMetrics(res.t,res.w,res.th,res.i,Va,mot.ref,mot.mode,lim,opts);
        if ~M.vAvail && isfield(res,'uNote') && ~isempty(res.uNote)
            M.notes{end+1} = res.uNote;
        end
    end

    function showMetrics(nm,mot,res,lim,originTxt)
        % Compute and display the metrics of the simulation that has JUST been run.
        try
            M = collectMetrics(res,mot,lim);
            ctx.name     = nm;
            if strcmp(mot.mode,'position'), ctx.modeName = 'Position Control';
            else,                           ctx.modeName = 'Speed Control'; end
            ctx.enfV     = isfield(res,'enforceV') && res.enforceV;
            ctx.enfI     = isfield(res,'enforceI') && res.enforceI;
            ctx.unstable = isfield(res,'stable') && ~res.stable;
            ctx.origin   = originTxt;
            S = metricsText(M,ctx);
            MH.valCase.Text = S.caseTxt;   MH.valSS.Text   = S.ssTxt;
            MH.valISE.Text  = S.iseTxt;    MH.valIpk.Text  = S.ipkTxt;
            MH.valJu.Text   = S.juTxt;     MH.valVrms.Text = S.vrmsTxt;
            MH.valVmax.Text = S.vmaxTxt;   MH.valImax.Text = S.imaxTxt;
            MH.valVviol.Text = S.vviolTxt; MH.valIviol.Text = S.iviolTxt;
            MH.valVviol.FontColor = violColor(S.vviolState);
            MH.valIviol.FontColor = violColor(S.iviolState);
            % Kalman-filter estimation metrics (n/a when this run has no estimator - never 0)
            E = computeEstimationMetrics(res);
            if E.avail
                MH.valRmseI.Text = sprintf('%.4e A',E.rmseI);
                MH.valRmseW.Text = sprintf('%.4e rad/s',E.rmseW);
                MH.valRmseT.Text = sprintf('%.4e rad',E.rmseT);
            else
                MH.valRmseI.Text = 'n/a (no estimator)';
                MH.valRmseW.Text = 'n/a (no estimator)';
                MH.valRmseT.Text = 'n/a (no estimator)';
            end
            if isfield(res,'fbText'), MH.valFb.Text = res.fbText; else, MH.valFb.Text = 'n/a (not MPC)'; end
            MH.valNote.Text = S.noteTxt;
            if S.noteWarn, MH.valNote.FontColor = COL.warn; else, MH.valNote.FontColor = COL.info; end
            metricsShown = true;
        catch MEm
            resetMetricsPanel(['Metrics could not be computed: ' MEm.message]);
        end
    end

    function c = violColor(state)
        % state: 0 = none, 1 = violations, -1 = unavailable
        if state > 0, c = COL.err; elseif state == 0, c = COL.ok; else, c = COL.info; end
    end

    function resetMetricsPanel(msg)
        % Clear every number (no stale metrics) but show the definition and units of the
        % CURRENT control mode and the currently selected controller.
        if nargin < 1 || isempty(msg)
            msg = 'No metrics for this selection yet - press "Run Selected Controller".';
        end
        if strcmp(modeKey(),'position')
            eDef = 'e = thref - theta';   eu = 'rad';    iu = 'rad^2.s';
        else
            eDef = 'e = wref - omega';    eu = 'rad/s';  iu = '(rad/s)^2.s';
        end
        MH.valCase.Text = sprintf('%s | %s | %s | no valid result',ddCtrl.Value,curMode(),eDef);
        MH.valSS.Text   = sprintf('-- %s',eu);
        MH.valISE.Text  = sprintf('-- %s',iu);
        MH.valIpk.Text  = '-- A';
        MH.valJu.Text   = '-- V^2.s';
        MH.valVrms.Text = '-- V';
        MH.valVmax.Text = '-- V';
        MH.valImax.Text = '-- A';
        MH.valVviol.Text = '--';
        MH.valIviol.Text = '--';
        MH.valRmseI.Text = '-- A';  MH.valRmseW.Text = '-- rad/s';  MH.valRmseT.Text = '-- rad';
        MH.valFb.Text = '--';
        MH.valVviol.FontColor = COL.info;
        MH.valIviol.FontColor = COL.info;
        MH.valNote.Text = msg;
        MH.valNote.FontColor = COL.info;
        metricsShown = false;
    end

    function markMetricsStale()
        % An input changed after the last simulation: keep the numbers (they are still the
        % true result of the previous run) but say clearly that they are out of date.
        if ~metricsShown, return; end
        MH.valNote.Text = ['OUT OF DATE: an input (motor, reference, gains, MPC, noise, Kalman filter, limits) changed after ' ...
                           'this simulation. The values above belong to the previous run - run again to update.'];
        MH.valNote.FontColor = COL.warn;
        MH.valCase.Text = [MH.valCase.Text '  [out of date]'];
        metricsShown = false;
    end

    function plotSingle(nm,mot,res)
        cc = ctrlColor(nm);
        isPos = strcmp(mot.mode,'position');
        hasEst = isfield(res,'hasEst') && res.hasEst;
        if res.stable
            suffix = sprintf(' - %s Controller%s',nm,modeTag());
            tc = [0 0 0];
        else
            suffix = sprintf(' - %s Controller%s (UNSTABLE)',nm,modeTag());
            tc = COL.err;
        end
        if hasEst, suffix = [suffix ' [actual vs KF estimate]']; end
        prepareAxes(axI, ['Armature Current i(t)' suffix],  'Current i(t) [A]',            mot.T,tc);
        prepareAxes(axTh,['Angular Position \theta(t)' suffix],'Position \theta(t) [rad]', mot.T,tc);
        prepareAxes(axW, ['Angular Speed \omega(t)' suffix],'Speed \omega(t) [rad/s]',      mot.T,tc);
        prepareAxes(axV, ['Applied Armature Voltage V_a(t)' suffix],'Voltage V_a(t) [V]',     mot.T,tc);
        cEst  = [0.000 0.447 0.741];   % Kalman-filter ESTIMATE (dashed)
        cMeas = [0.65 0.65 0.65];      % noisy MEASUREMENT (grey dots)
        if hasEst
            % noisy sensor readings first, so they stay underneath the other curves
            plot(axI,res.t,res.zi,'.','Color',cMeas,'MarkerSize',5,'DisplayName','measured i (noisy)');
            plot(axW,res.t,res.zw,'.','Color',cMeas,'MarkerSize',5,'DisplayName','measured \omega (noisy)');
            if isfield(res,'zth') && ~isempty(res.zth)
                plot(axTh,res.t,res.zth,'.','Color',cMeas,'MarkerSize',5,'DisplayName','measured \theta (noisy)');
            end
        end
        plot(axI, res.t,res.i, 'Color',cc,'LineWidth',1.8,'DisplayName',ternary(hasEst,[nm ' actual i'],nm));
        if hasEst
            plot(axI,res.t,res.ih,'--','Color',cEst,'LineWidth',1.6,'DisplayName','estimated i (KF)');
        end
        if isPos
            plot(axTh,res.t,res.th,'Color',cc,'LineWidth',1.8,'DisplayName',ternary(hasEst,[nm ' actual \theta'],[nm ' position \theta(t)']));
            if hasEst
                plot(axTh,res.t,res.thh,'--','Color',cEst,'LineWidth',1.6,'DisplayName','estimated \theta (KF)');
            end
            plot(axTh,[0 mot.T],[mot.thref mot.thref],'--','Color',[0.15 0.15 0.15], ...
                 'LineWidth',1.4,'DisplayName','Reference position');
            plot(axW, res.t,res.w, 'Color',cc,'LineWidth',1.8,'DisplayName',ternary(hasEst,[nm ' actual \omega'],[nm ' speed \omega(t)']));
            if hasEst
                plot(axW,res.t,res.wh,'--','Color',cEst,'LineWidth',1.6,'DisplayName','estimated \omega (KF)');
            end
            legend(axTh,'Location','best');
            legend(axW,'Location','best');
        else
            if hasEst
                plot(axTh,res.t,res.th,'Color',cc,'LineWidth',1.8,'DisplayName',[nm ' actual \theta']);
                plot(axTh,res.t,res.thh,'--','Color',cEst,'LineWidth',1.6,'DisplayName','estimated \theta (KF)');
                legend(axTh,'Location','best');
            else
                plot(axTh,res.t,res.th,'Color',cc,'LineWidth',1.8,'DisplayName',nm);
            end
            plot(axW, res.t,res.w, 'Color',cc,'LineWidth',1.8,'DisplayName',ternary(hasEst,[nm ' actual \omega'],[nm ' speed \omega(t)']));
            if hasEst
                plot(axW,res.t,res.wh,'--','Color',cEst,'LineWidth',1.6,'DisplayName','estimated \omega (KF)');
            end
            plot(axW,[0 mot.T],[mot.wref mot.wref],'--','Color',[0.15 0.15 0.15], ...
                 'LineWidth',1.4,'DisplayName','Reference speed');
            legend(axW,'Location','best');
        end
        if hasEst, legend(axI,'Location','best'); end
        % applied voltage = the voltage really applied to the plant (u(k) held over Ts for MPC)
        vName = 'applied V_a';
        if isfield(res,'uImpulse') && res.uImpulse ~= 0
            vName = 'applied V_a (t>0; Dirac impulse at t=0 not drawn)';
        end
        plotVoltage(axV,res,cc,vName);
        % limit lines only where MPC really ENFORCES them
        if isfield(res,'enforceV') && res.enforceV && isfinite(res.Vmin) && isfinite(res.Vmax)
            plot(axV,[0 mot.T],[res.Vmax res.Vmax],':','Color',[0.4 0.4 0.4],'LineWidth',1.3,'DisplayName','enforced V limits');
            plot(axV,[0 mot.T],[res.Vmin res.Vmin],':','Color',[0.4 0.4 0.4],'LineWidth',1.3,'HandleVisibility','off');
        end
        if isfield(res,'enforceI') && res.enforceI && isfinite(res.Imax)
            plot(axI,[0 mot.T],[res.Imax res.Imax],':','Color',[0.4 0.4 0.4],'LineWidth',1.3,'DisplayName','enforced I limit');
            plot(axI,[0 mot.T],[-res.Imax -res.Imax],':','Color',[0.4 0.4 0.4],'LineWidth',1.3,'HandleVisibility','off');
            legend(axI,'Location','best');
        end
        legend(axV,'Location','best');
    end

    function runExperiment()
        % Experiments 1-4: several MPC runs under identical motor / reference / MPC weights / limits /
        % seed, differing only in the quantity under study. Results: overlaid plots + table window.
        try
            setStatus('Running experiment...','info');
            drawnow;
            mot  = readMotorAndReference();
            lim  = readLimits();
            mpc0 = readMPC();
            expId = find(strcmp(DEF.expNames,ddExp.Value),1);
            if isempty(expId), expId = 1; end
            cases = struct('label',{},'mpc',{},'plot',{});
            switch expId
                case 1      % direct-state MPC vs KF-based MPC (same noise realisation, same everything)
                    a = mpc0;  a.fb = 'direct';
                    b = mpc0;  b.fb = 'kf';
                    cases(end+1) = makeCase('Direct-state MPC',a,true);
                    cases(end+1) = makeCase('KF-based MPC',b,true);
                    descr = 'Same motor, reference, MPC weights, constraints, process-noise realisation.';
                case 2      % measurement-noise variation (KF mode)
                    for L = 2:4
                        c = mpc0;  c.fb = 'kf';  c.noise.sigR = DEF.measPresets(L,:).';
                        cases(end+1) = makeCase(sprintf('Meas. noise %s',DEF.levelNames{L}),c,true); %#ok<AGROW>
                    end
                    descr = sprintf('KF-based MPC; process noise as set in G; Q_KF/R_KF scale as set in H. sigma_r presets: %s.', ...
                                    presetText(DEF.measPresets(2:4,:),DEF.levelNames(2:4)));
                case 3      % process-noise variation (KF mode)
                    for L = 2:4
                        c = mpc0;  c.fb = 'kf';  c.noise.sigQ = DEF.procPresets(L,:).';
                        cases(end+1) = makeCase(sprintf('Proc. noise %s',DEF.levelNames{L}),c,true); %#ok<AGROW>
                    end
                    descr = sprintf('KF-based MPC; measurement noise as set in G; Q_KF/R_KF scale as set in H. sigma_q presets: %s.', ...
                                    presetText(DEF.procPresets(2:4,:),DEF.levelNames(2:4)));
                otherwise   % combined 3 x 3 grid (KF mode); only the diagonal cases are drawn
                    for Lm = 2:4
                        for Lp = 2:4
                            c = mpc0;  c.fb = 'kf';
                            c.noise.sigR = DEF.measPresets(Lm,:).';
                            c.noise.sigQ = DEF.procPresets(Lp,:).';
                            lab = sprintf('Meas %s / Proc %s',DEF.levelNames{Lm},DEF.levelNames{Lp});
                            cases(end+1) = makeCase(lab,c,Lm == Lp); %#ok<AGROW>
                        end
                    end
                    descr = 'KF-based MPC; all 9 combinations are in the table, only the matched (diagonal) cases are plotted.';
            end
            nCase = numel(cases);
            results = cell(1,nCase);
            rows = cell(nCase,10);
            failed = {};
            for q = 1:nCase
                setStatus(sprintf('Experiment: case %d of %d (%s)...',q,nCase,cases(q).label),'info');
                try
                    r = simulateMPC(mot,cases(q).mpc,NPTS);
                    M = collectMetrics(r,mot,lim);
                    E = computeEstimationMetrics(r);
                    results{q} = r;
                    rows(q,:) = experimentRow(cases(q).label,M,E);
                catch ME2
                    failed{end+1} = sprintf('%s: %s',cases(q).label,ME2.message); %#ok<AGROW>
                    rows(q,:) = [{cases(q).label} repmat({'FAILED'},1,9)];
                end
            end
            if all(cellfun(@isempty,results))
                error('DCMotorGUI:SimFailed','Every case of the experiment failed. %s',strjoin(failed,' | '));
            end
            % ---- overlaid plots of the ACTUAL signals (and the KF estimate in experiment 1)
            ttl = sprintf(' - experiment %d%s',expId,modeTag());
            prepareAxes(axI, ['Armature Current i(t)' ttl],  'Current i(t) [A]',            mot.T,[0 0 0]);
            prepareAxes(axTh,['Angular Position \theta(t)' ttl],'Position \theta(t) [rad]', mot.T,[0 0 0]);
            prepareAxes(axW, ['Angular Speed \omega(t)' ttl],'Speed \omega(t) [rad/s]',      mot.T,[0 0 0]);
            prepareAxes(axV, ['Applied Armature Voltage V_a(t)' ttl],'Voltage V_a(t) [V]',     mot.T,[0 0 0]);
            pal = [0.635 0.078 0.184; 0.000 0.447 0.741; 0.466 0.674 0.188; 0.494 0.184 0.556; ...
                   0.850 0.325 0.098; 0.301 0.745 0.933; 0.929 0.694 0.125];
            kPlot = 0;
            for q = 1:nCase
                r = results{q};
                if isempty(r) || ~cases(q).plot, continue; end
                kPlot = kPlot + 1;
                col = pal(mod(kPlot-1,size(pal,1))+1,:);
                lab = cases(q).label;
                plot(axI, r.t,r.i, 'Color',col,'LineWidth',1.5,'DisplayName',[lab ' (actual)']);
                plot(axTh,r.t,r.th,'Color',col,'LineWidth',1.5,'DisplayName',[lab ' (actual)']);
                plot(axW, r.t,r.w, 'Color',col,'LineWidth',1.5,'DisplayName',[lab ' (actual)']);
                plotVoltage(axV,r,col,lab);
                if expId == 1 && isfield(r,'hasEst') && r.hasEst
                    plot(axI, r.t,r.ih, '--','Color',[0 0 0],'LineWidth',1.2,'DisplayName','KF estimate');
                    plot(axTh,r.t,r.thh,'--','Color',[0 0 0],'LineWidth',1.2,'DisplayName','KF estimate');
                    plot(axW, r.t,r.wh, '--','Color',[0 0 0],'LineWidth',1.2,'DisplayName','KF estimate');
                end
            end
            if strcmp(mot.mode,'position')
                plot(axTh,[0 mot.T],[mot.thref mot.thref],'--','Color',[0.15 0.15 0.15], ...
                     'LineWidth',1.4,'DisplayName','Reference position');
            else
                plot(axW,[0 mot.T],[mot.wref mot.wref],'--','Color',[0.15 0.15 0.15], ...
                     'LineWidth',1.4,'DisplayName','Reference speed');
            end
            legend(axI,'Location','best');   legend(axTh,'Location','best');
            legend(axW,'Location','best');   legend(axV,'Location','best');
            % ---- results table
            if strcmp(mot.mode,'position'), iseU = 'rad^2.s';  essU = 'rad';  trk = 'theta';
            else,                           iseU = '(rad/s)^2.s'; essU = 'rad/s'; trk = 'omega'; end
            cols = {'Case','RMSE ia [A]','RMSE omega [rad/s]','RMSE theta [rad]', ...
                    ['ISE ' trk ' [' iseU ']'],['e_ss ' trk ' [' essU ']'],'Peak |ia| [A]', ...
                    'V_RMS [V]','V violations','I violations'};
            showExperimentTable(expId,cols,rows,descr,mot,mpc0,failed);
            resetMetricsPanel(['Experiment finished - see the results table window (the metrics panel ' ...
                               'shows single runs only).']);
            if isempty(failed)
                setStatus(sprintf('Experiment %d finished (%d cases). Results table opened in a separate window.%s', ...
                                  expId,nCase,mpcEnforceNote()),'ok');
            else
                wTxt = ['Some cases FAILED and are not plotted: ' strjoin(failed,' | ')];
                setStatus(['WARNING: ' wTxt],'warn');
                uialert(fig,wTxt,'Experiment Warning','Icon','warning');
            end
        catch ME
            resetMetricsPanel('Experiment failed - there is no valid result, so no metrics are shown.');
            setStatus(['Error: ' ME.message],'error');
            uialert(fig,ME.message,'Experiment Error','Icon','error');
        end
        drawnow;
    end

    function showExperimentTable(expId,cols,rows,descr,mot,mpc0,failed)
        if ~isempty(expFig) && isvalid(expFig), delete(expFig); end
        expFig = uifigure('Name',sprintf('Experiment %d results',expId),'Position',[120 120 1180 400]);
        gl = uigridlayout(expFig,[2 1]);
        gl.RowHeight = {'1x',92};
        tb = uitable(gl,'Data',rows,'ColumnName',cols,'RowName',{},'FontSize',12);
        tb.ColumnWidth = {200,95,125,115,150,150,100,85,95,95};
        txt = sprintf(['%s Mode: %s. Np=%d, Nc=%d, Q_MPC=%.4g, R_MPC=%.4g, seed=%d. Noise values are ' ...
                       'standard deviations; RMSE = sqrt(mean((actual - estimated).^2)); "n/a" = not available ' ...
                       '(no estimator in that run), never zero. Process noise is injected in both MPC modes.'], ...
                      descr,curMode(),mpc0.Np,mpc0.Nc,mpc0.Q,mpc0.R,mpc0.noise.seed);
        if ~isempty(failed), txt = [txt ' FAILED: ' strjoin(failed,' | ')]; end
        uilabel(gl,'Text',txt,'WordWrap','on','VerticalAlignment','top','FontSize',12);
    end

    function highlightSelectedRow()
        % Visually mark the selected controller's row (cosmetic only).
        try
            tbs = {tblGains,tblGainsPos};
            for q = 1:2, removeStyle(tbs{q}); end
            selRow = find(strcmp(DEF.names,ddCtrl.Value),1);
            if isempty(selRow), return; end    % MPC has no row in the gain tables
            for q = 1:2
                addStyle(tbs{q},uistyle('BackgroundColor',[0.80 0.89 1.00], ...
                                        'FontWeight','bold'),'row',selRow);
            end
        catch
            % Styling is optional; ignore if unsupported.
        end
    end

    function setStatus(txt,kind)
        switch kind
            case 'ok',    statusLabel.FontColor = COL.ok;
            case 'warn',  statusLabel.FontColor = COL.warn;
            case 'error', statusLabel.FontColor = COL.err;
            otherwise,    statusLabel.FontColor = COL.info;
        end
        statusLabel.Text = txt;
        drawnow limitrate;
    end

end  % ======================= end of main function =======================


% ========================================================================
%                             LOCAL FUNCTIONS
% ========================================================================

function okFlag = hasControlToolbox()
% True if the Control System Toolbox functions we need are available.
okFlag = false;
try
    okFlag = license('test','Control_System_Toolbox') || license('test','Control_Toolbox');
    okFlag = okFlag && ~isempty(which('tf')) && ~isempty(which('lsim'));
catch
    okFlag = false;
end
end

function [g,hdr] = createSection(parent,row,titleText,rowHeights,colWidths,COL)
% Panel with a blue header bar; returns the inner grid (row 1 = header).
p = uipanel(parent,'BorderType','line','BackgroundColor',COL.panel);
p.Layout.Row = row;
nCols = numel(colWidths);
g = uigridlayout(p,[numel(rowHeights)+1, nCols]);
g.RowHeight       = [{26}, rowHeights];
g.ColumnWidth     = colWidths;
g.Padding         = [6 6 6 6];
g.RowSpacing      = 5;
g.ColumnSpacing   = 6;
g.BackgroundColor = COL.panel;
hdr = uilabel(g,'Text',['  ' titleText],'FontWeight','bold','FontSize',13, ...
              'FontColor',[1 1 1],'BackgroundColor',COL.blue, ...
              'HorizontalAlignment','left');
hdr.Layout.Row = 1;
% A one-column grid requires a scalar column index; [1 1] is invalid.
if nCols > 1
    hdr.Layout.Column = [1 nCols];
else
    hdr.Layout.Column = 1;
end
end

function h = addLabel(parent,txt,row,col)
h = uilabel(parent,'Text',txt,'FontSize',12,'HorizontalAlignment','left');
h.Layout.Row = row;
h.Layout.Column = col;
end

function h = addNumField(parent,val,row,col,lo,hi,loInclusive)
h = uieditfield(parent,'numeric','Value',val,'Limits',[lo hi], ...
                'LowerLimitInclusive',loInclusive, ...
                'ValueDisplayFormat','%.6g','HorizontalAlignment','right','FontSize',12);
h.Layout.Row = row;
h.Layout.Column = col;
end

function v = onoff(tf)
if tf, v = 'on'; else, v = 'off'; end
end

function T = makeGainTable(DEF,G)
% Table with columns Controller | Kp | Ki | Kd (one independent row per controller).
% G = gain matrix [Kp Ki Kd] per controller; defaults to the speed-mode gains.
if nargin < 2, G = DEF.gains; end
T = table(DEF.names(:), G(:,1), G(:,2), G(:,3), ...
          'VariableNames',{'Controller','Kp','Ki','Kd'});
end

function prepareAxes(ax,ttl,ylab,T,ttlColor)
% Clear an existing uiaxes and apply title / labels / grid (no new figures).
try, legend(ax,'off'); catch, end %#ok<NOCOM>
cla(ax);
hold(ax,'on');
grid(ax,'on');
box(ax,'on');
ax.FontSize = 12;
title(ax,ttl,'FontSize',13,'FontWeight','bold','Color',ttlColor);
xlabel(ax,'Time [s]','FontSize',12);
ylabel(ax,ylab,'FontSize',12);
xlim(ax,[0 T]);
ylim(ax,'auto');
end

function c = ctrlColor(nm)
switch nm
    case 'P',   c = [0.000 0.447 0.741];
    case 'PI',  c = [0.850 0.325 0.098];
    case 'PD',  c = [0.466 0.674 0.188];
    case 'PID', c = [0.494 0.184 0.556];
    case 'MPC', c = [0.635 0.078 0.184];
    otherwise,  c = [0 0 0];
end
end

function s = gainString(nm,kp,ki,kd)
switch nm
    case 'P',   s = sprintf('Kp=%.4g',kp);
    case 'PI',  s = sprintf('Kp=%.4g, Ki=%.4g',kp,ki);
    case 'PD',  s = sprintf('Kp=%.4g, Kd=%.4g',kp,kd);
    otherwise,  s = sprintf('Kp=%.4g, Ki=%.4g, Kd=%.4g',kp,ki,kd);
end
end

function s = mpcString(r)
% One-line description of the MPC configuration actually used.
s = sprintf('%s tracking: Q_MPC=%.4g, R_MPC=%.4g, Np=%d, Nc=%d, Ts=%.4g ms, %d samples, feedback: %s', ...
            r.mode,r.Q,r.R,r.Np,r.Nc,r.Ts*1e3,r.N,r.fbText);
end

function s = mpcKfText(r)
% Extra status sentence about the feedback configuration / Kalman filter of an MPC run.
if r.useKF
    E = computeEstimationMetrics(r);
    s = sprintf([' MPC received the Kalman-filter ESTIMATES (plant evolves with its own actual state). ' ...
                 'Estimation RMSE: ia = %.3g A, omega = %.3g rad/s, theta = %.3g rad. The pole check above ' ...
                 'refers to the state-feedback law.'],E.rmseI,E.rmseW,E.rmseT);
    if strcmp(r.mode,'position') && ~r.measTheta
        s = [s [' WARNING: theta is NOT measured and is unobservable from ia and omega, so the estimated ' ...
                'position keeps the initial error e0 and drifts; MPC is regulating the ESTIMATED position. ' ...
                'Tick "Position sensor" in section G to measure theta.']];
    end
    if r.kfSkipped > 0
        s = [s sprintf([' WARNING: %d Kalman correction step(s) were skipped (ill-conditioned innovation ' ...
                        'covariance); the filter used its prediction only at those samples.'],r.kfSkipped)];
    end
else
    s = ' MPC received the actual simulated states (no estimator).';
end
end

function s = mpcStabilityText(r)
if r.constrained
    parts = {};
    if r.enforceV
        parts{end+1} = sprintf('voltage limits [%.4g, %.4g] V',r.Vmin,r.Vmax);
    end
    if r.enforceI
        parts{end+1} = sprintf('armature-current limit |ia| <= %.4g A (at the sampling instants)',r.Imax);
    end
    s = sprintf(['Enforced inside the MPC optimisation with quadprog: %s (the classical controllers ' ...
                 'are NOT limited by these settings; their limits are only monitored). '], ...
                strjoin(parts,' and '));
    if r.enforceV
        s = [s 'The motor is open-loop stable and the input is bounded, so the response is bounded. '];
    else
        s = [s sprintf('Stability is judged from the unconstrained MPC law (largest discrete pole magnitude %.4f). ',r.rho)];
    end
    s = [s sprintf('Peak |Va| = %.4g V.',r.peakV)];
else
    s = sprintf(['Stable (largest discrete pole magnitude %.4f). ' ...
                 'No voltage or current limit enforced; peak |Va| = %.4g V.'],r.rho,r.peakV);
end
end

% ---------------------------- validation --------------------------------

function validateMPC(p,maxNp)
if ~(isnumeric(p.Q) && isscalar(p.Q) && isreal(p.Q) && isfinite(p.Q) && p.Q > 0)
    error('DCMotorGUI:InvalidInput','MPC coefficient Q must be a finite, real, strictly positive number.');
end
if ~(isnumeric(p.R) && isscalar(p.R) && isreal(p.R) && isfinite(p.R) && p.R > 0)
    error('DCMotorGUI:InvalidInput','MPC coefficient R must be a finite, real, strictly positive number.');
end
if ~(isnumeric(p.Np) && isscalar(p.Np) && isreal(p.Np) && isfinite(p.Np) && ...
     p.Np >= 1 && p.Np == round(p.Np) && p.Np <= maxNp)
    error('DCMotorGUI:InvalidInput', ...
          'MPC prediction horizon Np must be an integer between 1 and %d.',maxNp);
end
if ~(isnumeric(p.Nc) && isscalar(p.Nc) && isreal(p.Nc) && isfinite(p.Nc) && ...
     p.Nc >= 1 && p.Nc == round(p.Nc))
    error('DCMotorGUI:InvalidInput','MPC control horizon Nc must be a positive integer.');
end
if p.Nc > p.Np
    error('DCMotorGUI:InvalidInput','MPC control horizon Nc (%d) must not exceed the prediction horizon Np (%d).',p.Nc,p.Np);
end
if ~(p.Vmin < p.Vmax)   % also false for NaN
    error('DCMotorGUI:InvalidInput','Voltage limits must satisfy Vmin < Vmax.');
end
if ~(isnumeric(p.Imax) && isscalar(p.Imax) && ~isnan(p.Imax) && p.Imax > 0)   % Inf = not enforced
    error('DCMotorGUI:InvalidInput','The enforced current limit must be strictly positive.');
end
validateNoiseAndKF(p);
end

function validateNoiseAndKF(p)
% Feedback mode, noise standard deviations, seed and Kalman-filter settings.
if ~(ischar(p.fb) && any(strcmp(p.fb,{'direct','kf'})))
    error('DCMotorGUI:InvalidInput','Unknown MPC feedback configuration.');
end
okVec = @(v,n) isnumeric(v) && isreal(v) && numel(v) == n && all(isfinite(v(:)));
if ~(okVec(p.noise.sigR,3) && all(p.noise.sigR >= 0))
    error('DCMotorGUI:InvalidInput','Measurement-noise standard deviations (sigma_r) must be finite and non-negative.');
end
if ~(islogical(p.measTheta) && isscalar(p.measTheta))
    error('DCMotorGUI:InvalidInput','Invalid position-sensor setting.');
end
if ~(okVec(p.noise.sigQ,3) && all(p.noise.sigQ >= 0))
    error('DCMotorGUI:InvalidInput','Process-noise standard deviations (sigma_q) must be finite and non-negative.');
end
sd = p.noise.seed;
if ~(isnumeric(sd) && isscalar(sd) && isreal(sd) && isfinite(sd) && sd >= 0 && sd == round(sd) && sd <= 4294967295)
    error('DCMotorGUI:InvalidInput','The random seed must be an integer between 0 and 4294967295.');
end
if ~(okVec(p.kf.qScale,1) && p.kf.qScale > 0 && okVec(p.kf.rScale,1) && p.kf.rScale > 0)
    error('DCMotorGUI:InvalidInput','The Kalman-filter scales Q_KF scale and R_KF scale must be finite and strictly positive.');
end
if ~(okVec(p.kf.P0,3) && all(p.kf.P0 > 0))
    error('DCMotorGUI:InvalidInput','The initial error covariance P0 (diagonal) must be finite and strictly positive.');
end
if ~okVec(p.kf.e0,3)
    error('DCMotorGUI:InvalidInput','The initial estimate error e0 must be finite real numbers.');
end
end

function validateInputs(m)
posFields = {'La','Ra','J','Kt','Ke'};
for k = 1:numel(posFields)
    v = m.(posFields{k});
    if ~(isnumeric(v) && isscalar(v) && isreal(v) && isfinite(v) && v > 0)
        error('DCMotorGUI:InvalidInput', ...
              'Motor parameter %s must be a finite, real, strictly positive number.',posFields{k});
    end
end
if ~(isnumeric(m.B) && isscalar(m.B) && isreal(m.B) && isfinite(m.B) && m.B >= 0)
    error('DCMotorGUI:InvalidInput', ...
          'Motor parameter B must be a finite, real, non-negative number.');
end
if ~(isnumeric(m.wref) && isscalar(m.wref) && isreal(m.wref) && isfinite(m.wref))
    error('DCMotorGUI:InvalidInput','Reference speed must be a finite real number.');
end
if ~(isnumeric(m.thref) && isscalar(m.thref) && isreal(m.thref) && isfinite(m.thref))
    error('DCMotorGUI:InvalidInput','Reference position must be a finite real number.');
end
if ~(isnumeric(m.T) && isscalar(m.T) && isreal(m.T) && isfinite(m.T) && m.T > 0 && m.T <= 1000)
    error('DCMotorGUI:InvalidInput', ...
          'Simulation duration must be a finite number in the range (0, 1000] seconds.');
end
end

function validateGains(nm,kp,ki,kd)
gNames = {'Kp','Ki','Kd'};
vals   = {kp,ki,kd};
switch nm
    case 'P',   used = [true false false];
    case 'PI',  used = [true true  false];
    case 'PD',  used = [true false true];
    case 'PID', used = [true true  true];
    otherwise
        error('DCMotorGUI:InvalidInput','Unknown controller "%s".',nm);
end
for k = 1:3
    if used(k)
        v = vals{k};
        if ~(isnumeric(v) && isscalar(v) && isreal(v) && isfinite(v) && v >= 0)
            error('DCMotorGUI:InvalidInput', ...
                  '%s gain %s must be a finite, non-negative number.',nm,gNames{k});
        end
    end
end
end

% ------------------------- performance metrics ---------------------------

function M = computeMotorMetrics(t,w,th,ia,Va,ref,mode,lim,opts)
%COMPUTEMOTORMETRICS  Performance metrics from SIMULATED signals (no GUI dependency).
%
%   M = computeMotorMetrics(t,w,th,ia,Va,ref,mode,lim,opts)
%     t    simulation time vector [s] (finite, strictly increasing, may be nonuniform)
%     w    angular speed omega(t) [rad/s]        (same length as t)
%     th   angular position theta(t) [rad]       (same length as t)
%     ia   armature current [A]                  (same length as t)
%     Va   applied armature voltage [V]; [] if unavailable. Length numel(t)
%          (opts.vKind = 'sampled', trapezoidal integration) or numel(t)-1
%          (opts.vKind = 'zoh': Va(k) is held on [t(k),t(k+1)), exact integration)
%     ref  reference (scalar or vector) in the unit of the tracked output
%     mode 'speed' (e = ref - w, rad/s) or 'position' (e = ref - theta, rad)
%     lim  struct with finite Vmin, Vmax [V] and Imax [A] (violation counts are
%          reported as unavailable if a limit is missing / not finite)
%     opts optional struct: ssWindow (0.10), settleTol (1e-3), tolRel (0), vKind
%          ('sampled'), impulseV (0; weight [V.s] of a Dirac impulse in Va at t(1)),
%          minTau (NaN; fastest closed-loop time constant, only used for a note)
%
%   Unavailable metrics are NaN with the matching *Avail flag false (never zero).
%   A sample ON a limit is not a violation: Va > Vmax+tol, Va < Vmin-tol and
%   |ia| > Imax+tol are counted, with tol = opts.tolRel*max(1,|limit|).

if nargin < 9 || isempty(opts), opts = struct(); end
isModeOk = (ischar(mode) && isrow(mode)) || (isstring(mode) && isscalar(mode));
if ~isModeOk || ~any(strcmpi(char(mode),{'speed','position'}))
    M = emptyMetrics('speed');
    M.notes{end+1} = 'Unknown operating mode (expected speed or position): no metrics computed.';
    return;
end
mode  = lower(char(mode));
isPos = strcmp(mode,'position');
M = emptyMetrics(mode);
[opts,optNotes] = fillMetricOpts(opts);
M.notes = [M.notes optNotes];
M.ssFrac = opts.ssWindow;

% ---- time vector --------------------------------------------------------
if ~(isnumeric(t) && isreal(t) && numel(t) >= 2 && all(isfinite(t(:))))
    M.notes{end+1} = 'Invalid time vector (need >= 2 finite real samples): no metrics computed.';
    return;
end
t  = t(:);
dt = diff(t);
if ~all(dt > 0)
    M.notes{end+1} = 'Time vector is not strictly increasing: no metrics computed.';
    return;
end
n = numel(t);
M.nSamples = n;
M.T = t(end) - t(1);

% ---- tracking signal and reference ---------------------------------------
if isPos, y = th; else, y = w; end
yOK   = okSignal(y,n);
refOK = isnumeric(ref) && isreal(ref) && ~isempty(ref) && all(isfinite(ref(:))) && ...
        (isscalar(ref) || numel(ref) == n);

if yOK && refOK
    y = y(:);
    if isscalar(ref), r = ref*ones(n,1); else, r = ref(:); end
    e = r - y;                                    % tracking error, ref - output

    % Metric 2: ISE = integral of e^2 dt (trapezoidal rule on the actual time vector)
    ise = trapz(t,e.^2);
    if isfinite(ise)
        M.ISE = ise;  M.iseAvail = true;
    else
        M.notes{end+1} = 'ISE is not finite (overflow of a diverging response).';
    end

    % Metric 1: steady-state error from the final window (time-weighted mean)
    tStart = t(end) - opts.ssWindow*M.T;
    idx = find(t >= tStart - 1e-12*M.T);
    if numel(idx) < 3, idx = (max(1,n-2):n).'; end      % always >= 2-3 samples
    tw = t(idx);  yw = y(idx);  rw = r(idx);
    span = tw(end) - tw(1);
    yss  = trapz(tw,yw)/span;
    rss  = trapz(tw,rw)/span;
    M.yss = yss;  M.refEnd = rss;
    M.ess = rss - yss;                                  % signed  e_ss = r - y_ss
    M.essAbs = abs(M.ess);
    M.ssWindowStart = tw(1);
    M.ssPeakToPeak  = max(yw) - min(yw);
    scale = max([abs(rss), max(abs(yw)), realmin]);
    M.settled = M.ssPeakToPeak <= opts.settleTol*scale;
    M.essAvail = isfinite(M.ess);
else
    M.notes{end+1} = 'Tracking signal or reference is missing / invalid: steady-state error and ISE are unavailable.';
end

% ---- Metric 3: peak armature current -------------------------------------
if okSignal(ia,n)
    ia = ia(:);
    M.Ipeak = max(abs(ia));
    M.ipeakAvail = true;
else
    ia = [];
    M.notes{end+1} = 'Armature-current signal is missing / invalid: peak current and current violations are unavailable.';
end

% ---- Metric 4: control effort from the ACTUAL applied voltage ----------------
isZoh = strcmp(opts.vKind,'zoh');
if isZoh, nv = n-1; else, nv = n; end
impOK = isnumeric(opts.impulseV) && isscalar(opts.impulseV) && isreal(opts.impulseV) && isfinite(opts.impulseV);
v = [];
if isempty(Va)
    M.notes{end+1} = 'Applied-voltage signal is not available: control effort and voltage metrics are unavailable (not assumed).';
elseif ~okSignal(Va,nv)
    M.notes{end+1} = sprintf('Applied-voltage signal is invalid or has the wrong length (expected %d samples).',nv);
elseif ~impOK
    M.notes{end+1} = 'Impulse weight of the applied voltage is invalid: voltage metrics are unavailable.';
else
    v = Va(:);
    if isZoh
        JuReg = sum(v.^2 .* dt);                       % exact for a held (piecewise-constant) input
    else
        JuReg = trapz(t,v.^2);
    end
    if isfinite(JuReg)
        M.JuReg = JuReg;
        M.VrmsReg = sqrt(JuReg/M.T);
        M.VmaxAbsReg = max(abs(v));
        M.vImpulse = opts.impulseV;
        M.vHasImpulse = (opts.impulseV ~= 0);
        M.vZOH = isZoh;
        if M.vHasImpulse        % energy of a Dirac impulse is infinite
            M.Ju = Inf;  M.Vrms = Inf;  M.VmaxAbs = Inf;
        else
            M.Ju = JuReg;  M.Vrms = M.VrmsReg;  M.VmaxAbs = M.VmaxAbsReg;
        end
        M.vAvail = true;
    else
        v = [];
        M.notes{end+1} = 'Control effort is not finite (overflow of a diverging response).';
    end
end

% ---- Metric 5: limit violations (monitoring) -----------------------------------
if ~isempty(v) && limitOK(lim,'Vmin') && limitOK(lim,'Vmax') && lim.Vmin < lim.Vmax
    M.Vmin = lim.Vmin;  M.Vmax = lim.Vmax;
    tolHi = opts.tolRel*max(1,abs(lim.Vmax));
    tolLo = opts.tolRel*max(1,abs(lim.Vmin));
    M.nVHigh = nnz(v > lim.Vmax + tolHi);
    M.nVLow  = nnz(v < lim.Vmin - tolLo);
    M.nVViol = M.nVHigh + M.nVLow;
    M.nVSamples = numel(v);
    M.vImpulseViol = M.vHasImpulse;     % a Dirac impulse exceeds ANY finite limit
    M.vViolAvail = true;
elseif ~isempty(v)
    M.notes{end+1} = 'Voltage limits are not configured (finite Vmin < Vmax): voltage violations are unavailable.';
end
if ~isempty(ia) && limitOK(lim,'Imax') && lim.Imax > 0
    M.Imax = lim.Imax;
    tolI = opts.tolRel*max(1,abs(lim.Imax));
    M.nIViol = nnz(abs(ia) > lim.Imax + tolI);
    M.nISamples = numel(ia);
    M.iViolAvail = true;
elseif ~isempty(ia)
    M.notes{end+1} = 'Current limit is not configured (finite Imax > 0): current violations are unavailable.';
end

% ---- resolution warning ------------------------------------------------------
if isfinite(opts.minTau) && opts.minTau > 0 && max(dt) > 0.5*opts.minTau
    M.notes{end+1} = sprintf(['The time step (%.3g ms) is not small compared with the fastest closed-loop time ' ...
        'constant (%.3g ms): ISE, effort and sampled peaks are only approximate (errors of a few %%).'], ...
        1e3*max(dt),1e3*opts.minTau);
end

M.valid = M.essAvail || M.iseAvail || M.ipeakAvail || M.vAvail;
end

function M = emptyMetrics(mode)
% All metrics start as UNAVAILABLE (NaN + false flag), never as zero.
if strcmp(mode,'position')
    M.mode = 'position';  M.trackName = 'theta';  M.errUnit = 'rad';    M.iseUnit = 'rad^2.s';
else
    M.mode = 'speed';     M.trackName = 'omega';  M.errUnit = 'rad/s';  M.iseUnit = '(rad/s)^2.s';
end
M.valid = false;   M.notes = {};   M.nSamples = 0;   M.T = NaN;
M.ssFrac = NaN;
M.essAvail = false;  M.ess = NaN;  M.essAbs = NaN;  M.yss = NaN;  M.refEnd = NaN;
M.settled = false;   M.ssWindowStart = NaN;  M.ssPeakToPeak = NaN;
M.iseAvail = false;  M.ISE = NaN;
M.ipeakAvail = false;  M.Ipeak = NaN;
M.vAvail = false;    M.vZOH = false;  M.vHasImpulse = false;  M.vImpulse = 0;
M.Ju = NaN;  M.Vrms = NaN;  M.VmaxAbs = NaN;
M.JuReg = NaN;  M.VrmsReg = NaN;  M.VmaxAbsReg = NaN;
M.vViolAvail = false;  M.Vmin = NaN;  M.Vmax = NaN;
M.nVHigh = NaN;  M.nVLow = NaN;  M.nVViol = NaN;  M.nVSamples = NaN;  M.vImpulseViol = false;
M.iViolAvail = false;  M.Imax = NaN;  M.nIViol = NaN;  M.nISamples = NaN;
end

function [o,notes] = fillMetricOpts(o)
notes = {};
if ~isstruct(o), o = struct(); end
def = struct('ssWindow',0.10,'settleTol',1e-3,'tolRel',0,'vKind','sampled','impulseV',0,'minTau',NaN);
fn = fieldnames(def);
for k = 1:numel(fn)
    if ~isfield(o,fn{k}) || isempty(o.(fn{k})), o.(fn{k}) = def.(fn{k}); end
end
if ~(isnumeric(o.ssWindow) && isscalar(o.ssWindow) && isfinite(o.ssWindow) && o.ssWindow > 0 && o.ssWindow <= 1)
    o.ssWindow = def.ssWindow;
    notes{end+1} = 'Invalid steady-state window; the default 10 % of the run was used.';
end
if ~(isnumeric(o.settleTol) && isscalar(o.settleTol) && isfinite(o.settleTol) && o.settleTol > 0)
    o.settleTol = def.settleTol;
    notes{end+1} = 'Invalid settling tolerance; the default 0.1 % was used.';
end
if ~(isnumeric(o.tolRel) && isscalar(o.tolRel) && isfinite(o.tolRel) && o.tolRel >= 0)
    o.tolRel = def.tolRel;
    notes{end+1} = 'Invalid limit tolerance; strict comparison (0) was used.';
end
if ~(ischar(o.vKind) && any(strcmp(o.vKind,{'sampled','zoh'})))
    o.vKind = def.vKind;
    notes{end+1} = 'Invalid voltage-signal kind; "sampled" was assumed.';
end
end

function tf = okSignal(x,n)
% Numeric, real, finite vector of exactly n samples.
tf = isnumeric(x) && isreal(x) && isvector(x) && numel(x) == n && all(isfinite(x(:)));
end

function tf = limitOK(lim,f)
tf = isstruct(lim) && isfield(lim,f) && isnumeric(lim.(f)) && isscalar(lim.(f)) && ...
     isreal(lim.(f)) && isfinite(lim.(f));
end

function S = metricsText(M,ctx)
% Turn the metrics structure into the strings shown in the panel (unavailable -> "n/a").
if strcmp(M.mode,'position')
    eDef = 'e = thref - theta';
else
    eDef = 'e = wref - omega';
end
S.caseTxt = sprintf('%s | %s | %s | T = %.4g s',ctx.name,ctx.modeName,eDef,M.T);
if ~M.valid
    S.caseTxt = sprintf('%s | %s | %s | no valid result',ctx.name,ctx.modeName,eDef);
end
notes = {};
if ~isempty(ctx.origin), notes{end+1} = ['Metrics ' ctx.origin '.']; end
if ctx.unstable
    notes{end+1} = ['UNSTABLE closed loop: the metrics describe a diverging response, ' ...
                    'not a successful control result.'];
end

% 1. steady-state error
if M.essAvail
    if M.settled, tag = 'settled'; else, tag = 'NOT settled'; end
    S.ssTxt = sprintf('e_ss = %+.4e %s    |e_ss| = %.4e %s    [%s]',M.ess,M.errUnit,M.essAbs,M.errUnit,tag);
    if ~M.settled
        notes{end+1} = sprintf(['e_ss is the mean error over the last %.3g %% of the run, but the output ' ...
            'is still changing there (spread %.3g %s): it is not a true steady-state error.'], ...
            100*M.ssFrac,M.ssPeakToPeak,M.errUnit);
    end
else
    S.ssTxt = sprintf('n/a [%s]',M.errUnit);
end
% 2. ISE
if M.iseAvail, S.iseTxt = sprintf('%.4e %s',M.ISE,M.iseUnit); else, S.iseTxt = 'n/a'; end
% 3. peak current (also the maximum observed |ia|)
if M.ipeakAvail
    S.ipkTxt  = sprintf('%.4g A',M.Ipeak);
    S.imaxTxt = sprintf('%.4g A',M.Ipeak);
else
    S.ipkTxt = 'n/a';  S.imaxTxt = 'n/a';
end
% 4. control effort (and RMS voltage, max |Va|)
if M.vAvail
    if M.vHasImpulse
        S.juTxt   = 'Inf (impulse in Va)';
        S.vrmsTxt = 'Inf (impulse in Va)';
        S.vmaxTxt = 'Inf (impulse in Va)';
        notes{end+1} = sprintf(['Va contains a Dirac impulse of %.4g V.s at t = 0 (ideal derivative acting on ' ...
            'the reference step), so Ju, V_RMS and max|Va| are infinite. Excluding it (t > 0): ' ...
            'Ju = %.4g V^2.s, V_RMS = %.4g V, max|Va| = %.4g V.'], ...
            M.vImpulse,M.JuReg,M.VrmsReg,M.VmaxAbsReg);
    else
        S.juTxt   = sprintf('%.4e V^2.s',M.Ju);
        S.vrmsTxt = sprintf('%.4g V',M.Vrms);
        S.vmaxTxt = sprintf('%.4g V',M.VmaxAbs);
    end
else
    S.juTxt = 'n/a';  S.vrmsTxt = 'n/a';  S.vmaxTxt = 'n/a';
end
% 5. violations: counts + whether the limit is ENFORCED by the controller or only MONITORED
if M.vViolAvail
    if ctx.enfV
        how = 'ENFORCED by MPC (at sampling instants)';
    else
        how = 'MONITORED only (not enforced)';
    end
    S.vviolTxt = sprintf('%d of %d samples outside [%.4g, %.4g] V - %s',M.nVViol,M.nVSamples,M.Vmin,M.Vmax,how);
    if M.vImpulseViol
        S.vviolTxt = [S.vviolTxt '; + impulse at t = 0'];
    end
    if M.nVViol > 0 || M.vImpulseViol, S.vviolState = 1; else, S.vviolState = 0; end
    if ctx.enfV && M.nVViol > 0
        notes{end+1} = 'WARNING: the voltage limit is enforced by MPC but samples exceed it.';
    end
else
    S.vviolTxt = 'n/a'; S.vviolState = -1;
end
if M.iViolAvail
    if ctx.enfI
        how = 'ENFORCED by MPC (at sampling instants)';
    else
        how = 'MONITORED only (not enforced)';
    end
    S.iviolTxt = sprintf('%d of %d samples with |ia| > %.4g A - %s',M.nIViol,M.nISamples,M.Imax,how);
    if M.nIViol > 0, S.iviolState = 1; else, S.iviolState = 0; end
    if ctx.enfI && M.nIViol > 0
        notes{end+1} = 'WARNING: the current limit is enforced by MPC but samples exceed it.';
    end
else
    S.iviolTxt = 'n/a'; S.iviolState = -1;
end
notes = [notes M.notes];
if isempty(notes)
    S.noteTxt = ['Steady state = mean of the final ' sprintf('%.3g',100*M.ssFrac) ' % of the run. ' ...
                 'Peaks and integrals are evaluated on the simulation time grid. ' ...
                 'Limits are monitored for every controller; only MPC can enforce them.'];
    S.noteWarn = false;
else
    S.noteTxt = strjoin(notes,' ');
    S.noteWarn = true;
end
end

function validateLimits(lim)
% Finite, physically meaningful limits are required: the app never assumes "unlimited".
ok = @(v) isnumeric(v) && isscalar(v) && isreal(v) && isfinite(v);
if ~(ok(lim.Vmin) && ok(lim.Vmax) && ok(lim.Imax))
    error('DCMotorGUI:InvalidInput', ...
          ['The voltage limits (Vmin, Vmax) and the current limit must be finite numbers; ' ...
           'the application does not assume unlimited voltage or current.']);
end
if ~(lim.Vmin < lim.Vmax)
    error('DCMotorGUI:InvalidInput','Voltage limits must satisfy Vmin < Vmax.');
end
if ~(lim.Imax > 0)
    error('DCMotorGUI:InvalidInput','The maximum armature current must be strictly positive.');
end
end

% ------------------------- polynomial helpers ----------------------------

function p = trimPoly(p)
% Remove leading zeros (keeps at least one coefficient).
p = p(:).';
k = find(p ~= 0,1,'first');
if isempty(k), p = 0; else, p = p(k:end); end
end

function c = polyAdd(a,b)
% Add polynomials of different degree (descending-power coefficient vectors).
a = a(:).';  b = b(:).';
n = max(numel(a),numel(b));
c = [zeros(1,n-numel(a)) a] + [zeros(1,n-numel(b)) b];
end

function [Nc,Dc] = controllerPolynomials(nm,kp,ki,kd)
% C(s) = Nc(s)/Dc(s) exactly as specified:
%   P  : Kp            PI : Kp + Ki/s = (Kp*s + Ki)/s
%   PD : Kp + Kd*s     PID: Kp + Ki/s + Kd*s = (Kd*s^2 + Kp*s + Ki)/s
% If Ki = 0 the integral term is exactly zero, so the 1/s factor is dropped
% (this is the same transfer function, without a cancelling pole/zero at s = 0).
useInt = ki > 0;
switch nm
    case 'P'
        Nc = kp;            Dc = 1;
    case 'PI'
        if useInt, Nc = [kp ki];     Dc = [1 0]; else, Nc = kp;        Dc = 1; end
    case 'PD'
        Nc = [kd kp];       Dc = 1;
    case 'PID'
        if useInt, Nc = [kd kp ki];  Dc = [1 0]; else, Nc = [kd kp];   Dc = 1; end
    otherwise
        error('DCMotorGUI:InvalidInput','Unknown controller "%s".',nm);
end
Nc = trimPoly(Nc);
end

% ----------------------------- simulation --------------------------------

function sim = simulateController(nm,m,kp,ki,kd,nPts)
% Closed-loop simulation of current, position and speed for one controller.

[Nc,Dc] = controllerPolynomials(nm,kp,ki,kd);

% D(s) = (La*s + Ra)(J*s + B) + Kt*Ke
Dp = conv([m.La m.Ra],[m.J m.B]);
Dp(end) = Dp(end) + m.Kt*m.Ke;

% Closed-loop characteristic polynomial.
%   speed    : Dcl = D*Dc + Kt*Nc      (from 1 + P_w*C = 0)
%   position : Dcl = s*D*Dc + Kt*Nc    (from 1 + P_theta*C = 0, extra integrator in the plant)
isPos = isfield(m,'mode') && strcmp(m.mode,'position');
if isPos
    Dcl = trimPoly(polyAdd(conv(conv([1 0],Dp),Dc), m.Kt*Nc));
else
    Dcl = trimPoly(polyAdd(conv(Dp,Dc), m.Kt*Nc));
end

% Stability from closed-loop poles (all must have strictly negative real part).
p = roots(Dcl);
if isempty(p) || any(~isfinite(p))
    error('DCMotorGUI:Numerical','Could not compute closed-loop poles (non-finite values).');
end
tol    = 1e-9*max(1,max(abs(p)));
maxRe  = max(real(p));
stable = all(real(p) < -tol);

% Closed-loop transfer functions  X/Wref = P_X*C/(1 + P_w*C)  (D(s) cancelled exactly)
if isPos
    % Closed loop from the position reference (P_X*C/(1 + P_theta*C), D(s) cancelled exactly)
    Ti  = tf(trimPoly(conv(conv([1 0],[m.J m.B]),Nc)), Dcl);  % current  I/Thref
    Tw  = tf(trimPoly(conv([1 0],m.Kt*Nc)),            Dcl);  % speed    W/Thref
    Tth = tf(trimPoly(m.Kt*Nc),                        Dcl);  % position Th/Thref
else
    Ti  = tf(trimPoly(conv([m.J m.B],Nc)), Dcl);      % current  I/Wref
    Tw  = tf(trimPoly(m.Kt*Nc),            Dcl);      % speed    W/Wref
    Tth = tf(trimPoly(m.Kt*Nc),            [Dcl 0]);  % position Theta/Wref  (extra 1/s)
end

% Time grid and constant reference input (speed in rad/s, or position in rad)
t = linspace(0,m.T,nPts).';
u = m.ref*ones(size(t));

try
    yi  = lsim(Ti, u,t,'zoh');
    yw  = lsim(Tw, u,t,'zoh');
    yth = lsim(Tth,u,t,'zoh');
catch ME
    error('DCMotorGUI:SimFailed','lsim failed: %s',ME.message);
end

if ~(all(isfinite(yi)) && all(isfinite(yw)) && all(isfinite(yth)))
    error('DCMotorGUI:Numerical', ...
          ['The simulation produced non-finite values (overflow or numerical failure). ' ...
           'Check the parameters and reduce the controller gains.']);
end

% Applied armature voltage Va(t) = controller output, from the SAME closed loop
[Va,vImp,vNote] = classicalVoltage(m,Nc,Dp,Dcl,isPos,t,u);

sim.t = t;   sim.i = yi;   sim.w = yw;   sim.th = yth;
sim.u = Va;                 % applied voltage [V] on the same time grid (proper part, t > 0)
sim.uImpulse = vImp;        % weight [V.s] of a Dirac impulse in Va at t = 0 (PD / PID), else 0
sim.uHold = 'sampled';      % Va is a continuous-time signal sampled on t
sim.uNote = vNote;          % non-empty if Va could not be derived (then it is UNAVAILABLE)
sim.enforceV = false;       % classical controllers are not saturated: limits are only monitored
sim.enforceI = false;
sim.poles = p;   sim.maxRe = maxRe;   sim.stable = stable;
end

function [Va,vImp,note] = classicalVoltage(m,Nc,Dp,Dcl,isPos,t,uRef)
% Applied armature voltage of the classical loop for the constant reference m.ref:
%   speed    : Va/Wref  = C/(1 + P_w*C)     = Nc*D / Dcl
%   position : Va/Thref = C/(1 + P_theta*C) = s*D*Nc / Dcl
% (C = Nc/Dc, D(s) = Dp, Dcl = closed-loop characteristic polynomial used above.)
% For PD / PID the numerator is one degree higher than Dcl (ideal derivative on a reference
% step), so   Tu(s) = q1*s + q0 + r(s)/Dcl(s).   q1*s*(ref/s) is a Dirac impulse of weight
% q1*ref (= Kd*ref) at t = 0; the remainder is proper and is simulated with lsim.
Va = [];  vImp = 0;  note = '';
try
    if isPos
        NumU = trimPoly(conv(conv([1 0],Dp),Nc));
    else
        NumU = trimPoly(conv(Dp,Nc));
    end
    d = Dcl(:).';
    degN = numel(NumU) - 1;
    degD = numel(d) - 1;
    if degN <= degD
        Tu = tf(NumU,d);
    elseif degN == degD + 1
        % explicit polynomial division NumU = (q1*s + q0)*d + r, deg(r) < deg(d)
        nU = NumU(:).';
        q1 = nU(1)/d(1);
        R1 = nU - q1*[d 0];
        q0 = R1(2)/d(1);
        R2 = R1(2:end) - q0*d;
        r  = R2(2:end);
        Tu = tf(trimPoly(q0*d + [0 r]),d);
        vImp = q1*m.ref;
    else
        error('DCMotorGUI:Numerical','Controller output has a derivative of an impulse (not supported).');
    end
    Va = lsim(Tu,uRef,t,'zoh');
    if ~all(isfinite(Va))
        error('DCMotorGUI:Numerical','The applied-voltage simulation produced non-finite values.');
    end
catch MEv
    Va = [];  vImp = 0;
    note = ['Applied voltage of the classical controller could not be derived: ' MEv.message];
end
end

% ------------------------------- MPC -------------------------------------

function sim = simulateMPC(m,mpc,nPts)
% Receding-horizon MPC of the SAME linearized motor used by the classical controllers.
% State x = [ia; w; theta], input u = Va, zero initial state, constant reference m.wref.

% ---- sampling interval: Ts <= TsMax; for short runs Ts = T/(nPts-1) so the
%      MPC grid has the same nPts samples as the classical simulations.
TsNom = min(mpc.TsMax, m.T/(nPts-1));
N     = max(1,ceil(m.T/TsNom - 1e-9));          % number of sampling intervals
if N > mpc.MaxSteps
    error('DCMotorGUI:MPCTooLong', ...
          ['MPC would need %d samples (Ts <= %.4g ms, T = %.4g s), which exceeds the practical ' ...
           'limit of %d samples. Reduce the simulation duration to at most %.4g s.'], ...
          N, mpc.TsMax*1e3, m.T, mpc.MaxSteps, mpc.MaxSteps*mpc.TsMax);
end
Ts = m.T/N;                                     % exact grid ending at T

% ---- continuous model from the CURRENT motor parameters
A  = [-m.Ra/m.La, -m.Ke/m.La, 0;
       m.Kt/m.J,  -m.B/m.J,   0;
       0,          1,         0];
Bc = [1/m.La; 0; 0];
isPos = strcmp(m.mode,'position');
if isPos
    Cw = [0 0 1];                               % tracked output = position theta [rad]
    sim.mode = 'position';
else
    Cw = [0 1 0];                               % tracked output = speed w [rad/s]
    sim.mode = 'speed';
end

% ---- exact zero-order-hold discretization: expm([A Bc; 0 0]*Ts) = [Ad Bd; 0 1]
E  = expm([A Bc; zeros(1,4)]*Ts);
Ad = E(1:3,1:3);
Bd = E(1:3,4);
if ~(all(isfinite(Ad(:))) && all(isfinite(Bd(:))))
    error('DCMotorGUI:Numerical','Discretization of the motor model produced non-finite values.');
end

% ---- prediction matrices:  w_pred = Phi*x(k) + Gam*U,  U = [u(k) ... u(k+Nc-1)]'
Np = mpc.Np;  Nc = mpc.Nc;
[Phi,Gam] = buildPredictionMatrices(Ad,Bd,Cw,Np,Nc);

% ---- QP:  min_U  Q*||Phi*x - wref + Gam*U||^2 + R*||U||^2
%      = U'*H*U + 2*f'*U + const,   H = Q*Gam'*Gam + R*I,   f = Q*Gam'*(Phi*x - wref)
H = mpc.Q*(Gam'*Gam) + mpc.R*eye(Nc);
H = (H+H')/2;
[Lc,pd] = chol(H,'lower');
if pd ~= 0 || any(~isfinite(Lc(:)))
    error('DCMotorGUI:MPCQP','The MPC Hessian is not positive definite; check Q, R and Np.');
end
% Constraint ENFORCEMENT (inside the optimisation):  voltage  Vmin <= u(k+j) <= Vmax  and
% current  -Imax <= ia(k+j|k) <= Imax,  j = 1..Np  (current predicted at the sampling instants).
enforceV = isfinite(mpc.Vmin) || isfinite(mpc.Vmax);
enforceI = isfinite(mpc.Imax);
constrained = enforceV || enforceI;
if enforceI
    [PhiI,GamI] = buildPredictionMatrices(Ad,Bd,[1 0 0],Np,Nc);   % output = armature current ia
end
if constrained
    if isempty(which('quadprog'))
        error('DCMotorGUI:NoToolbox', ...
              ['A voltage / current limit is set to be ENFORCED, so MPC needs quadprog (Optimization ' ...
               'Toolbox), which is not available. Untick "MPC enforces V limits" / "MPC enforces I ' ...
               'limit" to run the unconstrained MPC (the limits are then only monitored).']);
    end
    lb = mpc.Vmin*ones(Nc,1);  ub = mpc.Vmax*ones(Nc,1);   % +-Inf where not enforced
    qpOpts = optimoptions('quadprog','Display','off');
end

% ---- noise realisation and Kalman-filter set-up ---------------------------------------
% Unit-normal sequences are drawn from a private stream in a FIXED order and with sizes that depend
% only on N, then scaled by the standard deviations: Mode A and Mode B (and different noise levels)
% therefore see the SAME underlying noise realisation. Process noise acts on the ACTUAL plant in both
% modes; measurement noise only matters when the Kalman filter is used.
useKF = strcmp(mpc.fb,'kf');
useTh = mpc.measTheta;                          % optional position sensor
if useTh
    mIdx = [1 2 3];   Cm = eye(3);              % measured: ia, w and theta
else
    mIdx = [1 2];     Cm = [1 0 0; 0 1 0];      % measured: ia and w (theta NOT measured, unobservable)
end
rs    = RandStream('mt19937ar','Seed',mpc.noise.seed);
Wn    = mpc.noise.sigQ(:) .* randn(rs,3,N);     % process noise  w(k), k = 0..N-1  (plant)
Vn    = mpc.noise.sigR(:) .* randn(rs,3,N+1);   % measurement noise v(k), k = 0..N (sensors); rows = ia, w, theta
kfSkipped = 0;
if useKF
    Qkf = diag(max(mpc.kf.qScale*mpc.noise.sigQ(:).^2, mpc.kf.floor));   % Q_KF (variances)
    Rkf = diag(max(mpc.kf.rScale*mpc.noise.sigR(mIdx).^2, mpc.kf.floor)); % R_KF (variances)
    P0m = diag(mpc.kf.P0(:));                                            % P0   (variances)
    XH  = zeros(3,N+1);                         % XH(:,k+1) = xhat(k|k)   (ESTIMATE)
    Zm  = zeros(numel(mIdx),N+1);               % Zm(:,k+1) = z(k)        (noisy MEASUREMENT)
    if any(~isfinite(Qkf(:))) || any(~isfinite(Rkf(:))) || any(eig(Qkf) <= 0) || any(eig(Rkf) <= 0)
        error('DCMotorGUI:KalmanCov','Q_KF or R_KF is not a valid (finite, positive definite) covariance matrix.');
    end
end

% ---- receding-horizon closed loop
X  = zeros(3,N+1);                              % ACTUAL state: X(:,k+1) = x(k)
Uv = zeros(1,N);                                % applied voltage u(k)
x  = zeros(3,1);                                % ACTUAL state; same zero initial conditions as classical runs
for k = 1:N+1
    % ---- state that MPC will use at this sampling instant
    if useKF
        z = Cm*x + Vn(mIdx,k);                  % noisy sensor reading of the ACTUAL state
        Zm(:,k) = z;
        if k == 1
            xPrior = x + mpc.kf.e0(:);          % prior xhat^-(0) = x(0) + e0 (initial estimate error)
            PPrior = P0m;                       % prior P^-(0)   = P0
        else
            xPrior = Ad*xh + Bd*Uv(k-1);        % predict with the PREVIOUS applied input u(k-1)
            PPrior = Ad*P*Ad' + Qkf;
            PPrior = (PPrior+PPrior')/2;
        end
        [xh,P,skipped] = kalmanCorrect(xPrior,PPrior,z,Cm,Rkf);
        if skipped, kfSkipped = kfSkipped + 1; end
        if any(~isfinite(xh)) || any(~isfinite(P(:)))
            error('DCMotorGUI:KalmanNumerical', ...
                  'The Kalman filter produced non-finite values at sample %d (t = %.4g s).',k,(k-1)*Ts);
        end
        XH(:,k) = xh;
        xm = xh;                                % Mode B: MPC receives the ESTIMATE
    else
        xm = x;                                 % Mode A: MPC receives the ACTUAL state
    end
    if k > N, break; end                        % last pass only produces the final estimate

    wPred0 = Phi*xm;                            % free response of the tracked output over the horizon
    f      = mpc.Q*(Gam'*(wPred0 - m.ref));
    if constrained
        if enforceI
            iPred0 = PhiI*xm;                   % free response of ia over the horizon
            Aineq  = [GamI; -GamI];
            bineq  = [mpc.Imax - iPred0; mpc.Imax + iPred0];
        else
            Aineq = [];  bineq = [];
        end
        [U,~,flag,out] = quadprog(2*H,2*f,Aineq,bineq,[],[],lb,ub,[],qpOpts);
        if flag <= 0 || isempty(U) || any(~isfinite(U))
            hint = '';
            if flag == -2
                hint = ' The enforced voltage / current limits are infeasible from this state - relax them.';
                if useKF
                    hint = [hint ' (Kalman mode: the QP uses the ESTIMATED state, which noise or the initial ' ...
                            'estimate error can push outside the limits.)'];
                end
            end
            error('DCMotorGUI:MPCQP', ...
                  'MPC optimization failed or is infeasible at sample %d (t = %.4g s): quadprog exitflag %d (%s).%s', ...
                  k, (k-1)*Ts, flag, strtrim(out.message), hint);
        end
    else
        U = -(Lc'\(Lc\f));                      % solves H*U = -f
        if any(~isfinite(U))
            error('DCMotorGUI:MPCQP', ...
                  'MPC optimization returned NaN/Inf at sample %d (t = %.4g s).',k,(k-1)*Ts);
        end
    end
    u = U(1);                                   % apply ONLY the first move
    Uv(k)  = u;
    x      = Ad*x + Bd*u + Wn(:,k);             % advance the ACTUAL plant (own state + process noise)
    if any(~isfinite(x))
        error('DCMotorGUI:Numerical', ...
              'MPC simulation produced non-finite states at sample %d (t = %.4g s).',k,k*Ts);
    end
    X(:,k+1) = x;
end

% ---- closed-loop pole check (unconstrained law). theta does not influence
%      (ia,w), so the (ia,w) block decides stability; theta is a free integrator,
%      exactly as in the classical speed-loop check.
Kx  = -(H\(mpc.Q*(Gam'*Phi)));                  % u(k) = e1'*U*(x) = Kx(1,:)*x + const
Acl = Ad + Bd*Kx(1,:);
if isPos
    rho = max(abs(eig(Acl)));                   % theta is fed back: all three poles count
else
    rho = max(abs(eig(Acl(1:2,1:2))));
end
if ~isfinite(rho), error('DCMotorGUI:Numerical','Could not compute the MPC closed-loop poles.'); end

sim.t  = (0:N).'*Ts;
sim.i  = X(1,:).';
sim.w  = X(2,:).';
sim.th = X(3,:).';
sim.u  = Uv(:);
sim.Q = mpc.Q;  sim.R = mpc.R;  sim.Np = Np;  sim.Nc = Nc;  sim.Ts = Ts;  sim.N = N;
sim.rho = rho;  sim.constrained = constrained;  sim.Vmin = mpc.Vmin;  sim.Vmax = mpc.Vmax;
sim.Imax = mpc.Imax;  sim.enforceV = enforceV;  sim.enforceI = enforceI;
sim.uHold = 'zoh';      % u(k) is held (zero-order hold) on [t(k), t(k+1)): numel(u) = numel(t)-1
sim.uImpulse = 0;       % MPC applies a bounded sampled voltage: no impulse
sim.uNote = '';
sim.peakV  = max(abs(Uv));
% bounded voltage + open-loop stable plant => bounded response (needs the VOLTAGE limit enforced)
sim.stable = enforceV || (rho < 1 - 1e-9);
sim.maxRe  = NaN;
% ---- feedback configuration and Kalman-filter results
sim.useKF = useKF;
sim.measTheta = useTh;
sim.hasEst = useKF;
sim.kfSkipped = kfSkipped;
sim.noise = mpc.noise;
if useKF
    sim.fbText = 'Kalman Filter estimates';
    sim.ih  = XH(1,:).';   sim.wh = XH(2,:).';   sim.thh = XH(3,:).';    % ESTIMATES
    sim.zi  = Zm(1,:).';   sim.zw = Zm(2,:).';                            % noisy MEASUREMENTS (ia, w)
    if useTh, sim.zth = Zm(3,:).'; else, sim.zth = []; end                % noisy theta reading (if sensor)
    sim.Qkf = Qkf;  sim.Rkf = Rkf;  sim.P0 = P0m;  sim.Pfinal = P;
else
    sim.fbText = 'Direct (actual) states';
end
end

function [xh,P,skipped] = kalmanCorrect(xm,Pm,z,Cm,Rkf)
% Kalman correction step.   K = P^- C' / (C P^- C' + R_KF),   xh = xh^- + K (z - C xh^-)
% Covariance in Joseph form, P = (I-KC) P^- (I-KC)' + K R_KF K'  (identical to (I-KC)P^- for the
% optimal K, but symmetric / positive semidefinite in finite precision). P is symmetrised and any
% negative eigenvalue (round-off only) is clamped to zero.
% If the innovation covariance S is not finite or ill-conditioned the correction is SKIPPED
% (prediction only, skipped = true); the caller counts and reports it.
skipped = false;
S = Cm*Pm*Cm' + Rkf;
S = (S+S')/2;
if any(~isfinite(S(:))) || rcond(S) < 1e-14
    xh = xm;  P = Pm;  skipped = true;
    return;
end
K  = (Pm*Cm')/S;
xh = xm + K*(z - Cm*xm);
IKC = eye(size(Pm,1)) - K*Cm;
P  = IKC*Pm*IKC' + K*Rkf*K';
P  = (P+P')/2;
[V,D] = eig(P);
d = diag(D);
if any(d < 0)
    P = V*diag(max(d,0))*V';
    P = (P+P')/2;
end
end

function E = computeEstimationMetrics(res)
% State-estimation metrics from the ACTUAL simulated state and its Kalman estimate:
%   RMSE = sqrt(mean((actual - estimated).^2)),  for ia, w and theta, over all samples.
% Unavailable (no estimator in the run) -> NaN with E.avail = false, never zero.
E.avail = false;  E.rmseI = NaN;  E.rmseW = NaN;  E.rmseT = NaN;
need = {'i','w','th','ih','wh','thh'};
if ~(isfield(res,'hasEst') && res.hasEst && all(isfield(res,need)))
    return;
end
n = numel(res.i);
sigs = {res.i,res.ih; res.w,res.wh; res.th,res.thh};
for q = 1:3
    if ~(numel(sigs{q,1}) == n && numel(sigs{q,2}) == n && all(isfinite(sigs{q,1}(:))) && all(isfinite(sigs{q,2}(:))))
        return;
    end
end
E.rmseI = sqrt(mean((res.i  - res.ih ).^2));
E.rmseW = sqrt(mean((res.w  - res.wh ).^2));
E.rmseT = sqrt(mean((res.th - res.thh).^2));
E.avail = true;
end

function plotVoltage(ax,r,cc,nm)
% Applied armature voltage actually fed to the plant: MPC u(k) is held over each Ts (staircase),
% classical Va(t) is a sampled continuous signal. Nothing is drawn if the voltage is unavailable.
if ~isfield(r,'u') || isempty(r.u), return; end
if isfield(r,'uHold') && strcmp(r.uHold,'zoh') && numel(r.u) == numel(r.t)-1
    stairs(ax,[r.t(1:end-1); r.t(end)],[r.u(:); r.u(end)],'Color',cc,'LineWidth',1.6,'DisplayName',nm);
elseif numel(r.u) == numel(r.t)
    plot(ax,r.t,r.u,'Color',cc,'LineWidth',1.6,'DisplayName',nm);
end
end

function c = makeCase(label,mpc,doPlot)
c = struct('label',label,'mpc',mpc,'plot',doPlot);
end

function row = experimentRow(label,M,E)
% One row of the experiment table; unavailable values are 'n/a' (never 0).
fN = @(av,v) ternary(av,sprintf('%.4g',v),'n/a');
fV = @(av,v) ternary(av,sprintf('%d',v),'n/a');
row = {label, fN(E.avail,E.rmseI), fN(E.avail,E.rmseW), fN(E.avail,E.rmseT), ...
       fN(M.iseAvail,M.ISE), fN(M.essAvail,M.ess), fN(M.ipeakAvail,M.Ipeak), ...
       fN(M.vAvail,M.Vrms), fV(M.vViolAvail,M.nVViol), fV(M.iViolAvail,M.nIViol)};
end

function out = ternary(cond,a,b)
if cond, out = a; else, out = b; end
end

function s = presetText(P,names)
parts = cell(1,numel(names));
for q = 1:numel(names)
    parts{q} = sprintf('%s = [%s]',names{q},strtrim(sprintf('%.3g ',P(q,:))));
end
s = strjoin(parts,', ');
end

function [Phi,Gam] = buildPredictionMatrices(Ad,Bd,Cout,Np,Nc)
% y_pred(j) = Cout*x(k+j|k):   y_pred = Phi*x(k) + Gam*U,   U = [u(k) ... u(k+Nc-1)]'
% (inputs after the control horizon Nc are held at u(k+Nc-1)).
Phi = zeros(Np,3);
Gam = zeros(Np,Nc);
Ap  = eye(3);                                   % Ad^(j-1)
h   = zeros(Np,1);                              % h(j) = Cout*Ad^(j-1)*Bd  (impulse response)
for j = 1:Np
    h(j)     = Cout*Ap*Bd;
    Ap       = Ap*Ad;                           % now Ad^j
    Phi(j,:) = Cout*Ap;
end
for j = 1:Np
    for mm = 0:j-1                              % input u(k+mm) affects y(k+j|k)
        col = min(mm,Nc-1)+1;                   % inputs beyond Nc-1 are held = u(k+Nc-1)
        Gam(j,col) = Gam(j,col) + h(j-mm);
    end
end
end
