close all; clear; clc;

%% defining digits
global digits rows cols numbers;
numbers = 10;
rows = 10;
cols = 6;
digits = zeros(rows, cols, numbers);
LL = 0.00; gg = 3.00;  % increases readability (contrast)
    
digits(:, :, 1) =  [ LL LL gg gg LL LL ;
                     LL gg gg gg LL LL ;
                     gg gg gg gg LL LL ;
                     LL LL gg gg LL LL ;
                     LL LL gg gg LL LL ;
                     LL LL gg gg LL LL ;
                     LL LL gg gg LL LL ;
                     LL LL gg gg LL LL ;
                     LL LL gg gg LL LL ;
                     gg gg gg gg gg gg ] ;
      
digits(:, :, 2) =  [ LL gg gg gg gg LL ;
                     gg gg LL LL gg gg ;
                     LL LL LL LL gg gg ;
                     LL LL LL LL gg gg ;
                     LL LL LL LL gg gg ;
                     LL LL LL gg gg LL ;
                     LL LL gg gg LL LL ;
                     LL gg gg LL LL LL ;
                     gg gg gg gg gg gg ;
                     gg gg gg gg gg gg ] ;
                
digits(:, :, 3) =  [ LL gg gg gg gg LL ;
                     gg gg LL LL gg gg ;
                     LL LL LL LL LL gg ;
                     LL LL LL LL gg gg ;
                     LL LL gg gg gg LL ;
                     LL LL gg gg gg LL ;
                     LL LL LL LL gg gg ;
                     LL LL LL LL LL gg ;
                     gg gg LL LL gg gg ;
                     LL gg gg gg gg LL ] ;
                
digits(:, :, 4) =  [ LL LL LL LL gg LL ;
                     LL LL LL gg gg LL ;
                     LL LL gg gg gg LL ;
                     LL gg LL gg gg LL ;
                     gg LL LL gg gg LL ;
                     gg gg gg gg gg gg ;
                     LL LL LL gg gg LL ;
                     LL LL LL gg gg LL ;
                     LL LL LL gg gg LL ;
                     LL LL gg gg gg gg ] ;
                   
digits(:, :, 5) =  [ gg gg gg gg gg gg ; 
                     gg gg LL LL LL LL ;
                     gg gg LL LL LL LL ;
                     gg gg LL LL LL LL ;
                     gg gg gg gg gg LL ;
                     LL LL LL LL gg gg ;
                     LL LL LL LL LL gg ;
                     gg LL LL LL LL gg ;
                     gg gg LL LL gg gg ;
                     LL gg gg gg gg LL ] ;
            
digits(:, :, 6) =  [ LL gg gg gg gg LL ; 
                     gg gg LL LL gg gg ;
                     gg gg LL LL LL LL ;
                     gg gg LL LL LL LL ;
                     gg gg gg gg gg LL ;
                     gg gg LL LL gg gg ;
                     gg LL LL LL LL gg ;
                     gg LL LL LL LL gg ;
                     gg gg LL LL gg gg ;
                     LL gg gg gg gg LL ] ;
                         
digits(:, :, 7) =  [ gg gg gg gg gg gg ;
                     LL LL LL LL LL gg ;
                     LL LL LL LL gg gg ;
                     LL LL LL LL gg LL ;
                     LL LL LL gg LL LL ;
                     LL LL LL gg LL LL ;
                     LL LL gg LL LL LL ;
                     LL LL gg LL LL LL ;
                     LL gg LL LL LL LL ;
                     LL gg LL LL LL LL ] ;
            
digits(:, :, 8) =  [ LL gg gg gg gg LL ; 
                     gg gg LL LL gg gg ;
                     gg LL LL LL LL gg ;
                     gg gg LL LL gg gg ;
                     LL gg gg gg gg LL ;
                     LL gg gg gg gg LL;
                     gg gg LL LL gg gg ;
                     gg LL LL LL LL gg ;
                     gg gg LL LL gg gg ;
                     LL gg gg gg gg LL ] ;
            
digits(:, :, 9) =  [ LL gg gg gg gg LL ;
                     gg gg LL LL gg gg ;
                     gg LL LL LL LL gg ;
                     gg LL LL LL LL gg ;
                     gg gg LL LL gg gg ;
                     LL gg gg gg gg gg ;
                     LL LL LL LL LL gg ;
                     gg LL LL LL LL gg ;
                     gg gg LL LL gg gg ;
                     LL gg gg gg gg LL ] ;
                 
digits(:, :, 10) = [ LL gg gg gg gg LL ; 
                     LL gg LL LL gg LL ;
                     gg gg LL LL gg gg ;
                     gg LL LL gg LL gg ;
                     gg LL LL gg LL gg ;
                     gg LL gg LL LL gg ;
                     gg LL gg LL LL gg ;
                     gg gg LL LL gg gg ;
                     LL gg LL LL gg LL ;
                     LL gg gg gg gg LL ] ;
     
%% defining parameters
global nodes omega_0 h C omega_d gamma_0 A phi_0 absPsi_0 dt t_axis epsilon;
nodes = 7;
omega_0 = [1.00, 1.40, 2.00, 0.50, 1.60, 0.30, 2.10];
h       = [0.10, 0.40, 0.80, 0.05, 1.00, 1.10, 0.35];
C       = [5.00, 15.0, 10.0, 12.0, 8.00, 4.70, 2.00];
omega_d = [1.00, 1.20, 1.50, 0.50, 0.30, 2.00, 0.80];
gamma_0 = ones(1, nodes) .* 0.70;
A = h .* (1 - 1i .* gamma_0);
phi_0 = [ 0:(nodes - 1) ] .* pi ./ 2;
absPsi_0 = [1 1 1 1 1 1 0] ./ sqrt(cols);
dt = 0.001;
t_axis = [0:5999] .* dt;
epsilon = 2;


%% function calls


% Back-Propagating First Level
isNoisy = true;
% collect targetJ
time = datestr(now)
targetJ = J(~isNoisy);
%collect noisyJ of P many
P = 5;
noisyJ = zeros(length(t_axis), rows, numbers, P);
constAB = zeros(2, rows, numbers, P);
for pp = 1:P
    noiseProd = pp
    time = datestr(now)
    noisyJ(:, :, :, pp) = J(isNoisy);
    for digiNum = 1:numbers
        for digiRow = 1:rows 
            constAB(:, digiRow, digiNum, pp) = (targetJ(:, digiRow, digiNum) \ [ones(length(t_axis), 1), noisyJ(:, digiRow, digiNum, pp)])';
        end
    end
end
time = datestr(now)

global constA constB;
constA = zeros(rows, numbers);
constB = zeros(rows, numbers);
for digiNum = 1:numbers
    noiseOpt = pp
    time = datestr(now)
    for digiRow = 1:rows 
            constA(digiRow, digiNum) = sum(constAB(1, digiRow, digiNum, :)) ./ P;
            constB(digiRow, digiNum) = sum(constAB(2, digiRow, digiNum, :)) ./ P;
    end
end
time = datestr(now)


% Back-Propagating Secong Level
constEF = zeros(2, numbers, P);
noisyDeep = zeros(rows, numbers, P);
targetDeep = zeros(rows, numbers);
for digiNum = 1:numbers
    targetDeep(:, digiNum) = L(targetJ, targetJ, digiNum, digiNum)';
end
for pp = 1:P
    EFiter = pp
    time = datestr(now)
    for digiNum = 1:numbers
        noisyDeep(:, digiNum, pp) = L(targetJ(:, :, :), noisyJ(:, :, :, pp), 1, digiNum)';
        constEF(:, digiNum, pp) = (targetDeep(:, digiNum) \ [ ones(rows, 1), noisyDeep(:, digiNum, pp)])'
    end
end

global constE constF;
constE = zeros(numbers, 1);
constF = zeros(numbers, 1);
for digiNum = 1:numbers
    noiseOpt = pp
    time = datestr(now)
    constE(digiNum) = sum(constAB(1, digiNum, :)) ./ P;
    constF(digiNum) = sum(constAB(2, digiNum, :)) ./ P;
end

%

%% repo for old-ish code
ranks = zeros(1,10);
doPlot = true;
if doPlot
        for digiNum = 1:numbers
            figure(digiNum); hold on;
            for digiRow = 1:rows
                subplot(3,4,digiRow)
                plot(t_axis, targetJ(:, digiRow, digiNum));
            end
            subplot(3,4,12)
            plot(t_axis, targetJ(:, :, digiNum));
            %legend("1","2","3","4","5","6","7","8","9","0");
            ranks(digiNum) = rank(targetJ(:,:,digiNum));
        end
        ranks
end


% repo end


%% function definitions

% returns L(targetJ_q1, noisyJ_q2)
function sum = L(tJ, nJ, tDigit, nDigit)
    global t_axis constA constB;
    sum = 0;
    for time_k = 1:length(t_axis)
       sum = sum + ( tJ(time_k, :, tDigit) - constA(:, tDigit)' - constB(:, tDigit)' .* nJ(time_k, :, nDigit) ); 
    end
end

% returns J_q(t) for q = 1,2,3,...,10
function arrJ = J(isNoisy)
    global t_axis rows numbers;
    arrJ = zeros(length(t_axis), rows, numbers);
    for digiNum = 1:numbers
    %for digiNum = 10
        arrJ(:, :, digiNum) = J_digit(digiNum, isNoisy);
    end
end

% returns J_g for q = digiNum
function j_sum = J_digit(digiNum, isNoisy)
    global digits phi_0 absPsi_0 dt t_axis rows cols nodes;
    j_sum = zeros(length(t_axis), rows);
    
    digiChemPot = [ digits(:, :, digiNum), zeros(rows, 1) ];
    if isNoisy
        for digiRow = 1:rows
            digiChemPot(digiRow, :) = digiChemPot(digiRow, :) + noiseRow(nodes);
        end
    end
    
    for digiRow = 1:rows
        mu_0 = digiChemPot(digiRow, :);
        psi = zeros(length(t_axis), nodes);
        j_m7 = zeros(length(t_axis), nodes);
        
        % time_k = 1
        psi(1, :) =  absPsi_0 .* exp(1i .* phi_0);
        psi(1, :) = psi(1, :) ./ sqrt(sum(abs(psi(1, :).^2)));
        % 7th column stores the sum over other 6
        for digiCol = 1:cols
            j_m7(1, digiCol) = parCur(digiCol, nodes, psi(1, :));
        end
        j_m7(1, nodes) = sum(j_m7(1, 1:cols));
        
        for time_k = 1:(length(t_axis) - 1)
            k1 = psiDot(t_axis(time_k)        , psi(time_k)        , mu_0);
            k2 = psiDot(t_axis(time_k) + dt./2, psi(time_k) + k1./2, mu_0);
            k3 = psiDot(t_axis(time_k) + dt./2, psi(time_k) + k2./2, mu_0);
            k4 = psiDot(t_axis(time_k) + dt   , psi(time_k) + k3   , mu_0);
            psi(time_k + 1, :) = psi(time_k, :) + ( dt./6 .* (k1 + 2.*k2 + 2.*k3 + k4) );
            % scaling psi back to sum-one status (cancels magnitude slingshot)
            psi(time_k + 1, :) = psi(time_k + 1, :) ./ sqrt(sum(abs(psi(time_k + 1, :)).^2));
            
            for digiCol = 1:cols
                j_m7(time_k + 1, digiCol) = parCur(digiCol, nodes, psi(time_k + 1, :));
            end
            j_m7(time_k + 1, nodes) = sum(j_m7(time_k + 1, 1:cols));
        end
        j_sum(:, digiRow) = j_m7(:, nodes);
        
        
        % checking whether probabilities are sum-one %
        %absPsiSq = abs(psi).^2;
        %for kk = 1:length(absPsiSq(:,1))
        %    absPsiSq(kk, 8) = sum(absPsiSq(kk, :));
        %end
        
        %j_stbl = zeros(1,7);
        %stable_range = (length(t_axis) .*9 ./10):(length(t_axis) - 1);
        %for stbl = stable_range
        %    j_stbl = j_stbl + j_m7(stbl, :);
        %end
        %j_stable = j_stbl ./ length(stable_range);
        %arrJ(digiRow) = sum(j_stable);
    end
end

% returns psiDot(t)
function psidot = psiDot(t, psi, mu_0)
    global omega_0 gamma_0 A;
    % only t-dependence is mu
    psidot = zeros(1, length(psi));
    chempot = mu(t, mu_0);
    for m = 1:length(psi)
        % psidot(m) = psi(m) .* ( 1i .* omega_0(m) - gamma_0(m) + chempot(m) ) + A(m).*sum(psi);   
        psidot(m) = psi(m) .* ( 1i .* omega_0(m) - (abs(psi(m)).^2 .* gamma_0(m)) + chempot(m) ) + A(m).*(sum(psi) - psi(m));  
    end
end

% returns mu(t)
function chempot = mu(t, mu_0)
    global C omega_d;
    chempot = mu_0 .* C .* sin(omega_d .* t);
end

% returns j_mn, m,n = 1,2,...,7
function j_mn = parCur(m, n, psi)
    global A;
    j_mn = 2.*imag( A(m) .* psi(m) .* conj(psi(n)) );
end

% returns std normal noise
function epsTheta = noise()
    global epsilon;
    epsTheta = epsilon .* normrnd(0,1);
end

% returns a row of 7 noises
function nR = noiseRow(len)
    nR = zeros(1, len);
    for m = 1:len
        nR(m) = noise();
    end
end