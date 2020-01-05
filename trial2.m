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
global nodes omega_0 h C omega_d gamma_0 A phi_0 absPsi_0 dt t_axis epsilon isNoisy P;
nodes = 7;
omega_0 = [1.00, 1.40, 2.00, 0.50, 1.60, 0.30, 2.10];
h       = [0.10, 0.40, 0.80, 0.05, 1.00, 1.10, 0.35];
C       = [5.00, 15.0, 10.0, 12.0, 8.00, 4.70, 2.00];
omega_d = [1.00, 1.20, 1.50, 0.50, 0.30, 2.00, 0.80];
gamma_0 = ones(1, nodes) .* 0.30;
A = h .* (1 - 1i .* gamma_0);
phi_0 = [ 0:(nodes - 1) ] .* pi ./ 2;
absPsi_0 = [1 1 1 1 1 1 0] ./ sqrt(cols);
dt = 0.001;
t_axis = [0:5999] .* dt;
epsilon = 0.0;
P = 3;


%% function calls
global constAB
%doTrain()
% Sou Testy-C Boi
testJ = J_digit(5, true);

plot(t_axis, testJ)

deepL = zeros(rows, numbers);
deepK = zeros(rows, 1);
for digiNum = 1:numbers
    deepL(:, digiNum) = L(testJ, digiNum);
    deepK(digiNum) = K(deepL);
end

%deepK;

%% repo for old-ish code
ranks = zeros(1,10);
doPlot = false;
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

% Call This FIRST
function doTrain()
    global rows numbers t_axis constAB constE constF isNoisy P targetJ targetL prodJ;

    % OUT First Level
    isNoisy = true;
    % collect targetJ
    targetJ = J(~isNoisy);
    %collect noisyJ of P many
    noisyJ = zeros(length(t_axis), rows, numbers, P);
    constAB = zeros(2, length(t_axis), rows, numbers);
    for pp = 1:P
        noisyJ(:, :, :, pp) = J(isNoisy);
        for digiNum = 1:numbers
            for digiRow = 1:rows
                for time_k = 1:length(t_axis)
                    constAB(:, time_k, digiRow, digiNum) = constAB(:, time_k, digiRow, digiNum) + (targetJ(time_k, digiRow, digiNum) \ [ 1, noisyJ(time_k, digiRow, digiNum, pp)])';

                end
            end
        end
    end
    constAB
    
    prodJ = zeros(rows, numbers, P);
    for pp = 1:P
        noisyJ(:, :, :, pp) = J(isNoisy);
        for digiNum = 1:numbers
            for digiRow = 1:rows
                for time_k = 1:length(t_axis)
                    prodJ(digiRow, digiNum, pp) = prodJ(digiRow, digiNum, pp) + ([ 1, noisyJ(time_k, digiRow, digiNum, pp)] * constAB(time_k, digiRow, digiNum, 2)).^2;
                end
            end
        end
    end
    prodJ = prodJ ./ length(t_axis)
    

    % OUT Secong Level
    constEF = zeros(2, numbers, P);
    noisyL = zeros(rows, numbers, P);
    targetL = zeros(rows, numbers);
    for digiNum = 1:numbers
        targetL(:, digiNum) = L(targetJ(:, :, digiNum), digiNum)';
    end
    for pp = 1:P
        for digiNum = 1:numbers
            noisyL(:, digiNum, pp) = L(noisyJ(:, :, digiNum, pp), digiNum)';
            constEF(:, digiNum, pp) = (targetL(:, digiNum) \ [ ones(rows, 1), noisyL(:, digiNum, pp)]);
        end
    end

    constE = zeros(1, numbers);
    constF = zeros(1, numbers);
    for digiNum = 1:numbers
        constE(digiNum) = sum(constEF(1, digiNum, :)) ./ P;
        constF(digiNum) = sum(constEF(2, digiNum, :)) ./ P;
    end

end

% returns K(targetK)
function sum2 = K(nL)
    global rows constE constF targetL;
    sum2 = 0;
    for digiRow = 1:rows
        %sum2 = sum2 + ( targetL(digiRow) - constE(digiRow) - constF(digiRow) .* nL(digiRow) ).^2;
        sum2 = sum2 + constE(digiRow) + constF(digiRow) .* nL(digiRow);
    end
end

% returns L(targetJ_q1, noisyJ_q2)
function sum = L(nJ, tDigit)
    global t_axis constA constB targetJ;
    sum = 0;
    for time_k = 1:length(t_axis)
       %sum = sum + ( targetJ(time_k, :, tDigit) - constA(:, tDigit)' - constB(:, tDigit)' .* nJ(time_k, :) ).^2;
       sum = sum + constA(:, tDigit)' + constB(:, tDigit)' .* nJ(time_k, :);
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
        
        
    end
end

% returns psiDot(t)
function psidot = psiDot(t, psi, mu_0)
    global omega_0 gamma_0 A;
    % only t-dependence is mu
    psidot = zeros(1, length(psi));
    chempot = mu(t, mu_0);
    for m = 1:length(psi)
        psidot(m) = psi(m) .* ( 1i .* omega_0(m) - gamma_0(m) + chempot(m) ) + A(m).*sum(psi);   
        % psidot(m) = psi(m) .* ( 1i .* omega_0(m) - (abs(psi(m)).^2 .* gamma_0(m)) + chempot(m) ) + A(m).*(sum(psi) - psi(m));  
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