
% Author: Enes Haxhimali
% Axon Type: unmyelinated 

% Parameters

%These simulations are all done with an 1 mew m thick axon%

d_x = 10; % compartment size 
d_t = 10; % step size 

axon_length = 400;
total_time = 500; 
J = axon_length/d_x;% number of compartment    %%%%%%%%%%%%%%%

%{
d_n = 2.5 * 10^(-6); % node length
d_b = 6 *10^(-6); % buoton diameter 
d_a = 1*10^(-6); % axon diameter
myl_d = 1.5*10^(-6); %myleinated axon diameter
myl_sheath = (myl_d-d_a)/2*10^(-6); % the thickness of the myelin sheath
a = d_a/2*10^(-6);  % radius of unmyelentiated node
a_my = myl_d/2*10^(-6); % radius of myelentiated node
r_b = d_b/2*10^(-6);   % radius of buoton
%}
%%{
d_n = 2.5 * 10^(-2); % node length
d_b = 6 *10^(-2); % buoton diameter 
d_a = 1*10^(-2); % axon diameter
myl_d = 1.5*10^(-2); %myleinated axon diameter
myl_sheath = (myl_d-d_a)/2*10^(-2); % the thickness of the myelin sheath
a = d_a/2*10^(-2);  % radius of unmyelentiated node
a_my = myl_d/2*10^(-2); % radius of myelentiated node
r_b = d_b/2*10^(-2);   % radius of buoton
%}
c_m = .1 * 10^(-2);%memebrane capacitance in mu F/cm^2
c_mye = 0.04*10^(-2);%myelin capacitance in mu F/cm^2

% calculating the capacitances
C_a= 2*pi*d_x*a*c_m;
C_my = 2*pi*d_x*a*c_mye;
C_b = c_m*(pi*d_b^2/2-2*r_b^2*acos(1-2*d_a^2/r_b^2));

% conductance
g_Na_n = 740 *10^(-2);% nodal sodium conductance mS/cm^2
g_leak_n = 47*10^(-2);% nodal leak conductance mS/cm^2
g_Na_a = 105*10^(-2);% unymel axon sodium conductance mS/cm^2
g_leak_a = 4.70*10^(-2);% unmyel axon leak conductance mS/cm^2
g_Na_m = 0;% myel axon sodium conductance mS/cm^2
g_leak_m = 0.188*10^(-2);% myel axon leak conductance mS/cm^2
g_Na_b = 0;% buoton sodium conductance mS/cm^2
g_leak_b = 4.7*10^(-2);% buoton leak conductance mS/cm^2

% resistance
R_a = 70*10^(-2);% specific restitivty of the axoplasm at 37 C in omega cm
R_b = 524*10^(-2); % Resistance of bouton in K omega

E_na = 51*10^(-3); % sodium reversal potential in mV
E_l = -80*10^(-3); % leak reversal potential in mV  80
E_cl = -40*10^(-3); % chloride reversal potential in mV

I_inj = 500*10^(-12); % injected current value in pA - picoamperes 

V_values = zeros(total_time/d_t,J) ; % this will hold all the potential values at each compartment and at each unit of time
R_values = zeros(J,1) ; % this will hold all the resistance values at each compartment and at each unit of time
C_values = zeros(J,1) ; % this will hold all the capacitance values at each compartment and at each unit of time
I_values = zeros(total_time/d_t,J) ; % this will hold all the injected current values at each compartment and at each unit of time
m_values = zeros(total_time/d_t,J) ;
h_values = zeros(total_time/d_t,J) ; 

%function [V_values] = potential_solver(axon_type,g_cl,)
% -- - - - - - - - - - -  - --  -%
% setting the V and R values, the buoton is in the middle
R_values(:)  = R_a*d_x/(pi*a^2) *10^-8;
R_values(ceil(J/2)) = R_b;

C_values(:)  = C_a;
C_values(ceil(J/2)) = C_b;

%Rate function values 
    
alpha_m = @(V) (0.029* V + 10.1)./(1+exp(-0.19*V - 9.31));
m_inf = @(V) 1./(1+exp(-0.24*V - 13.44));
beta_m = @(alpha_m,m_inf) alpha_m.*(1./m_inf - 1);
h_inf = @(V) 1./(1+exp(0.1775.*V + 13.26));
beta_h = @(V) 1.25./(1+exp(-0.1.*V - 5.6));
alpha_h = @(beta_h,h_inf) (h_inf.*beta_h)./(1-h_inf);

% before beginning the loop set the values of m and h and V etc for the
% first time step
V_values(1,:) = -80 ;
m_values(1,:) = (alpha_m(V_values(1,:) )*d_t)./(1+(alpha_m(V_values(1,:)) + beta_m(alpha_m(V_values(1,:)),m_inf(V_values(1,:)) ))*d_t);
h_values(1,:) = (alpha_h(beta_h(V_values(1,:) ),h_inf(V_values(1,:)))*d_t)./(1+(alpha_h(beta_h(V_values(1,:)),h_inf(V_values(1,:))) + beta_h(V_values(1,:)))*d_t);

V_0 = V_values(1,:);
m_0 = m_values(1,:);
h_0 = h_values(1,:);

V_1 = 10*V_values(1,:);
A = zeros(J,J);
    B = zeros(J,1);
for n  = 2:total_time/d_t - 1
    disp(['timestep: ' num2str(n)])
    % we need to solve linear system of equations A*V(n) = B at each time step, where A
    % is the matrix of coefficients for the trigonal system of equations.    
   
    %{
    if n == 2
        V_0 = V_values(1,:);
        m_0 = m_values(1,:);
    end
    %}

    count = 0;
   %while abs(V_0-V_1) > 10^(-4)
        count = count +1;
        %disp('here')
        %disp(abs(V_0-V_1))
        for j = 1:J
            
            % BC1: if j = 1, v_-1 = v_1 
            % BC2: if j = J+1, v_J+1 = v_J-1
            
            % also have a condition if we are the axon or the buoton -->
            % resitance and capatiance values will change 
            
            % calculating the rate functions for this specific V

            alpha_m_current = alpha_m(V_values(n,j));
            m_inf_current = m_inf(V_values(n,j));
            beta_m_current = beta_m(alpha_m_current,m_inf_current);
            h_inf_current = h_inf(V_values(n,j));
            beta_h_current = beta_h(V_values(n,j));
            alpha_h_current = alpha_h(beta_h_current,h_inf_current);
          
    
            % if j = ceil(J/2) then buoton else 
            %g_cl = 26.6*10^(-3);
            g_cl = 0;

            g_leak_a = g_leak_a + g_cl; % adding a cl component to the leak conductane
            if g_cl ~= 0 
                E_l = E_cl;
            end
             C_values(j) = C_values(j)
             R_values = R_values
            % first solve the triagonal system of equations for U
                if j ~= 1 && j ~= J
                    L  = -2*d_t/ (C_values(j)*(R_values(j-1) + R_values(j)));  
                    U  = -2*d_t/ (C_values(j)*(R_values(j) + R_values(j+1)));
                    % assigning the coefficients
                    A(j,j+1) = U;  
                    A(j,j-1) = L;
                elseif j == 1 % BC1
                    L  = 0;  
                    U  = - 4*d_t/ (C_values(1)*(R_values(1) + R_values(2)));         
                    A(j,j+1) = U;  % assigning the coefficients
                else  % BC2
                    L  = - 4*d_t/ (C_values(J)*(R_values(J-1) + R_values(J)));  
                    U = 0;
                    A(j,j-1) = .2;%L % assigning the coefficients
                end          

                disp(['L = ' num2str(L)])
                disp(['U = ' num2str(U)])
                disp(['m0 = ' num2str(m_0(j))])
                disp(['h0 = ' num2str(h_0(j))])
                disp(['rest =' num2str(d_t/C_values(j)*(g_Na_a*real(m_0(j))^3*real(h_0(j))+g_leak_a))])       
                %disp(['rest =' num2str(d_t/C_values(j)*(g_Na_a*m_0(j)^3*h_0(j)+g_leak_a))])
        
                D = 1 - L - U + d_t/C_values(j)*(g_Na_a*real(m_0(j))^3*real(h_0(j))+g_leak_a); %sum through conductances            
                A(j,j) =  D; % assigning the coefficients 
                disp(['D = ' num2str(D)])
                if j ~= ceil(J/2) % case when not buoton
                    if  (n == 1  ) && j == 1  % we have an injected current in the first time step in first compartment
                      
                        %disp(V_0 + d_t/C_values(j) * (g_Na_a*m_0(j)^3*h_0(j)*E_na+g_leak_a*E_l) + d_t/C_values(j) * I_inj)
                        B(j) = V_0(j) + d_t/C_values(j) * (g_Na_a*real(m_0(j))^3*real(h_0(j))*E_na+g_leak_a*E_l) + d_t/C_values(j) * I_inj;
                    else
                        B(j) = V_0(j) + d_t/C_values(j) * (g_Na_a*real(m_0(j))^3*real(h_0(j))*E_na+g_leak_a*E_l) + d_t/C_values(j);            
                    end    
                else
                    B(j) = V_0(j) + d_t/C_values(j) * (g_Na_b*real(m_0(j))^3*real(h_0(j))*E_na+g_leak_b*E_l) + d_t/C_values(j);        
                end                                        
                disp(['zb = ' num2str(B(j))])
        end 
        % updating V_0
        if  count ~= 1
            V_0 = V_1;
        end
        % Now solve AV = B to get new V_1:     
        V_1 = inv(A)*B;      

        %m_0= (m_0 + alpha_m(V_1)*d_t)./(1+(alpha_m(V_1) + beta_m(V_1))*d_t);
        %h_0 = (h_0 + alpha_h(V_1)*d_t)./(1+(alpha_h(V_1) +
        %beta_h(V_1))*d_t);
       if count ==20
        disp(m_0)
       end

       m_0 =  (m_0+d_t.*alpha_m(V_1))./(1+(alpha_m(V_1) + d_t*beta_m(alpha_m(V_1),m_inf(V_1) )));
      
       h_0= (h_0+alpha_h(beta_h(V_1 ),h_inf(V_1))*d_t)./(1+(alpha_h(beta_h(V_1),h_inf(V_1)) + beta_h(V_1))*d_t);

        %m_values(n,j) = (m_values(n-1,j) + alpha_m(V_values(n,j))*d_t)/(1+(alpha_m(V_values(n,j)) + beta_m(V_values(n,j)))*d_t);
        %h_values(n,j) = (h_values(n-1,j) + alpha_h(V_values(n,j))*d_t)/(1+(alpha_h(V_values(n,j)) + beta_h(V_values(n,j)))*d_t);
           
    %end
    V_values(n,:) = V_1;
    m_values(n,:) = m_0(:,1);
    h_values(n,:) = h_0(:,1);
end

V_values = V_values;  % converting from V to mV for graphs

time = 0:d_t:total_time-d_t;
time = time/1000;
figure()
plot(time,V_values(:,1))
hold on;
plot(time,V_values(:,100/d_x))
hold on;

plot(time,V_values(:,300/d_x))
hold on;
plot(time,V_values(:,J)) 
hold on;
plot(time,V_values(:,ceil(J/2)))


title("No Chloride Conductance")
xlabel("Time (msecs)")
ylabel("Membrane Potential (mV)")
legend('1 \mum - axon', '100 \mum - axon',  '300 \mum - axon', '400 \mum - axon','200 \mum - buoton')

