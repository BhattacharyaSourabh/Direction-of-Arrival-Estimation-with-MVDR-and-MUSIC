%Sourabh Bhattacharya
%MVDR algorithm
%21EE64R18
%Defining the constants
DOA = [-30 10 50];      
datapoint = 4500;      
num_ang = length(DOA);  
num_elem = 20;          
frequency=(8*(10^9));
wavelength=(3*(10^8)/frequency);      
separation = wavelength/2;    
SNR = 50;
%%
A_steer = zeros(num_elem,num_ang);   %Steering Matrix 
for iter =1:num_ang 
    A_steer(:,iter) = exp(-1j*2*pi*separation*sind(DOA(iter))*(0:num_elem-1)'/wavelength); %Assignment matrix 
end 
noise_diag = diag(sqrt((10.^(SNR/10))/2));
signal = noise_diag * ( randn(num_ang,datapoint) + 1j*randn(num_ang,datapoint) );
white_noise = sqrt(1/2)*(randn(num_elem,datapoint)+1j*randn(num_elem,datapoint));
X_orig = A_steer*signal; 
X_orig = X_orig+white_noise; 
%% MVDR (Capon)
tot_ang = -90:0.1:90;   
cov_mat = cov(X_orig');                     
inv_cov_mat = cov_mat^(-1); 
for count=1:length(tot_ang) 
    mat = zeros(num_elem,1); 
    mat = exp(-1j*2*pi*separation*(0:num_elem-1)'*sind(tot_ang(count))/wavelength);
    mat2 = mat'*inv_cov_mat*mat;
    iter_mvdr(count) = 1/ mat2; 
end
iter_mvdr = real(10*log10(iter_mvdr));
[pks,locs] = findpeaks(iter_mvdr,tot_ang,'SortStr','descend','Annotate','extents');
MVDR_net = sort(locs(1:num_ang))

%%
figure
plot(tot_ang,iter_mvdr,'-b',locs(1:num_ang),pks(1:num_ang),'r*'); hold on
xlim([min(tot_ang) max(tot_ang)])
xlabel('Angle \theta (degree)'); ylabel('Spatial Power Spectrum P(\theta) (dB)') 
title('DOA estimation based on MVDR algorithm ') 
xlim([min(tot_ang) max(tot_ang)])
grid on