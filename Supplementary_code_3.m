% function calAlphaBeta
%  clear

% load('E:\JJW/tif_data/RSSSM/rsssm_2003_1800_3600.mat')
load('E:\JJW/tif_data/RSSSM/rsssm_2018_1800_3600.mat')

str =["af","am","bf","bt","ta","tf","tg","trf","trg","min_all"];


for m=1: 10
    m
    load(strcat('E:\JJW\CH4_flux/alpha_value_', str(m),'.mat'));
    load(strcat('E:\JJW\CH4_flux/beta_value_', str(m),'.mat'));
    load(strcat('E:\JJW\CH4_flux/c_value_', str(m),'.mat'));

% pause 

  for p = 1:12
      p;
    c_value = eval(strcat('c_value_', str(m)));
    alpha_value = eval(strcat('alpha_value_', str(m)));
    beta_value = eval(strcat('beta_value_', str(m)));     
    ch4_flux_18(:,:,p) = real(c_value.*rsssm_2018_full(:,:,p).^alpha_value.*(1-rsssm_2018_full(:,:,p)).^beta_value);
    
    
    
  end
    ch4_flux_18(isnan(ch4_flux_18)) = 0;
    ch4_flux_18(isinf(ch4_flux_18)) = 0;
    ch4_flux_18(ch4_flux_18>=500)=500;
 
    ch4_flux = ch4_flux_18;
    sum_ch4_flux = sum(sum(ch4_flux_18))*11100^2*24*365/10^18;  %g/m2/yr-1
    mean_ch4_flux = mean(sum_ch4_flux);
    ch4_flux_18_mean = mean(ch4_flux,3);


    eval(['ch4_flux_RS18_',num2str(m),'=','ch4_flux',';']);
    save(strcat('ch4_flux_RS18_',num2str(m)),strcat('ch4_flux_RS18_',num2str(m)),'-v7.3') %variable name
    
    data(m,:) =sum_ch4_flux;

end
    xlswrite('ch4_flux_RS18',data)   


% figure(2)
% pcolor(flip(real(c_value_tf)));
% shading flat
% colorbar





