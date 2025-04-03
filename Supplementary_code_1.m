% function calAlphaBeta
clear
options = optimset('Display','off');
load('E:\JJW\CH4_flux/class_tag.mat')
load('../tif_data/soilgris_BDSSC/dry_bd_1800_3600.mat')
load('../tif_data/soilgris_BDSSC/correct_sand_1800_3600.mat')
load('../tif_data/soilgris_BDSSC/correct_silt_1800_3600.mat')
load('../tif_data/soilgris_BDSSC/correct_clay_1800_3600.mat')
load('../tif_data/soilgris_BDSSC/correct_soc_1800_3600.mat')
load('../tif_data/soilgrids_pHTN/ph_wgs_2021_1800_3600.mat')
load('../tif_data/soilgrids_pHTN/TN_wgs_2021_1800_3600.mat')

m = 1800;
n = 3600;
theta_opt = ones(m,n)*nan;
f_max = ones(m,n)*nan;
f_half = ones(m,n)*nan;
c_value_af = ones(m,n)*nan;%
alpha_value_af = ones(m,n)*nan;%
beta_value_af = ones(m,n)*nan;%

flag = ones(m,n);
flag(isnan(dry_bd)) = nan; 
flag(isnan(correct_sand)) = nan; 
flag(isnan(correct_silt)) = nan; 
flag(isnan(correct_clay)) = nan; 
flag(isnan(soc_full)) = nan; 
flag(isnan(ph_full)) = nan; 
flag(isnan(TN_full)) = nan; 

flag(class_tag ~= -1220) = nan;  %af
for i=1:m
    i
    for j=1:n
        j;
        if(~isnan(flag(i,j)))
                theta_opt(i,j) = -0.0239* soc_full(i,j) +0.3899;% theta_opt
                    if theta_opt(i,j)<0 
                        theta_opt(i,j) = 0;
                    end
                    if  theta_opt>1
                        theta_opt(i,j) = 1;
                    end
                f_max(i,j) = (-8.9256)* ph_full(i,j)+105.7;% Fmax
                    if f_max(i,j)<0 
                        f_max(i,j) = 0;
                    end
                f_half(i,j) = (-9.1964)* ph_full(i,j) +87.83;% F0.5
                    if f_half(i,j)<0 
                        f_half(i,j) = 0;
                    end        
                
                fun = @(x) [x(2)/(x(2)+x(3))-theta_opt(i,j);
                    x(1)*(x(2)/(x(2)+x(3)))^x(2)*(x(3)/(x(2)+x(3)))^x(3)-f_max(i,j);
                   x(1)*0.5^(x(2)+x(3))-f_half(i,j)];
                x = fsolve(fun,[1000,2,4],options);
            
                c_value_af(i,j) = x(1);%
                alpha_value_af(i,j) = x(2);%
                beta_value_af(i,j) = x(3);%
             end    
        end
    end


save('alpha_value_af','alpha_value_af','-v7.3') %variable name
save('beta_value_af','beta_value_af','-v7.3')  %variable name
save('c_value_af','c_value_af','-v7.3')  %variable name  
          
                
