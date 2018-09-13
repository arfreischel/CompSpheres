function Exp11_100_LogisticCurveFit()

%need to go through the red data step by step to get all 4 to plot
 
clear all
close all
clc

format long

data_11= xlsread('Exp11&12 Matlab.xlsx', 'exp 11');

[n,m]= size(data_11); 
time = data_11(1,5:end);

mean_vectorGREEN = zeros(6,4); 
mean_vectorRED = zeros(6,4);


num_data_chunks = 4;

c= ['k' 'b' 'r' 'g']
 A= zeros(8,3);
 B= zeros(8,3);
 C= zeros(8,3);
 D= zeros(8,3);
 
 aa= zeros(8,3);
 bb= zeros(8,3);
 cc= zeros(8,3);
 dd= zeros(8,3);
 
for i = 1:1:4
    %grab a data chunk
    first_row = (i-1)*25+1;
    second_row = first_row +24;
    
    data_chunkGREEN1 = getDataChunk(data_11,i); %grab specified rows
    data_chunkRED1 = getDataChunk(data_11,i+4);
    
    size(data_chunkGREEN1);
    size(data_chunkRED1);
    
    mean_vector_replicates_green = zeros(6,4);
    mean_vector_replicates_red = zeros(6,4);
      
    [a,b]= size(data_chunkGREEN1);
    
    %cc= c(i)
     
    
        
    for exp_num = 1:1:1
        [dataSeriesGREEN1, dataSeriesRED1] = getDataSeries1(data_chunkGREEN1, data_chunkRED1,exp_num)
        
%         mean_vectorGREEN1(exp_num,i) = getMeanLastTen(dataSeriesGREEN1)    
%         mean_vectorRED1(exp_num,i) = getMeanLastTen(dataSeriesRED1);   
%         
        for n= 1:1:4
        
           
            
            figure
                
        a= axes;
        t=title ({'experiment 11 (+glutamine)';'';''});
        a.Visible = 'off';
        t.Visible = 'on'; 
        
        xlim ([0 546]);
        xticks([0: 100: 546]);
        
        
 
        for q= 1:1:4       
%             
%             nan_indices = isnan(dataSeriesGREEN1(q,:)); %this part deals with the skipped data
%             dataSeriesGREEN1_clean = dataSeriesGREEN1(q,~nan_indices);
%             time_clean = time(~nan_indices); 
%                        
            
            y= @(x,time) ((exp(x(1).*time))*x(2)*x(3))./(x(3)-x(2)+(exp(x(1).*time)*x(2))); %x(1)=r ;x(2)=n0 ;x(3)=K
             
            %SSE= sum(((dataSeriesGREEN1(q,:)-y(x,time))).^2)+1/ epsi*sum(x>=[2 50])*sum(dataSeriesGREEN1(:,q))+1/ epsi*sum(x<=[0 1])*sum(dataSeriesGREEN1(:,q));            
                        
            %OLS= @(x) sum(((dataSeriesGREEN1(q,:)-y(x,time))).^2)+1/ epsi*sum(x>=[0.1 100])*sum(dataSeriesGREEN1(:,q))+1/ epsi*sum(x<=[0 0.001])*sum(dataSeriesGREEN1(:,q));  
            %OLS= @(x) sum((dataSeriesGREEN1(q,:)-y(x,time)).^2);
            OLS= @(x) sum((y(x,time)-dataSeriesGREEN1(q,:)).^2);
            %opts= optimset('MaxFunEvals', 1e25, 'MaxIter', 1e4);
            opts= optimset()
            
                      %calculate a GR ::            
            t= time(1,1:13)
            f= log(dataSeriesGREEN1(q,1:13))%72 hours
            pGreen= polyfit(t,f,1)                                      
            GR_green= pGreen(1)
%             v= polyval(pGreen, t)
%             plot(t,f,'o',t,v, '-')

switch n
       case 1
           %%read-in parameters
           X= fminsearch(OLS, [GR_green dataSeriesGREEN1(q,1)  mean(dataSeriesGREEN1(q,end-10:end))], opts)
           F= [GR_green dataSeriesGREEN1(q,1)  mean(dataSeriesGREEN1(q,end-10:end))]
       case 2
           %%initial guess
           X= fminsearch(OLS, [0.005 5 30], opts)
           F= [0.005 5 30]
       case 3
            %%r and N0 only
            X= fminsearch(OLS, [GR_green dataSeriesGREEN1(q,1) 30], opts)
            F= [GR_green dataSeriesGREEN1(q,1) 30]
       case 4
           %%K & N0
           X= fminsearch(OLS, [0.005 dataSeriesGREEN1(q,1) mean(dataSeriesGREEN1(q,end-10:end))], opts) 
           F= [0.005 dataSeriesGREEN1(q,1) mean(dataSeriesGREEN1(q,end-10:end))]
end
       %try printing them after the switch-->NOPE!
           %xlswrite ('SpheroidParameters.xlsx', X(n))
         
if n==1
    A(q,1:3)= X
    aa(q,1:3)= F
end
if n==2
    B(q,1:3)= X
    bb(q, 1:3)= F
end
if n==3
    C(q,1:3)= X
    cc(q, 1:3)= F
end
if n==4
    D(q,1:3)= X
    dd(q, 1:3)= F
end
           
           
           
           logiFit= y(X,time);  
                        
             R=corrcoef(dataSeriesGREEN1(q,:), y(X,time), 'Rows','complete')
             Rsq= (R(1,2).^2)
            
            
            subplot(1,2,1);%exp_num
            hh= plot(time, dataSeriesGREEN1(q,:),'.','color', c(q))%,'.', 'color', c(q)); %'LineWidth', 1);
            hold on
            
            plot(time,logiFit,'color',c(q)); 
            set(get(get(hh,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            
            
            hold on
  
            
            %string{q}=strcat('$y=',num2str(X(1),3),'{\left( \frac{x}{',num2str(X(1),3),'}\right)^{e^{-',num2str(X(3),3),'t}}}',  ';', 'R^2=',num2str(Rsq),'$');
            string{q}=strcat('$y=\frac{ \left(x',num2str(X(3),3),'e^{',num2str(X(1),3),'t}\right)}{\left(',num2str(X(3),3),'-x+ x e^{',num2str(X(1),3),'t}\right)}',';', 'R^2=',num2str(Rsq,4),'$');
            h = legend(string);set(h,'Interpreter','Latex');
          q
          string
            
            title('Logistic Curve, MCF-7')
            ylim([ 0 100])
            
            hold on
        end
        
        %
        
           %xlswrite ('SpheroidParameters.xlsx', X(1))
               
        xlim ([0 546]);
        xticks([0: 100: 546]);
         
        % once for every experiment; 4 experiments per plot
        string = {};
        for qq= 1:1:4
             
                     
            %y= @(r,time) ((exp(r.*time))*N_0*K)./(K-N_0+(exp(r.*time)*N_0)); %when fitting a system of ode's, create system of functions then call that function script instead of single ode
             y= @(x,time) ((exp(x(1).*time))*x(2)*x(3))./(x(3)-x(2)+(exp(x(1).*time)*x(2))); %x(1)=r ;x(2)=n0 ;x(3)=K
            %epsi= 1e-10    
            
            %OLS= @(x) sum(((dataSeriesRED1(q,:)-y(x,time))).^2)+1/ epsi*sum(x>=[1 30])*sum(dataSeriesRED1(:,q))+1/ epsi*sum(x<=[0 0.001])*sum(dataSeriesRED1(:,q));
            OLS= @(x) sum((y(x,time)-dataSeriesRED1(qq,:)).^2);
            
            %opts= optimset('MaxFunEvals', 1e25, 'MaxIter', 1e4);
            opts= optimset()
            
              %calculate a GR ::            
            t= time(1,1:13)
            f= log(dataSeriesRED1(qq,1:13))
            pRed= polyfit(t,f,1)                                      
            GR_Red= pRed(1)
%             v= polyval(pRed, t)
%              plot(t,f,'o',t,v, '-')
        

%for n=1:1:4
   switch n
       case 1
           %%read-in parameters
           X= fminsearch(OLS, [GR_Red dataSeriesRED1(qq,1) mean(dataSeriesRED1(qq, end-10:end))], opts)
           F= [GR_Red dataSeriesRED1(qq,1) mean(dataSeriesRED1(qq, end-10:end))]
       case 2
           %%initial guess
           X= fminsearch(OLS, [0.005 3 10], opts)
           F= [0.005 3 10]
       case 3
            %X
            %r and N0 only
            X= fminsearch(OLS, [GR_Red dataSeriesRED1(qq,1) 10], opts)
            F= [GR_Red dataSeriesRED1(qq,1) 10]
       case 4
           %%K & N0
            X= fminsearch(OLS, [0.0001 dataSeriesRED1(qq,1) mean(dataSeriesRED1(qq,end-10:end))], opts)
            F= [0.0001 dataSeriesRED1(qq,1) mean(dataSeriesRED1(qq,end-10:end))]
   end
 %  end

if n==1
    A(qq+4,1:3)= X
    aa(qq+4, 1:3)= F
end
if n==2
    B(qq+4,1:3)= X
    bb(qq+4, 1:3)= F
end
if n==3
    C(qq+4,1:3)= X
    cc(qq+4,1:3)= F
end
if n==4
    D(qq+4,1:3)= X
    dd(qq+4, 1:3)= F
end

     logiFit= y(X,time);  

               R=corrcoef(dataSeriesRED1(qq,:), y(X,time))
               Rsq= (R(1,2).^2)
            
            
            subplot(1,2,2);%exp_num
            hh= plot(time, dataSeriesRED1(qq,:),'.','color', c(qq))%,'.', 'color', c(q)); %'LineWidth', 1);
            hold on
            
            plot(time,logiFit,'color',c(qq)); 
            set(get(get(hh,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            
          
            string{qq}=strcat('$y=\frac{ \left(x',num2str(X(3),3),'e^{',num2str(X(1),3),'t}\right)}{\left(',num2str(X(3),3),'-x+ x e^{',num2str(X(1),3),'t}\right)}',';', 'R^2=',num2str(Rsq,4),'$');
            h = legend(string);set(h,'Interpreter','Latex');            
          
            title('Logistic Curve, MDA-MB-231')

            
            ylim([ -5 16])
            %hold off
            hold on    
        end
        
A
B



        end
                           
    end
    
if i==1
  xlswrite('SpheroidParameters', A, 'Logistic11', 'B2:D9')
  xlswrite('SpheroidParameters', aa, 'Logistic11', 'B11:D18')
  xlswrite('SpheroidParameters', B, 'Logistic11', 'B20:D27')
  xlswrite('SpheroidParameters', bb, 'Logistic11', 'B29:D36')
  xlswrite('SpheroidParameters', C, 'Logistic11', 'B38:D45')
  xlswrite('SpheroidParameters', cc, 'Logistic11', 'B47:D54')
  xlswrite('SpheroidParameters', D, 'Logistic11', 'B56:D63')
  xlswrite('SpheroidParameters', dd, 'Logistic11', 'B65:D72')
end
if i==2
 xlswrite('SpheroidParameters', A, 'Logistic11', 'F2:H9')
 xlswrite('SpheroidParameters', aa, 'Logistic11', 'F11:H18')
 xlswrite('SpheroidParameters', B, 'Logistic11', 'F20:H27')
 xlswrite('SpheroidParameters', bb, 'Logistic11', 'F29:H36')
 xlswrite('SpheroidParameters', C, 'Logistic11', 'F38:H45')
 xlswrite('SpheroidParameters', cc, 'Logistic11', 'F47:H54')
 xlswrite('SpheroidParameters', D, 'Logistic11', 'F56:H63')
 xlswrite('SpheroidParameters', dd, 'Logistic11', 'F65:H72')
end

if i==3
 xlswrite('SpheroidParameters', A, 'Logistic11', 'J2:L9')
 xlswrite('SpheroidParameters', aa, 'Logistic11', 'J11:L18')
 xlswrite('SpheroidParameters', B, 'Logistic11', 'J20:L27')
 xlswrite('SpheroidParameters', bb, 'Logistic11', 'J29:L36')
 xlswrite('SpheroidParameters', C, 'Logistic11', 'J38:L45')
 xlswrite('SpheroidParameters', cc, 'Logistic11', 'J47:L54')
 xlswrite('SpheroidParameters', D, 'Logistic11', 'J56:L63')
 xlswrite('SpheroidParameters', dd, 'Logistic11', 'J65:L72')
end

if i==4
 xlswrite('SpheroidParameters', A, 'Logistic11', 'N2:P9')
 xlswrite('SpheroidParameters', aa, 'Logistic11', 'N11:P18')
 xlswrite('SpheroidParameters', B, 'Logistic11', 'N20:P27')
 xlswrite('SpheroidParameters', bb, 'Logistic11', 'N29:P36')
 xlswrite('SpheroidParameters', C, 'Logistic11', 'N38:P45')
 xlswrite('SpheroidParameters', cc, 'Logistic11', 'N47:P54')
 xlswrite('SpheroidParameters', D, 'Logistic11', 'N56:P63')
 xlswrite('SpheroidParameters', dd, 'Logistic11', 'N65:P72')
end
    
    
end
       

end




function data_chunk1 = getDataChunk(data_11,i)
    %grab a data chunk
    first_row = (i-1)*25+1;
    second_row = first_row +24;
    
    size(data_11);
    
    data_chunk1 = data_11(first_row + 1:second_row, 5:end);    
    
    size(data_chunk1);
end

function [dataSeriesGREEN, dataSeriesRED] = getDataSeries1(dataGREEN, dataRED,exp_num)
    
    [n,m] = size(dataGREEN);
    dataSeriesGREEN = zeros(4,m);
    dataSeriesRED = zeros(4,m);

    num_replicates = 4;
    for k = 1:1:num_replicates
        exp_num;
        %tell periodic rows
        index = (k-1)*6+exp_num ;
        indexRed = (k-1)*6+ (6-exp_num+1);
        
        dataSeriesGREEN(k,:) = dataGREEN(index,:);
        dataSeriesRED(k,:) = dataRED(indexRed,:);
        
    end
   
end


% 
% function meanFirstThree= meanFirstThree(data)
%     meanFirstThree= mean(data(:,1:3));
% end
% 
% function value = getMeanLastTen(data)
%     value = mean(mean(data(:,(end-10):end)));
% end