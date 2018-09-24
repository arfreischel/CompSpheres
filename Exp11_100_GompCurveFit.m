function Exp11_100_GompCurveFit()

clear all
close all
clc

format long

data_11= xlsread('Exp11&12 Matlab.xlsx', 'exp 11');

[n,m]= size(data_11); 
time = data_11(1,5:end);

num_data_chunks = 4;

c= ['k' 'b' 'r' 'g']
 A= zeros(8,3);
 B= zeros(8,3);
 C= zeros(8,3);
 D= zeros(8,3);
 
 R_squaredA= zeros(8,1);
 R_squaredB= zeros(8,1);
 R_squaredC= zeros(8,1);
 R_squaredD= zeros(8,1);
 
 aa= zeros(8,3);
 bb= zeros(8,3);
 cc= zeros(8,3);
 dd= zeros(8,3);
 
for i = 1:1:4%num_data_chunks
    %grab a data chunk
    first_row = (i-1)*25+1;
    second_row = first_row +24;
    
    data_chunkGREEN1 = getDataChunk(data_11,i); %grab specified rows
    data_chunkRED1 = getDataChunk(data_11,i+4);
    
    size(data_chunkGREEN1);
    size(data_chunkRED1);
      
    [a,b]= size(data_chunkGREEN1);
           
    for exp_num = 1:1:1
        [dataSeriesGREEN1, dataSeriesRED1] = getDataSeries1(data_chunkGREEN1, data_chunkRED1,exp_num)
        
        dataSeriesGREEN1(dataSeriesGREEN1 < 0) = NaN;
        
        for n= 1:1:4
            
            figure
            
            a= axes;
            t=title ({'experiment 11 (+glutamine)';'';''});
            a.Visible = 'off';
            t.Visible = 'on'; 
        
            xlim ([0 546]);
            xticks([0: 100: 546]);
        
         for q= 1:1:4       
                     
            %N(t)= K+(N-K)exp(-r*t)
            y= @(x, time) x(1)*(x(2)/x(1)).^(exp(-x(3).*time)); %x(1)= K x(2)= N0 x(3)= r
            OLS= @(x) sum((y(x,time)-dataSeriesGREEN1(q,:)).^2);
            %OLS= @(x) sum(((dataSeriesGREEN1(q,:)-y(x,time))).^2)+1/ epsi*sum(x>=[0.1 100])*sum(dataSeriesGREEN1(:,q))+1/ epsi*sum(x<=[0 0.001])*sum(dataSeriesGREEN1(:,q));
            %opts= optimset('MaxFunEvals', 1e25, 'MaxIter', 1e4);
            opts= optimset();
            %R= fminsearch(OLS,rand(1,1),opts);
            
                     %calculate a GR ::            
            t= time(1,1:13)
            f= log(dataSeriesGREEN1(q,1:13))%72 hours
            pGreen= polyfit(t,f,1)                                      
            GR_green= pGreen(1)
%             v= polyval(pGreen, t)
%             plot(t,f,'o',t,v, '-')
            
            
            %%read-in parameters
                %X= fminsearch(OLS, [mean(dataSeriesGREEN1(q,end-10:end)) dataSeriesGREEN1(1) GR_green], opts)
            %%initial guess
                %X= fminsearch(OLS, [30 5 0.005], opts)
            %%r and N0 only
                %X= fminsearch(OLS, [30 dataSeriesGREEN1(1) GR_green], opts)
            %%K & N0
                X= fminsearch(OLS, [mean(dataSeriesGREEN1(q,end-10:end)) dataSeriesGREEN1(1) 0.005], opts)
            
            GompFit= y(X,time);  
            leastSquares= sum((y(X,time)-dataSeriesGREEN1(q,:)).^2);
            %rsqr= 1-((SSE)/(SST))  Rsq1 = 1 - sum((y - yCalc1).^2)/sum((y - mean(y)).^2)
            %need this Rsq on the chart, also why are these so weird? 
            %Rsq1 = 1 - sum((dataSeriesGREEN1(q,:) - expFit).^2)/sum((dataSeriesGREEN1(q,:) - mean(dataSeriesGREEN1)).^2)
            
             R=corrcoef(dataSeriesGREEN1(q,:), y(X,time), 'Rows','complete')
             Rsq= (R(1,2).^2)
                        
            subplot(1,2,1);%exp_num
            hh= plot(time, dataSeriesGREEN1(q,:),'.','color', c(q))%,'.', 'color', c(q)); %'LineWidth', 1);
            hold on
            
            plot(time,GompFit,'color',c(q)); 
            set(get(get(hh,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            
            
           
            %x(1)= K x(2)= N0 x(3)= r             
            %string{q}=strcat('$y= e^{',num2str(X(1)),'t}+',num2str(X(2)),';', 'R^2=',num2str(Rsq),'$');    
            %string{q}= strcat('$y=', num2str(X(1),3),'\left( \frac{x}{',num2str(X(1),3),'}\right)^e^{',num2str(X(3),3),'t}',';''R^2=',num2str(Rsq),'$'); 
            string{q}=strcat('$y=',num2str(X(1),3),'{\left( \frac{x}{',num2str(X(1),3),'}\right)^{e^{-',num2str(X(3),3),'t}}}',  ';', 'R^2=',num2str(Rsq,4),'$');
            h = legend(string);set(h,'Interpreter','Latex');
          q
          string

            title('Gompertz Fit, MCF-7')

            ylim([ 0 100])
            %hold off
            hold on
        end
        hold off
               
        xlim ([0 546]);
        xticks([0: 100: 546]);
         
        for qq= 1:1:4
             
            
            %N(t)= K+(N-K)exp(-r*t) 
            y= @(x, time) x(1)*((x(2)/x(1)).^(exp(-x(3).*time))); %x(1)= K x(2)= N0 x(3)= r
            
            %OLS= @(r) sum((y(r,time)-dataSeriesRED1(qq,:)).^2);
            OLS= @(x) sum((y(x,time)-dataSeriesRED1(qq,:)).^2);
            
            %opts= optimset('MaxFunEvals', 1e25, 'MaxIter', 1e4);
            opts= optimset();
            
            
              %calculate a GR ::            
            t= time(1,1:13)
            f= log(dataSeriesRED1(qq,1:13))
            pRed= polyfit(t,f,1)                                      
            GR_Red= pRed(1)
%             v= polyval(pRed, t)
%              plot(t,f,'o',t,v, '-')
        
            %%read-in parameters
                %X= fminsearch(OLS, [mean(dataSeriesRED1(q, end-10:end)) dataSeriesRED1(1) GR_Red], opts)
            %%initial guess
                %X= fminsearch(OLS, [10 3 0.005], opts)
            %%r and N0 only
                %X= fminsearch(OLS, [10 dataSeriesRED1(1) GR_Red], opts)
            %%K & N0
                X= fminsearch(OLS, [mean(dataSeriesRED1(q,end-10:end)) dataSeriesRED1(1) 0.005], opts)
            
            
            GompFit= y(X,time);  
            leastSquares= sum((y(X,time)-dataSeriesRED1(qq,:)).^2);
            
            
               R=corrcoef(dataSeriesRED1(qq,:), y(X,time))
               Rsq= (R(1,2).^2)
            
            
            subplot(1,2,2);%exp_num
            hh= plot(time, dataSeriesRED1(qq,:),'.','color', c(qq))%,'.', 'color', c(q)); %'LineWidth', 1);
            hold on
            
            plot(time,GompFit,'color',c(qq)); 
            set(get(get(hh,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            
            
            string{qq}=strcat('$y=',num2str(X(1),3),'{\left( \frac{x}{',num2str(X(1),3),'}\right)^{e^{-',num2str(X(3),3),'t}}}',  ';', 'R^2=',num2str(Rsq,4),'$');
            h = legend(string);set(h,'Interpreter','Latex');
          qq
          string

            title('Gompertz Fit, MDA-MB-231')

            ylim([ -5 16])
            %hold off
            hold on    
        end
        hold off
        
                      
        %calculating means etc--> this is K
        mean_vectorGREEN1(exp_num,i) = getMeanLastTen(dataSeriesGREEN1);    
        mean_vectorRED1(exp_num,i) = getMeanLastTen(dataSeriesRED1);   

        
    end
end

mean_vectorRED1
mean_vectorGREEN1

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

function meanFirstThree= meanFirstThree(data)
    meanFirstThree= mean(data(:,1:3));
end

function value = getMeanLastTen(data)
    value = mean(mean(data(:,(end-10):end)));
end