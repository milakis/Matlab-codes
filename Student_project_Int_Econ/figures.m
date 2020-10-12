function figures(data,data_dev,in,P,pigs)

    % these lines are a small hack so that figures aren't all on top of
    % each other when they are on screen
    top = get(0,'screensize');
    top = top(4);   
    menu_height = 80;
    next_top = top - menu_height;    

    % select variables to be plotted. if want to add/change any series,
    % have to do it in bigdata (names of vars in data struct), bigmodel
    % (names of vars in structural struct), labels and filename.
    bigdata = {repmat(P.BGoY,P.years,1).*(1+data_dev.bg_used), data_dev.n, data_dev.y, data_dev.rho};
    bigmodel = {repmat(P.BGoY,P.years,1).*(1+in.bg_nom), in.n, in.y_nom, in.rho};

    labels = {'govt debt','employment','output','rho'};     
    legends = [1,1,1,1];
    filename = {'bg','n','y','rho'};
    shock_string = P.shock_string;
    
    counter = max(size(bigdata));

    if strcmp(pigs,'pigs')
        pig_figs = 0;
        shock_string = [shock_string '_PIGS'];
        xx = 2;
        yy = 3;
        fig_set = [4,7,8,11,12];
    else
        pig_figs = 0;
        xx = 3;
        yy = 4;
        fig_set = 1:P.countries;        
    end
    
    if P.demean == 1
       shock_string = [shock_string '_dm'];
       N = length(bigdata);
       for i = 1:N
           bigmodel{i} = bigmodel{i} - repmat(mean(bigmodel{i}),P.years,1) + repmat(mean(bigdata{i}),P.years,1);
       end       
    end        
    
    ylims = zeros(counter,2);
    for jj = 1:counter
         ylims(jj,1) = max(max(max(bigdata{jj},bigmodel{jj}))) + 0.01;
         ylims(jj,2) = min(min(min(bigdata{jj},bigmodel{jj}))) - 0.01;
    end

    for jj = 1:counter            
        curr_data = bigdata{jj};
        curr_model = bigmodel{jj};
        figure('name',labels{jj})
        loc = get(gcf,'Position');
        loc(2) = next_top - loc(4);
        next_top = next_top - menu_height/3;
        set(gcf,'Position',loc);
        counter = 0;
        for ii = fig_set
            counter = counter + 1;
            hold on
            subplot(yy,xx,counter);        
            plot(P.start_year:1:P.end_year,curr_data(:,ii),'-',P.start_year:1:P.end_year,curr_model(:,ii),'--','LineWidth',1.5);
            xlim([P.start_year P.end_year])
            if P.same_scale == 1
                ylim([ylims(jj,2) ylims(jj,1)])                    
            else
                ylim('auto')
                if jj==3
                     ylim([0.89 1.11])
                end
                if jj==4
                     ylim([0.75 1.25])
                end                    
            end    
            set(gca,'XTick',P.start_year:1:P.end_year)          
            if P.end_year > 2012
                set(gca,'XTickLabel',{'2001','','','2004','','','','2008','','','','2012','','','','2016','','','','2020'})                
            else
                set(gca,'XTickLabel',{'2001','','','2004','','','','2008','','','','2012'})                
            end
            % Adjust figure positions and sizes for PIGS
            if pig_figs == 1
                pos = get(gca, 'Position');
                if ii == fig_set(1)
                    pos(1) = 0.065;
                    pos(2) = 0.58;
                end
                if ii == fig_set(2)
                    %pos(1) = 0.06;
                    pos(2) = 0.58;
                end                
                if ii == fig_set(3)
                    pos(1) = 0.065;
                    pos(2) = 0.08;
                end   
                if ii == fig_set(4)
                    %pos(1) = 0.06;
                    pos(2) = 0.08;
                end                   
                pos(3) = 0.395;
                pos(4) = 0.38;                
                %pos(3) = 0.9;
                set(gca, 'Position', pos)
                xlabel('year')
            end            
            title(P.names{ii});                    
            hold off
        end  
        hold on
        if pig_figs==0
            xlabel('year')
        end
        if legends(jj) == 1
            h = legend('data','benchmark model'); 
        elseif legends(jj) == 2
            h = legend('data','benchmark model');
        else
            h = legend('benchmark model');
        end
        
        if P.fiscal_counterfactual == 1 || P.mp_counterfactual == 1 || P.rho_counterfactual == 1
            if legends(jj) == 1
                h = legend('data','counterfactual model'); 
            elseif legends(jj) == 2
                h = legend('data','counterfactual model');
            else
                h = legend('counterfactual model');
            end
        end
        legend boxoff
        set(h,'Location','Best')
        if pig_figs == 0
            set(h, 'Position',[0.7,0.1,0.15,0.15]);
        end
        set(gcf, 'PaperPosition', [0 0 8 8*5.3/6])
        set(gcf, 'PaperSize', [8 8*5.3/6])        
        % Font size
        set(findall(gcf,'type','axes'),'fontsize',13)
        set(findall(gcf,'type','text'),'fontSize',14) 
        
    end
    
end

