function GC_Lin_LS
    GC_Lin_LS_Para
    
    while norm(x - x_star) > tol
        i = i + 1;
        x = x_list(:,i);
        W = W_list(:,:,MF_count);
        b = b_list(:,:,MF_count);
        delta_G = W*x + b;
        model_grad = d_hat_f(x);
        State = new_state;
        if State == 0
           %% Model-Based Method
            [eta, LS_MB_fail, count] = bt_LS_MB(f, x, model_grad + delta_G, rho, Beta, State);
            if LS_MB_fail == 0 
                new_state = 0;
                x_list(:,i+1) = x - eta * (model_grad + delta_G);
                disp('Model Base')
%                 keyboard
            else
                new_state = 1;
                x_list(:,i+1) = x;
            end
            last_count = total_count;
            total_count = total_count + count;         
            MB_Start_Index(MB_count) = last_count ;
            MB_End_Index(MB_count) = total_count;
            MB_count = MB_count + 1;
            x_list_count(:, last_count + 1: total_count - 1 ) = repmat(x,1,total_count - last_count - 1);
%             keyboard
            x_list_count(:, total_count) = x_list(:,i+1);
        else
           %% Model-Free Method
           disp('Model Free')
           true_grad = d_f(x);
%            keyboard
           last_count = total_count;
           total_count = total_count + 2 * n + 1;
           [eta,count] = bt_line_search_MF(f, x, true_grad, rho, Beta, d_f, x_star);
           total_count = total_count + count;
           MF_Start_Index(MF_count) = last_count ;
           MF_End_Index(MF_count) = total_count;

           x_list(:,i+1) = x - eta * true_grad;

           x_list_count(:, last_count + 1: total_count - 1 ) = repmat(x,1,total_count - last_count - 1);
           x_list_count(:, total_count) = x_list(:,i+1);

           [W,b] = Linear_Comp_Estimator(x,r0,d_hat_f,d_f);
           W_list(:,:,MF_count + 1) = W;
           b_list(:,:,MF_count + 1) = b;

           MF_count = MF_count + 1;
           new_state = 0;
        end
        
        if mod(i,5) == 0
            disp(['Iteration: ',num2str(i)])
        end    
    end
%     keyboard
    x_list_count = x_list_count(:,1:total_count);
    x_diff_count = x_list_count - x_star;
    for i = 1: MF_count - 1
        MF_steps = MF_Start_Index(i)  : MF_End_Index(i);
        MF_regime = x_diff_count(:, MF_steps);
        Model_free_progeress = norms(MF_regime,[], 1);
        semilogy( MF_steps, Model_free_progeress , 'b','LineWidth',1.5)
        hold on
    end

    for i = 1: MB_count - 1
        MB_steps = MB_Start_Index(i): MB_End_Index(i);
        MB_regime = x_diff_count(:, MB_steps);
        Model_based_progeress = norms(MB_regime,[], 1);
        semilogy( MB_steps, Model_based_progeress , 'r','LineWidth',1.5)
        hold on
    end
    xlabel('Number of Function Evaluations','FontSize',16)
    ylabel('$\|x - x^{*}\|$','Interpreter','latex','FontSize',16)
    
end