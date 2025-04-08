

evaluateModels <- function(data, iteration = 1, prop = 0.64, weighted = TRUE) {
  set.seed(iteration)
  data_split <- initial_split(data, strata = response, prop = prop)
  data_train <- training(data_split)
  data_test <- testing(data_split)
  
  train.y <-  data_train$response
  fraction_0 <- rep(1 - sum(train.y == 0) / length(train.y), sum(train.y == 0))
  fraction_1 <- rep(1 - sum(train.y == 1) / length(train.y), sum(train.y == 1))
  weights <- numeric(length(train.y))
  if (weighted == TRUE) {
    weights[train.y == 0] <- fraction_0
    weights[train.y == 1] <- fraction_1
  } else {
    weights <- rep(1, length(train.y))
  }
  data_train <-  data_train %>%
    mutate(
      case_wts = weights, 
      case_wts = importance_weights(case_wts) 
    )
  
  # Set up workflowset
  rec <- recipe(response ~ ., data =  data_train) %>%
    update_role(SUBJID, new_role = "ID") %>%
    update_role(Treatment, new_role = "Treatment") %>%
    step_impute_median(all_predictors()) %>%
    step_center(all_predictors()) %>% 
    step_scale(all_predictors())
  
  
  # Model specifications
  rf_spec <- rand_forest(trees = 1e3) %>% # untuned
    set_engine("ranger", importance = "impurity", seed = 2024) %>%
    set_mode("classification")
  
  
  lasso_reg_spec <-
    logistic_reg(mode = "classification", penalty = tune(), mixture = 1) %>%
    set_engine("glmnet")
  
  
  wfset <- 
    workflow_set(
      preproc = list(preprocess = rec), 
      models = list(lasso_reg_spec = lasso_reg_spec,
                    rf_spec = rf_spec
      ),
      case_weights = case_wts 
    )
  
  
  
  grid_ctrl <- control_grid(
    verbose = FALSE,
    save_pred = TRUE,
    save_workflow = TRUE,
    parallel_over = "everything"
  )
  
  data_folds <- vfold_cv(data_train, strata = response, repeats = 10, v = 5)
  
  suppressMessages(
    tune_results <-
      wfset %>%
      workflow_map(
        "tune_grid",
        seed = 2024,
        resamples = data_folds,
        grid = 100,
        control = grid_ctrl,
        verbose = FALSE,
        metrics = metric_set(roc_auc)
      )
  )
  
  
  models <- tune_results$wflow_id
  j <- 2
  coefs_tmp <- data.frame()
  preds_tmp <- data.frame()
  for(j in 1:length(models)) {
    
    # Fit the model
    tune_tmp <- tune_results %>%
      extract_workflow_set_result(models[j]) %>% 
      select_best(metric = "roc_auc")
    
    model_tmp <- tune_results %>% 
      extract_workflow(models[j]) %>%
      finalize_workflow(tune_tmp) %>%
      fit(data = data_train)
    
    # Predict 
    pred_tmp_mod <- cbind(data_test[,c(1:3)], predict(model_tmp, data_test, type = "prob")[,2]) %>%
      mutate(model = models[j],
             iteration = iteration)# Get prediction
    
    
    
    # Obtain variable importance scores
    if(models[j] %in% c("preprocess_lasso_reg_spec", "preprocess_elasticnet_reg_spec")) { # Get feature importance
      feature_tmp <- model_tmp %>%
        pull_workflow_fit() %>%
        tidy() %>%
        mutate(model = models[j],
               iteration = iteration)
      
    } else if (models[j] %in% c("preprocess_rf_spec", "preprocess_xgb_spec")) {
      feature_tmp <- model_tmp  %>% # For random forest and xgboost
        pull_workflow_fit() %>%
        vip(num_features = dim(data_train)[2]-4)
      feature_tmp <- feature_tmp$data %>%
        mutate(model = models[j],
               iteration = iteration) %>%
        dplyr::rename(term = Variable, 
                      estimate = Importance) %>%
        mutate(penalty = "NA")
    }
    
    
    preds_tmp <- rbind(preds_tmp, pred_tmp_mod)
    coefs_tmp <- rbind(coefs_tmp, feature_tmp)
    
  }
  
  
  
  
  ret <- list()
  ret$pred_mod <- preds_tmp
  ret$feature <- coefs_tmp
  return(ret)
  
  
}

safelyw <- safely(wilcox_test)


plotModels <- function(result) {
  
  # Gather predictions
  mod_predictions <- list()
  for (i in 1:length(result)) {
    mod_predictions[[i]] <- result[[i]]$pred_mod
  }
  mod_predictions <- mod_predictions %>%
    bind_rows() %>%
    group_by(SUBJID, response, model) %>%
    dplyr::summarise(pred = median(.pred_1)) %>%
    mutate(model = str_remove(model, pattern = "preprocess_")) %>%
    mutate(model = str_remove(model, pattern = "_spec")) %>%
    ungroup()
  

  # Gather features
  mod_features <- list()
  for (i in 1:length(result)) {
    mod_features[[i]] <- result[[i]]$feature
  }
  mod_features <- mod_features %>%
    bind_rows() %>%
    filter(term != "(Intercept)") %>%
    group_by(term, model) %>%
    dplyr::summarise(imp = median(estimate)) %>%
    arrange(desc(abs(imp))) %>%
    mutate(model = str_remove(model, pattern = "preprocess_")) %>%
    mutate(model = str_remove(model, pattern = "_spec")) %>%
    ungroup()
  
  
  # Multivariable models
  multivariable_plots <- mod_predictions %>%
    split(.$model) %>%
    purrr::map(., function(x) {
      
      model_i <- unique(x$model)

      p.value <- x %>%
        data.frame() %>%
        safelyw(pred ~ response)
      p.value <- p.value$result %>%
        mutate(p = paste0("p.value = ", p))
      
      max.pred <- max(x$pred) + 0.1
      
      ### Predictions
      preds.plot <- x %>% 
        ggplot(aes(x = response, y = pred)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.15, height = 0, shape = 21, size = 1, aes(fill = response)) +
        scale_fill_manual(values = c("0" = "#FFD729", "1" = "darkgrey")) +
        theme_bw() +
        theme(axis.title.x = element_text(size = 5),
              axis.title.y = element_text(size = 5),
              axis.text.x = element_text(size = 5),
              axis.text.y = element_text(size = 5),
              plot.title = element_text(size = 5),
              plot.subtitle = element_text(size = 5),
              legend.text = element_text(size = 5),
              legend.title = element_text(size = 5),
              strip.text.y = element_text(size = 5, angle = 0),
              strip.background = element_blank(),
              legend.key.width = unit(0.2, "cm"),
              legend.key.height = unit(0.2, "cm"),
              panel.grid = element_blank(),
              plot.caption = element_text(size = 4.5),
              axis.ticks = element_line(colour = "black", size = 0.2),
              panel.border = element_rect(fill=NA, colour = "black", size=0.2),
              legend.position = "none") +
        geom_text(data = p.value, aes(x = 1.5, y = max.pred, label = p), size = 5 / (14/5)) +
        labs(y = "prediction", x = "") 
      
      ### ROC curve
      roc_obj <- roc(x$response, x$pred, ci = TRUE, direction="<") # “<”: if the predictor values for the control group are lower or equal than the values of the case group
      r.stats <- data.frame(roc_obj$ci) %>%
        t() %>%
        data.frame() %>%
        dplyr::rename(ci.lower = X1,
                      auc = X2,
                      ci.upper = X3) %>%
        mutate(label = paste0("AUC = ", round(auc, 2), " (", round(ci.lower, 2), " - ", round(ci.upper, 2), ")"))
      
      roc_obj <- smooth(roc_obj, method = "density", bw = 0.001)
      r.data <- data.frame(x = 1- roc_obj$specificities,
                           y = roc_obj$sensitivities)
      r.ci <- ci.se(roc_obj, specificities = seq(0, 1, l = 25))
      r.ci <- data.frame(x = as.numeric(rownames(r.ci)),
                         lower = r.ci[, 1],
                         upper = r.ci[, 3])
      
      ROC <- r.data %>%
        ggplot(aes(x = x, y = y)) +
        geom_line(size = 0.5) +
        geom_ribbon(data = r.ci, aes(x = 1 - x, ymax = upper, ymin = lower), inherit.aes = F, alpha = 0.1) +
        geom_abline (intercept = 0, slope = 1, linetype = "dashed", size = 0.5) +
        theme_bw() +
        geom_text(data = r.stats, aes(x = 0.5, y = 1.1, label = label), size = 5 / (14/5)) +
        scale_x_continuous("False Positive Rate", limits = c(0,1.1), breaks = seq(0, 1, by = 0.2)) +
        scale_y_continuous("True Positive Rate", limits = c(0,1.1), breaks = seq(0, 1, by = 0.2)) +
        theme(axis.title.x = element_text(size = 5),
              axis.title.y = element_text(size = 5),
              axis.text.x = element_text(size = 5),
              axis.text.y = element_text(size = 5),
              plot.title = element_text(size = 5),
              legend.text = element_text(size = 5),
              legend.title = element_text(size = 5),
              strip.text.y = element_text(size = 5, angle = 0),
              strip.background = element_blank(),
              legend.key.width = unit(0.2, "cm"),
              legend.key.height = unit(0.2, "cm"),
              legend.position = "right",
              panel.grid = element_blank(),
              axis.ticks = element_line(colour = "black", size = 0.2),
              panel.border = element_rect(fill=NA, colour = "black", size=0.2),
              plot.subtitle = element_text(size = 5))
      
      ### Variable importance
      featureplot <- mod_features %>%
        filter(model == model_i) %>%
        mutate(imp = abs(imp)) %>%
        arrange(abs(imp)) %>%
        mutate(term = factor(term, levels = unique(.$term))) %>%
        ggplot(aes(y = imp, x = term)) +
        geom_col() +
        coord_flip() +
        theme_bw() +
        labs(x = "", y = "importance", title = "feature importance") +
        theme(axis.title.x = element_text(size = 5),
              axis.title.y = element_text(size = 5),
              axis.text.x = element_text(size = 5),
              axis.text.y = element_text(size = 5),
              plot.title = element_text(size = 5),
              legend.text = element_text(size = 5),
              legend.title = element_text(size = 5),
              strip.text.y = element_text(size = 5, angle = 0),
              strip.background = element_blank(),
              legend.key.width = unit(0.2, "cm"),
              legend.key.height = unit(0.2, "cm"),
              legend.position = "right",
              panel.grid = element_blank(),
              plot.subtitle = element_text(size = 5),
              axis.ticks = element_line(colour = "black", size = 0.2),
              panel.border = element_rect(fill=NA, colour = "black", size=0.2),
              plot.caption = element_text(size = 5)) 
      result <- list(preds.plot = preds.plot,
                     ROC = ROC,
                     featureplot = featureplot)
      
      
      result <- ggarrange(plotlist = result, nrow = 1, widths = c(1,1.5,1.5))
      result <- result %>%
        annotate_figure(., top = text_grob(model_i, face = "bold", size = 6.5))
      return(result)
      
    })
  
  
  
  
  result <- list(multivariable_plots = multivariable_plots)
  return(result)
}



