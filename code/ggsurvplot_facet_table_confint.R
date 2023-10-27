
#source("https://raw.githubusercontent.com/BingxinS/survminer-fix/master/ggsurvplot_facet_fix_12072018.R")

ggsurvplot_facet_fix <- function (fit, data, facet.by, color = NULL, palette = NULL, 
          legend.labs = NULL, pval = FALSE, pval.method = FALSE, pval.coord = NULL, 
          pval.method.coord = NULL, nrow = NULL, ncol = NULL, scales = "fixed", 
          short.panel.labs = TRUE, panel.labs = NULL, panel.labs.background = list(color = NULL, 
                                                                                    fill = NULL), panel.labs.font = list(face = NULL, color = NULL, 
                                                                                                                         size = NULL, angle = NULL), panel.labs.font.x = panel.labs.font, 
          panel.labs.font.y = panel.labs.font, ...) 
{
  if (length(facet.by) > 2) 
    stop("facet.by should be of length 1 or 2.")
  if (!is.null(panel.labs) & !is.list(panel.labs)) # changed from .is_list() to is.list()
    stop("Argument panel.labs should be a list. Read the documentation.")
  . <- NULL
  .labeller <- "label_value"
  if (short.panel.labs) 
    .labeller <- label_both
  .dots <- list(...)
  fit.ext <- .extract.survfit(fit, data)
  .formula <- fit.ext$formula
  surv.obj <- fit.ext$surv
  surv.vars <- fit.ext$variables
  all.variables <- c(surv.vars, facet.by) %>% unique()
  vars.notin.groupby <- setdiff(all.variables, facet.by)
  data <- fit.ext$data.all
  if (!is.null(panel.labs)) {
    for (.grouping.var in facet.by) {
      if (!is.null(panel.labs[[.grouping.var]])) {
        if (!is.factor(data[, .grouping.var])) 
          data[, .grouping.var] <- as.factor(data[, .grouping.var])
        levels(data[, .grouping.var]) <- panel.labs[[.grouping.var]]
      }
    }
  }
  if (!all(facet.by %in% surv.vars) | !is.null(panel.labs)) {
    fit <- .build_formula(surv.obj, all.variables) %>% surv_fit(., 
                                                                data = data)
  }
  if (length(vars.notin.groupby) == 1) {
    if (is.null(color)) 
      color <- vars.notin.groupby
    .survformula <- .build_formula(surv.obj, vars.notin.groupby)
  }
  else {
    new.strata <- .create_strata(data, vars.notin.groupby, 
                                 sep = "; ")
    strata.levels <- levels(new.strata)
    data <- data %>% dplyr::mutate(.strata. = new.strata)
    new.surv.vars <- paste(c(".strata.", facet.by), collapse = " + ")
    .new.formula <- .build_formula(surv.obj, new.surv.vars)
    fit <- surv_fit(formula = .new.formula, data = data)
    if (is.null(color)) 
      color <- ".strata."
    if (is.null(legend.labs)) 
      legend.labs <- gsub(";", ",", strata.levels)
    .survformula <- .build_formula(surv.obj, ".strata.")
  }
  ################# calling the ggsurvplot_fix function in this file  ####################
  #################    instead of ggsurvplot function in survminer    ####################
  ggsurv <- ggsurvplot_core_fix(fit, data = data, color = color, 
                                palette = palette, legend.labs = legend.labs, ...)
  p <- .facet(ggsurv$plot, facet.by, nrow = nrow, ncol = ncol, 
              scales = scales, short.panel.labs = short.panel.labs, 
              panel.labs.background = panel.labs.background, panel.labs.font = panel.labs.font, 
              panel.labs.font.x = panel.labs.font.x, panel.labs.font.y = panel.labs.font.y)
  if (pval) {
    grouped.d <- surv_group_by(data, grouping.vars = facet.by)
    sf <- surv_fit(.survformula, grouped.d$data, ...)
    grouped.d <- grouped.d %>% mutate(fit = sf)
    pvalue <- surv_pvalue(grouped.d$fit, grouped.d$data, 
                          pval.coord = pval.coord, pval.method.coord = pval.method.coord, 
                          ...) %>% dplyr::bind_rows() %>% tibble::as.tibble()
    pvals.df <- grouped.d %>% dplyr::select_(.dots = facet.by) %>% 
      dplyr::bind_cols(pvalue)
    pval.x <- pval.y <- pval.txt <- method.x <- method.y <- method <- NULL
    p <- p + geom_text(data = pvals.df, aes(x = pval.x, y = pval.y, 
                                            label = pval.txt), hjust = 0)
    if (pval.method) 
      p <- p + geom_text(data = pvals.df, aes(x = method.x, 
                                              y = method.y, label = method), hjust = 0)
  }
  p
}




ggsurvplot_facet_table_confint <- function (fit, data = NULL, 
                                            risktable = FALSE,
                                            risktable.title = NULL,
                                            risktable.ylab = NULL,
                                            conf.int = FALSE,
                                            facet.by = facet.by,
                                            text.size = 10, ...) 
{
  # get strata.by and facet.by from fit$call
  fit.str <- sub('.*~', '', as.character(fit$call)[2])
  # facet.by = 
  fit.factor <-  strsplit(fit.str, " ")
  fit.factor.length <- length(fit.factor[[1]])/2
  if (fit.factor.length > 3) 
    stop("facet.by has more than 2 factors. Please use 1 or 2 faceting factors.")
  if (fit.factor.length <= 1) 
    stop("facet.by missing. Please use 1 or 2 faceting factors.")
  strata.by <- fit.factor[[1]][2]
  # Generate risk table for each facet plot item
  if(fit.factor.length == 2) {
    facet.by <- fit.factor[[1]][4]
    
    #================get the plot  ================
    plt_fct <- ggsurvplot_facet_fix(fit, data, 
                                    facet.by = facet.by, short.panel.labs = TRUE)
    if (conf.int == TRUE){
      #================get the conf.int================
      ggsurv <- ggsurvplot(fit, data,
                           conf.int = conf.int, 
                           ggtheme = theme_bw())
      data.survplot <- .connect2origin_fix(ggsurv$data.survplot) 
      data.survplot <- .apply_surv_func(data.survplot)
      data.survplot.rect <- survdata2rect(data.survplot)
      # add conf.int to each facet plot
      facet.formula <- paste0("~", facet.by) %>% stats::as.formula()
      plt_fct <- plt_fct +
        #geom_ribbon(aes(ymin = lower, ymax = upper, fill = sex), alpha=0.3) +
        geom_rect(data = data.survplot.rect,
                  aes(xmin = time, xmax = timemax, ymin = lower, ymax = upper, fill = facet.by),
                  alpha = 0.3) +
        facet_wrap(facet.formula, labeller = label_both) + 
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              strip.background = element_rect(fill="white", size = 1),
              axis.text.x = element_text(color="black", size = text.size+2),
              axis.text.y = element_text(color="black", size = text.size),
              axis.title = element_text(color="black", size = text.size+4),
              strip.text = element_text(color="black", size = text.size))
    }
  }
  if(fit.factor.length == 3) {
    facet.by <- fit.factor[[1]][c(4,6)]
    
    #================get the plot  ================
    plt_fct <- ggsurvplot_facet_fix(fit, data, 
                                    facet.by = facet.by, short.panel.labs = TRUE)
    if (conf.int == FALSE){
      #================get the conf.int================
      ggsurv <- ggsurvplot(fit, data,
                           conf.int = conf.int, 
                           ggtheme = theme_bw())
      data.survplot <- .connect2origin_fix(ggsurv$data.survplot) 
      data.survplot <- .apply_surv_func(data.survplot)
      data.survplot.rect <- survdata2rect(data.survplot)
      # add conf.int to each facet plot
      facet.formula <- paste(facet.by, collapse = " ~ ") %>% stats::as.formula()
      plt_fct <- plt_fct +
        # geom_ribbon(aes(ymin = lower, ymax = upper, fill = sex), alpha=0.3) +
        geom_rect(data = data.survplot.rect,
                  aes(xmin = time, xmax = timemax, ymin = lower, ymax = upper, fill = facet.by),
                  alpha = 0.3) +
        facet_wrap(facet.formula, labeller = label_both) + 
        theme_bw() + 
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              strip.background = element_rect(fill="white", size = 1),
              axis.text.x = element_text(color="black", size = text.size+2),
              axis.text.y = element_text(color="black", size = text.size),
              axis.title = element_text(color="black", size = text.size+4),
              strip.text = element_text(color="black", size = text.size))
    }
  }
  if (risktable == FALSE){
    plt_fct
  } else {
    #================get the risk table================
    if(fit.factor.length == 2) {
      ggsurv <- ggsurvplot(fit, data,
                           risk.table = TRUE, 
                           ggtheme = theme_bw(), short.panel.labs = TRUE)
      tbl_fct <- 
        ggplot(ggsurv$table$data, ggplot2::aes_string("time", strata.by)) + 
        geom_text(aes(label = n.risk), size = max (2, (text.size-6)) ) +
        facet_wrap(facet.formula, labeller = label_both) + 
        theme_bw() + 
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              strip.background = element_rect(fill="white", size = 1),
              axis.text.x = element_text(color="black", size = text.size+2),
              axis.text.y = element_text(color="black", size = text.size),
              axis.title = element_text(color="black", size = text.size+4),
              strip.text = element_text(color="black", size = text.size))+
        ggtitle(risktable.title) +
        xlab("Time") + ylab(risktable.ylab) 
    }
    if(fit.factor.length == 3) {   
      ggsurv <- ggsurvplot(fit, data,
                           risk.table = TRUE, 
                           ggtheme = theme_bw(), short.panel.labs = TRUE)
      tbl_fct <- 
        ggplot(ggsurv$table$data, ggplot2::aes_string("time", strata.by)) + 
        geom_text(aes(label = n.risk), size = max (2, (text.size-6)) ) +
        facet_grid(facet.formula, labeller = label_both) + 
        theme_bw() + 
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              strip.background = element_rect(fill="white", size = 1),
              axis.text.x = element_text(color="black", size = text.size+2),
              axis.text.y = element_text(color="black", size = text.size),
              axis.title = element_text(color="black", size = text.size+4),
              strip.text = element_text(color="black", size = text.size))+
        ggtitle(risktable.title) +
        xlab("Time") + ylab(risktable.ylab) 
    }
    #========== combine plot and risk tables ===========
    g_plt_fct <- ggplotGrob(plt_fct)
    g_tbl_fct <- ggplotGrob(tbl_fct)
    min_ncol <- min(ncol(g_plt_fct), ncol(g_tbl_fct))
    g_plt_tbl <- gridExtra::gtable_rbind(g_plt_fct[, 1:min_ncol], g_tbl_fct[, 1:min_ncol], size="last")
    g_plt_tbl$widths <- grid::unit.pmax(g_plt_fct$widths, g_tbl_fct$widths)
    grid::grid.newpage() 
    grid::grid.draw(g_plt_tbl)
  }
}

survdata2rect <- function(data.survplot){
  data.survplot.rect = data.survplot[-nrow(data.survplot), ]
  data.survplot.rect$timemax = data.survplot$time[-1]
  strata.number <- length(levels(data.survplot$strata))
  ind <- which( (data.survplot.rect$timemax - data.survplot.rect$time)<0 )
  data.survplot.rect$timemax[ind] <- data.survplot.rect$time[ind]
  data.survplot.rect$timemax[c(1:strata.number)] <- data.survplot.rect$time[c(strata.number+1,ind+1)]
  data.survplot.rect
}