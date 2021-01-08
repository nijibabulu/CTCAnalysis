
lm_eqn <- function(df){
  m <- lm(y ~ x, df);
  f <- summary(m)$fstatistic
  pv <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(pv) <- NULL
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(p)~"<"~pv~~italic(r)^2~"="~r2,
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        pv = format(pv, digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}


plot_roc <- function(data, signal, classification,  response, color, title=signal, response_name=response,
                     regression = F, dot_panel = F, plot_optimal=F,
                     color_cutoffs=F, show_grid=F, fill_curve=F) {
  if(regression && !dot_panel) {
    stop("Regression is not possible without the dot panel display")
  }
  pred <- ROCR::prediction(dplyr::pull(data, signal), dplyr::pull(data,classification))
  perf <- ROCR::performance(pred, "tpr", "fpr")
  auc <- ROCR::performance(pred, "auc")@y.values[[1]]

  df <- tibble::tibble(`False Positive Rate`=perf@x.values[[1]],
                       `True Positive Rate`=perf@y.values[[1]],
                       Cutoff=perf@alpha.values[[1]])

  youden <- df[which.max(df$`True Positive Rate`-df$`False Positive Rate`),]

  line_aes <- if(color_cutoffs) ggplot2::aes(color=Cutoff) else ggplot2::aes()

  p <- ggplot2::ggplot(df , ggplot2::aes(x=`False Positive Rate`, ymin=0, ymax=`True Positive Rate`, y=`True Positive Rate`)) +
    ggplot2::geom_line(line_aes) +
    ggplot2::geom_abline(intercept = 0, slope = 1, lty=20) +
    ggplot2::annotate("text", x = .75, y=0.25, label=stringr::str_glue("AUC = {prettyNum(auc, digits=2)}")) +
    ggplot2::scale_color_distiller(palette = "Reds") +
    ggplot2::theme_bw() +
    ggplot2::coord_equal() +
    ggplot2::labs(title = title)
    ggplot2::theme()

  if(fill_curve)
    p <- p + ggplot2::geom_ribbon(fill="lightgrey")

  if(plot_optimal)
    p <- p + ggplot2::geom_point(data= youden) +
      ggrepel::geom_text_repel(data= youden, ggplot2::aes(label=stringr::str_glue("Optimal cutoff = {prettyNum(Cutoff, digits = 2)}")))

  responsemax <- max(dplyr::pull(data, response))
  eq <- lm_eqn(data %>% dplyr::select(y=signal, x=response))
  p1 <- ggplot2::ggplot(data, ggplot2::aes(x= !!rlang::ensym(response), y = !!rlang::ensym(signal), color=!!rlang::ensym(color))) +
    ggplot2::geom_point() +
    ggplot2::geom_hline(yintercept = youden$Cutoff, lty=20) +
    ggplot2::labs(x=response_name, y=stringr::str_glue("{title} Score")) +
    #labs(y="RB Loss Signature Score", x="PFS Time") +
    ggplot2::guides(color=ggplot2::guide_legend(title="Response Status")) +
    ggplot2::theme_bw()



  if(regression) {
    p1 <- p1 +  ggplot2::geom_smooth(mapping = ggplot2::aes(group = 1), color="black", method=lm) +
      ggplot2::annotate("text", x=responsemax/2, y=1.5, label=eq, parse=T)
  }

  if(dot_panel) {
    p <- patchwork::wrap_plots(p, p1, nrow=1, guides="collect")
  }
  if(!show_grid)
    p <- p & ggplot2::theme(panel.grid = ggplot2::element_blank())

  p
}
