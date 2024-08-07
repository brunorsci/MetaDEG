#' @title metaDEG
#' @description
#' A função realiza uma meta-análise nos resultados de expressão diferencial combinados de vários estudos.
#' @param DE_results_list Uma lista combinada dos resultados individuais de DEGs gerados pela função diff_expression.
#' @param sm Medida de efeito (padrão: "SMD" - Diferença de Média Padronizada).
#'        Outras opções: "MD" (Diferença de Médias), "OR" (Odds Ratio), "RR" (Risk Ratio).
#' @param method.tau Método para estimar tau² (padrão: "REML").
#'        Outras opções: "DL" (DerSimonian-Laird), "SJ" (Sidik-Jonkman), "EB" (Empirical Bayes).
#' @param hakn Ajuste de variância de Hartung-Knapp (padrão: TRUE).
#' @param method.ci Método para calcular intervalos de confiança (padrão: "z").
#'        Outras opções: "t" (intervalos de confiança baseados em t).
#' @param prediction Calcular intervalos de predição (padrão: TRUE).
#' @param fixed Incluir resultados de efeito fixo (padrão: FALSE).
#' @param random Incluir resultados de efeito aleatório (padrão: TRUE).
#' @return Um data frame com os resultados da meta-análise de diferentes experimentos de DEGs.
#' @author Bruno Rodrigo Assuncao(2024).
#' @export
#' @importFrom meta metagen
#' @importFrom stats p.adjust

metaDEG <- function(DE_results_list, sm = "SMD", method.tau = "REML", hakn = TRUE,
                    method.ci = "z", prediction = TRUE, fixed = FALSE, random = TRUE) {
### EXCEPTIONS ----####
  # Verificar e carregar os pacotes necessários
  required_packages <- c("meta", "stats")
  for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE)) {
      stop(paste("O pacote", pkg, "não está instalado. Por favor, instale-o antes de executar esta função."))
    }
  }
  # Checar valores dos parâmetros
  valid_sm <- c("SMD", "MD", "OR", "RR")
  if (!sm %in% valid_sm) {
    stop(paste("Valor inválido para 'sm'. Valores válidos são:", paste(valid_sm, collapse = ", ")))
  }
  
  valid_method_tau <- c("REML", "DL", "SJ", "EB")
  if (!method.tau %in% valid_method_tau) {
    stop(paste("Valor inválido para 'method.tau'. Valores válidos são:", paste(valid_method_tau, collapse = ", ")))
  }
  
  valid_method_ci <- c("z", "t")
  if (!method.ci %in% valid_method_ci) {
    stop(paste("Valor inválido para 'method.ci'. Valores válidos são:", paste(valid_method_ci, collapse = ", ")))
  }
  
# Checar se hakn, prediction, fixed, e random são booleanos
  if (!is.logical(hakn) || length(hakn) != 1) {
    stop("'hakn' deve ser um valor booleano (TRUE ou FALSE).")
  }
  if (!is.logical(prediction) || length(prediction) != 1) {
    stop("'prediction' deve ser um valor booleano (TRUE ou FALSE).")
  }
  if (!is.logical(fixed) || length(fixed) != 1) {
    stop("'fixed' deve ser um valor booleano (TRUE ou FALSE).")
  }
  if (!is.logical(random) || length(random) != 1) {
    stop("'random' deve ser um valor booleano (TRUE ou FALSE).")
  }

#### RUNNING FUNCTION ----####
  
  genes <- rownames(DE_results_list[[1]])
  meta_results <- data.frame(Gene = genes, EffectSize = NA, StdErr = NA, p.value = NA, stringsAsFactors = FALSE)
  
  for (gene in genes) {
    effect_sizes <- sapply(DE_results_list, function(res) res[gene, "logFC"]) # Extrai os valores de logFC(Log Fold Change-logFC)- tamanho do efeito.
    variances <- sapply(DE_results_list, function(res) res[gene, "t"]^2) # calcular a variância (valores de t ^2) para cada gene.
    
    # Verificar se todos os estudos possuem os coeficientes de interesse
    valid_studies <- !is.na(effect_sizes)
    
    if (sum(valid_studies) > 1) {
      # Realizar meta-análise
      meta_res <- metagen(TE = effect_sizes[valid_studies],
                          seTE = sqrt(variances[valid_studies]),
                          sm = sm,
                          method.tau = method.tau,
                          hakn = hakn,
                          method.ci = method.ci,
                          prediction = prediction,
                          fixed = fixed,
                          random = random)
      
      if (fixed && !is.null(meta_res$TE.fixed)) {
        meta_results[meta_results$Gene == gene, "EffectSize.fixed"] <- meta_res$TE.fixed
        meta_results[meta_results$Gene == gene, "StdErr.fixed"] <- meta_res$seTE.fixed
        meta_results[meta_results$Gene == gene, "p.value.fixed"] <- meta_res$pval.fixed
      }
      
      if (random && !is.null(meta_res$TE.random)) {
        meta_results[meta_results$Gene == gene, "EffectSize.random"] <- meta_res$TE.random
        meta_results[meta_results$Gene == gene, "StdErr.random"] <- meta_res$seTE.random
        meta_results[meta_results$Gene == gene, "p.value.random"] <- meta_res$pval.random
      }
    }
  }
  
  # Ajustar valores p
  if (fixed) {
    meta_results$p.adj.fixed <- p.adjust(meta_results$p.value.fixed, method = "fdr")
  }
  if (random) {
    meta_results$p.adj.random <- p.adjust(meta_results$p.value.random, method = "fdr")
  }
  
  return(meta_results)
}
