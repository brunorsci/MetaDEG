#'@title diff_expression
#'@description
#'Esta função realiza a análise de expressão diferencial usando contagens de RNA-Seq e metadados. Ela calcula as diferenças de expressão entre duas condições e realiza a análise diferencial usando o modelo linear.
#'@param counts Um data frame com os dados de contagem de expressão gênica. As colunas devem ser amostras e as linhas devem ser genes.
#'@param metadata Um data frame com informações sobre as amostras, incluindo as condições experimentais. Deve conter pelo menos duas colunas: uma para as amostras e outra para as condições.
#'@param sample_col Nome da coluna no data frame `metadata` que contém os identificadores das amostras. O padrão é `"sample"`.
#'@param condition_col Nome da coluna no data frame `metadata` que contém as condições experimentais (ex: caso vs controle). O padrão é `"condition"`.
#'@param condition1 Valor da condição de controle para a comparação. Por exemplo, `"controle"`.
#'@param condition2 Valor da condição experimental para a comparação. Por exemplo, `"caso"`.
#'@return Uma lista com dois data frames:
#'  - `DEG_Results`: Contém todos os resultados da expressão diferencial, incluindo log fold change, erro padrão, valor p ajustado, e comparação.
#'  - `Comparison_Info`: Contém informações sobre a comparação realizada.
#'@export
#'@importFrom edgeR DGEList calcNormFactors
#'@importFrom limma voom lmFit makeContrasts contrasts.fit eBayes
#'@importFrom stats p.adjust
diff_expression <- function(counts, metadata, sample_col = "sample", condition_col = "condition", condition1, condition2) {
  
  # Verificar e carregar os pacotes necessários
  required_packages <- c("edgeR", "limma", "stats")
  for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE)) {
      stop(paste("O pacote", pkg, "não está instalado. Por favor, instale-o antes de executar esta função."))
    }
  }
  
  # Verificar se os parâmetros de condição são válidos
  if (!condition1 %in% unique(metadata[[condition_col]])) {
    stop(paste("condition1:", condition1, "não encontrado em", condition_col, "no metadata."))
  }
  if (!condition2 %in% unique(metadata[[condition_col]])) {
    stop(paste("condition2:", condition2, "não encontrado em", condition_col, "no metadata."))
  }
  
  # Merge counts and metadata
  combined_data <- merge(metadata, counts, by.x = sample_col, by.y = sample_col)
  rownames(combined_data) <- combined_data[, sample_col]
  
  # Filtrar as condições de interesse
  combined_data <- combined_data[combined_data[, condition_col] %in% c(condition1, condition2), ]
  
  # Preparar os dados de expressão e a matriz de design
  expr_data <- combined_data[, !colnames(combined_data) %in% c(sample_col, condition_col)]
  
  # Remover quaisquer linhas ou colunas com valores NA
  expr_data <- expr_data[complete.cases(expr_data), ]
  
  # Certifique-se de que condition é um fator com níveis definidos
  combined_data[, condition_col] <- factor(combined_data[, condition_col], levels = c(condition1, condition2))
  
  # Criar a matriz de design
  design <- model.matrix(~0 + combined_data[, condition_col])
  colnames(design) <- levels(combined_data[, condition_col])
  
  # Criar o objeto DGEList
  dge <- DGEList(counts = t(expr_data))
  
  # Calcular os fatores de normalização usando TMM
  dge <- calcNormFactors(dge)
  
  # Transformação voom, aplicando os fatores de normalização calculados
  v <- voom(dge, design)
  
  # Ajuste do modelo linear
  fit <- lmFit(v, design)
  contrast <- paste(condition2, "-", condition1, sep="")
  contrast.matrix <- makeContrasts(contrasts = contrast, levels = design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  # Obter todos os resultados
  results <- as.data.frame(fit2$coefficients) # Coeficientes (logFC)
  results$StdErr <- fit2$sigma # Erro padrão dos coeficientes
  results$p.value <- fit2$p.value # Valor p
  results$logFC <- fit2$coefficients # Log fold change
  results$p.value.adj <- p.adjust(results$p.value, method = "fdr") # Valor p ajustado
  
  # Adicionar coluna de comparação aos resultados
  results$Comparison <- contrast
  
  # Exportar a comparação para um objeto separado
  comparison_info <- data.frame(
    Gene = rownames(results),
    Comparison = contrast,
    Condition1 = condition1,
    Condition2 = condition2
  )
  
  return(list(
    DEG_Results = results,
    Comparison_Info = comparison_info
  ))
}
