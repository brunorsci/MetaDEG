DEG Meta Analysis - RNASeq Data
================
Bruno Rodrigo Assunção
2024-07-23

WORKFLOW BIOINFORMÁTICO PARA META ANALISE DE EXPRESSÃO DIFERENCIAL EM
DADOS DE RNASEQ

**Preparação dos Dados:**

- Criamos dados fictícios de contagem e metadados para três estudos
  diferentes.

- Cada estudo tem amostras com dados de expressão gênica (`Gene1`,
  `Gene2`, `Gene3`), com amostras classificadas como “Case” ou
  “Control”.

### 1. Preparação dos Dados

``` r
# Carregar pacotes necessários
library(limma)
library(edgeR)
library(meta)
```

    ## Warning: package 'meta' was built under R version 4.3.3

    ## Carregando pacotes exigidos: metadat

    ## Warning: package 'metadat' was built under R version 4.3.3

    ## Loading 'meta' package (version 7.0-0).
    ## Type 'help(meta)' for a brief overview.
    ## Readers of 'Meta-Analysis with R (Use R!)' should install
    ## older version of 'meta' package: https://tinyurl.com/dt4y5drs

``` r
# Criar dados de contagem normalizados para três estudos
set.seed(123)
counts_study1 <- data.frame(
    SampleID = paste0("Sample", 1:8),
    Gene1 = rpois(8, lambda=10),
    Gene2 = rpois(8, lambda=15),
    Gene3 = rpois(8, lambda=20)
)

counts_study2 <- data.frame(
    SampleID = paste0("Sample", 9:16),
    Gene1 = rpois(8, lambda=11),
    Gene2 = rpois(8, lambda=16),
    Gene3 = rpois(8, lambda=21)
)

counts_study3 <- data.frame(
    SampleID = paste0("Sample", 17:24),
    Gene1 = rpois(8, lambda=9),
    Gene2 = rpois(8, lambda=14),
    Gene3 = rpois(8, lambda=19)
)
```

``` r
counts_study1
```

    ##   SampleID Gene1 Gene2 Gene3
    ## 1  Sample1     8     8    19
    ## 2  Sample2     9    19    16
    ## 3  Sample3    14    16    15
    ## 4  Sample4    10    16    18
    ## 5  Sample5    10    15    15
    ## 6  Sample6    15    12    17
    ## 7  Sample7    11    20    19
    ## 8  Sample8     5    18    14

``` r
# Criar metadados para cada estudo
metadata_study1 <- data.frame(
    SampleID = paste0("Sample", 1:8),
    Diagnosis = rep(c("Case", "Control"), each=4),
    Donor = rep(paste0("Donor", 1:4), times=2),
    Study = rep("Study1", 8)
)

metadata_study2 <- data.frame(
    SampleID = paste0("Sample", 9:16),
    Diagnosis = rep(c("Case", "Control"), each=4),
    Donor = rep(paste0("Donor", 5:8), times=2),
    Study = rep("Study2", 8)
)

metadata_study3 <- data.frame(
    SampleID = paste0("Sample", 17:24),
    Diagnosis = rep(c("Case", "Control"), each=4),
    Donor = rep(paste0("Donor", 9:12), times=2),
    Study = rep("Study3", 8)
)
```

``` r
metadata_study1
```

    ##   SampleID Diagnosis  Donor  Study
    ## 1  Sample1      Case Donor1 Study1
    ## 2  Sample2      Case Donor2 Study1
    ## 3  Sample3      Case Donor3 Study1
    ## 4  Sample4      Case Donor4 Study1
    ## 5  Sample5   Control Donor1 Study1
    ## 6  Sample6   Control Donor2 Study1
    ## 7  Sample7   Control Donor3 Study1
    ## 8  Sample8   Control Donor4 Study1

### 2. Análise de Expressão Diferencial Usando `limma-voom`

``` r
# Pacotes necessários
library(limma)
library(edgeR)
library(MASS)  # Para o pseudo count

# Função para realizar a análise de expressão diferencial
perform_DE_analysis <- function(counts, metadata) {
    # Merge counts and metadata
    combined_data <- merge(metadata, counts, by = "SampleID")
    rownames(combined_data) <- combined_data$SampleID
    
    # Preparar os dados de expressão e a matriz de design
    expr_data <- combined_data[, grepl("Gene", colnames(combined_data))]
    design <- model.matrix(~ Diagnosis, data = combined_data)
    
    # Criar o objeto DGEList
    dge <- DGEList(counts = t(expr_data))
    
    # Normalização usando TMM
    dge <- calcNormFactors(dge)
    
    # Transformação voom
    v <- voom(dge, design)
    
    # Ajuste do modelo linear
    fit <- lmFit(v, design)
    fit <- eBayes(fit)
    
    # Obter os resultados
    results <- topTable(fit, coef = "DiagnosisControl", adjust = "fdr", number = Inf)
    return(results)
}

# Realizar a análise de expressão diferencial para cada estudo
DE_results_study1 <- perform_DE_analysis(counts_study1, metadata_study1)
DE_results_study2 <- perform_DE_analysis(counts_study2, metadata_study2)
DE_results_study3 <- perform_DE_analysis(counts_study3, metadata_study3)
```

``` r
DE_results_study1
```

    ##             logFC  AveExpr          t   P.Value adj.P.Val         B
    ## Gene1 -0.07327334 17.86515 -0.2497132 0.8056367 0.8997907 -4.646540
    ## Gene2  0.04977992 18.45409  0.2121359 0.8343847 0.8997907 -4.680835
    ## Gene3  0.01905374 18.59079  0.1277148 0.8997907 0.8997907 -4.765059

### 3. Realizar a Meta-Análise

``` r
# Função para realizar meta-análise dos resultados de DE
meta_analysis <- function(DE_results_list) {
    genes <- rownames(DE_results_list[[1]])
    meta_results <- data.frame(Gene = genes, EffectSize = NA, StdErr = NA, p.value = NA)

    for (gene in genes) {
        effect_sizes <- sapply(DE_results_list, function(res) res[gene, "logFC"])
        variances <- sapply(DE_results_list, function(res) res[gene, "t"]^2)

        # Verificar se todos os estudos possuem os coeficientes de interesse
        valid_studies <- !is.na(effect_sizes)
        
        if (sum(valid_studies) > 1) {
            # Realizar meta-análise de efeitos aleatórios
            meta_res <- metagen(TE = effect_sizes[valid_studies], seTE = sqrt(variances[valid_studies]), sm = "SMD")
            meta_results[meta_results$Gene == gene, "EffectSize"] <- meta_res$TE.fixed
            meta_results[meta_results$Gene == gene, "StdErr"] <- meta_res$seTE.fixed
            meta_results[meta_results$Gene == gene, "p.value"] <- meta_res$pval.fixed
        }
    }

    meta_results$p.adj <- p.adjust(meta_results$p.value, method = "fdr")
    return(meta_results)
}

# Combinar resultados de DE para a meta-análise
DE_results_list <- list(DE_results_study1, DE_results_study2, DE_results_study3)
meta_results <- meta_analysis(DE_results_list)

# Visualizar os resultados da meta-análise
meta_results
```

    ##    Gene  EffectSize    StdErr   p.value     p.adj
    ## 1 Gene1 0.006095957 0.1715450 0.9716526 0.9716526
    ## 2 Gene2 0.027284262 0.1158034 0.8137366 0.9716526
    ## 3 Gene3 0.016970869 0.1266248 0.8933830 0.9716526

Interpretar os resultados da meta-análise de expressão diferencial
envolve entender como os tamanhos de efeito combinados e os valores de p
se traduzem em conclusões sobre a expressão gênica diferencial entre
casos e controles através de múltiplos estudos. Aqui está um guia
detalhado sobre como interpretar esses resultados:

### Resultados da Análise de Expressão Diferencial Individual

Para cada estudo individual (usando `limma-voom`), os resultados típicos
incluem:

- `logFC` (log fold-change): Representa a mudança logarítmica na
  expressão gênica entre os casos e controles. Valores positivos indicam
  aumento na expressão nos casos em comparação aos controles, e valores
  negativos indicam diminuição.

- `AveExpr` (expressão média): Expressão média do gene em todas as
  amostras.

- `t` (estatística t): Medida da relação entre o coeficiente estimado
  (logFC) e seu erro padrão.

- `P.Value` (valor p): Significância estatística do coeficiente
  estimado. Valores menores indicam maior evidência contra a hipótese
  nula.

- `adj.P.Val` (valor p ajustado): Valor p ajustado para múltiplos testes
  usando o método FDR (False Discovery Rate).

### Resultados da Meta-Análise

Após combinar os resultados dos estudos individuais, a tabela
`meta_results` inclui:

- `Gene`: Nome do gene.

- `EffectSize`: Tamanho do efeito combinado (geralmente uma média
  ponderada dos logFCs dos estudos individuais).

- `StdErr`: Erro padrão do tamanho do efeito combinado.

- `p.value`: Valor p da meta-análise, indicando a significância
  estatística do tamanho do efeito combinado.

- `p.adj`: Valor p ajustado para múltiplos testes na meta-análise.

  ``` r
  meta_results
  ```

      ##    Gene  EffectSize    StdErr   p.value     p.adj
      ## 1 Gene1 0.006095957 0.1715450 0.9716526 0.9716526
      ## 2 Gene2 0.027284262 0.1158034 0.8137366 0.9716526
      ## 3 Gene3 0.016970869 0.1266248 0.8933830 0.9716526

### Como Interpretar os Resultados da Meta-Análise

1.  **Tamanho do Efeito Combinado (`EffectSize`)**:

    - Um tamanho de efeito positivo indica que o gene é, em média, mais
      expresso nos casos do que nos controles através dos estudos.

    - Um tamanho de efeito negativo indica que o gene é, em média, menos
      expresso nos casos do que nos controles.

    - Quanto maior o valor absoluto do tamanho do efeito, mais forte é a
      diferença na expressão gênica entre os casos e controles.

2.  **Erro Padrão (`StdErr`)**:

    - Fornece uma medida da precisão do tamanho do efeito combinado.
      Valores menores indicam estimativas mais precisas.

3.  **Valor p (`p.value`)**:

    - Indica a significância estatística do tamanho do efeito combinado.

    - Valores p baixos (geralmente \< 0.05) sugerem que a diferença na
      expressão gênica é estatisticamente significativa.

4.  **Valor p Ajustado (`p.adj`)**:

    - Corrige os valores p para múltiplos testes, reduzindo a chance de
      falsos positivos.

    - Genes com valores `p.adj` abaixo de um limiar (geralmente 0.05)
      são considerados diferencialmente expressos de forma significativa
      em múltiplos estudos.

### Exemplo de Interpretação

Suponha que o resultado da meta-análise para um gene específico seja:

`Gene     EffectSize   StdErr   p.value   p.adj`

`Gene1    1.2          0.3      0.001     0.01`

- **Gene**: Gene1

- **EffectSize**: 1.2

  - Gene1 é mais expresso nos casos do que nos controles.

- **StdErr**: 0.3

  - A precisão da estimativa do tamanho do efeito é razoável.

- **p.value**: 0.001

  - A diferença na expressão é estatisticamente significativa.

- **p.adj**: 0.01

  - Mesmo após correção para múltiplos testes, a diferença na expressão
    de Gene1 entre casos e controles permanece estatisticamente
    significativa.

### Conclusão

Interpretar os resultados da meta-análise ajuda a identificar genes que
são consistentemente diferencialmente expressos entre casos e controles
através de múltiplos estudos. Esses genes podem ser candidatos para
estudos adicionais, validação experimental ou como biomarcadores
potenciais para a condição estudada.
