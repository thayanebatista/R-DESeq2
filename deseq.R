library(DESeq2)
library(TCGAbiolinks)
library(BiocParallel)
library(circlize)

data <-  #arquivo com dados geneticos (por exemplo baixados do banco TCGA)

type <- 'gene' # insira o tipo aqui

fileName <-  paste("/", type, ".rda", sep = '') # insira a pasta em '/' ou apenas mantenha como está e ficara salvo no seu Working Directory atual

# verifica se o arquivo com o deseq calculado já existe, caso não exista realiza o calculo
if (!file.exists(fileName) || reexecute) {
  
  # manipula os dados para testar
  if (F) {
    
    # mantem somente 10 linhas para teste
    data = data[1:10,]
    
    # nomes das amostras nao tumorais
    samplesNT = TCGAquery_SampleTypes(colnames(data), "NT")
    
    # nomes das amostras tumorais
    samplesTP = TCGAquery_SampleTypes(colnames(data), "TP")
    
    # mantem apenas 5 amostras tumorais e 5 nao tumorais
    samplesNT = samplesNT[1:5]
    samplesTP = samplesTP[1:5]
    
    data = data[, c(samplesNT, samplesTP)]
    
    data1 = data[, samplesNT]
    #data1 = apply(data1, 1:2, function (x) 500)
    
    data2 = data[, samplesTP]
    #data2 = apply(data2, 1:2, function (x) 1500)
    
    data = as.matrix(cbind(data1, data2))
    
    env = globalenv()
    
    # cria uma tabela com a condicao de cada amostra
    env$sampleCondition = data.frame(condition = c(rep("NT", length(samplesNT)), rep("TP", length(samplesTP))), case = colnames(data))
    rownames(env$sampleCondition) = colnames(data)
  }
  
  # prepara o deseq
  ddsHTSeq = DESeqDataSetFromMatrix(countData = data, 
                                    colData = sampleCondition, 
                                    design = ~condition)
  
  # mantem somente as linhas com leituras relevantes
  ddsHTSeq <- ddsHTSeq[rowSums(counts(ddsHTSeq)) >= 10,]
  
  # define qual é a condicao de controle
  ddsHTSeq$condition <- relevel(ddsHTSeq$condition, ref = "NT")
  
  # calcula o deseq
  dds = DESeq(ddsHTSeq, parallel = T)
  
  # calcula o vsd
  vsd = varianceStabilizingTransformation(dds, blind=TRUE)
  
  # resultados do deseq
  res = results(dds, parallel = T)
  
  # cria a variavel deseq com todos os valores calculados
  deseq = list(dds = dds, vsd = vsd, res = res)
  
  # persiste o deseq calculado
  save(deseq, file = fileName, compress = "xz")
}