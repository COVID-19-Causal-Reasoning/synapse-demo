##these are commands that looad the data from datahub to synapse

library(synapser)
library(readxl)
synLogin()

projectid='syn21788423'


storeDaProts<-function(fname,disease='MERS',cellType='CALU',study='MCL002'){
  library(dplyr)
    dat<-readxl::read_xlsx(fname,sheet='DA_Proteins_Only')
    cn<-c("Protein","REFERENCE","ACCESSION","SHORT_DESC","SHORT_DESCRIPTION","Num_Peptides","#TOTAL_PEPTIDES")
    mdat<-intersect(cn,colnames(dat))##get those colnams with metadata
  
    
      ##first try to auto-magically parse out column metadata
    updated.cn<-colnames(dplyr::select(dat,-mdat))%>%
      purrr::map_chr(~ stringr::str_replace(.x,"_RFP",""))%>%
      purrr::map_chr(~ stringr::str_replace(.x,'_vs_','vs'))
   # colnames(dat)<-c(mdat,updated.cn)
    updated.cn <- updated.cn%>%
      setNames(colnames(dplyr::select(dat,-mdat)))%>%
      purrr::map_dfr(~ stringr::str_split(.x,pattern='_',n=3,simplify=TRUE))

    tib.cn<-tibble::as_tibble(cbind(names(updated.cn),t(updated.cn)))%>%
      setNames(c('header','condition','time','type'))
    
    
    ##now tidy up original mdata, and join!
    tdf<-dat%>%tidyr::pivot_longer(cols=(-mdat),names_to="header",values_to = 'result',values_ptypes=list(result=character()))%>%
      left_join(tib.cn,by='header')
    tdf$result<-as.numeric(tdf$result)
    
    print(tdf)
    #now store table on synapse with metadata
    tab<-synapser::synBuildTable(paste(disease,cellType,study,'proteomics'),projectid,tdf)
    tab<-synapser::synStore(tab,used=list(disease=disease,cellType=cellType,dataType='proteomics',studyName=study))
    
    
    
}

storeDETranscripts<-function(fname,disease,cellType,study){

  library(dplyr)
  dat<-readxl::read_xlsx(fname,sheet='DE_Genes_Only')
  cn<-c("ProbeName","GeneSymbol","RefSeq","GeneName","ProbeID","Entrez Gene ID","Description")
  mdat<-intersect(cn,colnames(dat))##get those colnams with metadata
  
  
  ##first try to auto-magically parse out column metadata
  updated.cn<-colnames(dplyr::select(dat,-mdat))%>%
    purrr::map_chr(~ stringr::str_replace(.x,paste0(study,'_'),""))%>%
    purrr::map_chr(~ stringr::str_replace(.x,"_MERS","MERS"))%>%
    purrr::map_chr(~ stringr::str_replace(.x,"d3_5","d35"))%>%
    purrr::map_chr(~ stringr::str_replace(.x,'_vs_','vs'))
  # colnames(dat)<-c(mdat,updated.cn)
  updated.cn <- updated.cn%>%
    setNames(colnames(dplyr::select(dat,-mdat)))%>%
    purrr::map_dfr(~ stringr::str_split(.x,pattern='_',n=3,simplify=TRUE))
  
  tib.cn<-tibble::as_tibble(cbind(names(updated.cn),t(updated.cn)))%>%
    setNames(c('header','condition','time','type'))
  
  
  ##now tidy up original mdata, and join!
  tdf<-dat%>%tidyr::pivot_longer(cols=(-mdat),names_to="header",values_to = 'result',values_ptypes=list(result=character()))%>%
    left_join(tib.cn,by='header')
  tdf$result<-as.numeric(tdf$result)
  
  print(tdf)
  #now store table on synapse with metadata
  tab<-synapser::synBuildTable(paste(disease,cellType,study,'transcriptomics'),projectid,tdf)
  tab<-synapser::synStore(tab,used=list(disease=disease,cellType=cellType,dataType='rnaArray',studyName=study))
  
}

storeDaMetab<-function(fname,disease,cellType,study){
  library(dplyr)
  dat<-readxl::read_xlsx(fname,sheet='DA_Metabolites_Only')
  cn<-c("METABOLITE","KEGG","CAS","PubChem")
  mdat<-intersect(cn,colnames(dat))##get those colnams with metadata
  
  
  ##first try to auto-magically parse out column metadata
  updated.cn<-colnames(dplyr::select(dat,-mdat))%>%
    purrr::map_chr(~ stringr::str_replace(.x,"_RFP",""))%>%
    purrr::map_chr(~ stringr::str_replace(.x,'_vs_','vs'))
  # colnames(dat)<-c(mdat,updated.cn)
  updated.cn <- updated.cn%>%
    setNames(colnames(dplyr::select(dat,-mdat)))%>%
    purrr::map_dfr(~ stringr::str_split(.x,pattern='_',n=3,simplify=TRUE))
  
  tib.cn<-tibble::as_tibble(cbind(names(updated.cn),t(updated.cn)))%>%
    setNames(c('header','condition','time','type'))
  
  
  ##now tidy up original mdata, and join!
  tdf<-dat%>%tidyr::pivot_longer(cols=(-mdat),names_to="header",values_to = 'result',values_ptypes=list(result=character()))%>%
    left_join(tib.cn,by='header')
  tdf$result<-as.numeric(tdf$result)
  
  print(tdf)
  #now store table on synapse with metadata
  tab<-synapser::synBuildTable(paste(disease,cellType,study,'metabolomics'),projectid,tdf)
  tab<-synapser::synStore(tab,used=list(disease=disease,cellType=cellType,dataType='metabolomics',studyName=study))
  
  
}

#storeDaMetab('inst/MCL002_metab_stat.xlsx',disease='MERS',cellType='CALU',study="MCL002")
storeDaMetab('inst/ICL102_metab_stat.xlsx',disease='Influenza',cellType='CALU',study="ICL102")

##here are examples that have been run in the past.
#storeDETranscripts('inst/MCL001_mRNA_stat.xlsx',disease='MERS',cellType='CALU',study='MCL001')
#storeDETranscripts('inst/ICL102_mRNA_stat.xlsx',disease='Influenza',cellType='CALU',study='ICL102')

#storeDaProts('inst/MCL002_pro_stat.xlsx',disease='MERS',cellType='CALU',study='MCL002')
#storeDaProts('inst/ICL102_pro_stat.xlsx',disease='Influenza',cellType='CALU',study='ICL102')