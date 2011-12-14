funct.enrich <- function () {
  library(plyr)
  ################### Data section ###################
  # read pathway data file
    kegg.data = read.csv("Dropbox/dvh/dvu_kegg/dvu_pathway.tab.txt", header=F, sep="\t")
    colnames(kegg.data) = c("gene", "pathway")
  # read network siff file
    net.data =read.csv("~/Dropbox/dvh/DvHmotifNetwork.sif", sep="\t")
    colnames(net.data) = c("gene", "interaction", "bicluster")
  # read kegg annotations
    kegg.annot.data = read.csv("~/Dropbox/dvh/dvu_kegg/kegg_pathway_ids.csv", header=T, sep=",")
    colnames(kegg.annot.data) = c("pathway", "pathway.2", "pathway.1", "pathway.0")
  # read kegg orthologies 
    kegg.ko.data = read.csv("~/Dropbox/dvh/dvu_kegg/dvu_kegg_orthology_mapping.csv", header=T, sep=",")
  # read GO ontologies 
    go.data = read.csv("~/Dropbox/dvh/dvu_kegg/dvu_go_mapping.csv", header=T, sep=",")
    go.annot.data = read.csv("~/Dropbox/dvh/dvu_kegg/go_terms_ids.csv", header=T, sep=",")
    go.secondary = read.csv("Dropbox/dvh/dvu_kegg/go_secondary_ids.csv", header=T, sep=",")
  # read TIGR roles  
    tigr.data = read.csv("~/Dropbox/dvh/dvu_kegg/dvu_cog_tigr.csv", header=T, sep=",")
  
  ##### data merging
    #joined.data = join (net.data, kegg.data, by = "gene")
    merged.data = merge (net.data, kegg.data, by = "gene", all.x = T, all.y = T)
    merged.annotation = merge (merged.data, kegg.annot.data, by = "pathway", all.x = T, all.y = T, sort=F)
    merged.annotation.ko = merge (merged.data, kegg.ko.data, by = "gene", all.x = T, all.y = T, sort=T)
    merged.go = merge (go.data, go.annot.data, by = "go", all.x = T, all.y = F, sort=F)
    merged.go.secondary = merge (merged.go, go.secondary, by = "go", all.x = T, all.y = F, sort=F)
  
  ########### Kegg Functional Enrichment ####################
  ############                            ####################
  ############ by using hypergeometric    ####################
  # Collect KEGG targets
    kegg.targets = list()
    for(kegg.i in unique(kegg.data[,"pathway"])) {
      kegg.targets[[kegg.i]] = kegg.data[which(kegg.data[,"pathway"]==kegg.i),"gene"]
    }
  
  # Collect Bicluster members
    bicluster.members = list()
    for(bicluster.i in unique(net.data[,"bicluster"])) {
      bicluster.members[[bicluster.i]] = net.data[which(net.data[,"bicluster"]==bicluster.i),"gene"]
    }
      
  # Calculate hypergeometric p values
    pValues = matrix(ncol=length(kegg.targets),nrow=length(bicluster.members), dimnames= list(names(bicluster.members), names(kegg.targets)))
    for(kegg.j in names(kegg.targets)) {
      for(bicluster.j in names(bicluster.members)) {
            pValues[bicluster.j, kegg.j] = phyper(length(intersect(bicluster.members[[bicluster.j]],kegg.targets[[kegg.j]])), # q
                                              length(kegg.targets[[kegg.j]]), # m
                                              3490-length(kegg.targets[[kegg.j]]), # n
                                              length(bicluster.members[[bicluster.j]]), # k
                                              lower.tail=F)
    }
  }
  
  #Bonferroni correction 1
  #p.bonferroni = p.adjust(pValues, method='bonferroni')
  
  # Bonferroni correction
  unique.kegg.number = length(unique(kegg.data[,"pathway"]))
  unique.bicluster.number = length(unique(net.data[,"bicluster"]))
      
  L = pValues <= 0.05/(397*87)
  RN = matrix(rep(rownames(pValues), length(colnames(pValues))), nrow=unique.bicluster.number)
  CN = matrix(rep(colnames(pValues), length(rownames(pValues))), ncol=unique.kegg.number, byrow=T)
  filtered.results = data.frame(bicluster=RN[L], pathway=CN[L], p.value=pValues[L])
    
  annot.enrichment = merge (filtered.results, kegg.annot.data, by = "pathway", all.x = T, all.y = F)
  annot.enrichment2 = merge (annot.enrichment, merged.annotation.ko, by = "pathway", all.x = T, all.y = F)
  
  ############ TIGR Functional Enrichment ####################
  ############                            ####################
  ############ by using hypergeometric    ####################
      
  # ---------------------- TIGR Role.0 Enrichment ------------------------------------
    tigr.targets = list()
    for(tigr.i in unique(tigr.data[,"TIGRRoles.0"])) {
      if (! is.na(tigr.i)) {
        tigr.targets[[tigr.i]] = tigr.data[which(tigr.data[,"TIGRRoles.0"]==tigr.i),"gene"]
        }
      }
  
  # Bicluster members
    bicluster.members = list()
    for(bicluster.i in unique(net.data[,"bicluster"])) {
      bicluster.members[[bicluster.i]] = net.data[which(net.data[,"bicluster"]==bicluster.i),"gene"]
    }
  # hypergeometric p values
    pValues = matrix(ncol=length(tigr.targets),nrow=length(bicluster.members), dimnames= list(names(bicluster.members), names(tigr.targets)))
    for(tigr.j in names(tigr.targets)) {
      for(bicluster.j in names(bicluster.members)) {
            pValues[bicluster.j, tigr.j] = phyper(length(intersect(bicluster.members[[bicluster.j]],tigr.targets[[tigr.j]])), # q
                                              length(tigr.targets[[tigr.j]]), # m
                                              3490-length(tigr.targets[[tigr.j]]), # n
                                              length(bicluster.members[[bicluster.j]]), # k
                                              lower.tail=F)
    }
  }
  
  # Bonferroni correction
  unique.tigr.number = length(unique(tigr.data[,"TIGRRoles.0"]))
  unique.bicluster.number = length(unique(net.data[,"bicluster"]))
      
  L = pValues <= 0.05/(unique.bicluster.number*unique.tigr.number)
  RN = matrix(rep(rownames(pValues), length(colnames(pValues))), nrow=unique.bicluster.number)
  CN = matrix(rep(colnames(pValues), length(rownames(pValues))), ncol=unique.tigr.number, byrow=T)
  filtered.tigr.results = data.frame(bicluster=RN[L], TIGRRoles=CN[L], p.value=pValues[L])
    
  # ---------------------- TIGR Role.1 Enrichment ------------------------------------
    tigr.1.targets = list()
    for(tigr.1.i in unique(tigr.data[,"TIGRRoles.1"])) {
      if (! is.na(tigr.1.i)) {
        tigr.1.targets[[tigr.1.i]] = tigr.data[which(tigr.data[,"TIGRRoles.1"]==tigr.1.i),"gene"]
        }
      }
  # hypergeometric p values
    pValues = matrix(ncol=length(tigr.1.targets),nrow=length(bicluster.members), dimnames= list(names(bicluster.members), names(tigr.1.targets)))
    for(tigr.1.j in names(tigr.1.targets)) {
      for(bicluster.j in names(bicluster.members)) {
            pValues[bicluster.j, tigr.1.j] = phyper(length(intersect(bicluster.members[[bicluster.j]],tigr.1.targets[[tigr.1.j]])), # q
                                              length(tigr.1.targets[[tigr.1.j]]), # m
                                              3490-length(tigr.1.targets[[tigr.1.j]]), # n
                                              length(bicluster.members[[bicluster.j]]), # k
                                              lower.tail=F)
    }
  }
  
  # Bonferroni correction
  unique.tigr.1.number = length(unique(tigr.data[,"TIGRRoles.1"]))
  unique.bicluster.number = length(unique(net.data[,"bicluster"]))
    
  L = pValues <= 0.05/(unique.bicluster.number*unique.tigr.1.number)
  RN = matrix(rep(rownames(pValues), length(colnames(pValues))), nrow=unique.bicluster.number)
  CN = matrix(rep(colnames(pValues), length(rownames(pValues))), ncol=unique.tigr.1.number, byrow=T)
  filtered.tigr.1.results = data.frame(bicluster=RN[L], TIGRRoles.1=CN[L], p.value=pValues[L])
}