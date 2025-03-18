pacman::p_load(clusterProfiler)
library(GOSemSim)

if(FALSE){
  
  gene.df <- bitr(999, fromType = "ENTREZID",
                  toType = c("ENSEMBL", "SYMBOL"),
                  OrgDb = org.Hs.eg.db)
  
  ########################### GO
  # User can set computeIC=FALSE if they only want to use Wang’s method.
  mmGO = godata('org.Mm.eg.db', ont = "BP")
  hsGO = godata('org.Hs.eg.db', ont = "BP")
  
  # GO semantic similarity measurement
  goSim("GO:0004022", "GO:0005515", semData=hsGO, measure="Jiang")
  goSim("GO:0004022", "GO:0005515", semData=hsGO, measure="Wang")
  
  go1 = c("GO:0004022","GO:0004024","GO:0004174")
  go2 = c("GO:0009055","GO:0005515")
  mgoSim(go1, go2, semData=hsGO, measure="Wang", combine=NULL)
  
  
  # 2.4 Gene semantic similarity measurement
  #  GOSemSim implemented four combine methods, including max, avg, rcmax, and BMA
  GOSemSim::geneSim("2524", "10590", semData=hsGO, measure="Wang", combine="BMA")
  
  mgeneSim(genes=c("835", "5261","241", "994"),
           semData=hsGO, measure="Wang", combine="BMA", verbose=FALSE)
  
  # By default, godata function use ENTREZID as keytype, and the input ID type is ENTREZID. User can use other ID types such as ENSEMBL, UNIPROT, REFSEQ, ACCNUM, SYMBOL et al.
  
  # 2.5 Gene cluster semantic similarity measurement
  # GOSemSim also supports calculating semantic similarity between two gene clusters using clusterSim() function and measuring semantic similarity among multiple gene clusters using mclusterSim() function.
  gs1 <- c("835", "5261","241", "994", "514", "533")
  gs2 <- c("578","582", "400", "409", "411")
  clusterSim(gs1, gs2, semData=hsGO, measure="Wang", combine="BMA")
  
  x <- org.Hs.egGO
  hsEG <- mappedkeys(x)
  set.seed <- 123
  clusters <- list(a=sample(hsEG, 20), b=sample(hsEG, 20), c=sample(hsEG, 20))
  mclusterSim(clusters, semData=hsGO, measure="Wang", combine="BMA")
  
  
  ########################### DO
  # 3.1 DO (Disease Ontology (DO)) term semantic similarity measurement
  library(DOSE)
  a <- c("DOID:14095", "DOID:5844", "DOID:2044", "DOID:8432", "DOID:9146",
         "DOID:10588", "DOID:3209", "DOID:848", "DOID:3341", "DOID:252")
  b <- c("DOID:9409", "DOID:2491", "DOID:4467", "DOID:3498", "DOID:11256")
  
  doSim(a[1], b[1], measure="Wang")
  s = doSim(a, b, measure="Wang")
  s
  
  simplot(s,
          color.low="white", color.high="red",
          labs=TRUE, digits=2, labs.size=5,
          font.size=14, xlab="", ylab="")
  
  # 3.2 Gene semantic similarity measurement
  g1 <- c("84842", "2524", "10590", "3070", "91746")
  g2 <- c("84289", "6045", "56999", "9869")
  
  DOSE::geneSim(g1, g2, measure="Wang", combine="BMA")
  DOSE::geneSim("2524", "10590", measure="Wang", combine="BMA")
  
  # 3.3 Gene cluster semantic similarity measurement
  DOSE::clusterSim(g1, g2, measure="Wang", combine="BMA")
  
  g3 <- c("57491", "6296", "51438", "5504", "27319", "1643")
  clusters <- list(a=g1, b=g2, c=g3)
  DOSE::mclusterSim(clusters, measure="Wang", combine="BMA")
  
  
  ########################### MeSH
  # 4 MeSH semantic similarity analysis
  # MeSH (Medical Subject Headings) is the NLM (U.S. National Library of Medicine) controlled vocabulary used to manually index articles for MEDLINE/PubMed. MeSH is a comprehensive life science vocabulary. MeSH has 19 categories and MeSH.db contains 16 of them. That is:
  
  # Abbreviation	Category
  # A	Anatomy
  # B	Organisms
  # C	Diseases
  # D	Chemicals and Drugs
  # E	Analytical, Diagnostic and Therapeutic Techniques and Equipment
  # F	Psychiatry and Psychology
  # G	Phenomena and Processes
  # H	Disciplines and Occupations
  # I	Anthropology, Education, Sociology and Social Phenomena
  # J	Technology and Food and Beverages
  # K	Humanities
  # L	Information Science
  # M	Persons
  # N	Health Care
  # V	Publication Type
  # Z	Geographical Locations
  
  # MeSH terms were associated with Entrez Gene ID by three methods, gendoo, gene2pubmed and RBBH (Reciprocal Blast Best Hit).
  # Method	Way of corresponding Entrez Gene IDs and MeSH IDs
  # Gendoo	Text-mining
  # gene2pubmed	Manual curation by NCBI teams
  # RBBH	sequence homology with BLASTP search (E-value<10-50)
  
  # First, we need to load/fetch species-specific MeSH annotation database:
  library(AnnotationHub)
  library(MeSHDbi)
  ah <- AnnotationHub(localHub=TRUE)
  ah$species
  
  hsa <- query(ah, c("MeSHDb", "Homo sapiens"))
  file_hsa <- hsa[[1]]
  db <- MeSHDbi::MeSHDb(file_hsa)
  
  # The semantic data can be prepared by the meshdata() function:
  library(meshes)
  hsamd <- meshdata(db, category='A', computeIC=T, database="gendoo")
  # save(hsamd, file = "./sub_data/EGR1_study/SSA/hsamd.rda")
  
  # The meshSim() function is designed to measure semantic similarity between two MeSH term vectors.
  meshSim("D000009", "D009130", semData=hsamd, measure="Wang")
  meshSim(c("D001369", "D002462"), c("D017629", "D002890", "D008928"),
          semData=hsamd, measure="Wang")
  
  # 4.3 Gene semantic similarity measurement
  geneSim("241", "251", semData=hsamd, measure="Wang", combine="BMA")
  geneSim(c("241", "251"), c("835", "5261","241", "994"),
          semData=hsamd, measure="Wang", combine="BMA")
  
  
  
  
  
}
