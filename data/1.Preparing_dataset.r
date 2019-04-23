load("/home/lferrao/Dropbox/Doutorado_Tese/UChicago/Predicoes/pheno.geno.data/int/geno.pheno.int.fem.fes.Rdata")
snp.int = x.bayes.final
dim(snp.int)
int.fem = data.frame(population = "Int", Env = "fem",
                     production = prod.int.fem,
                     rust = rust.int.fem,
                     green = green.int.fem)
int.fes = data.frame(population = "Int", Env = "fes",
                     production = prod.int.fes,
                     rust = rust.int.fes,
                     green = green.int.fes)



head(snp.int)
snp.int[1:5,1:5]

write.csv(snp.int,
          file = "/home/lferrao/Dropbox/Meus_artigos_livros/(2017) Comp_Statistical_Methods_GS_CoffeaCanephora/Heredity/data_set/Pop_Int_SNPdata.csv")

write.csv(rbind(int.fem, int.fes),
          file = "/home/lferrao/Dropbox/Meus_artigos_livros/(2017) Comp_Statistical_Methods_GS_CoffeaCanephora/Heredity/data_set/Pop_Int_PHENOdata.csv")


load("/home/lferrao/Dropbox/Doutorado_Tese/UChicago/Predicoes/pheno.geno.data/prec/geno.pheno.prec.fem.fes.Rdata")

snp.prec = x.bayes.final
dim(snp.prec)
prec.fem = data.frame(population = "Prec", Env = "fem",
                     production = prod.prec.fem,
                     rust = rust.prec.fem,
                     green = green.prec.fem)
prec.fes = data.frame(population = "Prec", Env = "fes",
                     production = prod.prec.fes,
                     rust = rust.prec.fes,
                     green = green.prec.fes)
write.csv(snp.prec,
          file = "/home/lferrao/Dropbox/Meus_artigos_livros/(2017) Comp_Statistical_Methods_GS_CoffeaCanephora/Heredity/data_set/Pop_Prec_SNPdata.csv")

write.csv(rbind(prec.fem, prec.fes),
          file = "/home/lferrao/Dropbox/Meus_artigos_livros/(2017) Comp_Statistical_Methods_GS_CoffeaCanephora/Heredity/data_set/Pop_Prec_PHENOdata.csv")
