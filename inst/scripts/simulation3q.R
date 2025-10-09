library(dplyr)
library(magrittr)
library(stringr)
library(tidyr)
library(optparse)
library(readr)

option_list = list(
  make_option(c("--input"), type="character", default=NULL,
              help="input vcf file name", metavar="character"),

  make_option(c("--region"), type="character", default="chr1",
              help="chromosome selected for the UPD", metavar="character"),

  make_option(c("--pos1"), type="numeric", default=NULL,
              help="start of the interval", metavar="character"),

  make_option(c("--pos2"), type="numeric", default=NULL,
              help="end of the interval", metavar="character"),

  make_option(c("--out"), type="character", default="simulation.vcf",
              help="output vcf file name [default= %default]", metavar="character"),

  make_option(c("--type"), type="character", default="heterodisomy",
              help="type of UPD [default= %default]", metavar="character"),

  make_option(c("--parent"), type="character", default="father",
              help="select father or mother inheritance [default= %default]", metavar="character")
);





opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
tmp_vcf<-readLines(opt$input)
########################################################################
########################################################################
########################################################################



vcf<-as.data.frame(tmp_vcf)

header<-as.data.frame(vcf[grepl("^#", vcf$tmp_vcf), ])
body<-as.data.frame(vcf[!grepl("^#", vcf$tmp_vcf), ])
colnames(body)
colnames(body)<-c("one")
body_2<-body %>% separate(one, c("chr","pos", "omit","ref","alt","nothing","filter","info","format","geno_child","geno_father","geno_mother"), "\\t")
body_2$pos<-as.numeric(body_2$pos)
body_2$coordenadas<-paste0(body_2$chr,body_2$pos,body_2$ref,body_2$alt)
child<-body_2[,c(1,2,3,4,5,6,7,8,9,10,13)]
father<-body_2[,c(1,2,3,4,5,6,7,8,9,11,13)]
mother<-body_2[,c(1,2,3,4,5,6,7,8,9,12,13)]


print("todo ok")

subset_coord_arms<-data.frame(region=opt$region,start=opt$pos1,end=opt$pos2)


#######################################################
#######################################################
###########heterodisomy################################
#######################################################
#######################################################

subset_father<-subset(father,father$chr==opt$region)
subset_mother<-subset(mother,mother$chr==opt$region)


subset_father<-subset(subset_father,as.numeric(subset_father$pos)>=subset_coord_arms$start)
subset_father<-subset(subset_father,as.numeric(subset_father$pos)<=subset_coord_arms$end)




subset_mother<-subset(subset_mother,as.numeric(subset_mother$pos)>=subset_coord_arms$start)
subset_mother<-subset(subset_mother,as.numeric(subset_mother$pos)<=subset_coord_arms$end)


print("variantes:")
head(subset_father)



write.table(subset_father,file=paste0(opt$out,"_coordinates.csv"),quote = FALSE)

subset_2_father<-subset_father
subset_2_mother<-subset_mother


if (opt$type=="heterodisomy") {

  if (opt$parent=="father") {
    subset_2_child<-filter(child,coordenadas %in% subset_2_father$coordenadas)
    subset_2_child$geno_child[match(subset_2_father$coordenadas,subset_2_child$coordenadas)] <- subset_2_father$geno_father
    dat <- merge(body_2, subset_2_child,by = "coordenadas", all.x = TRUE)
    dat <- dat[match(body_2$coordenadas, dat$coordenadas),]
    dat_2<-transform(dat, geno_child.x = ifelse(is.na(geno_child.y),
    geno_child.x, as.character(geno_child.y)))
    body_def<-dat_2[,1:13]
  }

  if (opt$parent=="mother") {
    subset_2_child<-filter(child,coordenadas %in% subset_2_mother$coordenadas)
    subset_2_child$geno_child[match(subset_2_mother$coordenadas,subset_2_child$coordenadas)] <- subset_2_mother$geno_mother
    dat <- merge(body_2, subset_2_child, by = "coordenadas", all.x = TRUE)
    dat <- dat[match(body_2$coordenadas, dat$coordenadas),]
    dat_2<-transform(dat, geno_child.x = ifelse(is.na(geno_child.y),
    geno_child.x, as.character(geno_child.y)))
    body_def<-dat_2[,1:13]
  }

  rm(subset_2_child)
  rm(subset_2_father)
  rm(subset_father)
  rm(subset_mother)
  rm(dat)
  rm(dat_2)
  rm(mother)
  rm(father)
  rm(child)
  rm(body_2)


}



body_def$column_def<-paste(body_def$chr.x,body_def$pos.x,body_def$omit.x,body_def$ref.x,body_def$alt.x,body_def$nothing.x,body_def$filter.x,
                           body_def$info.x,body_def$format.x,body_def$geno_child.x,body_def$geno_father,body_def$geno_mother,sep = "\t")
body_def<-body_def[,14]
body_def<-as.data.frame(body_def)

names(header) <- c("body_def")
vcf_final<-rbind(header,body_def)
write_tsv(vcf_final,opt$out,col_names = FALSE)
exe<-paste0("sed -i \'s/\"\"/\"/g\' ",opt$out)
system2(exe)


