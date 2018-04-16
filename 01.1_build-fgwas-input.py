#!/usr/bin/python -O
# Jason Matthew Torres
'''
Input a gwas bed file and an annotation bed file
Usage: python JTbuild_fgwas_bedfile.py
'''
# libraries
import sys,os,gzip
import subprocess as sp

# globals
bedtools = "/apps/well/bedtools/2.24.0/bedtools"
cur_dir = "/well/mccarthy/users/jason/projects/ucsd_collaboration/"
input_dir = cur_dir + "fgwas_input/"
fgwas_head_list = ["CHR","POS0","POS","SNPID","F","Z","PVAL","NCASE","NCONTROL"] # fgwas field names for gwas columns
annot_bed_file = "/well/mccarthy/users/jason/datasets/bing/bing.bed"
gwas_bed_file = "/well/got2d/jason/reference/gwas/diamante-ukbb_hrc/ukbb_diamante-euro.bed"
outfile = input_dir+"ukbb_diamante-euro.fgwas"

print "Be sure to load R module before running script: module load R/3.3.1"


# functions

def get_annot_list(bed_file):
    '''
    Expecting no header
    '''
    sys.stdout.write("\nGetting list of annotations...\n")
    fin = open(bed_file,'r')
    annot_list = []
    count = 0
    for line in fin:
        count += 1
        sys.stdout.write("\r%d" % count)
        sys.stdout.flush()
        l = line.strip().split()
        annot_name = l[3]
        if annot_name not in annot_list:
            annot_list.append(annot_name)
    fin.close()
    sys.stdout.write("\nThere are %d annotations\n" % len(annot_list))
    return annot_list

def config_dic(bed_file):
    annot_list = get_annot_list(bed_file)
    sys.stdout.write("\nConfiguring dictionary...\n")
    dic = {}
    for annot in annot_list:
        dic[annot] = {}
    return annot_list, dic


def bed_intersect(gwas_bed,annot_bed):
    sys.stdout.write("\nBedtools intersect...\n")
    command_list1 = ["cat",gwas_bed,"|","cut","-f","1,2,3,4",">",input_dir+"temp.bed"]
    sp.check_call(" ".join(command_list1),shell=True)
    command_list2 = [bedtools, "intersect", "-a", input_dir+"temp.bed", "-b", annot_bed,
                    "-wb", ">", input_dir+"intersect.temp"]
    sp.check_call(" ".join(command_list2),shell=True)

def build_dic(dic):
    fin = open(input_dir+"intersect.temp",'r')
    for line in fin:
        l = line.strip().split()
        chrom,pos1,pos2,annot = l[0],l[1],l[2],l[7]
        snpid = chrom+":"+pos2
        dic[annot][snpid] = snpid
    fin.close()
    return dic

def build_matrix_file(gwas_bed, annot_list, dic, fgwas_head_list):
    '''
    Expecting four columns in gwas_bed, chrom, pos1, pos2, name
    '''
    sys.stdout.write("\nWriting Matrix file\n")
    fin = open(gwas_bed,'r')
    fout = open(input_dir+"intermediate.fgwas",'w')
    head_list = fgwas_head_list + annot_list
    fout.write("\t".join(head_list)+"\n")
    count = 0
    for line in fin:
        count += 1
        sys.stdout.write("\r%d" % count)
        sys.stdout.flush()
        l = line.strip().split()
        snpid = l[0]+":"+l[2]
        for annot in annot_list:
            try:
                len(dic[annot][snpid])
                l.append("1")
            except:
                l.append("0")
        fout.write("\t".join(l)+"\n")
    fin.close()
    fout.close()
    command_list = ["gzip", input_dir+"intermediate.fgwas"]
    sys.stdout.write("\nCompressing file...\n")
    sp.check_call(" ".join(command_list),shell=True)

def process_file(infile,outfile):
    '''
    Will use R to sort, remove duplicate SNPIDS, and determine distance to TSS
    '''
    r_file = cur_dir+"01.process.R"
    fout=open(r_file,'w')
    script='''
library("data.table")
library("dplyr")
library(GenomicRanges)

save.dir <- "%s"

command <- paste0("cat ","%s"," | zmore")
df <- fread(command)
print("Sorting file...")
chrom <- gsub("chr","",df$CHR)
df <- cbind(chrom,df)
df <- arrange(df,chrom,POS)
df <- dplyr::select(df,-chrom)
print(str(df))

print("Removing duplicates...")
df <- df[!duplicated(df$SNPID),]
print(str(df))

# TxDb for genomic features
library("GenomicFeatures")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")

print("Getting Distance to Nearest TSS for each SNP")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
trans <- transcripts(txdb)
snp.gr <- GRanges(df$CHR,IRanges(df$POS0, df$POS))
tss <- resize(trans, width=1, fix='start') # get TSS from transcripts
dist <- distanceToNearest(x=snp.gr,subject=tss)
distance_tss <- elementMetadata(dist)$distance
df$distance_tss <- distance_tss

write.table(df,"%s",sep="\t",row.names=F,quote=F)

    ''' % (cur_dir, infile, outfile)
    fout.write(script)
    fout.close()
    call = ["Rscript","--vanilla",r_file]
    sp.check_call(" ".join(call),shell=True)
    command_list = ["gzip", outfile]
    sys.stdout.write("\nCompressing file...\n")
    sp.check_call(" ".join(command_list),shell=True)


def main():
    bed_intersect(gwas_bed_file,annot_bed_file)
    mylist, mydic = config_dic(annot_bed_file)
    mydic = build_dic(mydic)
    build_matrix_file(gwas_bed_file,mylist,mydic,fgwas_head_list)
    process_file(input_dir+"intermediate.fgwas.gz",outfile)

if (__name__=="__main__"): main()
