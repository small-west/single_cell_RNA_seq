# 单细胞测序上游数据分析
## 概述
+ 单细胞测序技术主要包括以下流程：组织解离得到单细胞悬液，细胞裂解，RNA逆转录成cDNA，PCR扩增，高通量测序，数据分析等。  
+ 单细胞解离的三种方法：
  1. 人工显微操作
  2. LCM激光捕获显微切割
  3. 荧光激活细胞分选（FACS）  

+ 单细胞测序的数据分析主要包括以下三步：
1. 计算基因表达矩阵。根据测序reads上的barcode和UMI标签将reads比对到特定细胞的特定基因上并计数，以获得每个细胞中不同基因的表达量；
2. 质控。去除基因表达量很少和线粒体DNA含量较高的细胞；
3. 数据降维和聚类。通过主成分分析（Principal components analysis）及其他一些方法对基因表达数据进行降维，然后通过迭代性聚类分析对细胞进行分型。



## 预处理和可视化
原始测序数据经过预处理得到分子计数矩阵（count matrix）。这取决于单细胞文库构建方案中是否包含UMI  
这次的单细胞数据是直接由LCM显微切割使用samrt-seq2 protocol测得的5个细胞的数据，获得的是分开的5个细胞的fastq文件，而非10xgenomics的多个细胞放在一个fastq文件的情况，可以直接进行数据分析。  
与转录组的上游分析流程一样。

+ 用fastqc进行质量控制，multiqc统计结果。  
`fastqc seqfile1 seqfile2 -o fastqc_result`  
`multiqc fastqc_result_dir` 

+ 用trim-galore去接头
```
Trim_galore(){
    trim_galore -q 20 --phred33 --stringency 3 --length 20 -e 0.1 \
	--paired $1 $2 \
	--gzip -o trim_galore/
}
#$1、$2是双端测序的fastq文件
```

* 用STAR进行alignment，把reads比对到reference参考基因组序列上，对于每个read数据，STAR会找到和一个或多个参考基因组上的序列匹配的最长的序列。用这个方法STAR也能识别剪切事件。STAR的问题在于它需要大量的RAM，尤其当参考基因组很大的时候，比如小鼠和人。  
* 主要包括两个步骤，首先需要提供参考基因组（FASTA）和注释文件（GTF），STAR将创建一个基因组索引，其次，STAR会把提供的reads数据比对到基因组索引上。
* 由于STAR构建索引需要的内存很大，采用hisat2进行比对

### **STAR比对**
首先，我们需要建立索引，下载参考基因组序列，然后用STAR建立索引  

+ 下载参考基因组序列([UCSCmm10参考序列链接](ftp://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/))([Gencode参考序列链接](https://www.gencodegenes.org/mouse/))
```shell
wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/chromFa.tar.gz' -O chromFa.tar.gz
tar -zxvf chromFa.tar.gz
cat *.fa > mm10.fa
rm -rf chr*
#参考序列存放在/data/neurosys-svr2/zhouxueya/refer目录下
#从gencode下载的fasta参考序列也在这个目录下GRCm38.p6.genome.fa
```
+ 下载注释文件（gtf文件）([Gencode注释文件链接](https://www.gencodegenes.org/mouse/))
```
wget --timestamping 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz'
#注释文件放在/data/neurosys-svr2/zhouxueya/refer目录下
```
+ 用STAR建立索引
```
#workdir:/data/neurosys-svr2/zhouxueya/refer
mkdir index && cd index
mkdir star
STAR --runThreadN 6 --runMode genomeGenerate \
    --genomeDir /data/neurosys-svr2/zhouxueya/refer/index/star/ \
    --genomeFastaFiles /data/neurosys-svr2/zhouxueya/refer/GRCm38.p6.genome.fa \
    --sjdbGTFfile /data/neurosys-svr2/zhouxueya/refer/gencode.vM25.annotation.gtf.gz \
    --sjdbOverhang 150
#--sjdbOverhang这个参数最好是测序读长，从fastqc报告里可以得到
```
(需要的内存很大，建立索引中途屡次失败，换成hisat2试一下)

### **hisat2比对**
+ 建立索引 hisat2官网有小鼠和人的索引文件可以直接下载  
  ```
  wget --timestamping 'https://genome-idx.s3.amazonaws.com/hisat/mm10_genome.tar.gz'
  tar -zxvf *.tar.gz
  rm -rf *.tar.gz 
  ```
+ hisat2比对
  ```
  for i in {1..5}
  do
  echo "start mapping LCM-FF-${i}_R1_val_1.fq.gz LCM-FF-${i}_R2_val_2.fq.gz";
  hisat2 -t -x /data/neurosys-svr2/zhouxueya/refer/index/hisat2_index/mm10/genome \
  	-1 /home/zhouxueya/trim_galore/LCM-FF-${i}_R1_val_1.fq.gz \
  	-2 /home/zhouxueya/trim_galore/LCM-FF-${i}_R2_val_2.fq.gz \
  	-S /home/zhouxueya/alignment/LCM-FF-${i}_trimed.sam
  done
  #samtools transvert sam to bam and sorted and index
  #while read 按行读取文件
  ls *.sam | while read id; do (samtools sort -O bam -@ 5 -o $(basename $id ".sam").bam ${id});done
  rm *.sam
  #|的作用是把管道左侧的标准输出作为右侧命令的标准输入，但有些命令不接受标准输入作为参数
  #xargs的作用是将标准输入转化为命令行参数
  #xargs一般和|连用，xargs -i的作用是把前面的输出作为后面的参数，-i后面可以不指定名字。默认是{}，相当于变量替换，会把后面的{}替换成前面的每一个输出，-p可以输出最终会执行的命令，询问用户是否执行
  ls *.bam |xargs -i -p samtools index {}
  ```
+ bam文件质控
  ```
  geneBody_coverage.py -r /data/neurosys-svr2/zhouxueya/refer/mm10_Gencode_VM18.bed -i LCM-FF-1_trimed.bam,LCM-FF-2_trimed.bam,LCM-FF-3_trimed.bam,LCM-FF-4_trimed.bam,LCM-FF-5_trimed.bam -o bias_statustics/
  ```
### **featurecount定量**
  ```
  featureCounts -T 6 -p -g gene_id -a /data/neurosys-svr2/zhouxueya/refer/gencode.vM25.annotation.gtf -o count/expr_matrix *.bam
  ```
