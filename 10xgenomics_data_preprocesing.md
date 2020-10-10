## 10 x genomics 单细胞测序数据分析
10 x genomics的大体流程：样本->文库制备->测序->raw reads->alignment->count  
样本、细胞、转录本的区分：样本拥有样本barcode，细胞有细胞的barcode，转录本有独特的UMI
### 测序数据的理解和质控
+ 测序数据由两个fastq文件组成，文件名以R1和R2结尾
+ 参考10 x genomics官网10x chromium v3的帮助文档，R1由28bp的barcode+UMI组成，R2真正的3'端转录本
+ fastq对数据进行质控
  `fastqc -o fasqc_res/ /picb/neurosys/LJ/BigData/GLHu/20191218/*.fastq.gz`
  
### cell ranger处理数据
+ cell ranger的安装
+ 从download页面填写信息获取链接

+ cell ranger参考序列的准备
  人和小鼠的参考序列可以直接从官网下载得到
  ```
  #下载fasta文件和GTF文件，并解压
  wget ftp://ftp.ensembl.org/pub/release-101/fasta/macaca_fascicularis/dna/
  wget ftp://ftp.ensembl.org/pub/release-101/gtf/macaca_fascicularis/Macaca_fascicularis.Macaca_fascicularis_5.0.101.chr.gtf.gz
  gunzip *.gz
  #gtf文件过滤，原始GTF文件包含了所有的非poly-A的转录本，有很多重叠的注释，会造成reads被标记为multi-mapped
  cellranger mkgtf \
  Macaca_fascicularis.Macaca_fascicularis_5.0.101.chr.gtf \
  Macaca_fascicularis.Macaca_fascicularis_5.0.101.chr.filter.gtf \
  --attribute=gene_biotype:protein_coding
  #制作参考序列
  cellranger mkref \
  --genome=Macaca_fascicularis_genome \
  --fasta=../Macaca_fascicularis.Macaca_fascicularis_5.0.dna.toplevel.fa \
  --genes=Macaca_fascicularis.Macaca_fascicularis_5.0.101.chr.filter.gtf 
  #制作完成，可以在随后的cell-ranger使用中指定为
  #transcriptome参数cellranger --transcriptome=/data/neurosys-svr2/zhouxueya/macaca/cell_ranger_refer/Macaca_fascicularis_genome
  ```
  
+ cell ranger count
  注意修改fastq文件的名字
  参考官网修改  
  ```
  ls /data/neurosys-svr2/GH/20191218/
  19X118-2_S1_L005_R1_001.fastq.gz  19X118-2_S1_L005_R1_003.fastq.gz  19X118-2_S1_L005_R2_001.fastq.gz  19X118-2_S1_L005_R2_003.fastq.gz  md5_GLHu.txt
  19X118-2_S1_L005_R1_002.fastq.gz  19X118-2_S1_L005_R1_004.fastq.gz  19X118-2_S1_L005_R2_002.fastq.gz  19X118-2_S1_L005_R2_004.fastq.gz
  ```
 
  注：一个样本的fastq文件放在同一个目录下
  ```
  refer_path=/data/neurosys-svr2/zhouxueya/macaca/cell_ranger_refer/Macaca_fascicularis_genome/
  fastq_path=/data/neurosys-svr2/GH/20191218/
  cellranger count --id=19X118 \
                 --transcriptome=$refer_path \
                 --fastqs=$fastq_path \
                 --sample=19X118-2
  ```
