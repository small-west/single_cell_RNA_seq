## 绘制单细胞一个测序数据覆盖基因的深度图
+ 由于单细胞测序时采用的protocol差异，可能存在覆盖度偏好性
+ 从5'端到3'端，将每个基因上测序数据的深度图曲线画出来，等比缩放这些基因长度，在一个坐标轴内放下一个细胞内所有基因的深度曲线，最后将所有基因的深度曲线取均值，查看测序是否存在3'偏好性或5'偏好性

### 从gtf文件中提取基因位置的信息
```
grep -v ^# gencode.vM25.annotation.gtf |cut -f 1,4,5,9 |cut -f1 -d ";"|awk '{print $1,$2,$3,$5}'|sed -e 's/ /\t/g' |sed -e 's/\"//g'>gene_id.bed
sort -k1,1 -k2,2n gene_id.bed >sorted.bed
```

### 使用RseQC画出bias分布图
RseQC的geneBody_coverage.py模块能展示reads在转录本上的覆盖情况，用于评估测序数据是否存在5'或3' bias。需要的输入文件：
1. 经过sorted的bam文件；
2. 记录基因位置信息gene model的bed文件，可以从[RseQC官网](http://rseqc.sourceforge.net/#download-gene-models-update-on-08-07-2014)下载对应版本的bed文件。
```
geneBody_coverage.py -r /data/neurosys-svr2/zhouxueya/refer/mm10_Gencode_VM18.bed -i LCM-FF-1_trimed.bam,LCM-FF-2_trimed.bam,LCM-FF-3_trimed.bam,LCM-FF-4_trimed.bam,LCM-FF-5_trimed.bam -o bias_statustics/  
```

### 结果分析
+ 覆盖度曲线呈现两头高中间低的趋势

## bam文件coverage的统计
## coverage的两个层面
+ 覆盖度
+ 深度

## 统计Intronic reads, exonic reads, intergenic reads, ribosomal RNA reads的coverage
+ 使用qualimap对bam文件进行质控  
`qualimap rnaseq -bam LCM-FF-1_trimed.bam -gtf /data/neurosys-svr2/zhouxueya/refer/gencode.vM25.annotation.gtf -pe -s --java-mem-size=8G`
