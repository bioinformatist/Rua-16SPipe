# Rua-16SPipe

Rua-16SPipe是一个基于[Nextflow](https://www.nextflow.io/)的16S扩增子测序数据分析流程，尚在建设中。

## 示例数据

目前，Rua-16SPipe允许两种类型的输入文件，分别有对应的示例数据：

1. 加了Barcode之后混样测序，并未demux（最常见情况）：

使用了来自已发表paper的数据，经过[SeqKit](https://github.com/shenwei356/seqkit)作如下处理：

for x in *.fastq; do seqkit sample $x -p 0.0118121931188733 -o $x.0.02.fastq; done

2. Demuxed的数据：

[数据](data/demuxed)来源于(*acc number和ref待补充*)，经过[SeqKit](https://github.com/shenwei356/seqkit)作如下处理：

```shell
# 对所有fastq文件使用SeqKit抽样2%的reads
for x in *.fastq; do seqkit sample $x -p 0.02 -o $x.0.02.fastq; done
# 将全部结果文件重命名
rename fastq. "" *0.02*
```

由于SeqKit默认有固定的random seed，因此对于双端顺序一致的数据，是可以这样取得可用子集的。

## 安装

## 使用

使用`nextflow main.nf`来运行流程。第一次在你的账户下使用nextflow的时候，会见到如下提示：

```pre
Downloading nextflow dependencies. It may require a few seconds, please wait ..
```

这是nextflow在自动下载其所需要的依赖（某些JVM相关运行环境），耐心等待即可。

如果见到了如下奇怪的错误信息：

```pre
ERROR: Cannot download nextflow required file -- make sure you can connect to the internet
```

请去马路边强行扶老奶奶过几次马路（过来过去），并小心不要碰倒她。然后回来重试直到成功为止。