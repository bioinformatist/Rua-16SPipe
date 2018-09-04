# Rua-16SPipe

Rua-16SPipe是一个基于[Nextflow](https://www.nextflow.io/)的16S扩增子测序数据分析流程，尚在建设中。

> Rua！我是小僵尸我自豪！
> 别说了，枫哥不喜欢。

## 示例数据

目前，Rua-16SPipe允许两种类型的输入文件，分别有对应的[示例数据](data)：

1. [加了Barcode之后混样测序，并未demux（最常见情况）](data/mixed)：

使用了来自已发表paper的数据，经过[SeqKit](https://github.com/shenwei356/seqkit)作如下处理：

```shell
# 对所有fastq文件使用SeqKit抽样2%的reads
for x in *.fastq; do seqkit sample $x -p 0.0118121931188733 -o $x.0.02.fastq; done
# 将全部结果文件重命名
rename fastq. "" *0.02*
```

2. [Demuxed的数据](data/demuxed)：

数据来源于(*acc number和ref待补充*)，经过[SeqKit](https://github.com/shenwei356/seqkit)作如下处理：

```shell
# 对所有fastq文件使用SeqKit抽样2%的reads
for x in *.fastq; do seqkit sample $x -p 0.02 -o $x.0.02.fastq; done
# 将全部结果文件重命名
rename fastq. "" *0.02*
```

由于SeqKit默认有固定的random seed，因此对于双端顺序一致的数据，是可以这样取得序列名称配对完整的可用子集的。

> 插播一条广告：[SeqKit](https://github.com/shenwei356/seqkit)是由作者的[小伙伴](https://github.com/shenwei356)使用Go语言开发的，快到你根本感受不到时间的流逝，欢迎踊跃下载试用。

## 环境配置

关于运行流程所需的mothur等需要root账户安装的（或最好以root安装的）依赖程序，以及其他难以解决的依赖，请联系服务器管理员阅读[main.nf](main.nf)后处理。

> 作者将会在v0.1.0后开始提供docker image，你将可以借助docker实现完全零依赖的数据分析，敬请关注。

该流程不需要任何安装步骤:smile:。nextflow会自动从本仓库拉取对应的脚本文件等。

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

这说明你的人品出现了问题。**请去马路边强行扶老奶奶过几次马路（过来过去），并小心不要将她碰到在地。如果她强行摔倒，请转身快速跑掉。然后回来电脑旁重试，直到成功为止。**

> 开始使用之前请注意：暂时请使用相对路径（终端内使用`pwd`可以查看工作目录）指定输入文件（目录）的位置。作者在v0.1.0会解决这个问题。

### 获取最新的流程

Rua-16SPipe还在活跃开发中。请在每次使用之前运行`nextflow pull bioinformatist/Rua-16SPipe`来获取最新版本。

### 对于混样数据

命令行参数举例：

```shell
nextflow run bioinformatist/Rua-16SPipe --reads_mode mixed --ffastq data/mixed/SRR5834618_1.0.02.fastq --rfastq data/mixed/SRR5834618_2.0.02.fastq --oligo data/mixed/oligos_yu.txt
```

### 对于demuxed数据

命令行参数举例：

```shell
nextflow run bioinformatist/Rua-16SPipe --reads_mode demuxed --inputdir data/demuxed
```

## 运行结果

程序运行结束后，当前目录下会出现一个results目录，保存了所有结果。

## 建议与反馈

使用过程中遇到任何问题，或者要求增加功能特性等，请[留言](https://github.com/bioinformatist/Rua-16SPipe/issues) (English is better)。