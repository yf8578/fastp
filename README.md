# fastp
##相关参数描述
Fastp

filter out bad reads (too low quality, too short, or too many N…)

对每一个序列的头部或尾部,计算滑动窗内的质量均值，并将均值较低的子序列进行切除

trim all reads in front and tail

对所有reads 头部和尾部进行裁剪

cut adapters. Adapter sequences can be automatically detected, which means you don’t have to input the adapter sequences to trim them

自动去除adapter

correct mismatched base pairs in overlapped regions of paired end reads, if one base is with high quality while the other is with ultra low quality

当同一碱基一个高质量，一个低质量时，自动纠正低质量碱基

trim polyG in 3’ ends, which is commonly seen in NovaSeq/NextSeq data. Trim polyX in 3’ ends to remove unwanted polyX tailing (i.e. polyA tailing for mRNA-Seq data)
去除尾部的polyG

preprocess unique molecular identifier (UMI) enabled data, shift UMI to sequence name.

可以对带分子标签（UMI）的数据进行预处理，不管UMI在插入片段还是在index上

report JSON format result for further interpreting.
JSON格式的报告

visualize quality control and filtering results on a single HTML page (like FASTQC but faster and more informative).

HTML报告

split the output to multiple files (0001.R1.gz, 0002.R1.gz…) to support parallel processing. Two modes can be used, limiting the total split file number, or limitting the lines of each split file

拆分fq,使运行更快

support long reads (data from PacBio / Nanopore devices).

支持长reads

12 .support reading from STDIN and writing to STDOUT

参数

1.输入输出

fastp -i in.R1.fq.gz -o out.R1.fq.gz -I in.R2.fq.gz -O out.R2.fq.gz
1
2.输出unpaired reads

--unpaired1   # 输出在read1中未被过滤，read2中被过滤的reads，默认关闭**
--unpaired2   # 输出在read2中未被过滤，read1中被过滤的reads
1
2
3.输出过滤掉的reads

--failed_out  # 输出被过滤掉的reads，并给出过滤原因，例：failed_quality_filter,                                 failed_too_short。 如果没有指定（--unpaired1，--unpaired2），则unpaired reads                 过滤原因写为paired_read_is_failing 
1
4.提取部分reads进行分析

--reads_to_process # 设定需要分析的reads数，default=0,表示全部分析
1
5.是否覆盖已存在的文件

--dont_overwrite  # 加入该参数后 如果有已存在的任何fastp文件，则会中断，不会覆盖，默认关闭
1
filtering
1. quality filter （质量过滤）
从 N含量，低质量碱基比例，质量得分 三方面进行过滤

-Q, --disable_quality_filtering  # 关闭 quality filter，默认开启
-n, --n_base_limit               # 限制 N 个数
-q, --qualified_quality_phred    # 碱基质量值 default >=15
-u, --unqualified_percent_limit  # 允许的低质量碱基所占比例，default=40. 表示低质量碱基比例为40%                                    大于40则该条reads被丢弃
-e, --average_qual               # 平均质量得分阈值 如果reads 质量得分小于平均值则该reads被丢弃。                                          default=0 ，表示关闭



2. length filter （长度过滤）
从reads 最小，最大长度 两方面进行过滤

-L, --disable_length_filtering   # 关闭 length filter ，默认开启
-l, --length_required            # 设定最小reads长度
    --length_limit               # 设定最大reads长度，大于length_limit的会过滤掉。 0表示无限制



3. low complexity filter （低复杂度过滤）
complexity 定义：

当前碱基与下一个碱基不同（base[i] != base[i+1]），符合这种条件的碱基所占的比例 称为complexity 。

例：

# a 51-bp sequence, with 3 bases that is different from its next base
seq = 'AAAATTTTTTTTTTTTTTTTTTTTTGGGGGGGGGGGGGGGGGGGGGGCCCC'
complexity = 3/(51-1) = 6%

-y, --low_complexity_filter      # 开启过滤， 默认关闭
-Y， --complexity_threshold      # 设定complexity阈值（0-100），default=30

adapters
Adapter trimming 默认开启。通过检测reads间的overlap来检测adapter。

-A, --disable_adapter_trimming  # 关闭 Adapter trimming ，默认开启
--detect_adapter_for_pe         # 开启自动检测adapter，默认关闭。
--adapter_sequence              # 指定read1的adapter，如果fastp没有自动识别出adapter,则会使用指                                   定的adapter_sequence
--adapter_sequence_r2           # 指定read2的adapter，如果fastp没有自动识别出adapter,则会使用指                                   定的adapter_sequence_r2
--adapter_fasta                 # 或者提供FASTA格式文件，切除的特定序列

fastp首先裁剪自动检测的adapter，或者通过--adapter_sequence，--adapter_sequence_r2给定的序列，最后裁剪--adapter_fasta 给定的序列

#FASTA文件格式
>Illumina TruSeq Adapter Read 1
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
>Illumina TruSeq Adapter Read 2
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
>polyA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

per read cutting by quality score
通过质量值过滤每条reads

滑框方式，计算滑框平均质量得分，切除低于阈值的框

前端N也被裁掉

有三种切除方式，可单选或全选，默认关闭，这三种方式会影响最终的reads去重结果。

-5, --cut_front                    # 从5’端开始,若框内碱基平均质量得分小于阈值，则删掉
    --cut_front_window_size        # 设定窗口长度，default: 4
    --cut_front_mean_quality       # 设定裁剪阈值，default: 20 (Q20)
-3, --cut_tail                     # 从3’端开始,若框内碱基平均质量得分小于阈值，则删掉
    --cut_tail_window_size         # 设定窗口长度，default: 4
    --cut_tail_mean_quality        # 设定裁剪阈值，default: 20 (Q20)
-r, --cut_right                    # 从5’端开始，切除低于阈值的框，以及右侧所有序列
    --cut_right_window_size        # 设定窗口长度，default: 4
    --cut_right_mean_quality       # 设定裁剪阈值，default: 20 (Q20)
-W, --cut_window_size              # 滑框大小 Range: 1~1000, default: 4 
-M, --cut_mean_quality             # 平均质量得分 Range: 1~36 default: 20 (Q20)
#If you don't set window size and mean quality threshold for these function respectively, fastp will use the values from -W, --cut_window_size and -M, --cut_mean_quality

global trimming
max_len 是根据read最终长度计算，

-f, --trim_front1 and -t, --trim_tail1  # read1  front/tail trim 碱基数，default=0
-F, --trim_front2 and -T, --trim_tail2  # read2  front/tail trimming 若未设定，按照read1执行
-b, --max_len1                          # 按照read1最大长度trim
-B, --max_len2                          # 若未设定，按照read1执行
#如果 read1长度大于 --max_len1 ，则从tail端，按照--max_len1 进行trim.


–max_len 是根据read最终长度计算，可能影响read长度的参数：

1, UMI preprocessing (--umi)
2, global trimming at front (--trim_front)
3, global trimming at tail (--trim_tail)
4, quality pruning at 5' (--cut_front)
5, quality pruning by sliding window (--cut_right)
6, quality pruning at 3' (--cut_tail)
7, trim polyG (--trim_poly_g, enabled by default for NovaSeq/NextSeq data)
8, trim adapter by overlap analysis (enabled by default for PE data)
9, trim adapter by adapter sequence (--adapter_sequence, --adapter_sequence_r2. For PE data, this step is skipped if last step succeeded)
10, trim polyX (--trim_poly_x)
11, trim to max length (---max_len)

polyG tail trimming
去除read尾端由于光信号识别错误产生的G。NextSeq/NovaSeq 默认开启（通过FASTQ中machine ID识别）

-g or --trim_poly_g             # 开启
-G or --disable_trim_poly_g     # 关闭
--poly_g_min_len                # 设置最小长度，default=10

polyX tail trimming
默认关闭。polyx即polyA，通常在 mRNA-Seq reads 尾部找到。

当与polyG tail trimming 同时存在时，优先执行polyG tail trimming 。

-x or --polyX                  # 开启，默认关闭
--poly_x_min_len               # 设置最小长度，default=10
1
2
unique molecular identifier (UMI) processing（分子标签UMI处理）
对带UMI的FASTQ文件进行预处理

fastp会将标签提取出来，放在reads表头的第一部分，并且在sam/bam文件中也能体现。如果标签在reads中，会导致reads变短。

-U or --umi            # 开启
--umi_loc              # 指定标签位置（index1,index2,read1,read2,per_index,index1_index2
                          per_read）#If --umi_loc is specified with read1, read2 or                                per_read, the length of UMI should specified with --umi_len. 
--umi_prefix           # 指定输出标签的前缀，前缀与标签之间使用 _ 连接 例：UMI_AATTCCGG 
--umi_skip             # 如果标签在 read1/read2/per_read 上，fastp可以跳过标签后的碱基来trim标签
                       # --umi_skip 指定跳过的碱基数，默认关闭


output splitting
拆分输出文件，并行处理节约后续分析时间

两种拆分方式：文件个数；文件行数

--out1 or --out2                 # 拆分文件命名
-d or --split_prefix_digits      # 文件名前缀长度
-s or --split                    # 按照文件个数拆分
-S or --split_by_lines           # 按照文件行数拆分

overrepresented sequence analysis
过表达序列分析

-p, --overrepresentation_analysis  # 启用该分析，默认统计序列长度为10bp, 20bp, 40bp, 100bp或                                        cycle -2
-P, --overrepresentation_sampling  # 指定用于统计的reads数比例，默认20，即默认1/20的reads用于序列                                      统计. 1 表示将所有reads用于统计

fastp不仅给出了过表达序列的计数，还给出了它们如何在cycle中分布的信息。为每个检测到的过表达序列提供了一个数字，从中可以知道该序列在哪里最常见。

merge paired-end reads
合并PE序列overlap部分，overlap部分写入–merged_out。

--merged_out               # 指定输出文件
--out1 and --out2          # 表示没有合并成功的reads，并且reads1,reads2都filter
--unpaired1                # 未合并的reads，read1 passes filters but read2 doesn't.
--unpaired2                # 未合并的reads，read2 passes filters but read1 doesn't.
--include_unmerged         # 将所有不能合并的reads输出到一个文件中，默认关闭
--failed_out               # 存放所有没有filter的reads

其他参数

-6, --phred64              # 使用phred64质量体系（会被转换成33）
-z, --compression          # gzip压缩比例，1-9.1最快，9最慢，default=4
--stdin                    # 通过STDIN形式作为输入，如果为PE数据，需写入 --interleaved_in 参数
--interleaved_in           # 表示输入含有read1,read2
--stdout                   # 通过STDOUT形式输出，默认关闭

#filter reads with unwanted indexes (to remove possible contamination)
通过不需要的索引过滤reads(去除污染)
--filter_by_index1         # 指定包含指定index1的文件，一行一个index
--filter_by_index2         # 指定包含指定index2的文件，一行一个index
--filter_by_index_threshold # 允许的index错配碱基数， default 0

#base correction by overlap analysis options
碱基矫正
-c, --correction           # 开启碱基矫正（PE）
--overlap_len_require      # overlap的最小长度。 这将影响基于overlap的PE merge, adapter                                    trimming and correction. 30 by default.
--overlap_diff_limit       # 识别overlap时允许的最大错配数。这将影响基于overlap的PE merge,                                    adapter trimming and correction。 5 by default.
--overlap_diff_percent_limit  # 识别overlap时允许的最大错配碱基比例。 Default 20.

#reporting options
  -j, --json               # json格式报告
  -h, --html               # html格式报告
  -R, --report_title       # 报告标题
  
 #threading options
  -w, --thread             # 线程 , default is 2
