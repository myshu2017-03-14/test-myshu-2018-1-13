#----------------------------------------+
#        qiime启动与关闭                 |
#----------------------------------------+
# open qiime1
source activate qiime1
# close qiime1
source deactivate

source activate qiime2-2017.10
source deactivate

#----------------------------------------+
#  ete3包使用之前需要设定环境变量        |
#----------------------------------------+
# Activate the environment 
export PATH="/analysis/software_han/software/metagenome-analysis-software/Miniconda2/bin:$PATH"

#----------------------------------------+
#        pavian启动与关闭                |
#----------------------------------------+
sudo R
pavian::runApp(host="0.0.0.0",port=5000) # 最好加上host
# 访问地址 ：http://192.168.1.30:5000/

#----------------------------------------+
#        blast建库与比对命令             |
#----------------------------------------+
nohup makeblastdb -in ncbi_swissprot_filter.fasta -dbtype prot -input_type fasta -hash_index -title swissprot -out swissprot -parse_seqids &
nohup blastp -query AFprotein.fasta -db swissprot -out AFSwiss_Prot.out -outfmt 6 -evalue 1e-5 -num_threads 6 -max_target_seqs 5 &
nohup blastn -query AF.fasta -db swissprot -out AF_blastn.out -outfmt 6 -evalue 1e-5 -num_threads 6 -max_target_seqs 5 &
#BLASTN title
Query 	Query_length	Subject	Subject_length	%Identity Length	Mis-match	Gap	q.Start	q.End	s.Start	s.End 	E-value	Bit score	
Query	Subject	%Identity 	Length	Mis-match	Gap	q.Start	q.End	s.Start	s.End 	E-value	Bit score
Query	Subject	Annotation	%Identity 	Length	Mis-match	Gap	q.Start	q.End	s.Start	s.End 	E-value	Bit score
query id, subject id, % identity, alignment length, mismatches, gap opens, q.start, q.end, s.start, s.end, evalue, bit score



# EXCEL 两个单元格以“，”分隔合并
#HM.4.10,HM.4.1
=B60&","&C59

#14 barcode 
lbc41--lbc41 lbc9--lbc9 lbc89--lbc89 lbc1--lbc1 lbc2--lbc2 lbc49--lbc49 lbc57--lbc57 lbc10--lbc10 lbc33--lbc33 lbc25--lbc25 lbc17--lbc17 lbc81--lbc81 lbc65--lbc65 lbc73--lbc73

sudo apt-cache search <包名> #搜索查看当前系统中的安装包

less -SN file #按行显示文件，并显示行号

jobs -l #显示运行进程的ID
wget -c file #断点续传
#----------------------------------------+
#    shell编程技巧 （编辑文件）          |
#----------------------------------------+
1、在文件末尾插入换行符：
echo -e "\n" >>$GENE.1.gen.blastn_out
直接 echo "" >> text.txt 即可


2、sort排序（Shell）
sort -u file  #去除重复行
sort -r file  #降序输出
sort -r number.txt -o number.txt  #可以输出到源文件中，源文件直接被覆盖
sort -n file  #按照数值排序，默认是按照字符排序
sort -t : -k 2  #-t指定分隔符， 指定TAB分隔：sort -t $'\t' ； -k指定列
# 不对标题行（第一行）进行排序，从第二行开始按照第一列进行排序
head -1 $args && nawk 'NR>1' $args | sort -k 1 > $name\_plot.txt
# 根据第一列去重，相同的保留第二列值最大的那个
cat data.txt | sort -rnk2 | awk '{if (!keys[$1]) print $0; keys[$1] = 1;}'

3、批量替换文件名（Shell） 
for file in `ls *.ref`;do mv $file `echo $file|sed 's/multipleR-WHSW1700//g'`;done;

4、批量处理同一文件夹中的文件（shell）
#----------------------------------------
#!/bin/bash

IN=$1
for args in /analysis/software_han/1-data/HLA-datas/pacbio-data/HLA-raw-data/processed-raw-data/CCS-sequences/*
do
#        echo $args #输出文件目录下所有目录+文件
#	echo $(basename $args .fastq)  #提取文件名
 	#变量替换
 	sed "s/year/${year}/g"
 	#读取文件每一行
	cat $args | while read line
	do
	    echo "File:${line}"
	done
	# 如果行包括空白符，则读取的时候空白会变成换行，解决办法
	FS_old=$IFS
	FS=$'\n'
	or line in  `cat  rollback_config`
	o
	   echo "$line"
	one;
	FS=$IFS_old
	

done
#------ 两个文件 ----------------
#!/usr/bin/env bash  # 注意使用这个头才能使用(n++)；就不报错了
n=1
for i in S0671_06B_CHG028937-2017K14P150-1-LHT-CAGATC_L006_*
do
	n1=$(($n % 2))  #提交到SGE会报错 %2: syntax error: operand expected (error token is "%2")，更改方法是在“%”左右加上空格
	if [ $n1 -eq 0 ]
	then
		file2=$i
		name=$(basename $i _R2.fastq)
		fastqc -o . $name\_trim_R2.fq.gz
		((n++));
	else
		file1=$i
		((n++))
		#echo $n
	fi

done

#------ help ----------------
# echo the help if not input all the options
help()
{
#if [ $# -lt 3 ] ; then
cat <<HELP

USAGE: $0 otu_table_file map_file tree MAX_RARE_DEPTH COUNT output_dir
        or $0 -h  # show this message

EXAMPLE:
　$0 otu_table_no_singletons_even_3870.biom ../../../103_16S_samples_output/mapping_file.txt ../otus/rep_set.tre 3870 1:17,2:28,3:30,4:28 . 

HELP
        exit 0
#fi
}
[ -z "$1" ] && help
[ "$1" = "-h" ] && help


# 判断文件夹及文件是否存在
#如果文件夹不存在，创建文件夹
if [ ! -d "/myfolder" ]; then
  mkdir /myfolder
fi
# -f 参数判断 $file 是否存在
if [ ! -f "$file" ]; then
  touch "$file"
fi
#----------------------------------------

5、统计命令
#通过改变$1可以指定哪一列，默认$1为第一列
列求和： cat you.txt |awk '{a+=$1}END{print a}'
列求平均值：cat you.txt |awk '{a+=$1}END{print a/NR}'
列求最大值：cat you.txt |awk 'BEGIN{a=0}{if ($1>a) a=$1 fi}END{print a}'  

#统计fasta文件有多少条序列
grep ‘>’ S_pastorianus.fasta | wc –l
#统计fastq文件中reads数目
awk '{s++}END{print s/4}' file.fastq

#统计文件夹内文件个数以及文件夹个数
#统计某文件夹下文件的个数
ls -l |grep "^-"|wc -l
#统计某文件夹下目录的个数
ls -l |grep "^ｄ"|wc -l
#统计文件夹下文件的个数，包括子文件夹里的
ls -lR|grep "^-"|wc -l

#统计当前目录下所有文件夹的大小
du -ah --max-depth=1

#统计file2中有，file1中没有的行
grep -vFf file1 file2

# 统计某一列值重复值的个数
cut -f 3 blastn_97_otus_uniq_anno_cov.out | sort | uniq -c >count.txt

6、利用perl正则替换某一行的内容
#	-p 表示匹配每一行
# -e 表示后面程序写为一行
# -i 表示替换结果写回原文件
perl -p -i -e 's/^(>\S+?)\s.*/$1/g' combined_seqs.trim.fasta

7、根据ID提取FASTA文件中的seq（
#利用QIIME软件中的extract_seqs_by_sample_id.py脚本）
extract_seqs_by_sample_id.py -i polished_assembly_r.fasta -o Monascus_scf397.fasta -s scaffold397
# samtools faidx
samtools faidx filter.fasta 8e470bfd-0b02-45c6-8897-e1c51196d487 >one.fasta
# 根据file2中的id在file1中进行匹配，筛选出file1中包含file2中id的序列
grep -w -f file2.txt file1.txt > filtered.txt
# 根据ID提取fastq或者fasta文件中的序列
seqtk subseq in.fq name.lst > out.fq

8、提取文件
#提取文件的某行到某行区间
#比如从第3行到第10行
sed -n '3,10p' myfile > newfile
# 提取某个文件夹的前10个文件并cat
ls bacteria/*.fna | sed -n '1,10p' | xargs cat  >test


9、根据某一列的值提取某一范围的行($NF表示最后一行)
awk -F '\t' '$NF>=28510120 && $NF<=33480577{print}'
awk -F '\t' '{if($NF>=28510120 && $NF<=33480577)print $NF}' #输出最后一列
# cut取最后一列
cat S08-C-1-K_kraken.labels | rev |cut -f 1 -d ";" | rev 


10、删除以Love开头的行
sed -i '/^Love/d' 1.txt
sed -i '3d' 1.txt

11、搜索并列出路径
ll -R | grep "*.png"
find | grep "*.png"
find / -name gcc   #全局搜索

12、去重复
sort -n $file | uniq  
sort -n $file | awk '{if($0!=line)print; line=$0}'  
sort -n $file | sed '$!N; /^.?\n\1$/!P; D'  

#----------------------------------------+
#                 ubuntu相关命令         |
#----------------------------------------+
1、检查软件是否安装
dpkg -s gridengine-exec 
dpkg -L gridengine-exec  #列出安装路径




#----------------------------------------+
#               读取文件                 |
#----------------------------------------+
1、统计BAM文件
samtools view -c test.bam 
#统计所含的reads数目

~/software_han/samtools-1.4.1/bin/samtools stats test.bam
最后得到的文件中可以提取相关的信息
# This file was produced by samtools stats (1.4.1+htslib-1.4.1) and can be plotted using plot-bamstats
# This file contains statistics for all reads.
# The command line was:  stats 05-DPA1_02_02_02_gen_ref.bwa.sort.bam
# CHK, Checksum	[2]Read Names	[3]Sequences	[4]Qualities
# CHK, CRC32 of reads which passed filtering followed by addition (32bit overflow)
# Summary Numbers. Use `grep ^SN | cut -f 2-` to extract this part.
# First Fragment Qualitites. Use `grep ^FFQ | cut -f 2-` to extract this part.
# Columns correspond to qualities and rows to cycles. First column is the cycle number.
# Last Fragment Qualitites. Use `grep ^LFQ | cut -f 2-` to extract this part.
# Columns correspond to qualities and rows to cycles. First column is the cycle number.
# GC Content of first fragments. Use `grep ^GCF | cut -f 2-` to extract this part.
# GC Content of last fragments. Use `grep ^GCL | cut -f 2-` to extract this part.
# ACGT content per cycle. Use `grep ^GCC | cut -f 2-` to extract this part. The columns are: cycle; A,C,G,T base counts as a percentage of all A/C/G/T bases [%]; and N and O counts as a percentage of all A/C/G/T bases [%]
# Insert sizes. Use `grep ^IS | cut -f 2-` to extract this part. The columns are: insert size, pairs total, inward oriented pairs, outward oriented pairs, other pairs
# Read lengths. Use `grep ^RL | cut -f 2-` to extract this part. The columns are: read length, count
# Indel distribution. Use `grep ^ID | cut -f 2-` to extract this part. The columns are: length, number of insertions, number of deletions
# Indels per cycle. Use `grep ^IC | cut -f 2-` to extract this part. The columns are: cycle, number of insertions (fwd), .. (rev) , number of deletions (fwd), .. (rev)
# Coverage distribution. Use `grep ^COV | cut -f 2-` to extract this part.
# GC-depth. Use `grep ^GCD | cut -f 2-` to extract this part. The columns are: GC%, unique sequence percentiles, 10th, 25th, 50th, 75th and 90th depth percentile


#----------------------------------------+
#               文件转换                 |
#----------------------------------------+
1、fastq转fasta ：
cat file.fastq | paste - - - - | sed 's/^@/>/g'| cut -f1-2 | tr '\t' '\n' > file.fasta
#BAM＝＝>SAM
samtools view -h -o out.sam out.bam
#SAM＝＝>BAM
samtools view -bS out.sam >out.bam
-b 意思使输出使BAM format
-S 意思使输入使SAM，如果@SQ 缺剩， 要写-t

2、genebank to gff （利用bioperl中的脚本）
#可以直接利用在NCBI中的Genome Data Viewer中下载的GENBANK文件，*.flat格式进行转换
perl bp_genbank2gff3.pl NC_000006.11\[28477797..33448354\].flat
--------------------------------------------------------------------------------------------
Usage:
      bp_genbank2gff3.pl [options] filename(s)

      # process a directory containing GenBank flatfiles
      perl bp_genbank2gff3.pl --dir path_to_files --zip

      # process a single file, ignore explicit exons and introns
      perl bp_genbank2gff3.pl --filter exon --filter intron file.gbk.gz

      # process a list of files 
      perl bp_genbank2gff3.pl *gbk.gz

      # process data from URL, with Chado GFF model (-noCDS), and pipe to database loader
      curl ftp://ftp.ncbi.nih.gov/genomes/Saccharomyces_cerevisiae/CHR_X/NC_001142.gbk \
      | perl bp_genbank2gff3.pl -noCDS -in stdin -out stdout \
      | perl gmod_bulk_load_gff3.pl -dbname mychado -organism fromdata

        Options:
            --noinfer  -r  don't infer exon/mRNA subfeatures
            --conf     -i  path to the curation configuration file that contains user preferences
                           for Genbank entries (must be YAML format)
                           (if --manual is passed without --ini, user will be prompted to 
                            create the file if any manual input is saved)
            --sofile  -l  path to to the so.obo file to use for feature type mapping
                           (--sofile live will download the latest online revision)
            --manual   -m  when trying to guess the proper SO term, if more than
                           one option matches the primary tag, the converter will 
                           wait for user input to choose the correct one
                           (only works with --sofile)
            --dir      -d  path to a list of genbank flatfiles
            --outdir   -o  location to write GFF files (can be 'stdout' or '-' for pipe)
            --zip      -z  compress GFF3 output files with gzip
            --summary  -s  print a summary of the features in each contig
            --filter   -x  genbank feature type(s) to ignore
            --split    -y  split output to separate GFF and fasta files for
                           each genbank record
            --nolump   -n  separate file for each reference sequence
                           (default is to lump all records together into one 
                           output file for each input file)
            --ethresh  -e  error threshold for unflattener
                           set this high (>2) to ignore all unflattener errors
            --[no]CDS  -c  Keep CDS-exons, or convert to alternate gene-RNA-protein-exon 
                           model. --CDS is default. Use --CDS to keep default GFF gene model, 
                           use --noCDS to convert to g-r-p-e.
            --format   -f  Input format (SeqIO types): GenBank, Swiss or Uniprot, EMBL work
                           (GenBank is default)
            --GFF_VERSION  3 is default, 2 and 2.5 and other Bio::Tools::GFF versions available
            --quiet        don't talk about what is being processed 
            --typesource   SO sequence type for source (e.g. chromosome; region; contig)
            --help     -h  display this message
------------------------------------------------------------------------------------------------------------------
3、

#----------------------------------------+
#    50服务器上本地化运行antiSMASH       |
#----------------------------------------+

本地化运行antiSMASH（目前在50上安装了，30上还没有安装）
1、根据注释信息转换得到genbank格式的基因组文件
perl /home/myshu/output/program/gff_to_genbank.pl -i n -gff genome.gff3 –genome genome.fasta -protein protein.fasta -out genome.genbank
2、运行antiSMASH（具体的参数在50上运行run_antiSMASH.py -h）
nohup run_antiSMASH.py --cpus 8 --eukaryotic --smcogs --clusterblast -v --knownclusterblast --full-hmmer --asf --subclusterblast --all-orfs --input-type nucl --limit -1 --logfile P_new_out.log --outputfolder ./P_new_out ./PenAur.genbank &
 
