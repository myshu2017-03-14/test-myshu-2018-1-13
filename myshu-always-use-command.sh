#----------------------------------------+
#        qiime������ر�                 |
#----------------------------------------+
# open qiime1
source activate qiime1
# close qiime1
source deactivate

source activate qiime2-2017.10
source deactivate

#----------------------------------------+
#  ete3��ʹ��֮ǰ��Ҫ�趨��������        |
#----------------------------------------+
# Activate the environment 
export PATH="/analysis/software_han/software/metagenome-analysis-software/Miniconda2/bin:$PATH"

#----------------------------------------+
#        pavian������ر�                |
#----------------------------------------+
sudo R
pavian::runApp(host="0.0.0.0",port=5000) # ��ü���host
# ���ʵ�ַ ��http://192.168.1.30:5000/

#----------------------------------------+
#        blast������ȶ�����             |
#----------------------------------------+
nohup makeblastdb -in ncbi_swissprot_filter.fasta -dbtype prot -input_type fasta -hash_index -title swissprot -out swissprot -parse_seqids &
nohup blastp -query AFprotein.fasta -db swissprot -out AFSwiss_Prot.out -outfmt 6 -evalue 1e-5 -num_threads 6 -max_target_seqs 5 &
nohup blastn -query AF.fasta -db swissprot -out AF_blastn.out -outfmt 6 -evalue 1e-5 -num_threads 6 -max_target_seqs 5 &
#BLASTN title
Query 	Query_length	Subject	Subject_length	%Identity Length	Mis-match	Gap	q.Start	q.End	s.Start	s.End 	E-value	Bit score	
Query	Subject	%Identity 	Length	Mis-match	Gap	q.Start	q.End	s.Start	s.End 	E-value	Bit score
Query	Subject	Annotation	%Identity 	Length	Mis-match	Gap	q.Start	q.End	s.Start	s.End 	E-value	Bit score
query id, subject id, % identity, alignment length, mismatches, gap opens, q.start, q.end, s.start, s.end, evalue, bit score



# EXCEL ������Ԫ���ԡ������ָ��ϲ�
#HM.4.10,HM.4.1
=B60&","&C59

#14 barcode 
lbc41--lbc41 lbc9--lbc9 lbc89--lbc89 lbc1--lbc1 lbc2--lbc2 lbc49--lbc49 lbc57--lbc57 lbc10--lbc10 lbc33--lbc33 lbc25--lbc25 lbc17--lbc17 lbc81--lbc81 lbc65--lbc65 lbc73--lbc73

sudo apt-cache search <����> #�����鿴��ǰϵͳ�еİ�װ��

less -SN file #������ʾ�ļ�������ʾ�к�

jobs -l #��ʾ���н��̵�ID
wget -c file #�ϵ�����
#----------------------------------------+
#    shell��̼��� ���༭�ļ���          |
#----------------------------------------+
1�����ļ�ĩβ���뻻�з���
echo -e "\n" >>$GENE.1.gen.blastn_out
ֱ�� echo "" >> text.txt ����


2��sort����Shell��
sort -u file  #ȥ���ظ���
sort -r file  #�������
sort -r number.txt -o number.txt  #���������Դ�ļ��У�Դ�ļ�ֱ�ӱ�����
sort -n file  #������ֵ����Ĭ���ǰ����ַ�����
sort -t : -k 2  #-tָ���ָ����� ָ��TAB�ָ���sort -t $'\t' �� -kָ����
# ���Ա����У���һ�У��������򣬴ӵڶ��п�ʼ���յ�һ�н�������
head -1 $args && nawk 'NR>1' $args | sort -k 1 > $name\_plot.txt
# ���ݵ�һ��ȥ�أ���ͬ�ı����ڶ���ֵ�����Ǹ�
cat data.txt | sort -rnk2 | awk '{if (!keys[$1]) print $0; keys[$1] = 1;}'

3�������滻�ļ�����Shell�� 
for file in `ls *.ref`;do mv $file `echo $file|sed 's/multipleR-WHSW1700//g'`;done;

4����������ͬһ�ļ����е��ļ���shell��
#----------------------------------------
#!/bin/bash

IN=$1
for args in /analysis/software_han/1-data/HLA-datas/pacbio-data/HLA-raw-data/processed-raw-data/CCS-sequences/*
do
#        echo $args #����ļ�Ŀ¼������Ŀ¼+�ļ�
#	echo $(basename $args .fastq)  #��ȡ�ļ���
 	#�����滻
 	sed "s/year/${year}/g"
 	#��ȡ�ļ�ÿһ��
	cat $args | while read line
	do
	    echo "File:${line}"
	done
	# ����а����հ׷������ȡ��ʱ��հ׻��ɻ��У�����취
	FS_old=$IFS
	FS=$'\n'
	or line in  `cat  rollback_config`
	o
	   echo "$line"
	one;
	FS=$IFS_old
	

done
#------ �����ļ� ----------------
#!/usr/bin/env bash  # ע��ʹ�����ͷ����ʹ��(n++)���Ͳ�������
n=1
for i in S0671_06B_CHG028937-2017K14P150-1-LHT-CAGATC_L006_*
do
	n1=$(($n % 2))  #�ύ��SGE�ᱨ�� %2: syntax error: operand expected (error token is "%2")�����ķ������ڡ�%�����Ҽ��Ͽո�
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
��$0 otu_table_no_singletons_even_3870.biom ../../../103_16S_samples_output/mapping_file.txt ../otus/rep_set.tre 3870 1:17,2:28,3:30,4:28 . 

HELP
        exit 0
#fi
}
[ -z "$1" ] && help
[ "$1" = "-h" ] && help


# �ж��ļ��м��ļ��Ƿ����
#����ļ��в����ڣ������ļ���
if [ ! -d "/myfolder" ]; then
  mkdir /myfolder
fi
# -f �����ж� $file �Ƿ����
if [ ! -f "$file" ]; then
  touch "$file"
fi
#----------------------------------------

5��ͳ������
#ͨ���ı�$1����ָ����һ�У�Ĭ��$1Ϊ��һ��
����ͣ� cat you.txt |awk '{a+=$1}END{print a}'
����ƽ��ֵ��cat you.txt |awk '{a+=$1}END{print a/NR}'
�������ֵ��cat you.txt |awk 'BEGIN{a=0}{if ($1>a) a=$1 fi}END{print a}'  

#ͳ��fasta�ļ��ж���������
grep ��>�� S_pastorianus.fasta | wc �Cl
#ͳ��fastq�ļ���reads��Ŀ
awk '{s++}END{print s/4}' file.fastq

#ͳ���ļ������ļ������Լ��ļ��и���
#ͳ��ĳ�ļ������ļ��ĸ���
ls -l |grep "^-"|wc -l
#ͳ��ĳ�ļ�����Ŀ¼�ĸ���
ls -l |grep "^��"|wc -l
#ͳ���ļ������ļ��ĸ������������ļ������
ls -lR|grep "^-"|wc -l

#ͳ�Ƶ�ǰĿ¼�������ļ��еĴ�С
du -ah --max-depth=1

#ͳ��file2���У�file1��û�е���
grep -vFf file1 file2

# ͳ��ĳһ��ֵ�ظ�ֵ�ĸ���
cut -f 3 blastn_97_otus_uniq_anno_cov.out | sort | uniq -c >count.txt

6������perl�����滻ĳһ�е�����
#	-p ��ʾƥ��ÿһ��
# -e ��ʾ�������дΪһ��
# -i ��ʾ�滻���д��ԭ�ļ�
perl -p -i -e 's/^(>\S+?)\s.*/$1/g' combined_seqs.trim.fasta

7������ID��ȡFASTA�ļ��е�seq��
#����QIIME����е�extract_seqs_by_sample_id.py�ű���
extract_seqs_by_sample_id.py -i polished_assembly_r.fasta -o Monascus_scf397.fasta -s scaffold397
# samtools faidx
samtools faidx filter.fasta 8e470bfd-0b02-45c6-8897-e1c51196d487 >one.fasta
# ����file2�е�id��file1�н���ƥ�䣬ɸѡ��file1�а���file2��id������
grep -w -f file2.txt file1.txt > filtered.txt
# ����ID��ȡfastq����fasta�ļ��е�����
seqtk subseq in.fq name.lst > out.fq

8����ȡ�ļ�
#��ȡ�ļ���ĳ�е�ĳ������
#����ӵ�3�е���10��
sed -n '3,10p' myfile > newfile
# ��ȡĳ���ļ��е�ǰ10���ļ���cat
ls bacteria/*.fna | sed -n '1,10p' | xargs cat  >test


9������ĳһ�е�ֵ��ȡĳһ��Χ����($NF��ʾ���һ��)
awk -F '\t' '$NF>=28510120 && $NF<=33480577{print}'
awk -F '\t' '{if($NF>=28510120 && $NF<=33480577)print $NF}' #������һ��
# cutȡ���һ��
cat S08-C-1-K_kraken.labels | rev |cut -f 1 -d ";" | rev 


10��ɾ����Love��ͷ����
sed -i '/^Love/d' 1.txt
sed -i '3d' 1.txt

11���������г�·��
ll -R | grep "*.png"
find | grep "*.png"
find / -name gcc   #ȫ������

12��ȥ�ظ�
sort -n $file | uniq  
sort -n $file | awk '{if($0!=line)print; line=$0}'  
sort -n $file | sed '$!N; /^.?\n\1$/!P; D'  

#----------------------------------------+
#                 ubuntu�������         |
#----------------------------------------+
1���������Ƿ�װ
dpkg -s gridengine-exec 
dpkg -L gridengine-exec  #�г���װ·��




#----------------------------------------+
#               ��ȡ�ļ�                 |
#----------------------------------------+
1��ͳ��BAM�ļ�
samtools view -c test.bam 
#ͳ��������reads��Ŀ

~/software_han/samtools-1.4.1/bin/samtools stats test.bam
���õ����ļ��п�����ȡ��ص���Ϣ
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
#               �ļ�ת��                 |
#----------------------------------------+
1��fastqתfasta ��
cat file.fastq | paste - - - - | sed 's/^@/>/g'| cut -f1-2 | tr '\t' '\n' > file.fasta
#BAM����>SAM
samtools view -h -o out.sam out.bam
#SAM����>BAM
samtools view -bS out.sam >out.bam
-b ��˼ʹ���ʹBAM format
-S ��˼ʹ����ʹSAM�����@SQ ȱʣ�� Ҫд-t

2��genebank to gff ������bioperl�еĽű���
#����ֱ��������NCBI�е�Genome Data Viewer�����ص�GENBANK�ļ���*.flat��ʽ����ת��
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
3��

#----------------------------------------+
#    50�������ϱ��ػ�����antiSMASH       |
#----------------------------------------+

���ػ�����antiSMASH��Ŀǰ��50�ϰ�װ�ˣ�30�ϻ�û�а�װ��
1������ע����Ϣת���õ�genbank��ʽ�Ļ������ļ�
perl /home/myshu/output/program/gff_to_genbank.pl -i n -gff genome.gff3 �Cgenome genome.fasta -protein protein.fasta -out genome.genbank
2������antiSMASH������Ĳ�����50������run_antiSMASH.py -h��
nohup run_antiSMASH.py --cpus 8 --eukaryotic --smcogs --clusterblast -v --knownclusterblast --full-hmmer --asf --subclusterblast --all-orfs --input-type nucl --limit -1 --logfile P_new_out.log --outputfolder ./P_new_out ./PenAur.genbank &
 
