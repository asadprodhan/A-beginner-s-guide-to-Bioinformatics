<h1 align="center">A beginner's guide to Bioinformatics</h1>

<h3 align="center">M. Asaduzzaman Prodhan<sup>*</sup> </h3>

<div align="center"><b> School of Biological Sciences, The University of Western Australia </b></div>
<div align="center"><b> 35 Stirling Highway, Perth, WA 6009, Australia. <sup>*</sup>Correspondence: prodhan82@gmail.com </b></div>

<br />

<p align="center">
  <a href="https://github.com/asadprodhan/A-beginner-s-guide-to-Bioinformatics/tree/main?tab=GPL-3.0-1-ov-file#readme"><img src="https://img.shields.io/badge/License-GPL%203.0-yellow.svg" alt="License GPL 3.0"></a>
  <a href="https://orcid.org/0000-0002-1320-3486"><img src="https://img.shields.io/badge/ORCID-green?style=flat-square&logo=ORCID&logoColor=white" alt="ORCID"></a>
</p>


<br />

## **Introduction**


Welcome to bioinformatics! This hands-on guide introduces you to the world where biology meets computing — exploring DNA, proteins, and genomes using software tools.  

By the end of this session, you will be able to:

- Perform basic sequence searches  
- Extract and organize biological data  
- Share reproducible analyses  

We will start by installing essential software with **Conda** and perform a simple BLAST search to get you started on analyzing real biological data (Fig. 1).


<p align="center">
  <img 
    src="https://github.com/asadprodhan/A-beginner-s-guide-to-Bioinformatics/blob/main/A_beginner%E2%80%99s_guide_to_Bioinformatics.png"
  >
</p>
<p align = "center">
Figure 1. How you start a bioinformatics analysis.
</p>


## **Contents**  

- [What is conda?](https://github.com/asadprodhan/A-beginner-s-guide-to-Bioinformatics/edit/main/README.md#what-is-conda)  


- [Do I have ‘conda’ in my computer?](https://github.com/asadprodhan/A-beginner-s-guide-to-Bioinformatics/edit/main/README.md#do-i-have-conda-in-my-computer)  


- [How to install ‘conda’ in my computer?](https://github.com/asadprodhan/A-beginner-s-guide-to-Bioinformatics/edit/main/README.md#how-to-install-conda-in-my-computer)


- [How do I install softwares using ‘conda’?](https://github.com/asadprodhan/A-beginner-s-guide-to-Bioinformatics/edit/main/README.md#how-do-i-install-softwares-using-conda)


- [How do I find the installation command for a bioinformatics software?](https://github.com/asadprodhan/A-beginner-s-guide-to-Bioinformatics/edit/main/README.md#how-do-i-find-the-installation-command-for-a-bioinformatics-software)


- [Has my software installation been successful?](https://github.com/asadprodhan/A-beginner-s-guide-to-Bioinformatics/edit/main/README.md#has-my-software-installation-been-successful)


- [How do I share my ‘conda environment’?](https://github.com/asadprodhan/A-beginner-s-guide-to-Bioinformatics/edit/main/README.md#how-do-i-share-my-conda-environment)


- [How do I run a ‘conda environment’ that my colleague shared with me?](https://github.com/asadprodhan/A-beginner-s-guide-to-Bioinformatics/edit/main/README.md#how-do-i-run-a-conda-environment-that-my-colleague-shared-with-me)


- [Let’s carry out a simple bioinformatics task](https://github.com/asadprodhan/A-beginner-s-guide-to-Bioinformatics/edit/main/README.md#lets-carry-out-a-simple-bioinformatics-task)


- [How do I extract the sequences of the blast hits in fasta format?](https://github.com/asadprodhan/A-beginner-s-guide-to-Bioinformatics/edit/main/README.md#how-do-i-extract-the-sequences-of-the-blast-hits-in-fasta-format)


<a name="headers"/>

<br />

### **What is conda?**

Simply put, Conda is a tool that quickly installs, updates, and manages softwares. During installation, it will also automatically install any other required softwares (called dependencies). Furthermore, it can create ‘environment’ called conda environment. Conda environment is a directory that compartmentalises the installed softwares, thus does not allow the installed softwares to interfere with any other softwares that are outside this environment and vice versa.


### **Do I have ‘conda’ in my computer?**

Run *conda* in the terminal. If you already have conda, you will see instructions on how to use conda. If you don’t have it installed, you get an error saying ‘conda command not found’


### **How to install ‘conda’ in my computer?**

- visit the anaconda archive to find the latest installer of Anaconda (https://repo.anaconda.com/archive/Anaconda3-2020.02-Linux-x86_64.sh) 


- download using 'wget' 


  ```
  wget https://repo.anaconda.com/archive/Anaconda3-2020.02-Linux-x86_64.sh
  ```
  
  
- make the downloaded installer executable


  ```
  chmod +x Anaconda3-2020.02-Linux-x86_64.sh
  ```
  
  
- run it
  ```
  ./Anaconda3-2020.02-Linux-x86_64.sh
  ```
  
  
- this will install Anaconda3 in the home directory


- add conda to your PATH (How to add? See the commands at https://github.com/asadprodhan/About-the-PATH) 



### **How do I install softwares using ‘conda’?**

- create a conda environment as follows 


```
conda create -n bioinf
```

-n flag is to assign a name to the environment


- activate your conda environment


```
conda activate bioinf
```


- check that your conda environment has been activated. Two ways to check:


  - run the following command and the currently activated environment is indicated by an asterisk 
  
  
    ```
    conda env list
    ```
  
  
  - the currently activated environment is also indicated on the left-hand side of your prompt name (i.e., user name)
  
  
- now install the softwares in your environment as follows


```
conda install -c bioconda samtools
```

> ‘c’ stands for Chanell 


> conda will automatically install all dependent softwares for ‘samtools’ for example


```
conda install -c bioconda trimmomatic
```


```
conda install -c bioconda bedtools
```


```
conda install -c bioconda blast 
```


### **How do I find the installation command for a bioinformatics software?**


Google conda install <name of the software>, you will find something like below:
 
  
### **Has my software installation been successful?**  
  

Run *bedtools* in the terminal to see if ‘bedtools’ has successfully been installed
  
  
### **How do I share my ‘conda environment’?**  
  
  
```
conda env export > environment.yml
```
  
  
Conda environment might be asked by the manuscript reviewers
  

### **How do I run a ‘conda environment’ that my colleague shared with me?**   
  

```
conda env create -f my_environment.yml 
```

The above command will create a conda environment using someone’s shared environment called *my_environment.yml*


## **Let’s carry out a simple bioinformatics task**

  
Say you have a DNA sequence. For example, *AtNRT1.1* (an *Arabidopsis* *thaliana* gene). You want to find out its matches in Rice (*Oryza* *sativa*) and extract the matched sequences. 
  

**How do you do that?**
  
  
You can use the BLASTn tool. BLAST stands for Basic Local Alignment Search Tool. It will find the closest match of *AtNRT1.1* gene in rice. Once the closest matches are identified, their sequences can be extracted from the rice genome.
  
<br />
  
**Let’s run this task.**

<br />

**Step 1: Collect the sequences of the *AtNRT1.1* gene and rice genome from the NCBI database.**
  
  
*AtNRT1.1* gene sequence in fasta format can be collected as follows (Fig. 2):
  
  
<p align="center">
  <img 
    src="https://github.com/asadprodhan/A-beginner-s-guide-to-Bioinformatics/blob/main/AtNRT1.1_download.PNG"
  >
</p>
<p align = "center">
Figure 2. Downloading *AtNRT1.1* gene sequence from the NCBI database.
</p>
<br />

Rice genome sequence in fasta format can be collected as follows:
  
  
```  
sudo wget ftp://ftp.ensemblgenomes.org/pub/plants/release-47/fasta/oryza_sativa/dna/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa.gz
```
  

**Step 2: Turn the rice genome sequence into a searchable blast database**
  
  
Unzip the downloaded file:
  
  
```
gunzip Oryza_sativa.IRGSP-1.0.dna.toplevel.fa.gz
```
  
  
Make it executable by changing the file permission:
  
  
```
chmod +x Oryza_sativa.IRGSP-1.0.dna.toplevel.fa
```
  
  
Make the BLAST database:
  
  
```
makeblastdb -in Oryza_sativa.IRGSP-1.0.dna.toplevel.fa -out blastdb_rice -dbtype 'nucl' -hash_index
```
  

> -in: input fasta file that will be turned into a blast database  
  
  
> -out: name of the newly generated database
  
  
> -dbtype: type of the database is nucleotide (‘nucl’) 
  
  
> -hash_index: creates a hash table that is used to find sequence matches
  

Run the blastn command:
  
  
```
blastn -query AtNRT1.1.fasta -task blastn -db blastdb_rice -outfmt 6 -out blast_hits.tsv -evalue 1e-10 -num_threads 18 
```
  
  
> -query: query sequence i.e., the sequence that you are using to find matches in rice 
  
  
> -task: blastn, searching nucleotide against a nucleotide database 
  
  
> -db: database that you made earlier 
  
  
> -outfmt 6: it will produce the output in a tabular format
  
  
> -out: output file name
  
  
> -evalue: one of the blast hits parameters. Is the blast hit a by-chance match to your query sequence? evalue 1e-10 means that the chance is 1 in 10,000,000,000. Therefore, the lower the evalue is, the better it is
  
  
> -num_threads: number of computer threads to be used.
  
  

Now, the blastn output file has 12 columns (Fig. 3).


<p align="center">
  <img 
    src="https://github.com/asadprodhan/A-beginner-s-guide-to-Bioinformatics/blob/main/blast_hits.PNG"
  >
</p>
<p align = "center">
Figure 3. Blastn_hits.
</p>
  
  
The column headers are as follows:
  
  
> Column 1: qseqid, query sequence id. If you open you query file, you will find the id after ‘>’ sign
  
  
> Column 2: sseqid, subject i.e., reference sequence id
  
  
> Column 3: pident, percentage of identical matches
  
  
> Column 4: length, alignment length 
  
  
> Column 5: mismatch, number of mismatches
  
  
> Column 6: gapopen, number of gap openings
  
  
> Column 7: qstart, start of alignment in query
  
  
> Column 8: qend, end of alignment in query
  
  
> Column 9: sstart, start of alignment in subject
  
  
> Column 10: send, end of alignment in subject
  
  
> Column 11: evalue, expect value
  
  
> Column 12: bitscore, bit score

<br />
  
By default, you will get all 12 columns in your blast output file.
  
  
> You can customise it by replacing the ‘-outfmt 6’ flag by the -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"

  
The blast output file can be sorted in ascending or descending order using its column headers. For example, the following command will sort the blast hits in descending order using column 12 (the bitscore value). The higher the bitscore, the better it is.

<br />
  
```  
sort -n -r -k 12 blast_hits.tsv > blast_hits_sorted.tsv  
```
   
> sort -n means sorting numerically
  
  
> sort -r means reverse order
  
  
> sort -k 12 means sorting on column 12

<br />
  
See the sorted version of the above blastn output file (Fig. 4).


<p align="center">
  <img 
    src="https://github.com/asadprodhan/A-beginner-s-guide-to-Bioinformatics/blob/main/blast_hits_sorted.PNG"
  >
</p>
<p align = "center">
Figure 4. Blastn hits sorted in descending order of the bitscores (Column 12).
</p>

<br />
  
## **How do I extract the sequences of the blast hits in fasta format?** 
  
  
**Run the following commands step-by-step:**
  
  
```
awk '{print $2,$9,$10,"Rice_"$1"_"NR}' blast_hits_sorted.tsv > blast_hits.bed
```
  
  
This command line prints column 2, 9 and 10. It also assigns a name to each hit. Finally, saves the print out in ‘bed’ format, suitable for ‘bedtools’ 
  
  
```  
cat  blast_hits.bed | awk -v OFS='\t' '{ if ($2 < $3) {print $1,$2,$3,$4} else {print $1,$3,$2,$4} }' >  blast_hits_sorted.bed
```
  
  
This command line organises the chromosome positions from start to end
  
  
```  
bedtools getfasta -fi Brara_Chiifu_V3.5.fa -bed blast_hits_sorted.bed –name > blast_hits_seqs.fa
```
  
  
This command line extracts the sequences for each blast hit in fasta format (Fig. 5).

<br />
  
<p align="center">
  <img 
    src="https://github.com/asadprodhan/A-beginner-s-guide-to-Bioinformatics/blob/main/blast_hits_seqs.PNG"
  >
</p>
<p align = "center">
Figure 5. Sequences extracted for each Blastn hits (fasta format).
</p>


<br />
  
**Congratulations! You've just completed your first analysis in Bioinformatics. You've also learned how to install and manage your next set of softwares for your next analysis.**
  
