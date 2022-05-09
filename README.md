# Bloodhound

Source code for the Bloodhound project by Vaibhav Mangipudy, Linghao Kong, and Hunter King for the Spring 2022 Computational Genomics course.

File descriptions:
* bloodhound.py - Runnable program containing main method. Takes a csv file with genemark output names and sequences and the output directory of Jpred mass scheduler as input, prints tables of candidate aca and acr loci and produces figure showing candidate regions along the input range.
* acahmms.py - Contains methods for building the aca HMMS from aca secondary structure sequences,scoring input secondary structure sequences, and identifying candidate aca loci.
* acrhmm.py - Contains methods for training the acr HMM on known acr sequences, scoring input peptide sequences, and identifying candidate acr loci.
* acrproteins.csv - File containing the acr peptide sequences that the acr HMM is trained on
* PaeruginosaSample.csv - File containing the sample input peptide sequences (in genemark output format) from the region in Pseudomonas Aeruginosa genome containing a known aca protein
* PaeruginosaSample.faa_dir_output/ - Directory containing predicted secondary structures of sequences from PaeruginosaSample.csv outputted by Jpred mass scheduler 

Python Dependencies:
* Python 3.8.9
* pandas 1.4.1
* numpy 1.22.3
* sklearn 0.0
* hmmlearn 0.2.7
* matplotl 3.5.1

Instructions to run sample input:
1. Ensure all files and python dependencies are installed
2. Run "python3 bloodhound.py PaeruginosaSample.csv PaeruginosaSample.faa_dir_output"

Expected output from sample input:

![image](https://user-images.githubusercontent.com/87866533/167475918-4e06034e-72cc-4fe4-a665-d7c3db327d12.png)
![image](https://user-images.githubusercontent.com/87866533/167475944-c540e59c-9b93-464d-9e91-3dc860ce628c.png)

Instructions to run on new species:
1. Obtain species genome as a fasta file
2. Input genome into genemark via, *ensure "Protein sequence" checkbox is checked* http://exon.gatech.edu/GeneMark/gmhmmp.cgi
3. Concatenate genemark sequences to single line for jpred (following notebook code may help): https://colab.research.google.com/drive/1Z25WXRsxeRr4p67NTGlRpdXLpu-flI0i#scrollTo=QHwtVL82jRyd
4. Run jpred mass scheduler on formatted fasta file (tutorial here: http://www.compbio.dundee.ac.uk/jpred4/downloads/Jpred_4_Tutorial_on_Scheduled_Mass_Submission.pdf) on formatted peptide sequence fasta file (note: Jpred is extremely slow for many proteins; this step may be sped up by running jnet software on your own machine rather than using the jpred server)
5. Create a csv for the peptide sequences with two columns: 'seq_name' and 'seq_aa' (keep whole genemark output name in 'seq_name')
6. Run "python3 bloodhound.py <peptide_seqs.csv> <jpred_output.faa_dir_output>

Link to Colab notebook used to prototype and generate figures for data analysis: https://colab.research.google.com/drive/1F8mX9EzZoSRYDdkyBJJEhngCshYDkfW9#scrollTo=kD4yPkitxwsF
