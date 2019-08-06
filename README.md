# LPI_BLS
Predicting lncRNA-protein interactions with a broad learning system-based stacked ensemble classifier


Dependencies 
(1) numpy
(2) scikit-learn

Prerequisites:
(1) BroadLearningSystem
    We downloaded the code of Broad Learning System from http://www.broadlearning.ai/broad-learning-system-an-effective-and-efficient-incremental-learning-system-without-the-need-for-deep-architecture/
(2) Pse-in-One
    We downloaded the code of Pse-in-One from downloaded from http://bioinformatics.hitsz.edu.cn/Pse-in-One/

Usage:
(1) perform 5 fold cross validation on RPI488 or RPI7317.
$ python LPI_BLS.py -dataset RPI488/RPI7317

(2) predict a new lncRNA-protein pair.
$ python pse.py (fasta file of lncRNA sequence) (lncRNA_pse_feature) RNA PC-PseDNC-General -lamad=6 -w=0.9
$ python pse.py (fasta file of protein sequence) (protein_pse_feature) Protein PC-PseAAC-General -lamad=9 -w=0.5
$ python LPI_BLS.py -pair (lncRNA-protein pair needed to be predicted) -rf (fasta file of lncRNA sequence) -rp (fasta file of protein sequence) -rP (lncRNA_pse_feature) -pP (protein_pse_feature)

