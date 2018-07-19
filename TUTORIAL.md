# GraphProt Tutorial

## Installation

We recommend the installation of GraphProt via Bioconda. Once Bioconda is set up (see [https://bioconda.github.io/](https://bioconda.github.io/)), GraphProt including all it's dependencies can be installed easily.

```bash
conda install graphprot
```

In general, it is recommended to create a dedicated virtual environment for GrapProt. This ensures that upgrades or the installation of additional tools do not interfere with GraphProt or its dependencies.

```bash
conda create -n graphprot graphprot
```

After successfull creation of the environment, it can be activated via

```bash
conda activate graphprot
```

## Basic Usage

This section gives an example on

* how to create a GraphProt classification model based on tow small sets of bound and unbound sequences derived from a CLIP-Seq experiment and
* how to create a motif representation of the learned sequence-and-structure preferences.

### Preliminaries

GraphProt classification models are trained using sets of bound and unbound sequences in `fasta` format. Of importance here are the different meanings that GraphProt assigns to nucleotides represented by **uppercase** and **lowercase** IUPAC codes. Nucleotides set in **uppercase** letters will be used to  determine sequence and structure preferences and usually correspond to the bound regions as defined by CLIP-Seq peaks. Nucleotides set in **lowercase** letters are used for secondary structure prediction. In addition, these sequence parts will be incorporated into the model if they are located near the bound regions, measured as distance along the RNA backbone and bonds within the predicted secondary structures interactions.

### 1. Get training data

For the purpose of this tutorial, download a sample of 1000 positive and negative sequences from a TAF15 PAR-CLIP experiment (Tollervey, J. R. et al. Characterizing the RNA targets and position-dependent splicing regulation by TDP-43. Nat. Neurosci. 14, 452–458 (2011).). Please note that length of the sequences around the binding sites was reduced to allow quick calculation. This will reduce the quality of secondary struture predictions!

```
wget bioinf.uni-freiburg.de/~maticzkd/GraphProtTutorialData.tar.bz2 && \
tar xf GraphProtTutorialData.tar.bz2
cd GraphProtTutorialData
```

### Train a GraphProt model

These files can be directly used to train a sequence-and-structure GraphProt model.

```bash
GraphProt.pl \
--action train \
--fasta    TDP43_iCLIP_subset.train.positives.fa \
--negfasta TDP43_iCLIP_subset.train.negatives.fa \
--prefix   TDP43_iCLIP_subset_default
```

This will run GraphProt with fast default parameters and write a GraphProt model to file `TDP43_iCLIP_subset_default.model`.

### Motif detection

Now we use the model to extract a motif from the bound sites. Since the training data for this tutorial contains 1,000 sequences, we will lower the number of high-scoring sequences the motif is extracted from to 200 using option `--motif_top_n 200` an. The following code sets the motif length to 12 using option `--motif-len 12`, which is a good fit for the TDP43 motif.

```bash
GraphProt.pl \
--action motif \
--model  TDP43_iCLIP_subset_default.model \
--fasta  TDP43_iCLIP_subset.train.positives.fa \
--prefix TDP43_iCLIP_subset_default \
--motif_top_n 200 \
--motif_len 12
```

This code will produce three motif files:

* `TDP43_iCLIP_subset_default.sequence_motif.png`: GraphProt sequence motif showing the TDP43 consonsus motif "UGUGUGUGUGUG"
* `TDP43_iCLIP_subset_default.structure_motif.png`: GraphProt structure motif.
* `TDP43_iCLIP_subset_default.structure_motif_pairedunpaired.png`: GraphProt structure motif in simplified 'paired/unpaired' format.

## Extended Usage

### Parameter Optimization

GraphProt can be run with optimized hyperparameters. This will increase accuracy as well ass time and memory consumption. To determine optimized parameters, run

```bash
GraphProt.pl \
--action ls \
--fasta    TDP43_iCLIP_subset.ls.positives.fa \
--negfasta TDP43_iCLIP_subset.ls.negatives.fa \
--prefix   TDP43_iCLIP_subset_optim
```

This will create a file, ``, containing optimized parameters. These parameters can then be used for model training and application via option `--params`.

```bash
GraphProt.pl \
--action train \
--fasta    TDP43_iCLIP_subset.ls.positives.fa \
--negfasta TDP43_iCLIP_subset.ls.negatives.fa \
--prefix   TDP43_iCLIP_subset_optim \
--params   TDP43_iCLIP_subset_optim.params
```

Note that you will have to specify the appropriate optimized parameters every time the model is used, e.g. to determine a motif.

```bash
GraphProt.pl \
--action motif \
--model  TDP43_iCLIP_subset_optim.model \
--params TDP43_iCLIP_subset_optim.params \
--fasta  TDP43_iCLIP_subset.train.positives.fa \
--prefix TDP43_iCLIP_subset_optim \
--motif_top_n 200 \
--motif_len 12
```
