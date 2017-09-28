from scripttest import TestFileEnvironment
from filecmp import cmp

testdir = "tests/testenv_graphprot_classification_ls/"

env = TestFileEnvironment(testdir)


def test_classification_ls():
    "Run parameter linesearch."
    call = """../../GraphProt.pl -mode classification -action ls --onlyseq \
              -fasta ../testclip.train.positives.fa \
              -negfasta ../testclip.train.negatives.fa -prefix CL_ls --keep-tmp"""
    env.run(call)
    assert cmp("tests/CL_ls.params", testdir + "CL_ls.params")


def test_classification_train_from_lsparam():
    "Train a model using parameter file."
    call = """../../GraphProt.pl -mode classification -action train --onlyseq\
              -fasta ../testclip.train.positives.fa \
              -negfasta ../testclip.train.negatives.fa
              -params ../CL_ls.params \
              -prefix CL_train_from_ls --keep-tmp"""
    env.run(call)
    assert cmp("tests/CL_train_from_ls.model", testdir + "CL_train_from_ls.model")
