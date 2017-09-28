from scripttest import TestFileEnvironment
# import re
from filecmp import cmp

testdir = "tests/testenv_graphprot_classification/"

env = TestFileEnvironment(testdir)


def test_invocation_no_params():
    "Call without parameters should return usage information."
    call = "../../GraphProt.pl"
    run = env.run(
        call,
        expect_error=True)
    assert run.returncode == 2


def test_eden_accuracy1():
    "Train a model."
    call = "../../GraphProt.pl --action train --fasta ../testedenacc.train.positives.fa --negfasta ../testedenacc.train.negatives.fa -prefix test_edenacc_train --keep-tmp"
    env.run(call)


def test_eden_accuracy2():
    "Now test on train."
    call = "../../GraphProt.pl --action predict --fasta ../testedenacc.train.positives.fa --negfasta ../testedenacc.train.negatives.fa -model test_edenacc_train.model -prefix test_edenacc_manualtestontrain --keep-tmp"
    run = env.run(call)
    assert "test_edenacc_manualtestontrain.predictions" in run.files_created
    cmp("tests/test_edenacc_manualtestontrain.predictions", testdir + "test_edenacc_manualtestontrain.predictions")
