from scripttest import TestFileEnvironment
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
    call = """../../GraphProt.pl --action train \
              --fasta ../testedenacc.train.positives.fa \
              --negfasta ../testedenacc.train.negatives.fa \
              -prefix test_edenacc_train --keep-tmp"""
    env.run(call)


def test_eden_accuracy2():
    "Now test on train."
    call = """../../GraphProt.pl --action predict \
              --fasta ../testedenacc.train.positives.fa \
              --negfasta ../testedenacc.train.negatives.fa \
              -model ../test_edenacc_train.model \
              -prefix test_edenacc_manualtestontrain --keep-tmp"""
    run = env.run(call)
    assert "test_edenacc_manualtestontrain.predictions" in run.files_created
    cmp("tests/test_edenacc_manualtestontrain.predictions", testdir + "test_edenacc_manualtestontrain.predictions")


def test_classification_cv():
    "Test crossvalidation with classifiation."
    call = """../../GraphProt.pl -mode classification -action cv \
              -fasta ../testclip.train.positives.fa \
              -negfasta ../testclip.train.negatives.fa \
              -abstraction 1 -R 0 -D 0 -epochs 2 -lambda 0.1 -bitsize 10 -prefix CL_cv --keep-tmp"""
    env.run(call)
    cmp("tests/CL_cv.cv_results", testdir + "CL_cv.cv_results")


def test_classification_train():
    "Test model training with classification."
    call = """../../GraphProt.pl -mode classification -action train \
              -fasta ../testclip.train.positives.fa \
              -negfasta ../testclip.train.negatives.fa \
              -abstraction 1 -R 0 -D 0 -epochs 2 -lambda 0.1 -bitsize 10 \
              -prefix CL_train --keep-tmp"""
    env.run(call)
    cmp("tests/CL_train.model", testdir + "CL_train.model")


def test_classification_predict():
    "Test predicion with classification."
    call = """../../GraphProt.pl -mode classification -action predict \
              -fasta ../testclip.train.positives.fa \
              -negfasta ../testclip.train.negatives.fa \
              -model ../CL_train.model -abstraction 1 -R 0 -D 0 \
              -epochs 2 -lambda 0.1 -bitsize 10 -prefix CL_predict --keep-tmp"""
    env.run(call)
    cmp("tests/CL_predict.predictions", testdir + "CL_predict.predictions")


def test_classification_predict_only_positives():
    "Test prediction using only positives."
    call = """../../GraphProt.pl -mode classification -action predict \
              -fasta ../testclip.train.positives.fa \
              -model ../CL_train.model \
              -abstraction 1 -R 0 -D 0 -epochs 2 -lambda 0.1 -bitsize 10 \
              -prefix CL_predict_onlypos --keep-tmp"""
    env.run(call)
    cmp("tests/CL_predict_onlypos.predictions", testdir + "CL_predict_onlypos.predictions")


def test_classification_predict_ntmargins():
    "Test prediction of nucleotide-wise margins."
    call = """../../GraphProt.pl -mode classification -action predict_profile \
              --onlyseq \
              -fasta ../testclip.train.positives.fa \
              -model ../CL_train.model \
              -abstraction 1 -R 0 -D 0 -epochs 2 -lambda 0.1 -bitsize 10 \
              -prefix CL_ntmargins --keep-tmp"""
    env.run(call)
    cmp("tests/CL_ntmargins.profile", testdir + "CL_ntmargins.profile")


def test_classification_predict_has():
    "Test prediction of high affinity sites."
    call = """../../GraphProt.pl -mode classification -action predict_has --onlyseq \
              -fasta ../testclip.train.positives.fa \
              -model ../CL_train.model \
              -abstraction 1 -R 0 -D 0 -epochs 2 -lambda 0.1 -bitsize 10 -percentile 55 \
              -prefix CL_has --keep-tmp"""
    env.run(call)
    cmp("tests/CL_has.has", testdir + "CL_has.has")


def test_classification_motif():
    "Test motif calculation."
    call = """../../GraphProt.pl -mode classification -action motif \
              -fasta ../testclip.train.positives.fa \
              -model ../CL_train.model \
              -abstraction 1 -R 0 -D 0 -epochs 2 -lambda 0.1 -bitsize 10 \
              -prefix CL_motif --keep-tmp"""
    run = env.run(call)
    assert "CL_motif.structure_motif.png" in run.files_created
    assert "CL_motif.structure_motif_pairedunpaired.png" in run.files_created
    assert "CL_motif.sequence_motif.png" in run.files_created
    cmp("tests/CL_motif.sequence_motif", testdir + "CL_motif.sequence_motif")
    cmp("tests/CL_motif.structure_motif", testdir + "CL_motif.structure_motif")


def test_classification_motif_onlyseq():
    "Test motif calculation."
    call = """../../GraphProt.pl -mode classification -action motif \
              --onlyseq \
              -fasta ../testclip.train.positives.fa \
              -model ../CL_train.model \
              -R 1 -D 0 -bitsize 10 -prefix CL_onlyseq --keep-tmp"""
    run = env.run(call)
    assert "CL_onlyseq.sequence_motif.png" in run.files_created
    cmp("tests/CL_onlyseq.sequence_motif", testdir + "CL_onlyseq.sequence_motif")


def test_classification_profile_all_vp():
    "Calculate profiles for pure viewpoint sequences."
    call = """../../GraphProt.pl -mode classification -action predict_profile \
              -fasta ../testclip.twoseqs300nt_allcaps.fa \
              -model ../CL_train.model \
              -abstraction 1 -R 0 -D 0 -epochs 2 -lambda 0.1 -bitsize 10 \
              -prefix testclip.twoseqs300nt_allcaps.margins --keep-tmp"""
    env.run(call)
    cmp("tests/testclip.twoseqs300nt_allcaps.margins.profile",
        testdir + "testclip.twoseqs300nt_allcaps.margins.profile")


def test_classification_profile_center_vp():
    "Calculate profiles for center viewpoint seqeunces."
    call = """../../GraphProt.pl -mode classification -action predict_profile \
              -fasta ../testclip.twoseqs300nt.fa \
              -model ../CL_train.model \
              -abstraction 1 -R 0 -D 0 -epochs 2 -lambda 0.1 -bitsize 10 \
              -prefix testclip.twoseqs300nt.margins --keep-tmp"""
    env.run(call)
    cmp("tests/testclip.twoseqs300nt.margins.profile",
        testdir + "testclip.twoseqs300nt.margins.profile")


def test_classification_motif_center_vp():
    "Calculate motif for center viewpoint jsequences."
    call = "../../GraphProt.pl -mode classification -action motif \
            -fasta ../testclip.twoseqs300nt.fa \
            -model ../CL_train.model \
            -abstraction 1 -R 0 -D 0 -epochs 2 -lambda 0.1 -bitsize 10 \
            -prefix testclip.twoseqs300nt.motif --keep-tmp"""
    run = env.run(call)
    assert "testclip.twoseqs300nt.motif.structure_motif.png" in run.files_created
    assert "testclip.twoseqs300nt.motif.structure_motif_pairedunpaired.png" in run.files_created
    assert "testclip.twoseqs300nt.motif.sequence_motif.png" in run.files_created
    cmp("tests/testclip.twoseqs300nt.motif.sequence_motif", testdir + "testclip.twoseqs300nt.motif.sequence_motif")
    cmp("tests/testclip.twoseqs300nt.motif.structure_motif", testdir + "testclip.twoseqs300nt.motif.structure_motif")
