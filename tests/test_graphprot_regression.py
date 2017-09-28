from scripttest import TestFileEnvironment
from filecmp import cmp

testdir = "tests/testenv_graphprot_regression/"

env = TestFileEnvironment(testdir)


def test_regression_cv():
    "Crossvalidation with regression."
    call = """../../GraphProt.pl -mode regression -action cv \
              -fasta ../test_data_full_A.train.fa
              -affinities ../test_data_full_A.train.affys -prefix REG_cv \
              -abstraction 1 -R 0 -D 0 -epsilon 0.11 -c 11 -bitsize 10 --keep-tmp"""
    env.run(call)
    assert cmp("tests/REG_cv.cv_results", testdir + "REG_cv.cv_results")


def test_regression_train():
    "Train model with regression."
    call = """../../GraphProt.pl -mode regression -action train \
              -fasta ../test_data_full_A.train.fa \
              -affinities ../test_data_full_A.train.affys -prefix REG_train \
              -abstraction 1 -R 0 -D 0 -epsilon 0.11 -c 11 -bitsize 10 --keep-tmp"""
    env.run(call)
    assert cmp("tests/REG_train.model", testdir + "REG_train.model")


def test_regression_predict():
    "Predict from model with regression."
    call = """../../GraphProt.pl -mode regression -action predict \
              -fasta ../test_data_full_A.train.fa \
              -affinities ../test_data_full_A.train.affys \
              -model ../REG_train.model -prefix REG_predict \
              -abstraction 1 -R 0 -D 0 -epsilon 0.11 -c 11 -bitsize 10 --keep-tmp"""
    env.run(call)
    assert cmp("tests/REG_predict.predictions", testdir + "REG_predict.predictions")


def test_regression_predict_noaffy():
    "Predict from model with regression without affinity input."
    call = """../../GraphProt.pl -mode regression -action predict \
            -fasta ../test_data_full_A.train.fa \
            -model ../REG_train.model \
            -prefix REG_predict_noaffys \
            -abstraction 1 -R 0 -D 0 -epsilon 0.11 -c 11 -bitsize 10 --keep-tmp"""
    env.run(call)
    assert cmp("tests/REG_predict_noaffys.predictions", testdir + "REG_predict_noaffys.predictions")
