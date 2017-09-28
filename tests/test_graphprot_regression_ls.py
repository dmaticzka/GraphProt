from scripttest import TestFileEnvironment
from filecmp import cmp

testdir = "tests/testenv_graphprot_regression_ls/"

env = TestFileEnvironment(testdir)


def test_regression_cv():
    "Crossvalidation with regression."
    call = """../../GraphProt.pl -mode regression -action cv \
              -fasta ../test_data_full_A.train.fa
              -affinities ../test_data_full_A.train.affys -prefix REG_cv \
              -abstraction 1 -R 0 -D 0 -epsilon 0.11 -c 11 -bitsize 10 --keep-tmp"""
    env.run(call)
    cmp("tests/REG_cv.cv_results", testdir + "REG_cv.cv_results")


def test_regression_train():
    "Train model with regression."
    call = """../../GraphProt.pl -mode regression -action train \
              -fasta ../test_data_full_A.train.fa \
              -affinities ../test_data_full_A.train.affys -prefix REG_train \
              -abstraction 1 -R 0 -D 0 -epsilon 0.11 -c 11 -bitsize 10 --keep-tmp"""
    env.run(call)
    cmp("tests/REG_train.model", testdir + "REG_cv.cv_results")


def test_regression_predict():
    "Predict from model with regression."
    call = """../../GraphProt.pl -mode regression -action predict \
              -fasta ../test_data_full_A.train.fa \
              -affinities ../test_data_full_A.train.affys \
              -model ../REG_train.model -prefix REG_predict \
              -abstraction 1 -R 0 -D 0 -epsilon 0.11 -c 11 -bitsize 10 --keep-tmp"""
    env.run(call)
    cmp("tests/REG_predict.predictions", testdir + "REG_predict.predictions")


def test_regression_predict_noaffy():
    "Predict from model with regression without affinity input."
    call = """../../GraphProt.pl -mode regression -action predict \
            -fasta ../test_data_full_A.train.fa \
            -model ../REG_train.model \
            -prefix REG_predict_noaffys \
            -abstraction 1 -R 0 -D 0 -epsilon 0.11 -c 11 -bitsize 10 --keep-tmp"""
    env.run(call)
    cmp("tests/REG_predict_noaffys.predictions", testdir + "REG_predict_noaffys.predictions")


def test_regression_ls():
    "Optimize regression parameters using linesearch."
    call = """../../GraphProt.pl -mode regression -action ls \
              -fasta ../test_data_full_A.train.fa
              -affinities ../test_data_full_A.train.affys \
              -prefix REG_ls -abstraction 1 -R 0 -D 0 -epsilon 0.11 -c 11 -bitsize 10 --keep-tmp"""
    env.run(call)
    cmp("tests/REG_ls.params", testdir + "REG_ls.params")


def test_regression_train_from_ls():
    "Train a regression model using parameters from file."
    call = """../../GraphProt.pl -mode regression -action train \
              -fasta ../test_data_full_A.train.fa \
              -affinities ../test_data_full_A.train.affys
              -params ../REG_ls.params -prefix REG_train_from_ls --keep-tmp"""
    env.run(call)
    cmp("tests/REG_train_from_ls.model", testdir + "REG_train_from_ls.model")
