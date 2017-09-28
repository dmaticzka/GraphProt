from scripttest import TestFileEnvironment
# import re
# from filecmp import cmp

testdir = "tests/testenv_graphprot_classification"

env = TestFileEnvironment(testdir)


def test_invocation_no_params():
    "Call without parameters should return usage information."
    call = "GraphProt.pl"
    run = env.run(
        call,
        expect_error=True)
    assert run.returncode == 2
