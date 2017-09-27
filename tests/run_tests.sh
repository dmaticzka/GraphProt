#export PATH=/usr/local/perl/bin/:$PATH
#export PATH=/home/maticzkd/src/libsvm-3.12/:$PATH
#export PATH=/usr/local/R/2.15.1-lx/bin:$PATH

bash graphprot_test_misc.sh 2>&1 | grep --line-buffered -v "^Elapsed time" | sed -e 's/tmp-\S\{6\}/tmp-XXXXX/g' #&> graphprot_test_misc.log

bash graphprot_test_classification_misc.sh 2>&1 | grep --line-buffered -v "^Elapsed time" | sed -e 's/tmp-\S\{6\}/tmp-XXXXX/g' #&> graphprot_test_classification_misc.log

bash graphprot_test_regression_misc.sh 2>&1 | grep --line-buffered -v "^Elapsed time" | sed -e 's/tmp-\S\{6\}/tmp-XXXXX/g' #&> graphprot_test_regression_misc.log

bash graphprot_test_regression_ls.sh 2>&1 | grep --line-buffered -v "^Elapsed time" | sed -e 's/tmp-\S\{6\}/tmp-XXXXX/g' #&> graphprot_test_regression_ls.log

bash graphprot_test_classification_ls.sh 2>&1 | grep --line-buffered -v "^Elapsed time" | sed -e 's/tmp-\S\{6\}/tmp-XXXXX/g' #&> graphprot_test_classification_ls.log

gunzip -f *.gz
