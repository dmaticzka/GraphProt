for a in SGD SVR; do for b in CLIP RNACOMPETE; do touch *.fa; make distclean; make test -e EVAL_TYPE=RNACOMPETE -e SVM=SVR; done; done
