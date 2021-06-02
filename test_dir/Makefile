.PHONY:test
test:test-ntm test-ctsm test-pgs
	@echo "Congratulations. Tests passed"

.PHONY:data-ntm
data-ntm:
	@../prophylaxisEstimator.py ntm --inputfile ../PK_default.csv --timesteps 10 10 1
	@./help.py ntm.csv

.PHONY:test-ntm
test-ntm:data-ntm
	diff ntm.csv ntm_test.csv

.PHONY:data-ctsm
data-ctsm:
	@../prophylaxisEstimator.py ctsm --inputfile ../PK_default.csv
	@./help.py ctsm.csv

.PHONY:test-ctsm
test-ctsm:data-ctsm
	diff ctsm.csv ctsm_test.csv

.PHONY:data-pgs
data-pgs:
	@../prophylaxisEstimator.py pgs --inputfile ../PK_default.csv
	@./help.py pgs.csv
.PHONY:test-pgs
test-pgs:data-pgs
	diff pgs.csv pgs_test.csv

.PHONY:clean
clean:
	${RM} ntm.csv ctsm.csv pgs.csv
