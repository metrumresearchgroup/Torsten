all:
	dune build src/stanc/stanc.exe

TEST_DUNES := $(foreach d,$(shell find test/integration -type d),$(d)/dune)
TEST_TORSTEN := $(foreach d,$(shell find test/integration/good/torsten -type d),$(d)/dune)
test: $(TEST_DUNES)

.PHONY: doc test

test:
	dune runtest

test_torsten: $(TEST_TORSTEN)
	dune runtest

testcoverage:
	@find . -name '*.coverage' | xargs rm -f
	dune clean
	BISECT_FILE=`pwd`/bisect dune runtest --instrument-with bisect_ppx --force
	bisect-ppx-report html --expect src/ --do-not-expect src/stancjs/
	bisect-ppx-report summary --expect src/ --do-not-expect src/stancjs/
	@rm *.coverage

format:
	dune build @fmt

cross:
	dune build src/stanc/stanc.exe -x windows

static:
	dune build src/stanc/stanc.exe --profile static

clean:
	dune clean

doc:
	dune build @doc

re: clean all
