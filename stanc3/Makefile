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
