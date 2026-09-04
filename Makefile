# Dark Calculator & Scientific Engine Automation Makefile
.PHONY: all install build start test test-coverage clean

all: install build test

install:
	@echo "Installing project dependencies..."
	npm install --silent

build:
	@echo "Building application assets and validating modules..."
	node build.js

start:
	@echo "Starting Dark Calculator server on port 3000..."
	node src/server.js

test:
	@echo "Running automated verification test suite..."
	node tests/runner.js

test-coverage:
	@echo "Generating test coverage report..."
	node tests/runner.js --coverage

clean:
	@echo "Cleaning transient build artifacts..."
	rm -rf coverage .nyc_output
