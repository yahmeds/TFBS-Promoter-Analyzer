all: test doc

# Chemins et variables
PYTHON = python3
DOC_DIR = docs
MODULES = pwm.py utils.py scan_pwm.py putative_TFBS

test:
	python3 -m unittest discover -v src

# Nettoyer les fichiers générés
clean:
	rm -rf $(DOC_DIR)/_build
	find . -name "*.pyc" -delete
	find . -name "__pycache__" -delete

.PHONY: all doc test clean
