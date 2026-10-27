.PHONY: install test reproduce
install:
	python -m pip install -e '.[dev]'
test:
	pytest -q
reproduce:
	python -m kc200gt.analysis
