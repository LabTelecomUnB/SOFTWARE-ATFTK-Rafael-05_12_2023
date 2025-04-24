all: build

build:
	python setup.py build_ext --build-temp ./__build --build-lib ./__build --inplace

dist:
	#python -m build --outdir ./__dist
	python setup.py sdist bdist_wheel

to_pypi:
	# python -m twine upload --repository pypi dist/*
	python -m twine upload --repository pypi dist/aftk-0.0.0.dev2.tar.gz

install_requirements:
	pip install --upgrade -r requirements.txt

clean:
	rm -rf build dist *.egg-info wheelhouse __build

clean_pip:
	pip uninstall -y -r <(pip freeze)

.PHONY: all build clean dist