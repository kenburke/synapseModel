# Synapse Model

[![Build
Status](https://travis-ci.org/kenburke/synapseModel.svg?branch=master)](https://travis-ci.org/kenburke/synapseModel)

Reduced model of synaptic transmission, with continuous integration testing.

## purpose

1. Develop reduced model of synaptic transmission
2. Model presynaptic calcium dynamics and stochastic vesicular release
3. Fit model parameters to experimental data
4. Modulate calcium parameters to predict effect on synaptic transmission
5. Compare results to observed biological effects to predict biophysical mechanisms

In English, I'm trying to model a lot of neuronal synapses onto a single neuron. This model can be used to change specific biophysical properties of the synapses and simulate the effect on the post-synaptic neuron's electrical currents.


## structure


```
.
├── README.md
│   ...
├── input
│   ├── ...
│   └── <possible input files>
├── output
│   ├── ...
│   └── <output files>
├── synapse
│   ├── __init__.py
│   ├── __main__.py
│   ├── model.py
│   ├── io.py
│   └── utils.py
└── test
    ├── test_model.py
    ├── test_io.py
    └── test_utils.py
    
```

## usage

To use the package, first create a conda environment

```
conda env create
```
to install all dependencies outlined in `environment.yml` by default. **

Then activate the env

```
source activate synapseModel
```

Then the package's main function (located in `synapse/__main__.py`) 
can be run as follows

```
python -m synapse [-L| -M] <output file>
```

where ``-L`` and ``-M`` allow you to load an input file from the ``input/`` folder or manually enter model parameters, respectively.
The output will be saved into ``output/<output file>``

**
If this command results in a permission error, try

```
sudo conda env create
```

followed by your password, depending on whether you have permissions or own this folder.


## testing

Testing is as simple as running

```
python -m pytest
```

from the root directory of this project.

## interactive

This project comes with an interactive Jupyter notebook in ``synapseModel.ipynb``.


## contributors

Original implementation and theoretical design by Ken Burke.