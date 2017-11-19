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
│   └── <input files>
├── session
│   ├── session1_plots/
│   │   ├── ...
│   │   └── ...
│   ├── <session1.pkl>
│   ├── session2_plots/
│   │   ├── ...
│   │   └── ...
│   ├── <session2.pkl>
│   └── ...
├── synapse
│   ├── __init__.py
│   ├── __main__.py
│   ├── i_o.py
│   ├── utils.py
│   ├── model.py
│   ├── check.py
│   └── math_funcs.py
└── test
    ├── test_i_o.py
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
can be run from the `synapseModel/` folder as follows

```
python -m synapse [-L| -M] <session_name>
```

where ``-L`` and ``-M`` allow you to load an input file from the ``input/`` folder or manually enter model parameters, respectively.

The program will then prompt you to specify what parameter you would like to modulate over repeated runs, given these baseline parameters.

The output will be saved into ``session/<session_name>.pkl``, along with standard plots in ``session/<session_name>_plots/<plot_names>.svg``


**
If this command results in a permission error, try

```
sudo conda env create
```

followed by your password, depending on whether you have permissions or own this folder.


## advanced usage

If you would like to work more with a specified simulation, including varying more than one parameter at a time, simply start python then import the ``i_o`` and ``utils`` modules:

```
python
from synapse import i_o, utils
```

You can then initialize a new Simulation as follows:

```
SIM = utils.Simulation()
```
with the following optional arguments:
```
name=None                   -- requests user text input by default
params=None                 -- pass parameter dictionary directly (NOT RECOMMENDED, used with i_o.load_input_pickle)
params_from_file=False      -- pass pickle filename in input/, otherwise default.pkl
params_from_user=False      -- requests user text input for each model parameter
```

You can load an existing session with
```
SIM = i_o.load_session_pickle(<filename>)
```
where `<filename>` is a .pkl from a previous run of the script (or `SIM.save()` call)

If you print the Simulation, you can see the current state of the simulation:

```
>>> print(SIM)

NAME : test

PARAMS :
num_cav		    :	1
ca_coop		    :	3.72
cav_p_open	    :	0.83
ca_decay	    :	0.05
num_syn		    :	100
ca_ec50		    :	0.7
quantal_size        :   -10
depletion_on        :	False
vesicle_prox        :	0.25
stim1_time	    :	0.1
a_tau		    :	(0.001, 0.005)
num_trials	    :	300
pool_tau	    :	1.0
sweep_length        :   1
num_stim	    :	2
fs              :   20000.0
num_cav_ratio   :   1
cav_i		    :	1
stim_int	    :	0.05

Runs stored in DEFAULT_RUNS = 1

Runs stored in MOD_RUNS = 0

```

When working with a Simulation object, you have the following methods available to you:

```
save                        -- saves Simulation as pickle
replicate                   -- given a 'simulation_run' object, rerun the model with those parameters
run_default                 -- run the model with the current stored params, save output as simulation_run object appended to self.default_runs
run_modulation              -- run the model many times, multiplying the specified parameter's default value by the value in mod_range each time. Saves output in tuple of ( [list of sim_runs] , {mod_param_name : [mod_range]} ) appended to self.mod_runs
run_analysis                -- analyze the data in given list of simulation_runs (defaults to self.default_runs) and saves plots
plot_epsc_trace             -- for given simulation_run, plot EPSC traces (default is to average them, but can plot all trials with average=False)
plot_hill_func              -- plot where the calcium influx falls on the hill function for a certain synapse and trial in a simulation_run (default is trace=0,synapse=0)
plot_I_ca_trace             -- plot presynaptic calcium influx vs time for a synapse and trial (default is synapse=0, trial=0, average=False. first two can be tuples for more than one synapse or trace)
_runModel                   -- run the model for a given param set (def = self.params), return simulation_run object (def text_display=False can turn on to see progress in terminal)
```


## interactive

This project comes with an interactive Jupyter notebook in ``synapseModel.ipynb`` to explore the advanced applications.


## testing

Testing is as simple as running

```
python -m pytest
```

from the root directory of this project.



## contributors

Original implementation and theoretical design by Ken Burke, Nov 2017.
