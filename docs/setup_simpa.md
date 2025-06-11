# SIMPA Installation Guide

This document describes how to install the [SIMPA](https://github.com/IMSY-DKFZ/simpa) package and external dependencies used for optical and acoustic forward modelling.

## Install SIMPA from GitHub

```
pip install git+https://github.com/IMSY-DKFZ/simpa.git@main
```

Check that the installation was successful:

```
python -c "import simpa; print('SIMPA version:', simpa.__version__)"
```

## Install MCX

Download the latest nightly build from [mcx.space](http://mcx.space/nightly/github/) and extract it. In this repository we installed it under `/tmp/mcx`.

```
wget https://mcx.space/nightly/github/mcx-linux-x64-github-latest.zip -O mcx.zip
unzip mcx.zip -d /tmp/mcx
```

## Install k-Wave

Download the latest k-Wave release and extract it. If MATLAB is available you only need the MATLAB toolbox. In this repository the toolbox is stored under `/tmp/kWave`.

```
# Example (download may require manual steps)
mkdir /tmp/kWave
# extract files here
```

## Configure `path_config.env`

Create a file `path_config.env` and specify the paths to the forward model binaries and the location where simulation output should be stored:

```
SIMPA_SAVE_DIRECTORY=/workspace/Faster_Convex_Cone/data
MCX_BINARY_PATH=/tmp/mcx/mcx/bin/mcx
MATLAB_BINARY_PATH=/usr/local/MATLAB/bin/matlab
```

## Using the PathManager

SIMPA reads these variables via the `PathManager` class:

```python
from simpa.utils import PathManager
pm = PathManager(environment_path="./path_config.env")
print(pm.get_mcx_binary_path())
print(pm.get_matlab_binary_path())
```

