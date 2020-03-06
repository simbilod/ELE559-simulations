# ELE559-simulations
Workshop on photonic device and circuit simulation using open-source tools. Created for Princeton's Spring 2020 ELE559 : Photonic Systems class.

## Log in to Adroit

https://myadroit.princeton.edu/

## Checkout environment

```
git clone https://github.com/simbilod/ELE559-simulations
```

## Create and enable isolated conda environment

Change directory (cd) into the course folder, then type:

```
# cd ELE559-simulations
conda env create -f environment.yml
conda activate mp
```

## Install JupyterLab Ipykernel

```
conda activate mp
conda install ipykernel
ipython kernel install --user --name=PyMeep
conda deactivate
```

## (For 3D plotting : optional) Make sure extension is enabled in JupyterLab

conda install -c conda-forge nodejs
jupyter labextension install @jupyter-widgets/jupyterlab-manager
jupyter labextension install ipyvolume
jupyter labextension install jupyter-threejs

