# ELE559-simulations
Workshop on photonic device and circuit simulation using open-source tools. Originally created for Princeton's Spring 2020 ELE559 : Photonic Systems class.

## Log in to Adroit

https://myadroit.princeton.edu/
Requesting access (Princeton students) : https://researchcomputing.princeton.edu/access

## Checkout environment

Navigate to the shell

![checkout](images/shell_access.png)

Run this command :

```
git clone https://github.com/simbilod/ELE559-simulations
```
![clone](images/cloning.PNG)

## Start up the Jupyter server

Under "My interactive sessions", click on Jupyter for classes :

![jupyterlab](images/lab_access.png)

Optionally, check the box for Jupyter Lab. Fill in your requested time and cores (be considerate to other users). Then click launch. When the environment is ready, launch it.

![jupyterlab](images/jupyterlab.png)

## Opening a Notebook

You will see an interface like this :

![interface](images/interface.png)

You should see multiple Python environments when you log in (mp, pmp, layout, etc.). Make sure you select the "mp" environments when you boot the Notebooks. You can change what environment you are running on by clicking the text besides the dot on the top right :

![env](images/env.png)

## Using a Notebook

Execute a cell by pressing `shift+enter`. The output will be displayed right underneath. The added text should provide enough context to follow!

More on Jupyter Notebooks : https://www.dataquest.io/blog/jupyter-notebook-tutorial/

## Manually making an environment visible

Jupyter should auto-detect the environments, but if for some reason it does not you can manually add them. 

Go to the shell (see "Checkout environment") and type

```
conda info --envs
```
to list the environments. Activate the environment you want to add by typing
```
conda activate environment_name
```
Then link the environment "environment_name" to your Jupyter kernel with label "my_environment_name" by running
```
python -m ipykernel install --user --name=my_environment_name
```
You should get the following response :
![add_env](images/manual_add_env.PNG)
Restart the Jupyter for Classes server and you should see the environment under the name you assigned.

## Accessing KLayout GUI from the myAdroit Desktop 

To see the GDS files you have designed on the server itself, you can use the preinstalled KLayout on the Desktop app.

First, (optionally) create an alias to open KLayout in your .bashrc. In the shell (see "Checkout environment") type

```
cd
nano .bashrc
```
and add the lines
```
# ELE559 initialization script
source /home/ELE559/ELE559.bashrc
```
Then, open the desktop app like you would a Jupyter for Classes session. Open a terminal and type
```klayout```
(or full link from the alias) to open the program.

![desktop_klayout](images/desktop_klayout.png)

## Acknowledgements

This material is based in part upon work supported by the National Science Foundation under Grant Number E2CDA-1740262. Any opinions, findings, and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Science Foundation.
