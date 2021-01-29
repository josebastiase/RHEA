<p align="center">
 <img src="images/RHEA1.png" width="1200" height="300">
</p>

## A verified numerical simulator for hydro-geomechanical heterogeneity based on [MOOSE](https://mooseframework.inl.gov/)


## About

RHEA (Real HEterogeneity App), is an open-source fully coupled finite element application capable of including element-resolution hydro-geomechanical properties in coupled simulations. The code was developed by [José Bastías](https://ingeo.agw.kit.edu/21_99.php) and [Andy Wilkins](https://research.csiro.au/mgt/andy-wilkins/) in collaboration between Karlsruhe Institute of Technology (KIT) and Commonwealth Scientific and Industrial Research Organisatio (CSIRO).

RHEA is a MOOSE based application, for more information visit [MOOSE Framework](https://mooseframework.inl.gov/). To use RHEA, you will need to have [installed MOOSE](https://mooseframework.inl.gov/getting_started/installation/index.html), which can take some time.

RHEA is released under the GNU Lesser General Public License Version 2.1.  See the [LICENSE](LICENSE) file for details.

## Getting started

For system requirements and MOOSE installation visit [Getting Started](https://mooseframework.inl.gov/getting_started/installation/index.html) page of the MOOSE framework.

### Clone RHEA

The stable version of RHEA can be cloned directly for the GitHub repository. As usual in any MOOSE based app, RHEA will be located by default in the `cd ~/projects` foulder generated during the installation of MOOSE.

```
cd ~/projects
git clone https://github.com/josebastiase/RHEA.git
cd ~/projects/RHEA
git checkout master

```

### Compile and test RHEA

As any MOOSE based app, you first need to compile RHEA 

```
cd ~/projects/RHEA
make -j4

```

and then test RHEA

```
cd ~/projects/RHEA
./run_tests

```
If RHEA has compiled successfully, you should see various output, ending with the line

```
2 passed, 0 skipped, 0 pending, 0 failed
```

## Examples

### Terzaghi consolidation

Consolidation of a soil sample due to an external load was studied by Terzaghi.  JOSE TO FILL IN DETAILS.

The RHEA files for this scenario are found in `test/tests/terzaghi/`.  There are three important files:

- `test/tests/terzaghi/Workflow_TerzaghiImportData.ipynb`.  This is a [Jupyter notebook](https://jupyter.org/) that creates files that define the hydraulic conductivity, porosity, bulk modulus and shear modulus throughout the Terzaghi soil sample.  In this case, these properties are homogeneous.  The files created are `K.data`, `p.data`, `L.data` and `G.data` (which are part of this repository, so you don't need to create them yourself).
- `test/tests/terzaghi/TerzaghiImportData.i`.  This is the RHEA input file.  Run it using the `rhea-opt` executable you created during compilation: `rhea-opt -i TerzaghiImportData.i`.
- `test/tests/terzaghi/plot_results.py`.  This is a python file that plots the results, demonstrating agreement between RHEA and the analytical formulae derived by Terzaghi:

![Image](test/tests/terzaghi/terzaghi_p.png)

PERHAPS JOSE WANTS TO REPLACE THIS FIG WITH THE NICE ONE IN THE PAPER.  IN THIS CASE, YOU SHOULD MODIFY plot_results.py accordingly, i think

### Consolidation of a heterogeneous sample

The RHEA files for this scenario are found in `test/tests/terzaghi_layers/`.

JOSE


