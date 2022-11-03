
Table of contents
=================

<!--ts-->
   * [How to cite](#how-to-cite)
   * [General idea](#general-idea)
   * [How to Install the Project](#how-to-install-the-project)
      * [Running in a local installation](#running-in-a-local-installation)
      * [Running on Google drive](#running-on-google-drive)
   * [How to Run the Project](#how-to-run-the-project)
   * [Output data](#output-data)
   * [Examples used as an input](#examples-used-as-an-input)
   * [Acknowledgements](#acknowledgements)
   * [License](#license)
<!--te-->

# How to cite

This code is the active version of the code archived here and documented in de Almeida Ribeiro, Meister & Molinero, "Optimization of the underlying distribution of heterogeneous ice nucleation temperatures from drop-freezing experiments", (submitted). 

You are welcome to use and distribute this code as you see fit, but it remains the intellectual property of the authors and must be cited appropriately (please cite the paper). Direct any questions about this code to: Ingrid de A. Ribeiro (ingrid.ribeiro@utah.edu), or create an issue in this Github repository.

# General idea

We present the “Heterogeneous Underlying-based” (HUB) method and associated Python codes that model (HUB-forward code) and interpret (HUB-backward code) the results of droplet freezing experiments under the assumptions that each ice nucleating site in the sample has a characteristic nucleation temperature that is time independent. 

HUB-forward predicts $f_{ice} (T)$ and $N_m (T)$ from a proposed distribution of IN temperatures.
The underlying distribution $P_u(T)$ or the differential freezing spectrum $n_m (T)$ is a mixture of normalized probability distribution functions: 

$P_u(T) = c_1 * P_1(T) + c_2 * P_2(T)  + c_3 * P_3(T) +  ... +  c_p * P_p(T)$ 

where $p$ is the **number of sub populations** or classes, $c_1, c_2, c_3, \ ..., c_p$ are the **weights** such that $c_1 + c_2 + c_3 \ ... + \ c_p = 1$, and $P_1(T), P_2(T), P_3(T), \ ..., P_p(T)$ represent normalized Gaussian distributions of each sub population of ice nucleators (INs). 

Each distribution $P_i(T)$ is defined by a set of parameters: 
- the **mode** $T_{mode,i}$ which is the  most commonly observed value in a set of data, i.e. the x-axis location of the highest peak in that sub population; 
- the **scaling factor** $s_{i}$ which controls the spread of the distribution.

HUB-backward uses a non-linear optimization method to find the differential freezing spectrum $n_m (T)$ that best represents the experimental target cumulative freezing spectrum $N_m (T)$ or fraction of frozen droplets $f_{ice} (T)$ in the experiments. 

-----
# How to Install the Project
## Running in a local installation

Launch with:
```
python HUB-forward.py
```
or
```
python HUB-backward.py
```

## Running on Google drive

Colaboratory, or [Colab](https://colab.research.google.com/?utm_source=scs-index) for short, is a product from Google Research. It allows anybody to write and execute arbitrary Python code through the browser (Firefox, Safari, Chrome etc).

To use it, follow the steps:

- download the directory HUB-forward (or HUB-backward), and then upload it in a **Google drive** directory;

- inside the directory, click on "add more apps",

![image](https://user-images.githubusercontent.com/60315074/199114124-7ffb328d-dd1d-44f4-9d04-246bba6f7538.png)

- type *colaboratory* and install it

<img width="969" alt="image" src="https://user-images.githubusercontent.com/60315074/184701819-8e4baaf3-f2b7-47d2-b067-bb3947151fba.png">

- open the Python notebook HUB-forward.ipynb (or HUB-backward.ipynb)

<img width="848" alt="image" src="https://user-images.githubusercontent.com/60315074/199114807-10b687f2-4221-444b-a48f-0b299eb307fb.png">

# How to Run the Project

- run it by clicking on the play icon and allow the notebook to access your Google Drive files. While the code is executing, a series of questions will appear,

<img width="613" alt="image" src="https://user-images.githubusercontent.com/60315074/199116848-1bd67052-2f33-456f-8005-598beb18a1c0.png">

# Output data

After answering all of these questions, the HUB-forward code will generate in the same directory a few files in .pdf and .txt extensions: the underlying and the concentration-dependent distributions, the fraction of frozen droplets $f_{ice} (T)$ and the cumulative freezing spectrum $N_m (T)$.

After answering a series of questions, HUB-backward will generate a few files in .pdf and .txt extensions: a spline fit of the input data, the optimized cumulative freezing spectrum $N_m (T)$ data, and the optimized differential freezing spectrum $n_m (T)$.

# Examples used as an input

HUB-backward requires an input file located in the directory called input that has to be located in the same directory as the code. The file has to be a text format .txt, the first column is temperature, and the second column is $N_m (T)$ or $f_{ice} (T)$. 

Examples are located inside the input directory:\
Nm_bacteria.txt (Lukas et al., 2022; Schwidetzky et al., 2021)\
Nm_fusarium_kunert2019_strain_3-68.txt (Kunert et al., 2019)\
Nm_pollen_thesis.txt (Dreischmeier, 2019)\
Nm_pH_6p2.txt, Nm_pH_5p6.txt and Nm_pH_4p4.txt (Lukas et al., 2020)\
fractionofice_cholesterol_fig10.txt and fractionofice_cholesterol_fig13.txt (Zhang and Maeda, 2022)

# Acknowledgements 

I. de A. R and V. M. gratefully acknowledge support by AFOSR through MURI Award No. FA9550-20-1-0351. K. M. acknowledges support by the National Science Foundation under Grant No. (NSF 2116528).

# License

The details on the license are shown in [LICENSE](https://github.com/Molinero-Group/underlying-distribution/blob/main/LICENSE).
