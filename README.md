# underlying-distribution
 
We present the “Heterogeneous Underlying-based” (HUB) method and associated Python codes that model (HUB-forward code) and interpret (HUB-backward code) the results of droplet freezing experiments under the assumptions that each ice nucleating site in the sample has a characteristic nucleation temperature that is time independent. 

-----

HUB-forward predicts $f_{ice} (T)$ and $N_m (T)$ from a proposed distribution of IN temperatures.
The underlying distribution $P_u(T)$ or the differential freezing spectrum $n_m (T)$ is a mixture of normalized probability distribution functions: 

$P_u(T) = c_1 * P_1(T) + c_2 * P_2(T)  + c_3 * P_3(T) +  ... +  c_p * P_p(T)$ 

where $p$ is the **number of sub populations** or classes, $c_1, c_2, c_3, \ ..., c_p$ are the **weights** such that $c_1 + c_2 + c_3 \ ... + \ c_p = 1$, and $P_1(T), P_2(T), P_3(T), \ ..., P_p(T)$ represent normalized Gaussian distributions of each sub population of ice nucleators (INs). 

Each distribution $P_i(T)$ is defined by a set of parameters: 
- the **mode** $T_{mode,i}$ which is the  most commonly observed value in a set of data, i.e. the x-axis location of the highest peak in that sub population; 
- the **scaling factor** $s_{i}$ which controls the spread of the distribution.

-----

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

## Running on Colab

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

After answering all of these questions, the HUB-forward code will generate in the same directory a few files in .pdf and .txt extensions: the underlying and the concentration-dependent distributions, the fraction of frozen droplets $f_{ice} (T)$ and the cumulative freezing spectrum $N_m (T)$.

HUB-backward requires an input file located in the directory called input. The file has to be a text format .txt, the first column is temperature, and the second column is $N_m (T)$ or $f_{ice} (T)$. After answering a series of questions, HUB-backward will generate a few files in .pdf and .txt extensions: a spline fit of the input data, the optimized cumulative freezing spectrum $N_m (T)$ data, and the optimized differential freezing spectrum $n_m (T)$.

# Acknowledgements 

I. de A. R and V. M. gratefully acknowledge support by AFOSR through MURI Award No. FA9550-20-1-0351. K. M. acknowledges support by the National Science Foundation under Grant No. (NSF 2116528).

# License

This code is the active version of the code archived in XXX at (doi:?) and documented in de Almeida Ribeiro, Meister & Molinero (submitted). 

You are welcome to use and distribute this code as you see fit, but it remains the intellectual property of the authors and must be cited appropriately (please cite both the paper and the code). Direct any questions about this code to: Ingrid Ribeiro (ingrid.ribeiro@utah.edu), or create an issue in this Github repository.

-----

MIT License

Copyright (c) 2022 Ingrid Ribeiro

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
