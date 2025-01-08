# Global plate model choice impacts reconstructions of the latitudinal biodiversity gradient

Author(s): [Lewis A. Jones](mailto:LewisA.Jones@outlook.com), [William Gearty](willgearty@gmail.com), [Lucas Buffan](mailto:lucas.buffan@ens-lyon.fr), and [Bethany J. Allen](mailto:Bethany.Allen@bsse.ethz.ch).

This repository contains the data and code required to run the analyses of the article, "Global plate model choice impacts reconstructions of the latitudinal biodiversity gradient" (Jones et al. 2025). 

To cite the paper: 
> Lewis A. Jones, William Gearty, Lucas Buffan, and Bethany J. Allen. 2025. Global plate model choice impacts reconstructions of the latitudinal biodiversity gradient.

-------

## Study details

In this study, we evaluate the influence of Global Plate Model choice on reconstructions of latitudinal biodiversity gradients in deep time. Specifically, we test whether different types of gradient (unimodal, flat, bimodal) might emerge based on Global Plate Model choice. This has implications for our understanding of deep time macroecolgical patterns and their drivers. Our study focused on three widely used open-source models, which are available via the [GPlates Web Service](https://gwsdoc.gplates.org/reconstruction-models):

* MERDITH2021 model (Merdith et al., 2021) - GPlates ID = MERDITH2021
* PALEOMAP model (Scotese & Wright, 2018) - GPlates ID = PALEOMAP
* GOLONKA model (Wright et al., 2013) - GPlates ID = GOLONKA
* TorsvikCocks2017 model (Torsvik and Cocks, 2016) - GPlates ID = TorsvikCocks2017

-------
## Repository structure

In this repository, files and code are organised as:

* **Data** files are stored in the `/data/` folder
* **Analysis** code in the `/R/` folder
* **Results** in the `/results/` folder
* **Figures** in the `/figures/` folder
* **Manuscript** files in the `/manuscript/` folder

-------

## Workflow

The workflow and documentation for data analysis can be found in: `/R/run_analysis.R`.

Documentation and comments relating to the workflow can be found within the aforementioned scripts, as well as the relevant subscripts.
