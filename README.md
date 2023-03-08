# Global plate model choice impacts reconstructions of the latitudinal biodiversity gradient

Author(s): [Lewis A. Jones](mailto:LewisA.Jones@outlook.com), [Bethany J. Allen](mailto:Bethany.Allen@bsse.ethz.ch), [William Gearty](willgearty@gmail.com) and [Lucas Buffan](mailto:lucas.buffan@ens-lyon.fr).

This repository contains the data and code required to run the analyses of the article, "Global plate model choice impacts reconstructions of the latitudinal biodiversity gradient" (Jones et al. 2023). 

To cite the paper: 
> Lewis A. Jones, Bethany J. Allen, William Gearty, Lucas Buffan. 2023. Global plate model choice impacts reconstructions of the latitudinal biodiversity gradient.

-------

## Study details

In this study, we evaluate the influence of Global Plate Model choice on reconstructions of latitudinal biodiversity gradients in deep time. Specifically, we test whether different types of gradient (unimodal, flat, bimodal) might emerge based on Global Plate Model choice. This has implications for our understanding of deep time macroecolgical patterns and their drivers. Our study focused on three widely used open-source models, which are available via the [GPlates Web Service](https://gwsdoc.gplates.org/reconstruction-models):

* MERDITH model (Merdith et al., 2021) - GPlates ID = MERDITH2021
* PALEOMAP model (Scotese & Wright, 2018) - GPlates ID = PALEOMAP
* GOLONKA model (Wright et al., 2013) - GPlates ID = GOLONKA

-------
## Repository structure

In this repository, files and code are organised as:

* **Data** files are stored in the `/data/` folder
* **Analysis** code in the `/R/` folder
* **Results** in the `/results/` folder
* **Figures** in the `/figures/` folder

-------

## Workflow

The workflow and documentation for data analysis can be found in: `/R/run_analysis.R`.

The workflow for data visualisation can be found in: `/R/generate_figures.R`.

Documentation and comments relating to the workflow can be found within the aforementioned scripts, as well as the relevant subscripts.
