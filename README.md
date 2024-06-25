# FluModel_ParamSearch
Parameterizing a three-subtype influenza model for tropical influenza, based on Vietnam

This repository contains the results of a parameter-space search to produce acyclic influenza dynamics as seen in tropical world regions. There are 18 forms of the model, changing immunity waning dynamics, types of subpopulations, and presence of case importation. The document 'Methods_Description.pdf' provides an explanation of how the model was forumalated and parameterized. Within each folder are the model code and the parameterizations that meet a set of criteria to produce desired influenza dynamics. The following naming conventions are used to differentiate the model forms:

### Number of stages in the recovered class
i) RR: 2 stages <br>
ii) RRRR: 4 stages <br>
iii) RRRRRR: 6 stages

### Subpopulations
i) L: one population <br>
ii) 2L: two populations of equal size <br>
iii) Ll: two populations of unequal size (80/20 ratio)

### Importation
i) nI: no importation <br>
ii) I: importation at regular intervals
