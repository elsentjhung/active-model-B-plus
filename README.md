## Numerical Implemantation of Active Model B+

This contains the source code used in the numerical simulation of Active Model B+, as described in this paper [Tjhung, Nardini, and Cates, PRX (2018)].

List of files:
* `./source_code/active_model_B_plus.c` = main source file
* `./source_code/system_parameters.h` = this file contains all the parameters in the model
* `./source_code/mathematical_functions.h` = this file contains implementation of numerical derivatives
* `./source_code/animate.sh` = script file to generate a movie (requires gnuplot and fffmpeg)
* `./output/output.mp4` = example of the output

Also check this link [Discretization in Active Model B+] for the numerical discretization used in our codes.

[Discretization in Active Model B+]: https://elsentjhung.github.io/2020/12/26/discretization.html
[Tjhung, Nardini, and Cates, PRX (2018)]: https://journals.aps.org/prx/abstract/10.1103/PhysRevX.8.031080
