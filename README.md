# pythondice2013implementation
This is a port of the 2013 dice model by William Nordhaus from GAMS to python using pyomo.  

This project runs the 2013 version of the DICE model obtaining the same result as if using the GAMS model available from http://www.econ.yale.edu/~nordhaus/homepage/Web-DICE-2013-April.htm.  This model still requires a non-linear solver such as ipopt, which is available at https://projects.coin-or.org/Ipopt or in binary form at http://ampl.com/products/solvers/open-source/.   All the parameters are in the diceParameters.csv file and are named the same way as in the original GAMS model.  All parameters are imported into a pandas dataframe.  The results are saved to a dataframe and exported to a csv file.
