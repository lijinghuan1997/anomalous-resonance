# anomalous-resonance
calculate single particle trajectory

  In the program of fortran code, we can change the single particle energy/ pitch angle/ phase angle in the spacecraft rest frame (FAC coordinate), which is next adjusted into plasma rest frame. 
  
  Subroutines Bfield and Efield give the wave fields. 
  
  Subroutine outf outputs the particle information and saves it into the FSD.txt file.

  The plot_gyrophase.m code gives an example to load data and plot the gyro-phase distribution, which can also generate one file in the '.mat' format to save the data.
  
  One example is provided in the repository.
  
  The initial particle distribution function is in the plot_gyrophase.m file, including the temperature and number density parameters.  
