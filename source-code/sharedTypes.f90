module sharedTypes
!
! DESCRIPTION
!   specify global variables
!
! REVISIONS
!  06/07/09 wittenburg
!
  use kindnumber
  implicit none

!!! parameters that are read in via parameter file 'input.par'

  character(len=100) :: snpfile,& ! genotype coding without header
       &datafile,& ! datafile with header
       &weightFile,& !file with weights per SNP
       &corrfile='phenotype_corr.txt',& ! copy of datafile modified in iteration steps
       &outfile,& ! output file with locus number and estimate of genetic effect
       &resfile,& ! output file with residuals=ycorr
       &nuisanceFile,& ! output file for nuisance effects estimated via BLUE
       &logFile !output file with information about SNP with allele frequency<maf

  character(len=50) :: solveMethod,& ! solving strategy, either
       !  'fast'Bayes or 'svs'
       &NuisanceMethod ! use ASReml ('asreml') or Gibbs sampler
       ! ('gibbs') for solving nuisance effects or 'blup' for BLUP/BLUE
 
  integer(ik4) :: posobs,& ! column of phenotypic value in datafile
       &skiplines,& !number of lines to skip in datafile
       &maxit  ! number of maximum iteration in 'ice' or Gibbs sampling rounds
       
  integer(ik2) :: reparamMethod,& ! method to reparametrize genotypes
       &stepnum ! number of steps for interval of credibility

  logical :: solveNuisance,& ! additional nuisance effects have to be solved
       &solveDom,& ! solve for dominance effects
       &solveEpi,& ! solve for epistatic effects
       &writeFreq,& ! write estimated allele frequencies
       &gammaFix,& ! don't update gamma iteratively
       &lambdaFix ! don't update lambda iteratively

  real(rk8) :: maf,& !lower limit for Minor Allele Frequency
       &nonZeroMain,& !proportion of assumed main effects (a-priori)
       &nonZeroEpi,& !proportion of assumed epistatic effects
       &errorlevel !type-I error for genetic effects

end module sharedTypes
