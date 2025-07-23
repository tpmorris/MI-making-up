version 18
clear all

* Control panel: number of reps, full-data sample size, number of imputations
local run     0
local reps  800
local n     500
local nimps 100

set seed 1

* Generate n x reps dataset with id and rep_id
if `n'>`reps' set obs `n' // Note: half is the complete data (n rows); the other half is the missing
else set obs `reps'
gen int rep_id = _n
gen int id = _n
fillin rep_id id
    drop _fillin
    drop if rep_id > `reps'
    drop if id > `n'

* Generate complete data for (all reps)
matrix c = (1, .5 \ .5, 1)
drawnorm y x , corr(c)
gen byte z = rbinomial(1,.5)

save sim_data, replace


frame create estimates int(rep_id) str4(method) int(nimps) str3(target) float(est se df rvi fmi)

quietly {
if `run' {
set coeftabresults off
noi _dots 0, title("Simulation running...")
forval i = 1/`reps' {
    noi _dots `i' 0
    
    use sim_data if rep_id==`i', clear
    expand 2 , gen(miss)
    replace y = . if miss
    replace x = . if miss
    sort rep_id miss id

    * Full data analyses
    regress y x
        frame post estimates (`i') ("Full") (0) ("x") (_b[x]) (_se[x]) (`e(df_r)') (0) (0)
    regress y x z
        frame post estimates (`i') ("Full") (0) ("z") (_b[z]) (_se[z]) (`e(df_r)') (0) (0)
    
    * Multiple imputation
    mi set wide
    mi register imputed y x
    mi impute monotone (regress) y x = i.z, add(`nimps')

    * Use first 5 imputations
    mi estimate, post vartable nimp(5): regress y x
        frame post estimates (`i') ("MI") (5) ("x") (_b[x]) (_se[x]) (e(df_mi)[1,1]) (e(rvi_mi)[1,1]) (e(fmi_mi)[1,1])
    mi estimate, post vartable nimp(5): regress y x z
        frame post estimates (`i') ("MI") (5) ("z") (_b[z]) (_se[z]) (e(df_mi)[1,2]) (e(rvi_mi)[1,2]) (e(fmi_mi)[1,2])
    * Use all 100 imputations
    mi estimate, post vartable: regress y x
        frame post estimates (`i') ("MI") (`nimps') ("x") (_b[x]) (_se[x]) (e(df_mi)[1,1]) (e(rvi_mi)[1,1]) (e(fmi_mi)[1,1])
    mi estimate, post vartable: regress y x z
        frame post estimates (`i') ("MI") (`nimps') ("z") (_b[z]) (_se[z]) (e(df_mi)[1,2]) (e(rvi_mi)[1,2]) (e(fmi_mi)[1,2])
}
set coeftabresults on
}
}

frame change estimates
label data "Estimates dataset from simulation study to see whether multiple imputation is making up information"
gen float   true = .5 if target == "x" // In regression of y|x: true beta_1 = rho(sigma_y / sigma_x) = .5(1/1)
    replace true =  0 if target == "z"

save estimates.dta, replace

***
* Analysis: estimates -> performance
use estimates , clear
simsum est , id(rep_id) method(nimps) by(target) se(se) df(df) true(true) relprec mcse saving(perf_relprec, replace)
simsum fmi , id(rep_id) method(nimps) by(target) true(0) bias mcse saving(perf_fmi, replace)

* Reshape, recode, clean and label estimates dataset for figure
use perf_fmi, clear
rename fmi* est*
    append using perf_relprec
gen str9 estimand = "β{sub:1}" if target=="x"
    replace estimand = "β{sub:2}" if target=="z"
rename est*_mcse est_mcse*
reshape long est est_mcse , i(perfmeasnum target) j(method)
recode method (5=1) (100=2)
    lab def method 0 "Full data" 1 "5 imputations" 2 "100 imputations" 3 "{bf:Estimand:} {&beta}_Z" 4"5 imputations" 5 "100 imputations" 6 "{bf:Estimand:} {&beta}_X", modify
    lab val method method
recode method (1=4) (2=5) if target=="x"
gen float ll = est - est_mcse
gen float ul = est + est_mcse
rename est estimate
gen byte pmeas = perfmeasnum==3
    lab def pmeas 0 "Relative precision" 1 "FMI" , modify
    lab val pmeas pmeas

* Figure
#delimit ;
twoway 
    (rspike ll ul method if method!=0, lc(maroon) lw(*9) hor)
    (scatter method estimate if method!=0, msym(|) mc(white) msize(medlarge))
    ,
    by(pmeas , note("Note: white vertical lines are performance point estimates;" "maroon bars are {&plusmn} 1×Monte Carlo SE") iscale(*1.8) legend(off) xrescale)
    ytitle("") ylab(1 2 3 4 5 6 , nogrid tl(0) valuelabels) ysca(range(.5 2.5) noline) xsca(noline)
    ytic(0.5 3.4 6.5, grid glcolor(gs10) glwidth(medium) glpattern(solid) tl(0))
    xlab(, tl(0) nogextend)
    xtitle("") ysize(4)
;
graph export performance.png, replace
graph export performance.pdf, replace
