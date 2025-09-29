/******************************************************************************* 
    Applied Economics — Synthetic Control & Power Simulations
    Robles Urquiza · Burgermeister · Huerta ·  Vitaliani — UdeSA (2025)
*******************************************************************************/

clear all
set more off

*--- Paths (usa rutas relativas del repo)
global ROOT "."
global INPUT "$ROOT/input"
global OUTPUT "$ROOT/output"
global TEMP "$ROOT/temp"

cap mkdir "$OUTPUT"
cap mkdir "$TEMP"

********************************************************************************
*                                S Y N T H                                     *
*    São Paulo: synthetic control, robustness (placebo, LOO, permutations)
********************************************************************************

* Datos panel (df.csv debe estar en /input)
import delimited "$INPUT/df.csv", clear
tsset code year

*----------------------------- Figure 1 ---------------------------------------
* São Paulo vs. promedio Brasil (sin SP)
quietly {
    forval y = 1990/2009 {
        egen avg`y' = mean(homiciderates) if year==`y' & code!=35
    }
}
egen Brazil_average = rowmin(avg1990-avg2009)
set scheme s1color

twoway (tsline homiciderates if code==35, lcolor(black)) ///
       (tsline Brazil_average, lcolor(gs8) lp(dash)), ///
       xscale(range(1990 2009)) xlabel(1990 1995 2000 2005) ///
       ytitle("Homicide Rate") xtitle("Year") ///
       legend(order(1 "São Paulo" 2 "Brazil (avg.)") pos(7) cols(1) ring(0) size(*0.8)) ///
       xline(1998, lcolor(black) lp(dot)) ///
       note("Figure 1. Source: IPEA.")
graph export "$OUTPUT/Fig1.png", replace

*------------------------ Prep + synth baseline --------------------------------
rename proportionextremepoverty propexpov
save "$INPUT/df.dta", replace

cap which synth
if _rc ssc install synth, replace

synth homiciderates stategdpcapita stategdpgrowthpercent populationprojectionln ///
      yearsschoolingimp homiciderates(1990(1)1998) propexpov(1990(1)1998) ///
      giniimp(1990(1)1998), trunit(35) trperiod(1999) nested

matrix gaps        = e(Y_treated) - e(Y_synthetic)
matrix Y_treated   = e(Y_treated)
matrix Y_synth     = e(Y_synthetic)
matrix RMSPE       = e(RMSPE)
svmat gaps
svmat Y_treated
svmat Y_synth

*----------------------------- Figure 2 ---------------------------------------
twoway (line Y_treated year, lcolor(black)) ///
       (line Y_synth year,   lcolor(black) lp(dash)), ///
       xscale(range(1990 2009)) xlabel(1990 1995 2000 2005) ///
       yscale(range(10 60)) ///
       legend(order(1 "São Paulo" 2 "Synthetic São Paulo") pos(7) cols(1) ring(0) size(*0.7)) ///
       xline(1999, lcolor(black) lp(dot)) ///
       xtitle("Year") ytitle("Homicide Rates") ///
       note("Figure 2: São Paulo vs. Synthetic São Paulo")
graph export "$OUTPUT/Fig2.png", replace

*----------------------------- Figure 3 ---------------------------------------
twoway (line gaps year, lcolor(black)), ///
       xscale(range(1990 2009)) xlabel(1990 1995 2000 2005) ///
       yscale(range(-30 39)) ylabel(-30 -20 -10 0 10 20 30) ///
       xline(1999, lcolor(black) lp(dot)) yline(0, lcolor(black) lp(dot)) ///
       xtitle("Year") ytitle("Gap in Homicide Rates") ///
       text(27.5 1994.5 "Policy change →", place(e)) ///
       note("Figure 3: Gap = Treated − Synthetic")
graph export "$OUTPUT/Fig3.png", replace

*----------------------------- Figure 4 (placebo timing 1994) -----------------
preserve
synth homiciderates stategdpcapita stategdpgrowthpercent populationprojectionln ///
      yearsschoolingimp homiciderates(1990(1)1994) propexpov(1990(1)1994) ///
      giniimp(1990(1)1994), trunit(35) trperiod(1994) xperiod(1990(1)1994) nested

matrix Yt94 = e(Y_treated)
matrix Ys94 = e(Y_synthetic)
svmat Yt94
svmat Ys94

twoway (line Yt94 year, lcolor(black)) ///
       (line Ys94 year, lcolor(black) lp(dash)) if year<=1998, ///
       xscale(range(1990 1998)) xlabel(1990 1992 1994 1996 1998) ///
       yscale(range(0 50)) ylabel(0 10 20 30 40 50) ///
       legend(order(1 "São Paulo" 2 "Synthetic São Paulo") pos(7) cols(1) ring(0) size(*0.7)) ///
       xline(1995, lcolor(black) lp(dot)) ///
       xtitle("Year") ytitle("Homicide Rates") ///
       text(40 1992.1 "Placebo policy →", place(e)) ///
       note("Figure 4: Placebo implementation in 1994")
graph export "$OUTPUT/Fig4.png", replace
restore

*----------------------------- Figure 5 (leave-one-out) -----------------------
use "$INPUT/df.dta", clear
tsset code year

* Resultados base (SP tratado)
quietly synth homiciderates stategdpcapita stategdpgrowthpercent populationprojectionln ///
               yearsschoolingimp homiciderates(1990(1)1998) propexpov(1990(1)1998) ///
               giniimp(1990(1)1998), trunit(35) trperiod(1999) nested ///
               xperiod(1990(1)1998) keep("$OUTPUT/loo-resout35", replace)

* LOO: sacar cada estado del donor pool
forvalues j = 11/53 {
    if `j'==35 continue
    use "$INPUT/df.dta", clear
    tsset code year
    drop if code==`j'
    quietly synth homiciderates stategdpcapita stategdpgrowthpercent populationprojectionln ///
                   yearsschoolingimp homiciderates(1990(1)1998) propexpov(1990(1)1998) ///
                   giniimp(1990(1)1998), trunit(35) trperiod(1999) nested ///
                   xperiod(1990(1)1998) keep("$OUTPUT/loo-resout`j'", replace)
}

* Armar panel con múltiples sintéticos
use "$OUTPUT/loo-resout11.dta", clear
forvalues i = 12/53 {
    cap confirm file "$OUTPUT/loo-resout`i'.dta"
    if _rc==0 merge 1:1 _Co_Number _time using "$OUTPUT/loo-resout`i'.dta", nogen
}

* Renombrar sintéticos y gaps (ejemplo para varios ids + SP)
foreach i in 14 32 42 50 53 35 {
    cap ren _Y_synthetic _Y_synthetic_`i'
    cap ren _Y_treated   _Y_treated_`i'
    cap gen _Y_gap_`i' = _Y_treated_`i' - _Y_synthetic_`i'
}

twoway (line _Y_synthetic_14 _time, lcolor(gs14)) ///
       (line _Y_synthetic_32 _time, lcolor(gs14)) ///
       (line _Y_synthetic_42 _time, lcolor(gs14)) ///
       (line _Y_synthetic_50 _time, lcolor(gs14)) ///
       (line _Y_synthetic_53 _time, lcolor(gs14)) ///
       (line _Y_treated_35   _time, lcolor(black) lwidth(medthick)) ///
       (line _Y_synthetic_35 _time, lcolor(black) lp(dash)), ///
       xscale(range(1990 1998)) xlabel(1990 1995 2000 2005) ///
       yscale(range(0 50)) ylabel(0 10 20 30 40 50) ///
       legend(order(6 "São Paulo" 7 "Synthetic São Paulo" 3 "Synthetic SP (leave-one-out)") ///
              pos(7) cols(1) size(*0.7) ring(0)) ///
       xline(1999, lcolor(black) lp(dot)) ///
       xtitle("Year") ytitle("Homicide Rates") ///
       text(45 1994 "Policy change →", place(e)) ///
       note("Figure 5: Leave-one-out synthetic distribution")
graph export "$OUTPUT/Fig5.png", replace

*----------------------------- Figure 6 (permutation test: gaps) --------------
use "$INPUT/df.dta", clear
tsset code year

tempname RMS
forvalues i = 11/17 21/29 31/33 41/43 50/53 {
    quietly synth homiciderates stategdpcapita stategdpgrowthpercent populationprojectionln ///
                   yearsschoolingimp homiciderates(1990(1)1998) propexpov(1990(1)1998) ///
                   giniimp(1990(1)1998), trunit(`i') trperiod(1999) nested
    matrix gaps`i' = e(Y_treated) - e(Y_synthetic)
    svmat gaps`i'
    matrix RMSPE`i' = e(RMSPE)
    svmat RMSPE`i'
}

quietly synth homiciderates stategdpcapita stategdpgrowthpercent populationprojectionln ///
               yearsschoolingimp homiciderates(1990(1)1998) propexpov(1990(1)1998) ///
               giniimp(1990(1)1998), trunit(35) trperiod(1999) nested
matrix gaps35 = e(Y_treated) - e(Y_synthetic)
svmat gaps35
matrix RMSPE35 = e(RMSPE)
svmat RMSPE35

twoway ///
(line gaps111 year, lcolor(gs14)) (line gaps121 year, lcolor(gs14)) (line gaps131 year, lcolor(gs14)) ///
(line gaps141 year, lcolor(gs14)) (line gaps151 year, lcolor(gs14)) (line gaps161 year, lcolor(gs14)) (line gaps171 year, lcolor(gs14)) ///
(line gaps211 year, lcolor(gs14)) (line gaps221 year, lcolor(gs14)) (line gaps231 year, lcolor(gs14)) (line gaps241 year, lcolor(gs14)) (line gaps251 year, lcolor(gs14)) (line gaps261 year, lcolor(gs14)) (line gaps271 year, lcolor(gs14)) (line gaps281 year, lcolor(gs14)) (line gaps291 year, lcolor(gs14)) ///
(line gaps311 year, lcolor(gs14)) (line gaps321 year, lcolor(gs14)) (line gaps331 year, lcolor(gs14)) ///
(line gaps411 year, lcolor(gs14)) (line gaps421 year, lcolor(gs14)) (line gaps431 year, lcolor(gs14)) ///
(line gaps501 year, lcolor(gs14)) (line gaps511 year, lcolor(gs14)) (line gaps521 year, lcolor(gs14)) (line gaps531 year, lcolor(gs14)) ///
(line gaps351 year, lcolor(black)), ///
xscale(range(1990 2009)) xlabel(1990 1995 2000 2005) ///
yscale(range(-30 30)) ylabel(-30 -20 -10 0 10 20 30) ///
xline(1998, lcolor(black) lp(dash)) yline(0, lcolor(black)) ///
legend(order(28 "São Paulo" 1 "Control states") pos(7) cols(1) size(*0.7) ring(0)) ///
xtitle("Year") ytitle("Gap in Homicide Rates") ///
text(25 1993.3 "Policy change →", place(e)) ///
note("Figure 6: Permutation test (gaps)")
graph export "$OUTPUT/Fig6.png", replace

*----------------------------- Figure 7 (subset by MSPE threshold) ------------
gen threshold = 2*1.258143   // umbral usado en el TP

twoway (line gaps131 year, lcolor(gs14)) (line gaps151 year, lcolor(gs14)) (line gaps171 year, lcolor(gs14)) ///
       (line gaps211 year, lcolor(gs14)) (line gaps231 year, lcolor(gs14)) (line gaps241 year, lcolor(gs14)) ///
       (line gaps251 year, lcolor(gs14)) (line gaps311 year, lcolor(gs14)) (line gaps411 year, lcolor(gs14)) ///
       (line gaps421 year, lcolor(gs14)) (line gaps431 year, lcolor(gs14)) (line gaps531 year, lcolor(gs14)) ///
       (line gaps351 year, lcolor(black)), ///
       xscale(range(1990 2009)) xlabel(1990 1995 2000 2005) ///
       yscale(range(-30 30)) ylabel(-30 -20 -10 0 10 20 30) ///
       xline(1998, lcolor(black) lp(dash)) yline(0, lcolor(black)) ///
       legend(order(13 "São Paulo" 2 "Controls (MSPE < 2× SP)") pos(7) cols(1) ring(0) size(*0.7)) ///
       xtitle("Year") ytitle("Gap in Homicide Rates") ///
       note("Figure 7: Permutation test (selected controls)")
graph export "$OUTPUT/Fig7.png", replace


********************************************************************************
*                          P O W E R   S I M S                                 *
*   Monte Carlo: tamaños de muestra × tamaños de efecto × varianza × p(T=1)
********************************************************************************

* Helper: rutina de simulación (evita repetir bloques)
program drop _all
program define _power_sim
    * args: filename sigma treatshare controlflag
    args basefile sigma ptreated controlflag

    clear
    set seed 123
    set obs 17000
    gen ganancias_estimadas = rnormal(10000,2000)
    drop if ganancias_estimadas<0
    gen impuestos_pagados = 0.2*ganancias_estimadas + rnormal(0,`sigma')
    drop if impuestos_pagados<0
    save "$INPUT/`basefile'.dta", replace

    local efectos = 5
    local sizes   = 10
    clear
    set obs `=`efectos'*`sizes''
    gen sample_size = 1000
    forval i=2/10 {
        replace sample_size = 1000*`i' if _n<=`efectos'*`i' & sample_size==.
    }
    bysort sample_size: gen efecto = .
    foreach k in 0.01 0.025 0.05 0.075 0.1 {
        bysort sample_size: replace efecto = `k' if missing(efecto)
        bysort sample_size: replace efecto = . if _n>1 & missing(efecto)==0
    }
    * Reasignar efectos (1..5 por bloque)
    by sample_size: replace efecto = 0.01 if _n==1
    by sample_size: replace efecto = 0.025 if _n==2
    by sample_size: replace efecto = 0.05 if _n==3
    by sample_size: replace efecto = 0.075 if _n==4
    by sample_size: replace efecto = 0.1 if _n==5

    tempfile results
    save "`results'", replace

    * Loop por celdas: 500 reps por celda
    preserve
    use "`results'", clear
    gen st_power = .
    quietly {
        forval j=1/_N {
            su sample_size if _n==`j'
            local n = r(mean)
            su efecto if _n==`j'
            local eff = r(mean)

            tempname R
            matrix `R' = J(500,1,.)
            forvalues x=1/500 {
                preserve
                use "$INPUT/`basefile'.dta", clear
                sample `n', count
                gen T = (runiform()<`ptreated')
                replace impuestos_pagados = impuestos_pagados*(1+`eff') if T==1
                if `controlflag'==1 {
                    reg impuestos_pagados T ganancias_estimadas, robust
                }
                else {
                    reg impuestos_pagados T, robust
                }
                matrix `R'[`x',1] = _b[T]/_se[T]
                restore
            }
            clear
            svmat `R'
            gen reject = (R1>1.65)
            quietly su reject
            scalar pow = r(mean)
            restore
            replace st_power = pow in `j'
        }
    }
    save "$OUTPUT/power.dta", replace
end

*--- Escenarios (2.1 a 2.4)
* 2.1 Baseline: sigma=500, p(T)=0.5, sin controles
_power_sim base_sim_21 500 0.5 0

* 2.2 Más ruido: sigma=5000, p(T)=0.5
_power_sim base_sim_22 5000 0.5 0

* 2.3 p(T)=0.2 y 0.8 (mismo sigma alto)
_power_sim base_sim_23A 5000 0.2 0
_power_sim base_sim_23B 5000 0.8 0

* 2.4 Menos ruido + control por ganancias
_power_sim base_sim_24 500 0.5 1

*--- Gráfico de potencia (toma el último power.dta generado)
use "$OUTPUT/power.dta", clear
sort efecto sample_size
replace st_power = round(st_power,.01)
separate st_power, by(efecto)

set scheme s1color
graph set window fontface "Arial"
twoway (connected st_power1 sample_size, lcolor(black) mcolor(black)) ///
       (connected st_power2 sample_size, lcolor(blue)  mcolor(blue))  ///
       (connected st_power3 sample_size, lcolor(red)   mcolor(red))   ///
       (connected st_power4 sample_size, lcolor(green) mcolor(green)) ///
       (connected st_power5 sample_size, lcolor(purple) mcolor(purple)), ///
       xtitle("Number of observations") ytitle("Power") ///
       legend(rows(1) title("Effect size") ring(0) pos(6) ///
              label(1 "1%") label(2 "2.5%") label(3 "5%") label(4 "7.5%") label(5 "10%")) ///
       plotregion(color(white)) graphregion(color(white))
graph export "$OUTPUT/power_graph.png", replace width(2000) height(1200)

display as text "Done. Figures in /output."
