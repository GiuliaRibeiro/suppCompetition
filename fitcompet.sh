nohup ssh abacus0001 "cd ~/lotka_volterra_tecamebas; R CMD BATCH compet1.R &" &
nohup ssh abacus0001 "cd ~/lotka_volterra_tecamebas; R CMD BATCH compet2.R &" &
nohup ssh abacus0002 "cd ~/lotka_volterra_tecamebas; R CMD BATCH compet3.R &" &
nohup ssh abacus0002 "cd ~/lotka_volterra_tecamebas; R CMD BATCH compet4.R &" &
nohup ssh abacus0003 "cd ~/lotka_volterra_tecamebas; R CMD BATCH compet5.R &" &
nohup ssh abacus0003 "cd ~/lotka_volterra_tecamebas; R CMD BATCH compet6.R &" &
nohup ssh abacus0005 "cd ~/lotka_volterra_tecamebas; R CMD BATCH compet7.R &" &
