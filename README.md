# Code for simulation study on “Is multiple imputation making up information?” blog post

Simulation code to accompany my [blog post on whether multiple imputation invents information](https://open.substack.com/pub/tpmorris/p/is-multiple-imputation-making-up). The description of the simulation study is in that blog.

Note: code by me is now under "Stata" because @Elessenne contributed R code under "R". He used 10,000 repetitions, which was wise: although I said I was comfortable with the Monte Carlo error after 800 reps, this was WRT the reassurance about MI, and I subsequently undermined myself by saying that I was not sure if differences were just Monte Carlo error.

A note on the Stata code: I simulated all the full data in one go and saved it as `sim_data.dta`. For each rep, I then used the relevant chunk of data. Initially I tried to do `mi impute … , by(rep_id)` on `sim_data.dta`, and was pleased with that idea, but Stata hung when it was trying to append all the imputed datasets – shame! So the full dataset is saved in the repo but the imputed datasets are not.

The community-contributed `{simsum}` Stata command is needed to run this code [available here](https://github.com/UCL/simsum).

Caution: do not use the Stata code in this do-file as a model of a good simulation study. I did it offhand for a blog post and would be more structured and careful for a paper. See [here](github.com/tpmorris/TheRightWay) for how to do it better.
