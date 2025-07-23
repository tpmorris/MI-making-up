# Code for simulation study on “Is multiple imputation making up information?” blog post

Simulation code to accompany my [blog post on whether multiple imputation invents information](https://open.substack.com/pub/tpmorris/p/is-multiple-imputation-making-up). The description of the simulation study is in that blog.

A note on the code: I simulated all the full data in one go and saved it as `sim_data.dta`. For each rep, I then used the relevant chunk of data. Initially I tried to do `mi impute … , by(rep_id)` on `sim_data.dta`, and was pleased with that idea, but Stata hung when it was trying to append all the imputed datasets – shame! So the full dataset is saved in the repo but the imputed datasets are not.

The community-contributed `[simsum]` command is needed to run this code [available here](https://github.com/UCL/simsum).

Caution: do not use the code in this do-file as a model of a good simulation study in Stata. I did it offhand for a blog post and would be more structured and careful for a paper. See [here](github.com/tpmorris/TheRightWay) for how to do it better.
