# LCsim
Lung Cancer Risk Factor simulator.

Lung cancer causes significant mortality in the US population, with an estimated 156,000 deaths predicted for 2017. More people die from lung cancer per year than from colon, breast, and prostate cancers combined. Current national guidelines for lung screening are based on smoking and age; however, recent studies suggest that incorporating various other risk factors for lung cancer (e.g.: family history, COPD, etc) into selection criteria for screening can improve rates of lung cancer detection. These studies were based on specific trial populations (such as the National Lung Screening Trial) that may not represent the general population though.

The Cancer Intervention and Surveillance Modeling Network (CISNET) developed microsimulation models for identifying optimal lung screening strategies in the general poulation using the Smoking History Generator (SHG) algorithm. Unfortunately, there is currently no available population-level data of other risk factors for lung cancer - besides smoking history - with which to evaluate different screening strategies.

The LCRFsim package was developed to generate these data of other risk factors for lung cancer. LCRFsim simulates a set of correlated risk factors for lung cancer in the general U.S. population, exploiting (i) smoking history data obtained from SHG, (ii) correlation structures among risk factors as observed in the data available, and (iii) incorporating the prevalence information of risk factors in a given U.S. birth cohort.

To start off, I'd recommend checking out example.R! That code is a walkthrough that will run the SHG, process the SHG outputs, feed them into the lung cancer risk factor generator, and plot some diagnostic plots.
