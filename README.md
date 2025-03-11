# virtual-machine-cellGMM
This repository is related to the paper _Cellwise outlier detection in heterogeneous populations_ by Giorgia Zaccaria, Luis A. García-Escudero, Francesca Greselin, and Agustín Mayo-Íscar. It contains the codes for reproducing some of the analyses presented in the paper.

Specifically:
-  "Code for Reproducibility Simulation and Figure Technometrics_EntireScen1" (available in two formats: .Rmd and .html) $\rightarrow$ This script provides the complete simulation study for Scenario 1 with cellwise contamination levels of $0\%$, $5\%$, and $10\%$. It is used to reproduce Figures 1 and 2 of the Main Article for Scenario 1 and runs in approximately 1 hour.

-  "Code for Reproducibility Simulation Technometrics_Scen1_5out_0mis.R" $\rightarrow$ This is a sub-example of the simulation study for Scenario 1 with $5\%$ of cellwise contamination. At the end of the script, there is a code for reproducing the part of Table 1 in the Main Article corresponding to Scenario 1 with $5\%$ outlying values.

- "Code for Reproducibility Simulation Figure Technometrics_Fast.R" $\rightarrow$ This script generates Figures 1 and 2 of the Main Article for all three scenarios. The code starts by loading the file "Data for Figure Reproducibility.RData", which contains the results for the three scenarios presented in the simulation study of the Main Article.

We have also prepared a **Github Codespac**e that can be used to automatically set-up an Rstudio Server and explore the codes in an interactive way. Instructions to access the Github Codespace are reported below.

> [!WARNING]
> Please note that you need a Github account to access the Github Codespace.

1. **Click on** <https://codespaces.new/giorgiazaccaria/virtual-machine-cellGMM> to set-up your own version of the Codespace; **set** the Machine type option equal to **4-core**, and **click on** the green button named `Create codespace`.
   
![image](https://github.com/user-attachments/assets/6b2c4137-4ede-4950-ae28-14e7c89a6d83)

2. Wait for the Codespace to be created. This operation takes few minutes. At the end you should see something like the following screen.
   
4. Click on the `Ports` tab (which is highlighted in the previous image). Drag the mouse over the Forwarded Addresses field in the 8787 port and click on the globe icon (as displayed in the previous image). This will open a new tab in your browser with the Rstudio Server interface.
   
5. Login into the Rstudio server using the following credentials:
Username: `rstudio`
Password: `rstudio`
These operations will create an Rstudio session in your browser. Now you can explore the codes interactively!


> [!NOTE]
> GitHub Codespaces is paid for either by an organization, an enterprise, or a personal account. The Free and Pro plans for personal accounts include free use of GitHub Codespaces up to a fixed amount of usage every month.
 We refer to the <a href="https://docs.github.com/en/billing/managing-billing-for-your-products/managing-billing-for-github-codespaces/about-billing-for-github-codespaces">official docs</a> for more details.
