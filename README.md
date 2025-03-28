# virtual-machine-cellGMM
This repository is related to the paper _Cellwise outlier detection in heterogeneous populations_ by Giorgia Zaccaria, Luis A. García-Escudero, Francesca Greselin, and Agustín Mayo-Íscar. It contains the codes for reproducing some of the analyses presented in the paper via a Github Codespace.

Specifically, the relevant codes for reproducing the analyses are the following:

- "Code for Reproducibility Simulation Figure Technometrics_Fast.R" $\rightarrow$ This code generates Figures 1 and 2 of the Main Article for all three scenarios. The code starts by loading the file "Data for Figure Reproducibility.RData", which contains the results for the three scenarios presented in the simulation study of the Main Article.

-  "Code for Reproducibility Simulation Technometrics_Scen1_5out_0mis.R" $\rightarrow$ This is a sub-example of the simulation study for Scenario 1 with 5% of cellwise contamination. At the end of the script, there is a code for reproducing the part of Table 1 in the Main Article corresponding to Scenario 1 with $5\%$ outlying values.

> [!IMPORTANT]
> The number of samples is currently set to $10$ in the R script, although the simulations in the paper were run on $100$ samples. Increasing the number of samples increases the computational time required to execute the code. It is worth noting that the simulations were run on virtual machines with $12$ cores, whereas the one created below has only $4$ cores. The user can change the number of samples up to $100$. With $10$ samples, it runs on the GitHub Codespace in approximately 16 minutes.

-  "Code for Reproducibility Simulation and Figure Technometrics_EntireScen1.html" $\rightarrow$ This .html file provides the complete simulation study for Scenario 1 with cellwise contamination levels of $0\%$, $5\%$, and $10\%$. It is used to reproduce Figures 1 and 2 of the Main Article for Scenario 1 and runs in approximately 1 hour on a laptop with $12$ CPU cores ($10$ of which were used for the implementation of the code).

We have prepared a **Github Codespace** that can be used to automatically set-up an Rstudio Server and explore the codes in an interactive way. Instructions to access the Github Codespace are reported below.

> [!WARNING]
> Please note that you need a Github account to access the Github Codespace.

1. **Click** the following link using `Ctrl + Click` (Windows) or `Cmd + Click` (Mac) to set-up your own version of the Codespace in a new tab. 

<a href="https://codespaces.new/giorgiazaccaria/virtual-machine-cellGMM" target="_blank" rel="noopener noreferrer">https://codespaces.new/giorgiazaccaria/virtual-machine-cellGMM</a>

2.  **Set** the Machine type option equal to **4-core**, and **click on** the green button named `Create codespace` (see the pictures below).

![image](https://github.com/user-attachments/assets/6b2c4137-4ede-4950-ae28-14e7c89a6d83)

3. Wait for the Codespace to be created. This process may take a few minutes. Wait until the number 2 (or greater than 2) appears next to `Ports`. You should see a screen similar to the one below.
   
![image](https://github.com/user-attachments/assets/f97af82d-1ac5-4307-bb67-0b54c795cb7e)
   
5. **Click on** the `Ports` tab, which is highlighted in the previous image. Drag the mouse over the Forwarded Addresses field in the **8787 port** and click on the **globe icon**, as displayed in the following image.
   
![image](https://github.com/user-attachments/assets/7dd1b898-e1e3-4a62-9db5-acdfc25f4610)

6. The previous action will open a new tab in your browser with the Rstudio Server interface, as follows.

![image](https://github.com/user-attachments/assets/759db5a3-dbcc-4533-b63b-80fbfe79effb)
   
7. **Login** into the Rstudio server using the following credentials:
Username: `rstudio`
Password: `rstudio`
These operations will create an Rstudio session in your browser. When the Rstudio session is ready, **click on** the file "virtual-machine-cellGMM" (bottom-right panel in R Studio). Now you can explore the codes interactively!

> [!NOTE]
> GitHub Codespaces is paid for either by an organization, an enterprise, or a personal account. The Free and Pro plans for personal accounts include free use of GitHub Codespaces up to a fixed amount of usage every month.
 We refer to the <a href="https://docs.github.com/en/billing/managing-billing-for-your-products/managing-billing-for-github-codespaces/about-billing-for-github-codespaces">official docs</a> for more details.
