# virtual-machine-cellGMM
This repository is related to the paper _Cellwise outlier detection in heterogeneous populations_ by Giorgia Zaccaria, Luis A. García-Escudero, Francesca Greselin, and Agustín Mayo-Íscar available at https://www.tandfonline.com/doi/full/10.1080/00401706.2025.2497822. It contains the codes for reproducing some of the analyses presented in the paper via a Github Codespace. Specifically, the codes are the following:

- **Codes for running cellGMM:** "cellGMM.R" (main script). 
   This script automatically calls the following files:  
    - "InitializationFunctions_cellGMM.R"; 
    - "InternalFunctions_cellGMM.R".
     
- **Code for reproducing Figures 1 and 2 of the paper:** "Tech-Figs1-2.Rmd". 
     This R Markdown file automatically uploads "Data-Figs1-2-Complete.RData", which contains the results of the analyses reported in Section 3. 
    
- **Code for reproducing the analyses in Section 3 of the paper, Figures 1 and 2 and Table 1 for Scenario 1:** "Tech-Figs1-2-Tab1-Scenario1.Rmd". 
   This R Markdown file uses "snipEM_1.0.1.tar.gz" and "MixtureMissing_1.0.2.tar.gz" files if the user has not these packages already installed or a different version for *MixtureMissing*.

- **Code for reproducing the analysis in Section 4.2 of the paper and Figure 5:** "Tech-Fig5.Rmd". 
   This R Markdown file compares `cellGMM` with `cellMCD` and `DI`. "Carina Nebula_Reduced.png" is uploaded.

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

> [!IMPORTANT]
> Please note that the code was run on machines with 12 cores, whereas the GitHub Codespace has only 4 cores. Therefore, you should expect longer runtimes when using the Github Codespace.

> [!NOTE]
> GitHub Codespaces is paid for either by an organization, an enterprise, or a personal account. The Free and Pro plans for personal accounts include free use of GitHub Codespaces up to a fixed amount of usage every month.
 We refer to the <a href="https://docs.github.com/en/billing/managing-billing-for-your-products/managing-billing-for-github-codespaces/about-billing-for-github-codespaces">official docs</a> for more details.
