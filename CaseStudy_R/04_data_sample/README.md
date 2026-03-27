# Data Sample: MoBPS Simulation & GBLUP Results

This folder contains sample datasets generated using the **MoBPS** package.

## Population Abbreviations

| ID | Description |
| :--- | :--- |
| **Np1** | Parent population; phenotypically high-ranked trees from the foundation population. |
| **No** | Offspring; Candidate population for selection. |
| **Np2** | Pollen contaminants (External genetic contribution). |
| **Nr** | Randomly sampled trees from the foundation population. |
| **rep0** | Iteration index of the dataset. |

---

## File Descriptions

### 1. `rep0_GBLUP_results`
Contains the GBLUP estimation results for the populations.
* **Rows 2 – 101**: Np1 (Parents)
* **Rows 102 – 601**: No (Offspring)
* **Rows 602 – 1501**: Nr (Random samples)

### 2. `rep0_Gmatrix_tuned_??`
The Genomic Relationship Matrix (G-matrix) for all individuals and **Np2**. 
> *Note: These parameters were tuned for optimization compatibility.*

### 3. `rep0_??_data`
Contains the reference values for validation:
* **Pheno**: True phenotypes.
* **BV**: True breeding values.

---

## Simulation Settings
* **Heritability ($h^2$)**: 0.2
* **Model**: GBLUP implemented with the founder population as the reference.
