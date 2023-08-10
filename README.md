# SURF: a simple user-friendly reaction format

This repository containes example files and code to read, write and transform **SURF files** containing **chemical reaction data**.

## Installation
First, clone this repository and then install the necessary requirements by running the following command (assuming you have `pip` installed):
```bash
pip install -r requirements.txt
```
Now you should be ready to go.

## SURF Structure
SURF is a tabular file structure that is both human- and machine-readable. In a SURF spreadsheet, each row stores data of one reaction. The column headers structure the data and are split into constant (CC) and flexible (FC) categories. CCs never change and should be always present, independent of the number of reaction components. They capture the identifiers and provenance of the reaction as well as basic characteristics (reaction type, named reaction, reaction technology) and conditions (temperature, time, atmosphere, scale, concentration, stirring/shaking). Add-ons, such as the procedure or comments, belong to the CCs as well. The FCs describe the more variable part of a reaction, the different starting material(s), solvent(s), reagent(s) and product(s). Each reaction component is described by at an identifier such as the CAS number or molecule name, and a SMILES or InChI string storing the chemcial structure. While the SMILES/InChI string is available for every compound and can also serve as structural input for machine learning models, the CAS number, even though not always available, can be handy for chemists in the lab to order, itemize and find chemicals. For the starting material(s) and reagent(s), e.g. catalyst, ligand, additive, in addition to the two identifiers, a third column is added to cover the stochiometric amount (equivalents). For products, the respective yield and yield type is referenced. The flexibility of SURF allows capturing multiple starting materials and reagents, as these can be accommodated by adding three additional columns with (CAS, SMILES/InChI, and equivalents). If desired, further columns for additional identifiers like names or lot numbers can be added.

The file `data/surf_template.txt` provides an example of the SURF structure with five minisci-type reactions from literature.

## Combining Multiple SURF Files
To concatenate multiple SURF files into one larger file, put all SURF files to be combined into one folder and run the following:
```bash
python concat_surf.py <folder path> <output file>
```

#### Example:
```bash
ls data/hte-data
    hte-1.txt hte-2.txt hte-3.txt hte-4.txt ...
python concat_surf.py data/hte-data data/hte-data.txt
```

## Transformations for Interoperability
The idea of SURF is not to replace existing reaction file structures, but have a human- and machine-readable format that is interoperable. Below we describe examples of how to transform SURF into other existing formats and back.

### SURF to ORD
The [Open Reaction Database (ORD)](https://pubs.acs.org/doi/10.1021/jacs.1c09820) is an open-access schema and infrastructure for structuring and sharing organic reaction data. Translating SURF files into the protocol buffers format used by the open reaction database, run the following:
```bash
python surf2ord.py <input SURF file> <output ORD file>
```
To get all the options of the script, run:
```bash
python surf2ord.py --help
```

If the SURF file does not contain any or only partial provenance information, the user can provide personal information with the `--username`, `--email`, `--orcid` and `--organization` options.

#### Example:
```bash
python surf2ord.py data/surf_template.txt data/surf_template.pbtxt --username "Alex Mueller" --email "alex.mueller@roche.com"
```

### ORD to SURF
To translate protocol buffers files used by the open reaction database back into the SURF format, run the following:
```bash
python ord2surf.py <input ORD file> <output SURF file>
```
To get all the options of the script, run:
```bash
python ord2surf.py --help
```

#### Example:
```bash
python ord2surf.py data/ord_search_result.pb data/ord_search_result.txt --validate
```

### SURF to Reaction SMILES
Reaction SMILES are a frequently used representation for chemical reactions. However, they just represent the molecular structures involved in the reaction and lack detailed information on conditions, equivalents and analytics. We therefore only provide a way to translate SURF files into Reaction SMILES but not vice versa:
```bash
python surf2rxnsmiles.py <input SURF file> <output RXNSMILES file>
```
To get all the options of the script, run:
```bash
python surf2rxnsmiles.py --help
```

#### Example:
```bash
python surf2rxnsmiles.py data/surf_template.txt data/surf_template.rxnsmi
```

### SURF to UDM
The [Unified Data Model (UDM)](https://doi.org/10.1515/pac-2021-3013) is an open, extendable and freely available data format for the exchange of experimental information about compound synthesis and testing, developed by the [Pistoia Alliance](https://www.pistoiaalliance.org/projects/current-projects/unified-data-model/). To translate SURF files into UDM XML files, run the following: 

```bash
python surf2udm.py <input SURF file> <output UDM file>
```
To get all the options of the script, run:
```bash
python surf2udm.py --help
```

#### Example:
```bash
python surf2udm.py data/surf_template.txt data/surf_template.xml
```

### UDM to SURF
To translate UDM XML files into SURF files, run the following: 
```bash
python udm2surf.py <input UDM file> <output SURF file>
```
To get all the options of the script, run:
```bash
python udm2surf.py --help
```

#### Example:
```bash
python udm2surf.py data/udm_file.xml data/surf_file.txt
```

## Citation
If you are using SURF in your project, please cite the following reference:
```
@article{nippa_mueller2023surf,
  title={SURF - a simple, user-friendly reaction format},
  author={Nippa, David F. and M{\"u}ller, Alex T. and Atz, Kenneth and Konrad, David B. and Grether, Uwe and Martin, Rainer E. and Schneider, Gisbert},
  journal={},
  year={}
}

# ToDo: To be updated with DOI and journal.
```
