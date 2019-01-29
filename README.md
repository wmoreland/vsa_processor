# Vesicle Size Analysis Processor

This program will take 2D vesicle area measurements (e.g. output from [ImageJ](https://imagej.nih.gov/ij/)) and calculate 3D vesicle number densities, an important characteristic of multi-phase magma.

## Usage

Currently you must manually enter the paths to your input files, sample name, and clast vesicularity on lines 341 to 344 before running the program.

## Requirements
* [Python 3.7](https://www.python.org/)
* [Matplotlib](https://matplotlib.org/users/installing.html)
* [pandas](https://pandas.pydata.org/getpandas.html)

## Input format

The required inputs are two CSV files of image data and vesicle areas, clast vesicularity, and sample number. The column headings for the CSV files are described below and example data is provided with the program.

#### Vesicle areas file
This must be a CSV file with one column per image and the measured vesicle areas forming the rows. The column headings are the names of each image and must be without spaces. For the scanned images only the titles "scan" or "billet" are accepted. The remainder take the form like: 25a, 25b, 100aa, 100ab, 250aac, 250bba, etc. As long as the magnification goes first it does not matter what follows.

#### Image data file

This must be a CSV file with one column per parameter as follows:

| Column heading | Description                                                 | Data type |
|----------------|-------------------------------------------------------------|-----------|
| image          | name of image (see areas file description)                  | string    |
| scale_factor   | scale in pixels per millimetre                              | float     |
| width          | width of image in pixels                                    | integer   |
| height         | height of image in pixels                                   | integer   |
| edge_grey      | edge vesicles mean grey-scale value from ImageJ<sup>*</sup> | integer   |
| crystals_grey  | crystals mean grey-scale value from ImageJ<sup>*</sup>      | integer   |
| vesicles_grey  | vesicles mean grey-scale value from ImageJ<sup>*</sup>      | integer   |
<sup>* </sup>Documentation for this is in preparation

## References

This method is based on the following research which is cited in relevant sections of the code:

Sahagian DL, Proussevitch AA (1998) 3D particle size distributions from 2D observations: stereology for natural applications. Journal of Volcanology and Geothermal Research 84:173–196. doi: [10/fwnv97](https://doi.org/10/fwnv97)

Shea T, Houghton BF, Gurioli L, Cashman KV, Hammer JE, Hobden BJ (2010) Textural studies of vesicles in volcanic rocks: An integrated methodology. Journal of Volcanology and Geothermal Research 190:271–289. doi: [10/bmn9xt](https://doi.org/10/bmn9xt)
