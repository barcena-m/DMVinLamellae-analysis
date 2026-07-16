# DMV in lamellae analysis
Python script for the analysis of double-membrane vesicles (DMVs) and DMV-pore complexes in cryo-lamellae.

The script:
- Calculates DMV diameter from pore and DMV center coordinates.
- Fits upper and lower lamella planes from user-defined points.
- Calculates lamella thickness.
- Estimates the fraction of DMV surface embedded within the lamella.
- Corrects pore counts based on the embedded DMV surface fraction.
- Calculates pore density.
- Exports results to text and Excel formats.


## Requirements
Python 3.x

Required packages:
- numpy
- openpyxl

## Running the script
Run:
```bash
python dmv_in_lamellae_analysis.py
```

The program will prompt for:
- the input list filename
- the output filename

## Input file formats

The analysis requires:

1. An input list file containing calibration parameters and the names of the input files.
2. Vesicle coordinate files.
3. Lamella plane coordinate files.

Example files are provided in the `example_data` directory.

**Vesicle coordinate file:** DMV ID, xyz coordinates.
For each DMV, all coordinates except the final one are interpreted as pore coordinates, while the final coordinate is interpreted as the DMV centre coordinate. These files can be generated from an IMOD model using `model2point -contour`.

**Lamella plane coordinate file:** plane ID, point ID, xyz coordinates.
Points assigned to plane ID 1 are used to calculate the upper lamella plane, while points assigned to plane ID 2 are used to calculate the lower lamella
plane. At least three non-collinear points are required for each plane; additional points may be provided to improve the plane fit. The planes are
calculated using a best-fit approach based on all points assigned to each plane.
These files can be generated from an IMOD model using `model2point -contour -object`.

## Output

The script generates three output files:
- A text summary file (`.txt`) containing all measurements in tabulated format.
- A log file (`.log`) containing processing information and calculations.
- An Excel workbook (`.xlsx`) containing the same measurements in spreadsheet format.
